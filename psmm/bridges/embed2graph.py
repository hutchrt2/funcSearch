import argparse
import os
import glob
import json
import torch
import faiss
import numpy as np
import pandas as pd
from datetime import datetime
from Bio import SeqIO
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

def get_device(device_override=None):
    if device_override:
        return device_override
    if torch.cuda.is_available():
        return "cuda"
    elif torch.backends.mps.is_available():
        return "mps"
    return "cpu"

def load_model(model_name="esmc_300m", device="cpu"):
    print(f"Loading model {model_name} on {device}...")
    try:
        model = ESMC.from_pretrained(model_name, device=torch.device(device))
        model.eval()
        return None, model
    except Exception as e:
        raise RuntimeError(f"Failed to load model {model_name}: {e}")

def embed_sequences(sequences, tokenizer, model, device="cpu", batch_size=8):
    embeddings = []
    total_seqs = len(sequences)
    total_batches = (total_seqs + batch_size - 1) // batch_size
    
    for batch_idx, i in enumerate(range(0, total_seqs, batch_size)):
        batch_seqs = sequences[i:i + batch_size]
        
        batch_embeddings = []
        for seq in batch_seqs:
            protein = ESMProtein(sequence=seq)
            protein_tensor = model.encode(protein)
            
            with torch.no_grad():
                logits_output = model.logits(
                    protein_tensor, 
                    LogitsConfig(sequence=True, return_embeddings=True)
                )
                # logits_output.embeddings has shape (1, seq_len, hidden_size)
                # We perform mean pooling over the sequence dimension (dim=1)
                token_embeddings = logits_output.embeddings
                mean_pooled = token_embeddings.mean(dim=1) # Shape: (1, hidden_size)
                
                # L2 Normalize for cosine similarity calculation
                mean_pooled = torch.nn.functional.normalize(mean_pooled, p=2, dim=1)
                
                # Convert back to float32 to prevent numeric type mismatch downstream (e.g., in FAISS)
                batch_embeddings.append(mean_pooled.to(torch.float32).cpu().numpy())
                
        if batch_embeddings:
            embeddings.append(np.vstack(batch_embeddings))
            
        if (batch_idx + 1) % 5 == 0 or batch_idx + 1 == total_batches:
            print(f"  Processed batch {batch_idx + 1}/{total_batches} ({min((batch_idx + 1) * batch_size, total_seqs)}/{total_seqs} sequences)...")
            
    if not embeddings:
        return np.array([])
    return np.vstack(embeddings)

def init_database(args):
    print("Phase A: Initializing database...")
    fasta_path = args.db_fasta
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA not found at {fasta_path}")
        
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No valid sequences found in reference FASTA at {fasta_path}")
        
    sequences = [str(r.seq) for r in records]
    uniprot_ids = [r.id.split('|')[1] if '|' in r.id else r.id for r in records]
    
    # Define cache paths
    os.makedirs(args.db_dir, exist_ok=True)
    index_path = os.path.join(args.db_dir, "faiss_index.bin")
    emb_path = os.path.join(args.db_dir, "esmc_embeddings.npy")
    ids_path = os.path.join(args.db_dir, "index_uniprot_ids.txt")
    meta_path = os.path.join(args.db_dir, "index_metadata.json")
    
    # Caching check
    use_cache = False
    existing_embeddings = None
    existing_ids = []
    
    if getattr(args, "clean", False):
        print("Cleaning old database files for a fresh start...")
        for p in [index_path, emb_path, ids_path, meta_path]:
            if os.path.exists(p):
                try:
                    os.remove(p)
                except Exception as e:
                    print(f"Warning: Failed to remove cache file {p}: {e}")
                    
    if not args.force:
        cache_files_exist = all(os.path.exists(p) for p in [index_path, emb_path, ids_path, meta_path])
        if cache_files_exist:
            try:
                with open(meta_path, "r") as f:
                    cache_metadata = json.load(f)
                
                # Verify that the model name matches
                if cache_metadata.get("model_name") == args.model:
                    existing_embeddings = np.load(emb_path)
                    with open(ids_path, "r") as f:
                        existing_ids = [line.strip() for line in f]
                    
                    # Ensure alignment check
                    if len(existing_ids) == existing_embeddings.shape[0]:
                        use_cache = True
                        print(f"Found existing database cache with {len(existing_ids)} embedded sequences (Model: {args.model}).")
                    else:
                        print("Warning: Cache file mismatch (IDs count vs embeddings count). Rebuilding index.")
                else:
                    print(f"Notice: Model changed from {cache_metadata.get('model_name')} to {args.model}. Rebuilding database from scratch.")
            except Exception as e:
                print(f"Warning: Failed to load cache metadata: {e}. Rebuilding database from scratch.")
    else:
        print("Force rebuild requested. Bypassing database cache.")
        
    # Process incrementally or full build
    if use_cache:
        existing_id_set = set(existing_ids)
        new_indices = []
        for i, uid in enumerate(uniprot_ids):
            if uid not in existing_id_set:
                new_indices.append(i)
                
        if not new_indices:
            print("All sequences in reference FASTA are already embedded and indexed. Skipping database initialization.")
            return
            
        print(f"Incremental update: {len(new_indices)} new sequences detected to embed.")
        new_seqs = [sequences[i] for i in new_indices]
        new_uids = [uniprot_ids[i] for i in new_indices]
        
        device = get_device(args.device)
        tokenizer, model = load_model(args.model, device=device)
        print(f"Embedding {len(new_seqs)} new sequences...")
        new_embeddings = embed_sequences(new_seqs, tokenizer, model, device=device, batch_size=args.batch_size)
        
        # Combine
        embeddings = np.vstack([existing_embeddings, new_embeddings])
        combined_ids = existing_ids + new_uids
    else:
        device = get_device(args.device)
        tokenizer, model = load_model(args.model, device=device)
        print(f"Embedding all {len(sequences)} reference sequences...")
        embeddings = embed_sequences(sequences, tokenizer, model, device=device, batch_size=args.batch_size)
        combined_ids = uniprot_ids

    # Save outputs and build FAISS index
    print(f"Building FAISS index with {embeddings.shape[0]} sequences...")
    d = embeddings.shape[1]
    index = faiss.IndexFlatIP(d)
    
    try:
        index.add(embeddings)
        faiss.write_index(index, index_path)
        np.save(emb_path, embeddings)
        with open(ids_path, "w") as f:
            for uid in combined_ids:
                f.write(f"{uid}\n")
                
        # Write metadata JSON
        metadata = {
            "model_name": args.model,
            "dimension": d,
            "last_updated": datetime.now().isoformat(),
            "sequence_count": len(combined_ids)
        }
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2)
            
    except Exception as e:
        raise RuntimeError(f"Failed to create or save FAISS database index: {e}")
        
    print(f"Initialization complete. Total database size: {len(combined_ids)} sequences. Files saved in '{args.db_dir}'.")

def query_database(args):
    print("Phase B: Query Embedding & Vector Search...")
    index_path = os.path.join(args.db_dir, "faiss_index.bin")
    id_path = os.path.join(args.db_dir, "index_uniprot_ids.txt")
    if not os.path.exists(index_path) or not os.path.exists(id_path):
        raise FileNotFoundError(f"Database files not found in '{args.db_dir}'. Run with --init first.")
        
    try:
        index = faiss.read_index(index_path)
    except Exception as e:
        raise RuntimeError(f"Failed to load FAISS index: {e}")
        
    with open(id_path, "r") as f:
        uniprot_ids = [line.strip() for line in f]
        
    device = get_device(args.device)
    tokenizer, model = load_model(args.model, device=device)
    
    meta_path = args.metadata
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata CSV not found at {meta_path}")
    metadata_df = pd.read_csv(meta_path)
    if "target_accession" in metadata_df.columns:
        metadata_df = metadata_df.rename(columns={"target_accession": "uniprot_id"})
    if "uniprot_id" not in metadata_df.columns:
        raise ValueError(f"Metadata CSV at {meta_path} is missing the required 'uniprot_id' column (or 'target_accession' to map to it).")
    
    query_files = []
    if os.path.isdir(args.query):
        query_files = sorted(glob.glob(os.path.join(args.query, "*.fasta")) + glob.glob(os.path.join(args.query, "*.fa")))
    else:
        if not os.path.exists(args.query):
            raise FileNotFoundError(f"Query path not found: {args.query}")
        query_files = [args.query]
        
    all_results = []
    
    for q_file in query_files:
        print(f"Processing query file: {q_file}")
        try:
            records = list(SeqIO.parse(q_file, "fasta"))
        except Exception as e:
            print(f"Warning: Failed to parse {q_file} as FASTA: {e}. Skipping.")
            continue
            
        if not records:
            print(f"Warning: No sequences found in {q_file}. Skipping.")
            continue
            
        q_sequences = [str(r.seq) for r in records]
        q_ids = [r.id for r in records]
        
        q_embeddings = embed_sequences(q_sequences, tokenizer, model, device=device, batch_size=args.batch_size)
        
        k = min(args.k, index.ntotal)
        distances, indices = index.search(q_embeddings, k)
        
        # Phase C: Vector Similarity Results
        for i, (dists, idxs) in enumerate(zip(distances, indices)):
            q_id = q_ids[i]
            for d, idx in zip(dists, idxs):
                if idx == -1: continue
                t_id = uniprot_ids[idx]
                cosine_sim = float(d)
                
                all_results.append({
                    "query": q_id,
                    "target": t_id,
                    "score": cosine_sim,
                    "score_type": "cosine_similarity"
                })

                
    if not all_results:
        print("No search results found.")
        return
        
    results_df = pd.DataFrame(all_results)
    
    # Phase D: The Relational Handshake (The Join)
    print("Phase D: Joining with knowledge graph metadata...")
    
    # Check if target corresponds to global_node_id or uniprot_id
    sample_targets = results_df["target"].dropna().unique()
    is_global_node_id = any(val in metadata_df["global_node_id"].values for val in sample_targets)

    if is_global_node_id:
        print("Detected global_node_id targets. Performing join on global_node_id...")
        results_df = results_df.rename(columns={"target": "global_node_id"})
        joined_df = pd.merge(results_df, metadata_df, on="global_node_id", how="inner")
    else:
        print("Detected standard target IDs. Performing join on uniprot_id...")
        joined_df = pd.merge(results_df, metadata_df, left_on="target", right_on="uniprot_id", how="inner")
    
    if joined_df.empty:
        print("=" * 60)
        print("WARNING: The join with metadata resulted in an empty DataFrame.")
        print("Troubleshooting checks:")
        print(" 1. Ensure reference FASTA headers yield target IDs matching metadata 'uniprot_id' or 'global_node_id'.")
        sample_targets = list(sample_targets[:5])
        sample_meta = list(metadata_df["uniprot_id"].unique()[:5]) if "uniprot_id" in metadata_df.columns else []
        sample_nodes = list(metadata_df["global_node_id"].unique()[:5]) if "global_node_id" in metadata_df.columns else []
        print(f"    Sample target IDs extracted: {sample_targets}")
        print(f"    Sample metadata uniprot_ids: {sample_meta}")
        print(f"    Sample metadata global_node_ids: {sample_nodes}")
        print(" 2. Make sure they match exactly (case-sensitive, no extra whitespace).")
        print("=" * 60)
        
    # Phase E: Output Serialization
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    joined_df.to_csv(args.output, index=False)
    print(f"Phase E: Results saved to {args.output} ({len(joined_df)} entries joined)")

def main():
    parser = argparse.ArgumentParser(description="Embed2Graph Bridge - Protein Sequence Vector Search")
    parser.add_argument("--init", action="store_true", help="Initialize the database by embedding reference sequences")
    parser.add_argument("--query", type=str, help="Path to a query FASTA file or directory containing FASTA files")
    parser.add_argument("--output", type=str, default="output/vector_query_results.csv", help="Path to output CSV")
    parser.add_argument("--model", type=str, default="esmc_300m", help="ESM-C model name (e.g., esmc_300m, esmc_600m)")
    parser.add_argument("--k", type=int, default=5, help="Number of nearest neighbors to retrieve")
    parser.add_argument("--batch-size", type=int, default=8, help="Batch size for model inference")
    parser.add_argument("--device", type=str, default=None, choices=["cuda", "mps", "cpu"], help="Specify compute device")
    parser.add_argument("--force", action="store_true", help="Force rebuild of reference database, bypassing cache")
    parser.add_argument("--clean", action="store_true", help="Delete old database files to force a clean start")
    
    project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    # Configurable database/metadata input options
    parser.add_argument("--db-fasta", type=str, default=os.path.join(project_dir, "data", "build", "psfd_sequences.fasta"), help="Path to reference FASTA database")
    parser.add_argument("--metadata", type=str, default=os.path.join(project_dir, "data", "build", "sequence_metadata.csv"), help="Path to sequence metadata CSV")
    parser.add_argument("--db-dir", type=str, default=os.path.join(project_dir, "data", "embeddb"), help="Directory to save/load vector database index")
    
    args = parser.parse_args()
    
    if args.init:
        init_database(args)
    if args.query:
        query_database(args)
    if not args.init and not args.query:
        parser.print_help()
        
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
