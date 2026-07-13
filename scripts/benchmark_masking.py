#!/usr/bin/env python3
import os
import sys
import glob
import shutil
import hashlib
import argparse
import subprocess
import pandas as pd
import numpy as np

try:
    from Bio import SeqIO
except ImportError:
    print("Biopython not found. Please install it using: pip install biopython")
    sys.exit(1)

try:
    import faiss
except ImportError:
    print("FAISS not found. Assuming you don't need vector masking, or please install it: pip install faiss-cpu")

def get_sequence_hash(sequence: str) -> str:
    import re
    clean_seq = re.sub(r'\s+', '', sequence.upper())
    return hashlib.sha256(clean_seq.encode('utf-8')).hexdigest()

def backup_file_or_dir(path: str):
    if not os.path.exists(path):
        return
    backup_path = path + ".backup"
    if os.path.exists(backup_path):
        print(f"Backup already exists: {backup_path}")
        return
    print(f"Backing up {path} to {backup_path}...")
    if os.path.isdir(path):
        shutil.copytree(path, backup_path)
    else:
        shutil.copy2(path, backup_path)

def restore_file_or_dir(path: str):
    backup_path = path + ".backup"
    if not os.path.exists(backup_path):
        print(f"No backup found for {path}")
        return
    print(f"Restoring {path} from {backup_path}...")
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    if os.path.isdir(backup_path):
        shutil.copytree(backup_path, path)
    else:
        shutil.copy2(backup_path, path)

def get_masked_ids_and_hashes(target_dir: str):
    fasta_files = glob.glob(os.path.join(target_dir, "*.fasta")) + glob.glob(os.path.join(target_dir, "*.fa"))
    if not fasta_files:
        print(f"No FASTA files found in {target_dir}")
        return set(), set()
        
    masked_ids = set()
    masked_hashes = set()
    
    for f_path in fasta_files:
        for record in SeqIO.parse(f_path, "fasta"):
            # Try to extract UniProt ID
            header = record.id
            uid = header
            if "|" in header:
                parts = header.split("|")
                if len(parts) >= 2:
                    uid = parts[1]
            masked_ids.add(uid.upper())
            masked_hashes.add(get_sequence_hash(str(record.seq)))
            
    return masked_ids, masked_hashes

def mask_database(project_dir: str, target_dir: str, queries_dir: str):
    masked_ids, masked_hashes = get_masked_ids_and_hashes(target_dir)
    if not masked_ids:
        print("Nothing to mask.")
        return
        
    print(f"Identified {len(masked_ids)} sequences to mask.")
    print("Masked IDs:", masked_ids)
    
    # 1. Backups
    fasta_path = os.path.join(project_dir, "data", "build", "psfd_sequences.fasta")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    blastdb_dir = os.path.join(project_dir, "data", "blastdb")
    embeddb_dir = os.path.join(project_dir, "data", "embeddb")
    
    for path in [fasta_path, meta_path, blastdb_dir, embeddb_dir]:
        backup_file_or_dir(path)
        
    # 2. Mask FASTA
    print("Masking psfd_sequences.fasta...")
    if os.path.exists(fasta_path):
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            
        kept_lines = []
        skip_current = False
        current_seq_buffer = []
        current_header = ""
        
        for line in lines:
            if line.startswith(">"):
                if current_header:
                    # check if previous was masked
                    h_hash = get_sequence_hash("".join(current_seq_buffer))
                    if h_hash not in masked_hashes and not skip_current:
                        kept_lines.append(current_header)
                        kept_lines.extend(current_seq_buffer)
                
                current_header = line
                current_seq_buffer = []
                skip_current = False
                
                # Check header for ID
                uid = line.split()[0].replace(">", "")
                if "|" in line:
                    parts = line.split("|")
                    if len(parts) >= 2:
                        db_acc = parts[1].strip()
                        if ":" in db_acc:
                            uid = db_acc.split(":")[1].strip()
                        else:
                            uid = db_acc.strip()
                            
                if uid.upper() in masked_ids:
                    skip_current = True
            else:
                current_seq_buffer.append(line)
                
        # Last sequence
        if current_header:
            h_hash = get_sequence_hash("".join(current_seq_buffer))
            if h_hash not in masked_hashes and not skip_current:
                kept_lines.append(current_header)
                kept_lines.extend(current_seq_buffer)
                
        with open(fasta_path, "w") as f:
            f.writelines(kept_lines)
            
    # 3. Mask Metadata CSV
    print("Masking sequence_metadata.csv...")
    if os.path.exists(meta_path):
        df = pd.read_csv(meta_path)
        if "target_accession" in df.columns:
            initial_len = len(df)
            df = df[~df["target_accession"].str.upper().isin(masked_ids)]
            print(f"Removed {initial_len - len(df)} rows from metadata.")
            df.to_csv(meta_path, index=False)
            
    # 4. Mask Embeddings (esmc_embeddings.npy)
    print("Masking FAISS Vector Embeddings...")
    npy_path = os.path.join(embeddb_dir, "esmc_embeddings.npy")
    ids_path = os.path.join(embeddb_dir, "index_uniprot_ids.txt")
    idx_path = os.path.join(embeddb_dir, "faiss_index.bin")
    
    if os.path.exists(npy_path) and os.path.exists(ids_path):
        embeddings = np.load(npy_path)
        with open(ids_path, "r") as f:
            uids = [line.strip() for line in f]
            
        if len(uids) == embeddings.shape[0]:
            keep_indices = [i for i, uid in enumerate(uids) if uid.upper() not in masked_ids]
            
            print(f"Removing {len(uids) - len(keep_indices)} vectors from embeddings...")
            filtered_embeddings = embeddings[keep_indices]
            filtered_uids = [uids[i] for i in keep_indices]
            
            # Save new arrays
            np.save(npy_path, filtered_embeddings)
            with open(ids_path, "w") as f:
                for u in filtered_uids:
                    f.write(f"{u}\n")
                    
            # Rebuild FAISS index
            if "faiss" in sys.modules:
                d = filtered_embeddings.shape[1]
                index = faiss.IndexFlatIP(d)
                index.add(filtered_embeddings)
                faiss.write_index(index, idx_path)
                
                # Update json count if exists
                import json
                meta_json = os.path.join(embeddb_dir, "index_metadata.json")
                if os.path.exists(meta_json):
                    with open(meta_json, "r") as f:
                        j = json.load(f)
                    j["sequence_count"] = len(filtered_uids)
                    with open(meta_json, "w") as f:
                        json.dump(j, f, indent=2)
                print("FAISS index rebuilt successfully.")
            else:
                print("FAISS not imported, skipped rebuilding .bin index.")
                
    # 5. Anonymize Queries
    print("Generating anonymized query files...")
    os.makedirs(queries_dir, exist_ok=True)
    fasta_files = glob.glob(os.path.join(target_dir, "*.fasta")) + glob.glob(os.path.join(target_dir, "*.fa"))
    counter = 1
    for f_path in fasta_files:
        for record in SeqIO.parse(f_path, "fasta"):
            out_path = os.path.join(queries_dir, f"Query_{counter}.fasta")
            with open(out_path, "w") as f:
                f.write(f">Benchmark_Test_Protein_{counter}\n")
                f.write(str(record.seq) + "\n")
            counter += 1
            
    # 6. Rebuild MMseqs index
    print("Rebuilding MMseqs database index...")
    script_path = os.path.join(project_dir, "psmm", "bridges", "seq2graph.py")
    if os.path.exists(script_path):
        subprocess.run([sys.executable, script_path, "--init", "--clean"], cwd=os.path.dirname(script_path))
    else:
        print(f"Warning: seq2graph.py not found at {script_path}")
        
    print("\nMasking complete! You can now run your anonymized queries from the 'data/benchmark_queries' folder.")
def unmask_database(project_dir: str):
    print("Restoring original unmasked database state...")
    fasta_path = os.path.join(project_dir, "data", "build", "psfd_sequences.fasta")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    blastdb_dir = os.path.join(project_dir, "data", "blastdb")
    embeddb_dir = os.path.join(project_dir, "data", "embeddb")
    
    for path in [fasta_path, meta_path, blastdb_dir, embeddb_dir]:
        restore_file_or_dir(path)
        
    print("Unmasking complete! Normal database restored.")

def mask_single(project_dir: str, uid: str, h_hash: str):
    fasta_path = os.path.join(project_dir, "data", "build", "psfd_sequences.fasta")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    blastdb_dir = os.path.join(project_dir, "data", "blastdb")
    embeddb_dir = os.path.join(project_dir, "data", "embeddb")
    
    fasta_backup = fasta_path + ".backup"
    meta_backup = meta_path + ".backup"
    
    # 1. Mask FASTA (from pristine backup)
    if os.path.exists(fasta_backup):
        with open(fasta_backup, "r") as f:
            lines = f.readlines()
            
        kept_lines = []
        skip_current = False
        current_seq_buffer = []
        current_header = ""
        
        for line in lines:
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq_buffer)
                    curr_hash = get_sequence_hash(seq_str)
                    if curr_hash != h_hash and not skip_current:
                        kept_lines.append(current_header)
                        kept_lines.extend(current_seq_buffer)
                
                current_header = line
                current_seq_buffer = []
                skip_current = False
                
                # Check header for ID
                curr_uid = line.split()[0].replace(">", "")
                if "|" in line:
                    parts = line.split("|")
                    if len(parts) >= 2:
                        db_acc = parts[1].strip()
                        if ":" in db_acc:
                            curr_uid = db_acc.split(":")[1].strip()
                        else:
                            curr_uid = db_acc.strip()
                            
                if curr_uid.upper() == uid.upper():
                    skip_current = True
            else:
                current_seq_buffer.append(line)
                
        # Last sequence
        if current_header:
            seq_str = "".join(current_seq_buffer)
            curr_hash = get_sequence_hash(seq_str)
            if curr_hash != h_hash and not skip_current:
                kept_lines.append(current_header)
                kept_lines.extend(current_seq_buffer)
                
        with open(fasta_path, "w") as f:
            f.writelines(kept_lines)

    # 2. Mask Metadata CSV (from pristine backup)
    if os.path.exists(meta_backup):
        df = pd.read_csv(meta_backup)
        if "target_accession" in df.columns:
            df = df[df["target_accession"].str.upper() != uid.upper()]
            df.to_csv(meta_path, index=False)
            
    # 3. Mask FAISS Vector Embeddings (from pristine backup)
    npy_backup = os.path.join(embeddb_dir + ".backup", "esmc_embeddings.npy")
    ids_backup = os.path.join(embeddb_dir + ".backup", "index_uniprot_ids.txt")
    meta_json_backup = os.path.join(embeddb_dir + ".backup", "index_metadata.json")
    
    npy_path = os.path.join(embeddb_dir, "esmc_embeddings.npy")
    ids_path = os.path.join(embeddb_dir, "index_uniprot_ids.txt")
    idx_path = os.path.join(embeddb_dir, "faiss_index.bin")
    meta_json = os.path.join(embeddb_dir, "index_metadata.json")
    
    os.makedirs(embeddb_dir, exist_ok=True)
    
    if os.path.exists(npy_backup) and os.path.exists(ids_backup):
        embeddings = np.load(npy_backup)
        with open(ids_backup, "r") as f:
            uids = [line.strip() for line in f]
            
        if len(uids) == embeddings.shape[0]:
            keep_indices = [i for i, u in enumerate(uids) if u.upper() != uid.upper()]
            filtered_embeddings = embeddings[keep_indices]
            filtered_uids = [uids[i] for i in keep_indices]
            
            np.save(npy_path, filtered_embeddings)
            with open(ids_path, "w") as f:
                for u in filtered_uids:
                    f.write(f"{u}\n")
                    
            if "faiss" in sys.modules:
                d = filtered_embeddings.shape[1]
                index = faiss.IndexFlatIP(d)
                index.add(filtered_embeddings)
                faiss.write_index(index, idx_path)
                
                if os.path.exists(meta_json_backup):
                    import json
                    with open(meta_json_backup, "r") as f:
                        j = json.load(f)
                    j["sequence_count"] = len(filtered_uids)
                    with open(meta_json, "w") as f:
                        json.dump(j, f, indent=2)
                        
    # 4. Rebuild MMseqs index
    script_path = os.path.join(project_dir, "psmm", "bridges", "seq2graph.py")
    if os.path.exists(script_path):
        subprocess.run([sys.executable, script_path, "--init", "--clean"], cwd=os.path.dirname(script_path), capture_output=True)

def benchmark_query(project_dir: str, method: str, query_header: str, query_seq: str) -> pd.DataFrame:
    tmp_dir = os.path.join(project_dir, "data", "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    query_file = os.path.join(tmp_dir, "query_temp.fasta")
    output_csv = os.path.join(tmp_dir, "output_temp.csv")
    
    if os.path.exists(output_csv):
        os.remove(output_csv)
        
    with open(query_file, "w") as f:
        f.write(f">{query_header}\n")
        f.write(query_seq + "\n")
        
    if method == "embed2graph":
        module_name = "psmm.bridges.embed2graph"
    elif method == "seq2graph":
        module_name = "psmm.bridges.seq2graph"
    else:
        return pd.DataFrame()
        
    subprocess.run(
        [sys.executable, "-m", module_name, "--query", query_file, "--output", output_csv],
        cwd=project_dir,
        capture_output=True,
        text=True
    )
    
    if os.path.exists(query_file):
        os.remove(query_file)
        
    if os.path.exists(output_csv):
        df = pd.read_csv(output_csv)
        os.remove(output_csv)
        return df
    return pd.DataFrame()

def run_benchmark(project_dir: str, folder_path: str, output_path: str, methods: list):
    fasta_files = glob.glob(os.path.join(folder_path, "*.fasta")) + glob.glob(os.path.join(folder_path, "*.fa")) + glob.glob(os.path.join(folder_path, "*.faa"))
    if not fasta_files:
        print(f"No FASTA files found in {folder_path}")
        return
        
    all_sequences = []
    for f_path in fasta_files:
        filename = os.path.basename(f_path)
        for record in SeqIO.parse(f_path, "fasta"):
            header = record.description
            seq_str = str(record.seq)
            uid = record.id
            if "|" in record.id:
                parts = record.id.split("|")
                if len(parts) >= 2:
                    uid = parts[1]
            all_sequences.append({
                "file_name": filename,
                "header": header,
                "sequence": seq_str,
                "uid": uid.upper(),
                "hash": get_sequence_hash(seq_str)
            })
            
    print(f"Loaded {len(all_sequences)} test sequences from {len(fasta_files)} files.")
    
    fasta_path = os.path.join(project_dir, "data", "build", "psfd_sequences.fasta")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    blastdb_dir = os.path.join(project_dir, "data", "blastdb")
    embeddb_dir = os.path.join(project_dir, "data", "embeddb")
    
    print("Creating pristine database backups before benchmarking...")
    for path in [fasta_path, meta_path, blastdb_dir, embeddb_dir]:
        backup_file_or_dir(path)
        
    for path in [fasta_path, meta_path, blastdb_dir, embeddb_dir]:
        if not os.path.exists(path + ".backup"):
            print(f"Error: Backup failed for {path}. Aborting.")
            return
            
    all_results = []
    
    try:
        for idx, seq_info in enumerate(all_sequences):
            print(f"\n[Progress] Sequence {idx+1}/{len(all_sequences)}: {seq_info['uid']} (from {seq_info['file_name']})")
            print(f"Masking sequence {seq_info['uid']} / hash {seq_info['hash'][:8]}...")
            mask_single(project_dir, seq_info["uid"], seq_info["hash"])
            
            for method in methods:
                print(f"Benchmarking query using {method}...")
                res_df = benchmark_query(project_dir, method, seq_info["header"], seq_info["sequence"])
                
                if not res_df.empty:
                    res_df.insert(0, "query_file", seq_info["file_name"])
                    res_df.insert(1, "query_header", seq_info["header"])
                    res_df.insert(2, "query_uniprot_id", seq_info["uid"])
                    res_df.insert(3, "query_seq_hash", seq_info["hash"])
                    res_df.insert(4, "search_method", method)
                    all_results.append(res_df)
                else:
                    empty_row = pd.DataFrame([{
                        "query_file": seq_info["file_name"],
                        "query_header": seq_info["header"],
                        "query_uniprot_id": seq_info["uid"],
                        "query_seq_hash": seq_info["hash"],
                        "search_method": method,
                        "uniprot_id": None,
                        "global_node_id": None,
                        "score": None,
                        "evalue": None,
                        "selected_protein_name": "No hits found"
                    }])
                    all_results.append(empty_row)
    finally:
        print("\nRestoring database to original unmasked state...")
        for path in [fasta_path, meta_path, blastdb_dir, embeddb_dir]:
            restore_file_or_dir(path)
        print("Database restored successfully.")
        
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        final_df.to_csv(output_path, index=False)
        print(f"\n[Success] Benchmark results saved to {output_path}")
    else:
        print("\n[Warning] No benchmark results generated.")

def main():
    parser = argparse.ArgumentParser(description="Benchmark Masking Utility for Leave-One-Out Validation")
    parser.add_argument("command", choices=["mask", "unmask", "benchmark"], help="Action to perform")
    parser.add_argument("--folder", type=str, default=None, help="Folder containing validation FASTA files (defaults to data/benchmark_targets)")
    parser.add_argument("--output", type=str, default=None, help="Output file path (defaults to data/benchmark_results.csv)")
    parser.add_argument("--methods", type=str, default="seq2graph,embed2graph", help="Comma-separated search methods (seq2graph, embed2graph)")
    
    args = parser.parse_args()
    
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    if args.command == "mask":
        target_dir = os.path.join(project_dir, "data", "benchmark_targets")
        queries_dir = os.path.join(project_dir, "data", "benchmark_queries")
        mask_database(project_dir, target_dir, queries_dir)
    elif args.command == "unmask":
        unmask_database(project_dir)
    elif args.command == "benchmark":
        folder = args.folder if args.folder else os.path.join(project_dir, "data", "benchmark_targets")
        output = args.output if args.output else os.path.join(project_dir, "data", "benchmark_results.csv")
        methods = [m.strip() for m in args.methods.split(",") if m.strip()]
        
        run_benchmark(project_dir, folder, output, methods)

if __name__ == "__main__":
    main()
