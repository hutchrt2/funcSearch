import os
import csv
import json
import pandas as pd
from collections import defaultdict

def main():
    # Dynamically find the root directory so it works on both local and HPC
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    results_path = os.path.join(project_dir, "data", "benchmark_results.csv")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    db_path = os.path.join(project_dir, "data", "global_path_index.json")
    
    if not os.path.exists(results_path) or not os.path.exists(db_path):
        print(f"Error: Missing required files in {project_dir}/data")
        return
        
    print("Loading metadata...")
    meta_df = pd.read_csv(meta_path)
    
    acc_to_node = {}
    for _, row in meta_df.iterrows():
        uid = str(row['target_accession']).upper()
        nid = str(row['global_node_id'])
        acc_to_node[uid] = nid
        
    print("Loading graph index...")
    with open(db_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    node_to_concept = {}
    node_to_label = {}
    global_node_to_concept = {}
    for e in data.get("entities", []):
        cid = e.get("ontology_id")
        if cid:
            node_to_concept[e["id"]] = cid
            node_to_label[cid] = e.get("selected_label") or e.get("name") or cid
            g_node = "global.entity." + e["id"].split(".entity.")[-1] if ".entity." in e["id"] else e["id"]
            global_node_to_concept[g_node] = cid
            
    # Build separate indices for direct semantic triples (relations) vs event contexts
    print("Building relation and context indices...")
    gene_triples = defaultdict(set)
    gene_contexts = defaultdict(set)
    
    for rel in data.get("relations", []):
        sub_id = rel.get("subject_entity_id") or rel.get("subject_node_id")
        obj_id = rel.get("object_entity_id") or rel.get("object_node_id")
        ctx_ids = rel.get("context_entity_ids") or rel.get("context_node_ids") or []
        
        rel_nids = [n for n in [sub_id, obj_id] if n]
        rel_concepts = {node_to_concept[n] for n in rel_nids if n in node_to_concept}
        
        ctx_nids = [n for n in ctx_ids if n]
        ctx_concepts = {node_to_concept[n] for n in ctx_nids if n in node_to_concept}
        
        for n in rel_nids:
            g_node = "global.entity." + n.split(".entity.")[-1] if ".entity." in n else n
            gene_triples[g_node].update(rel_concepts)
            gene_contexts[g_node].update(ctx_concepts)
            
    for ev in data.get("events", []):
        participants = [p.get("node_id") for p in ev.get("participants", []) if isinstance(p, dict)]
        ev_concepts = {node_to_concept[n] for n in participants if n in node_to_concept}
        for n in participants:
            if n:
                g_node = "global.entity." + n.split(".entity.")[-1] if ".entity." in n else n
                gene_contexts[g_node].update(ev_concepts)
            
    # Evaluate benchmark results
    print("Evaluating results...")
    df = pd.read_csv(results_path)
    results_eval = []
    
    for _, row in df.iterrows():
        query_uid = str(row['query_uniprot_id']).upper()
        hit_uid = str(row['uniprot_id']).upper() if pd.notna(row['uniprot_id']) else "NAN"
        method = row['search_method']
        
        q_node = acc_to_node.get(query_uid)
        h_node = acc_to_node.get(hit_uid)
        
        q_triple_concepts = set(gene_triples.get(q_node, set()))
        h_triple_concepts = set(gene_triples.get(h_node, set()))
        
        q_ctx_concepts = set(gene_contexts.get(q_node, set()))
        h_ctx_concepts = set(gene_contexts.get(h_node, set()))
        
        # Remove self-references
        q_concept = global_node_to_concept.get(q_node)
        if q_concept:
            q_triple_concepts.discard(q_concept)
            q_ctx_concepts.discard(q_concept)
            
        h_concept = global_node_to_concept.get(h_node)
        if h_concept:
            h_triple_concepts.discard(h_concept)
            h_ctx_concepts.discard(h_concept)
        
        def jaccard(s1, s2):
            if not s1 or not s2:
                return 0.0
            union = s1.union(s2)
            return len(s1.intersection(s2)) / len(union) if len(union) > 0 else 0.0

        triple_overlap = jaccard(q_triple_concepts, h_triple_concepts)
        context_overlap = jaccard(q_ctx_concepts, h_ctx_concepts)
        combined_overlap = jaccard(q_triple_concepts | q_ctx_concepts, h_triple_concepts | h_ctx_concepts)
            
        is_hit_found = hit_uid != "NAN" and str(row.get("selected_protein_name", "")) != "No hits found"
        
        shared_triples = [str(node_to_label.get(c, c)) for c in q_triple_concepts.intersection(h_triple_concepts)]
        shared_contexts = [str(node_to_label.get(c, c)) for c in q_ctx_concepts.intersection(h_ctx_concepts)]
        
        results_eval.append({
            "Query UID": query_uid,
            "Method": method,
            "Hit UID": hit_uid,
            "Triple Overlap Score": round(triple_overlap, 3),
            "Context Overlap Score": round(context_overlap, 3),
            "Combined Overlap Score": round(combined_overlap, 3),
            "Shared Triples": " | ".join(shared_triples),
            "Shared Contexts": " | ".join(shared_contexts),
            "Hit Found": is_hit_found
        })
        
    eval_df = pd.DataFrame(results_eval)
    
    print("\n--- Benchmark Accuracy Summary (Triples vs. Contexts) ---")
    
    methods = eval_df['Method'].unique()
    for method in methods:
        method_df = eval_df[eval_df['Method'] == method]
        total = len(method_df)
        hits_found = method_df['Hit Found'].sum()
        avg_triple = method_df[method_df['Hit Found']]['Triple Overlap Score'].mean()
        avg_context = method_df[method_df['Hit Found']]['Context Overlap Score'].mean()
        avg_combined = method_df[method_df['Hit Found']]['Combined Overlap Score'].mean()
        
        print(f"\nMethod: {method}")
        print(f"  Total Queries: {total}")
        print(f"  Queries with a Hit: {hits_found} ({(hits_found/total)*100:.1f}%)")
        print(f"  Average Triple Overlap (Semantic Relations): {avg_triple:.3f}")
        print(f"  Average Context Overlap (Events/Conditions): {avg_context:.3f}")
        print(f"  Average Combined Overlap:                    {avg_combined:.3f}")
        
    out_eval_path = os.path.join(project_dir, "data", "benchmark_evaluation.csv")
    eval_df.to_csv(out_eval_path, index=False)
    print(f"\nDetailed evaluation saved to {out_eval_path}")

if __name__ == "__main__":
    main()
