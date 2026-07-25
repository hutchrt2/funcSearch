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
            
    # Build context mapping for all genes
    print("Building context index...")
    gene_contexts = defaultdict(set)
    for rel in data.get("relations", []):
        nids = [rel.get("subject_entity_id"), rel.get("object_entity_id")] + rel.get("context_entity_ids", [])
        nids = [n for n in nids if n]
        
        # Get all valid concepts in this relation
        concepts = {node_to_concept[n] for n in nids if n in node_to_concept}
        
        # Add all concepts to each node's context pool
        for n in nids:
            g_node = "global.entity." + n.split(".entity.")[-1] if ".entity." in n else n
            gene_contexts[g_node].update(concepts)
            
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
        
        q_ctx = gene_contexts.get(q_node, set())
        h_ctx = gene_contexts.get(h_node, set())
        
        # Remove self-references
        q_concept = global_node_to_concept.get(q_node)
        if q_concept in q_ctx:
            q_ctx.remove(q_concept)
            
        h_concept = global_node_to_concept.get(h_node)
        if h_concept in h_ctx:
            h_ctx.remove(h_concept)
        
        if not q_ctx or not h_ctx:
            overlap_score = 0.0
        else:
            intersection = q_ctx.intersection(h_ctx)
            union = q_ctx.union(h_ctx)
            overlap_score = len(intersection) / len(union) if len(union) > 0 else 0.0
            
        is_hit_found = hit_uid != "NAN" and str(row.get("selected_protein_name", "")) != "No hits found"
        
        shared_labels = [str(node_to_label.get(c, c)) for c in q_ctx.intersection(h_ctx)]
        missed_labels = [str(node_to_label.get(c, c)) for c in q_ctx - h_ctx]
        
        results_eval.append({
            "Query UID": query_uid,
            "Method": method,
            "Hit UID": hit_uid,
            "Context Overlap Score": round(overlap_score, 3),
            "Shared Contexts": " | ".join(shared_labels),
            "Missed Contexts": " | ".join(missed_labels),
            "Hit Found": is_hit_found
        })
        
    eval_df = pd.DataFrame(results_eval)
    
    print("\n--- Benchmark Context Accuracy Summary ---")
    
    methods = eval_df['Method'].unique()
    for method in methods:
        method_df = eval_df[eval_df['Method'] == method]
        total = len(method_df)
        hits_found = method_df['Hit Found'].sum()
        avg_overlap = method_df[method_df['Hit Found']]['Context Overlap Score'].mean()
        
        print(f"\nMethod: {method}")
        print(f"  Total Queries: {total}")
        print(f"  Queries with a Hit: {hits_found} ({(hits_found/total)*100:.1f}%)")
        print(f"  Average Context Overlap (when hit found): {avg_overlap:.3f}")
        
    out_eval_path = os.path.join(project_dir, "data", "benchmark_evaluation.csv")
    eval_df.to_csv(out_eval_path, index=False)
    print(f"\nDetailed evaluation saved to {out_eval_path}")

if __name__ == "__main__":
    main()
