import os
import csv
import pandas as pd
from collections import defaultdict
import re

def get_words(text):
    if not isinstance(text, str):
        return set()
    # Simple tokenization
    text = text.lower()
    text = re.sub(r'[^a-z0-9\s]', ' ', text)
    words = set(text.split())
    # Remove common stop words
    stop_words = {'protein', 'domain', 'containing', 'family', 'the', 'of', 'and', 'a', 'to', 'in', 'is', 'you', 'that', 'it', 'he', 'was', 'for', 'on', 'are', 'as', 'with', 'his', 'they', 'I', 'at', 'be', 'this', 'have', 'from', 'or', 'one', 'had', 'by', 'word', 'but', 'not', 'what', 'all', 'were', 'we', 'when', 'your', 'can', 'said', 'there', 'use', 'an', 'each', 'which', 'she', 'do', 'how', 'their', 'if', 'will', 'up', 'other', 'about', 'out', 'many', 'then', 'them', 'these', 'so', 'some', 'her', 'would', 'make', 'like', 'him', 'into', 'time', 'has', 'look', 'two', 'more', 'write', 'go', 'see', 'number', 'no', 'way', 'could', 'people', 'my', 'than', 'first', 'water', 'been', 'call', 'who', 'oil', 'its', 'now', 'find', 'long', 'down', 'day', 'did', 'get', 'come', 'made', 'may', 'part'}
    return words - stop_words

def main():
    project_dir = "/home/thomas/Projects/PlantStress-MechanismMap/BioHPC_Mount/5_PSMM"
    results_path = os.path.join(project_dir, "data", "benchmark_results.csv")
    meta_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    
    if not os.path.exists(results_path):
        print(f"Error: {results_path} not found.")
        return
        
    print("Loading metadata to find ground truth labels...")
    meta_df = pd.DataFrame()
    if os.path.exists(meta_path):
        meta_df = pd.read_csv(meta_path)
    
    # Create a mapping of uniprot_id to ground truth protein name
    ground_truth = {}
    if not meta_df.empty and 'target_accession' in meta_df.columns and 'protein_name' in meta_df.columns:
        for _, row in meta_df.iterrows():
            uid = str(row['target_accession']).upper()
            ground_truth[uid] = str(row['protein_name'])
            
    print(f"Loaded {len(ground_truth)} ground truth labels.")
    
    # Read the benchmark results
    df = pd.read_csv(results_path)
    
    # Calculate simple word overlap heuristic to judge "Correctness"
    results_eval = []
    
    for _, row in df.iterrows():
        query_uid = str(row['query_uniprot_id']).upper()
        hit_uid = str(row['uniprot_id'])
        method = row['search_method']
        hit_name = str(row['selected_protein_name'])
        
        # Get ground truth query name
        query_name = ground_truth.get(query_uid, str(row['query_header']))
        
        # Heuristic accuracy score (Jaccard similarity of words)
        q_words = get_words(query_name)
        h_words = get_words(hit_name)
        
        if not q_words or not h_words:
            overlap_score = 0.0
        else:
            intersection = q_words.intersection(h_words)
            union = q_words.union(h_words)
            overlap_score = len(intersection) / len(union) if len(union) > 0 else 0.0
            
        is_hit_found = hit_name != "No hits found" and hit_uid != "nan"
        
        results_eval.append({
            "Query UID": query_uid,
            "Method": method,
            "Query Ground Truth Name": query_name,
            "Hit UID": hit_uid,
            "Hit Name": hit_name,
            "Word Overlap Score": round(overlap_score, 3),
            "Hit Found": is_hit_found
        })
        
    eval_df = pd.DataFrame(results_eval)
    
    print("\n--- Benchmark Accuracy Summary ---")
    
    methods = eval_df['Method'].unique()
    for method in methods:
        method_df = eval_df[eval_df['Method'] == method]
        total = len(method_df)
        hits_found = method_df['Hit Found'].sum()
        avg_overlap = method_df[method_df['Hit Found']]['Word Overlap Score'].mean()
        
        print(f"\nMethod: {method}")
        print(f"  Total Queries: {total}")
        print(f"  Queries with a Hit: {hits_found} ({(hits_found/total)*100:.1f}%)")
        print(f"  Average Annotation Overlap (when hit found): {avg_overlap:.3f}")
        
    out_eval_path = os.path.join(project_dir, "data", "benchmark_evaluation.csv")
    eval_df.to_csv(out_eval_path, index=False)
    print(f"\nDetailed evaluation saved to {out_eval_path}")
    print("Review this CSV to manually grade if the retrieved hit is functionally equivalent to the query.")

if __name__ == "__main__":
    main()
