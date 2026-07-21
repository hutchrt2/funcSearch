import json
import scipy.stats as stats
from collections import Counter, defaultdict
from typing import Any

def calculate_enrichments(entities: list[dict[str, Any]], relations: list[dict[str, Any]]) -> None:
    def is_gene(t: str) -> bool:
        t = t.lower()
        return "gene" in t or "protein" in t
    
    def is_trait(t: str) -> bool:
        t = t.lower()
        return "trait" in t or "pathway" in t or "phenotype" in t or "process" in t or "disease" in t

    gene_entities = [e for e in entities if is_gene(e.get("type", ""))]
    trait_entities = [e for e in entities if is_trait(e.get("type", ""))]

    node_to_concept = {}
    for e in entities:
        concept_id = e.get("ontology_id")
        if concept_id:
            node_to_concept[e["id"]] = concept_id

    gene_concepts = {node_to_concept[e["id"]] for e in gene_entities if e["id"] in node_to_concept}
    trait_concepts = {node_to_concept[e["id"]] for e in trait_entities if e["id"] in node_to_concept}

    if not gene_concepts or not trait_concepts:
        print("No gene or trait concepts found.")
        return

    total_relations = len(relations)

    gene_counts = Counter()
    trait_counts = Counter()
    co_counts = Counter()

    for rel in relations:
        node_ids = [rel.get("subject_entity_id"), rel.get("object_entity_id")] + rel.get("context_entity_ids", [])
        node_ids = [n for n in node_ids if n]
        
        rel_concepts = {node_to_concept[n] for n in node_ids if n in node_to_concept}
        g_c = rel_concepts & gene_concepts
        t_c = rel_concepts & trait_concepts
        
        for g in g_c:
            gene_counts[g] += 1
        for t in t_c:
            trait_counts[t] += 1
            
        for g in g_c:
            for t in t_c:
                co_counts[(g, t)] += 1

    p_values = []
    tests = []

    for (g, t), v1 in co_counts.items():
        v2 = gene_counts[g] - v1
        v3 = trait_counts[t] - v1
        v4 = total_relations - v1 - v2 - v3
        
        _, p = stats.fisher_exact([[v1, v2], [v3, v4]], alternative="greater")
        
        tests.append((g, t))
        p_values.append(p)
    
    if not p_values:
        return

    fdr_adjusted = stats.false_discovery_control(p_values)
    
    significant_pairs = defaultdict(list)
    concept_to_label = {e["ontology_id"]: e.get("selected_label") or e["ontology_id"] for e in trait_entities if "ontology_id" in e}

    added_count = 0
    for i, (g, t) in enumerate(tests):
        if fdr_adjusted[i] < 0.05:
            significant_pairs[g].append({
                "trait_concept": t,
                "trait_label": concept_to_label.get(t, t),
                "p_value": p_values[i],
                "fdr": fdr_adjusted[i]
            })
            added_count += 1
            
    print(f"Found {added_count} significant enrichments.")
            
    for e in gene_entities:
        concept = e.get("ontology_id")
        if concept and concept in significant_pairs:
            sorted_enrichments = sorted(significant_pairs[concept], key=lambda x: x["fdr"])
            e["enrichments"] = sorted_enrichments
        else:
            e["enrichments"] = []


def main():
    # Update this path to the location on the BioHPC cluster
    db_path = "/local/storage/thomas/5_PSMM/data/global_path_index.json"
    
    print(f"Loading {db_path}...")
    with open(db_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    print("Calculating enrichments...")
    calculate_enrichments(data["entities"], data["relations"])
    
    print("Saving...")
    with open(db_path, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"))
    print("Done!")
    
if __name__ == "__main__":
    main()
