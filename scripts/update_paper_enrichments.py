import json
import glob
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Update paper JSONs with enrichments")
    parser.add_argument("--db", required=True, help="Path to global_path_index.json")
    parser.add_argument("--papers", required=True, help="Path to papers directory")
    args = parser.parse_args()

    print("Loading global_path_index.json...")
    with open(args.db, "r", encoding="utf-8") as f:
        data = json.load(f)

    enrichment_map = {}
    for e in data.get("entities", []):
        if "enrichments" in e:
            enrichment_map[e["id"]] = e["enrichments"]

    print(f"Loaded {len(enrichment_map)} entities with enrichments.")

    paper_files = glob.glob(os.path.join(args.papers, "*.json"))
    print(f"Updating {len(paper_files)} paper files...")

    updated = 0
    for file in paper_files:
        with open(file, "r", encoding="utf-8") as f:
            paper_data = json.load(f)
        
        changed = False
        for e in paper_data.get("entities", []):
            if e["id"] in enrichment_map:
                e["enrichments"] = enrichment_map[e["id"]]
                changed = True
        
        if changed:
            with open(file, "w", encoding="utf-8") as f:
                json.dump(paper_data, f, separators=(",", ":"))
            updated += 1

    print(f"Done! Updated {updated} paper files.")

if __name__ == "__main__":
    main()
