import os
import sys
import time

# Ensure we can import api_server
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from fastapi.testclient import TestClient
from funcSearch.api.server import app, db

client = TestClient(app)

def run_tests():
    # Use context manager to trigger @app.on_event("startup")
    print("Initializing test client (triggering startup event)...")
    with TestClient(app) as test_client:
        print("\n--- 1. Testing /api/stats endpoint ---")
        response = test_client.get("/api/stats")
        print("Status Code:", response.status_code)
        print("Response:", response.json())
        assert response.status_code == 200
        stats = response.json()
        assert stats["entities"] > 0
        assert stats["concepts"] > 0
        assert stats["relations"] > 0
        print("Stats endpoint test passed!")

        print("\n--- 1b. Testing /api/resolve_entities endpoint ---")
        resolve_res = test_client.post("/api/resolve_entities", json={"terms": ["GABA"], "category": "compound"})
        print("Status Code:", resolve_res.status_code)
        print("Response:", resolve_res.json())
        assert resolve_res.status_code == 200
        assert len(resolve_res.json()) > 0
        print("Resolve entities endpoint test passed!")

        print("\n--- 2. Testing /api/extract endpoint (compounds only) ---")
        # Let's extract relations for compound "GABA"
        payload = {
            "compounds": "GABA",
            "fasta": "",
            "attributes": {
                "genes": True,
                "metabolites": True,
                "pathways": True,
                "tissues": True,
                "species": True,
                "experimental_conditions": True,
                "plant_traits": True,
                "molecular_traits": True,
                "human_traits": True,
            },
            "method": "embed2graph"
        }
        
        t0 = time.time()
        res = test_client.post("/api/extract", json=payload)
        t1 = time.time()
        
        print("Status Code:", res.status_code)
        print(f"Response (first 2 rows of {len(res.json())}):", res.json()[:2])
        assert res.status_code == 200
        first_extract_time = t1 - t0
        print(f"First extraction time: {first_extract_time:.4f}s")
        assert len(res.json()) > 0
        
        print("\n--- 3. Testing /api/extract endpoint Caching (duplicate request) ---")
        t0 = time.time()
        res_cached = test_client.post("/api/extract", json=payload)
        t1 = time.time()
        cached_extract_time = t1 - t0
        print(f"Cached extraction time: {cached_extract_time:.4f}s")
        assert res_cached.status_code == 200
        assert len(res_cached.json()) == len(res.json())
        # Caching should be near instant (< 5ms)
        print(f"Speedup: {first_extract_time / max(1e-6, cached_extract_time):.1f}x")
        assert cached_extract_time < 0.050, "Cache did not return response within 50ms"
        print("Extract caching test passed!")

        print("\n--- 4. Testing /api/extract with FASTA logic (mocked search) ---")
        # Let's query with a sequence of a known entity and mock its search result
        # To do this, let's look at sequence_index.json accession A0A022PNA9, which links to entity: PMC4125134.entity.f9ce18b23f49 (elongation factor Tu)
        # Suffix is f9ce18b23f49.
        # We can mock perform_search_internal by putting search results directly into _search_cache to avoid executing bridge script subprocess
        from funcSearch.api.server import _search_cache
        sequence = "MGRAPCCDKASVKRGPWSPEEDEQLRSYVQSHGIGGNWIALPQKAGLNRCGKSCRLRWLNYLRPDIKHGGYTEQEDHIICSLYNSIGSRWSIIASKLPGRTDNDVKNYWNTKLKKKAMGAVQPRAAASAPSQCTSSAMAPALSPASSSVTSSSGDACFAAAATTTTTMYPPPTTPPQQQFIRFDAPPAAAAAASPTDLAPVPPPATVTADGDGGWASDALSLDDVFLGELTAGEPLFPYAELFSGFAGAAPDSKATLELSACYFPNMAEMWAASDHAYAKPQGLCNTLT"
        mock_results = [{
            "query": "query_sequence",
            "uniprot_id": "A0A0P0WQQ7",
            "global_node_id": "global.entity.10e2e5352583",
            "selected_protein_name": "A0A0P0WQQ7_ORYSJ",
            "selected_gene_name": "OsMYB55",
            "selected_organism": "Oryza sativa subsp. japonica"
        }]
        
        # Insert mock search result into cache
        _search_cache[(sequence.strip().upper(), "embed2graph", None, None, None)] = mock_results
        
        fasta_payload = {
            "compounds": "",
            "fasta": f">test_seq\n{sequence}",
            "attributes": None, # Should default to all
            "method": "embed2graph"
        }
        
        res_fasta = test_client.post("/api/extract", json=fasta_payload)
        print("Status Code:", res_fasta.status_code)
        print(f"Response (first 2 rows of {len(res_fasta.json())}):", res_fasta.json()[:2])
        assert res_fasta.status_code == 200
        assert len(res_fasta.json()) > 0
        print("FASTA resolution extraction test passed!")

    print("\nAll verification tests passed successfully!")

if __name__ == "__main__":
    run_tests()
