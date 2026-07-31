# funcSearch

funcSearch is an integrated system for mapping, analyzing, and querying plant stress biological entities at the protein level. It bridges structural language models with traditional sequence alignment tools to translate raw molecular sequence data into a navigable relational knowledge graph.

This repository unifies four foundational components into a single deployable service:

1. **Sequence Fetcher** — Aggregates and curates protein sequences from UniProt, NCBI, Ensembl, and Phytozome.
2. **Seq2Graph Bridge** — Ultra-fast sequence homology matching via MMseqs2, linking query proteins to knowledge graph entities.
3. **Embed2Graph Bridge** — Deep-learning protein embeddings via ESM-C (EvolutionaryScale), mapped through FAISS cosine similarity search.
4. **API Dispatcher** — A FastAPI server that routes biochemical queries to the appropriate bridge and returns structured knowledge graph relationships.

---

## Architecture

```text
1_funcSearch/
├── funcSearch/                         # Core Python package
│   ├── api/
│   │   ├── server.py             # FastAPI application (endpoints, CORS, caching)
│   │   └── verify.py             # Database integrity verification utilities
│   ├── bridges/
│   │   ├── seq2graph.py          # MMseqs2 sequence alignment bridge
│   │   └── embed2graph.py        # ESM-C embedding + FAISS vector search bridge
│   └── fetcher/
│       └── pipeline.py           # Multi-DB sequence retrieval engine
├── scripts/
│   ├── start_api.sh              # Launch the FastAPI server
│   ├── full_rebuild_pipeline.sh  # End-to-end database rebuild
│   ├── load_normalization_data.sh# Copy & decompress funcMap normalization data
│   └── benchmark_masking.py      # Benchmark masking utilities for evaluation
├── data/                         # All runtime data (gitignored)
│   ├── input/                    # Normalization CSVs from funcMap
│   ├── build/                    # Intermediate: sequence FASTA + metadata
│   ├── blastdb/                  # MMseqs2 database files (symlink)
│   ├── embeddb/                  # FAISS index + ESM-C embeddings (symlink)
│   ├── global_path_index.json    # Knowledge graph database (symlink)
│   ├── manifest.json             # Paper metadata manifest (symlink)
│   └── papers/                   # Individual paper JSON files (symlink)
├── requirements.txt
└── README.md
```

---

## Key Features

### Dual Search Engines
- **MMseqs2 (Seq2Graph)**: Alignment-based homology search with configurable E-value threshold and minimum sequence identity. Returns all hits meeting the threshold.
- **ESM-C + FAISS (Embed2Graph)**: Embedding-based semantic similarity search with configurable k-nearest-neighbors and minimum cosine similarity threshold.

### In-Process Execution (No Subprocesses)
Both bridge engines execute directly within the API server process:
- ESM-C model and FAISS index are lazy-loaded on the first embedding query, then held in memory.
- MMseqs2 searches invoke the core library in-process rather than spawning a child Python interpreter.

### Performance Optimizations
- **Binary Pickle Cache**: The ~5 GB `global_path_index.json` knowledge graph is compiled to a pickle cache on first load, reducing subsequent startup from minutes to seconds.
- **Search Result Caching**: In-memory LRU caches for both `/search` and `/api/extract` endpoints with compound cache keys (sequence + method + all parameters).
- **Extraction Result Caching**: Full relation extraction queries are cached by their compound input, FASTA input, attribute filters, and search parameters.

### Advanced Search Parameters
Users can configure search behavior through the API and UI:

| Parameter | Engine | Description | Default |
| --------- | ------ | ----------- | ------- |
| `evalue` | MMseqs2 | E-value significance threshold | `0.001` |
| `min_seq_id` | MMseqs2 | Minimum sequence identity (0–1) | None |
| `k` | ESM-C | Maximum nearest neighbors returned | `5` |
| `min_similarity` | ESM-C | Minimum cosine similarity threshold | None |

### Knowledge Graph Integration
- Queries are resolved against a structured knowledge graph derived from the funcMap corpus.
- Matched proteins are linked to their biological relationships (genes, metabolites, pathways, tissues, species, experimental conditions).
- Results include paper provenance (PMCIDs), ontology annotations, and entity context.

---

## System Requirements

- **Python** ≥ 3.10 (required for ESM-C / EvolutionaryScale package)
- **MMseqs2** — must be installed and available on `PATH`
- **CUDA** (optional) — GPU acceleration for ESM-C inference
- Dependencies in `requirements.txt`: PyTorch, FAISS, Transformers, ESM, FastAPI, Uvicorn, BioPython, Pandas, NumPy

---

## CLI Reference

The `funcSearch.py` script provides a unified interface for all pipeline and server operations.

### 1. Server Commands
Starts the FastAPI backend for handling sequence alignments and ontology relationships.

```bash
./funcSearch.py serve [OPTIONS]
```
| Option | Description |
|--------|-------------|
| `--port` | Port to run the server on (default: 8999) |
| `--host` | Host address to bind (default: 0.0.0.0) |
| `--dev` | Run in development mode with automatic code reloading |
| `--tunnel` | Publicly expose the API using an Ngrok tunnel |
| `--ngrok-token` | Auth token for the Ngrok tunnel (if required) |

### 2. Pipeline Rebuild Commands
Manages the end-to-end dataset pipeline, including raw knowledge graph construction and sequence index building.

**Full Pipeline:**
```bash
./funcSearch.py rebuild full [--sync-norm]
```
*(Runs the entire process sequentially: graph compilation, enrichment calculations, sequence fetching, and database indexing. Pass `--sync-norm` to download the latest normalization dataset before starting).*

**Knowledge Graph Pipeline:**
```bash
./funcSearch.py rebuild graph [all | compile | enrich]
```
- `all`: (Default) Compiles the graph and calculates enrichments.
- `compile`: Compiles the UI database from raw pipeline outputs (saves to `data/global_path_index.json`).
- `enrich`: Only calculates the Fisher/Chi-squared pathway enrichments and updates paper bundles.

**Sequence Search Pipeline:**
```bash
./funcSearch.py rebuild search [all | fetch | index] [--method METHOD]
```
- `all`: (Default) Fetches sequences and builds the indices.
- `fetch`: Aggregates physical FASTA sequences from upstream sources.
- `index`: Builds the `seq2graph` (MMseqs2) and `embed2graph` (FAISS) search indices.
- `--method`: Specify `seq2graph`, `embed2graph`, or `all` when building indices.

### 3. Data & Enrichments
Standalone commands for manipulating the compiled datasets and enrichments.

```bash
# Calculate enrichment stats on the graph
./funcSearch.py enrich calculate

# Export the enrichments as a flat CSV file
./funcSearch.py enrich export --output enrichments_export.csv

# Distribute the global enrichments into individual paper JSON bundles
./funcSearch.py enrich update-papers

# Compile raw graph data without running full rebuild
./funcSearch.py build-graph --outdir data
```

### 4. Benchmarking & Evaluation
Tools for validating the search performance of the two engines using sequence masking.

```bash
# Run the evaluation benchmark
./funcSearch.py benchmark run [--query FILE] [--output FILE] [--method seq2graph|embed2graph|both]

# Evaluate the generated benchmark results
./funcSearch.py benchmark evaluate
```

### 5. Utilities
```bash
# Verify the integrity of the built databases
./funcSearch.py verify

# Standalone sequence fetcher
./funcSearch.py fetch [--input FILE] [--output DIR] [--no-cache]

# Standalone DB builder
./funcSearch.py db [--init] [--clean] [--method METHOD]
```

---

## API Usage Examples

### Sequence Search (ESM-C)

```bash
curl -X POST http://localhost:8999/search \
  -H "Content-Type: application/json" \
  -d '{"sequence": "MKTLVL...", "method": "embed2graph", "k": 10, "min_similarity": 0.7}'
```

### Sequence Search (MMseqs2)

```bash
curl -X POST http://localhost:8999/search \
  -H "Content-Type: application/json" \
  -d '{"sequence": "MKTLVL...", "method": "seq2graph", "evalue": 0.001}'
```

### Relation Extraction

```bash
curl -X POST http://localhost:8999/api/extract \
  -H "Content-Type: application/json" \
  -d '{
    "fasta": ">query\nMKTLVL...",
    "method": "embed2graph",
    "k": 5,
    "attributes": {"genes": true, "pathways": true, "metabolites": true}
  }'
```

---

## Data Privacy

Raw datasets, sequence FASTAs, model weights, vector indexes, and the knowledge graph database are **strictly excluded** from version control. Only source code, pipeline scripts, and configuration are tracked. The `data/` directory is fully gitignored.

Additionally, the `input/` directory is provided as a blank template. Its contents are gitignored except for a `.gitkeep` file to prevent accidental uploads of confidential sequences or datasets.

---

## Required Input Files

To build the knowledge graph correctly, you must provide the resolved entity normalizations in the `input/` directory. The pipeline will accept these files as uncompressed `.csv` files or compressed `.tar.gz` archives, as the script will automatically uncompress them if needed.

Place the following three files in `input/` (or in a sub-folder matching the expected tarball logic):

1. **`gene_protein_entity_summary.csv`**
   - **Role:** Centralized summary of all resolved gene/protein entities mapped to global node IDs.
   - **Key Columns:** `pmcids`, `entity_instance_ids` (semicolon separated), `best_decision`, `selected_uniprot_accession`, `selected_phytozome_gene_id`.

2. **`gene_protein_taxon_summary.csv`**
   - **Role:** Species-specific sequence assignments and taxonomy contexts for ambiguous or multi-species genes.
   - **Key Columns:** Same as the entity summary (`pmcids`, `entity_instance_ids`, `best_decision`, etc.).

3. **`gene_protein_normalizations.csv`**
   - **Role:** The raw, per-paper normalization assignments containing full alias lists and row-level metadata.
   - **Key Columns:** `pmcid` and `entity_instance_id` (singular identifiers), `decision`, `canonical_form`, `aliases`.

