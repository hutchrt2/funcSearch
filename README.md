# PlantStress MechanismMap (PSMM)

The Plant Stress Mechanism Map (PSMM) is an integrated system for mapping, analyzing, and querying plant stress biological entities at the protein level. It bridges structural language models with traditional sequence alignment tools to translate raw molecular sequence data into a navigable relational knowledge graph.

This repository unifies four foundational components into a single deployable service:

1. **Sequence Fetcher** — Aggregates and curates protein sequences from UniProt, NCBI, Ensembl, and Phytozome.
2. **Seq2Graph Bridge** — Ultra-fast sequence homology matching via MMseqs2, linking query proteins to knowledge graph entities.
3. **Embed2Graph Bridge** — Deep-learning protein embeddings via ESM-C (EvolutionaryScale), mapped through FAISS cosine similarity search.
4. **API Dispatcher** — A FastAPI server that routes biochemical queries to the appropriate bridge and returns structured knowledge graph relationships.

---

## Architecture

```text
5_PSMM/
├── psmm/                         # Core Python package
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
│   ├── load_normalization_data.sh# Copy & decompress PSFD normalization data
│   └── benchmark_masking.py      # Benchmark masking utilities for evaluation
├── data/                         # All runtime data (gitignored)
│   ├── input/                    # Normalization CSVs from PSFD
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
- Queries are resolved against a structured knowledge graph derived from the PSFD corpus.
- Matched proteins are linked to their biological relationships (genes, metabolites, pathways, tissues, species, experimental conditions).
- Results include paper provenance (PMCIDs), ontology annotations, and entity context.

---

## System Requirements

- **Python** ≥ 3.10 (required for ESM-C / EvolutionaryScale package)
- **MMseqs2** — must be installed and available on `PATH`
- **CUDA** (optional) — GPU acceleration for ESM-C inference
- Dependencies in `requirements.txt`: PyTorch, FAISS, Transformers, ESM, FastAPI, Uvicorn, BioPython, Pandas, NumPy

---

## Quick Start

### 1. Virtual Environment Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 2. Building the Knowledge Graph

If you have raw pipeline outputs from PSFD (either as extracted folders or `.tar.gz` files), drop them directly into the `input/` folder. The system will automatically extract tarballs and resolve the required datasets.

```bash
./psmm.py build-graph
```

This will compile the knowledge graph and write it to `data/global_path_index.json`, along with the `manifest.json` and individual paper records.

### 3. Database Rebuild (Full Pipeline)

Sequentially loads normalization data, fetches sequences, and builds both MMseqs2 and FAISS indexes:

```bash
./psmm.py rebuild
```

Pipeline stages:
1. Load normalization data from PSFD (`load_normalization_data.sh`)
2. Run the sequence fetcher (`psmm.fetcher.pipeline`)
3. Build MMseqs2 database (`psmm.bridges.seq2graph --init --clean`)
4. Build FAISS index with ESM-C embeddings (`psmm.bridges.embed2graph --init --clean`)

### 3. Launch the API Server

```bash
./psmm.py serve --port 8999
```

The server starts on port **8999** by default and exposes the following endpoints:

| Endpoint | Method | Description |
| -------- | ------ | ----------- |
| `/search` | POST | Sequence similarity search (MMseqs2 or ESM-C) |
| `/api/extract` | POST | Full relation extraction with attribute filtering |
| `/api/stats` | GET | Database statistics (entity/concept/relation counts) |
| `/api/resolve_entities` | POST | Term-to-entity resolution |
| `/api/ontology_count` | GET | Unique ontology term count |
| `/api/data/*` | GET | Static data files (manifest, papers, metadata) |

### 4. Benchmarking

The system supports controlled evaluation by masking known entities from the databases and measuring retrieval accuracy:

```bash
./psmm.py benchmark run --help
./psmm.py benchmark evaluate
```

### 5. Enrichments

Run pathway enrichment tests directly on the knowledge graph using the unified CLI:

```bash
# Calculate and add enrichments to global_path_index.json
./psmm.py enrich calculate --db data/global_path_index.json

# Export global enrichments to CSV
./psmm.py enrich export --output enrichments_export.csv

# Update individual paper JSON files with enrichments
./psmm.py enrich update-papers
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

