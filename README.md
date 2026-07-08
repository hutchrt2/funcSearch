# Plant Stress Mechanism Map (PSMM) System

## Overview
The Plant Stress Mechanism Map (PSMM) is an integrated, high-performance system for mapping, analyzing, and retrieving complex plant stress biological entities. By bridging large-scale structural language models with traditional sequence alignment tools, PSMM translates raw molecular sequence data into a cohesive, relational knowledge graph framework.

This repository serves as the central hub of the PSMM infrastructure, unifying four foundational components:
1. **Sequence Fetcher Engine**: Aggregates, normalizes, and curates sequence data from multiple databases (UniProt, NCBI, etc.).
2. **Seq2Graph Bridge**: Provides ultra-fast sequence homology searching using MMseqs2 to link query sequences to established knowledge graph entities.
3. **Embed2Graph Bridge**: Employs deep-learning Evolutionary Scale Modeling (ESM-C) to generate protein language embeddings, mapping complex protein structures and semantic similarities through FAISS vector indexing.
4. **API Dispatcher**: A high-throughput FastAPI service that seamlessly routes incoming biochemical queries to the appropriate mapping bridge and serves knowledge graph relationships.

## Architecture

```text
5_PSMM/
├── data/                       # Local datasets, indexes, and models (Ignored by Git)
│   ├── input/                  # Base metadata and normalization mappings
│   ├── build/                  # Intermediate sequence databases
│   ├── blastdb/                # MMseqs2 indices for seq2graph mapping
│   └── embeddb/                # FAISS vector database and embeddings for embed2graph mapping
├── psmm/                       # Core python modules
│   ├── fetcher/                # Pipeline for dataset aggregation
│   ├── bridges/                # Core sequence alignment & embedding logic
│   └── api/                    # FastAPI server & endpoints
└── scripts/                    # Shell utilities for pipeline execution
```

## System Requirements
- Python 3.9+
- MMseqs2 (must be installed and available in the system PATH)
- Dependencies listed in `requirements.txt` (includes PyTorch, FAISS, Transformers, ESM)

## Quick Start & Pipeline Execution

The PSMM system relies on a sequential build process that processes normalization data, fetches sequences, and compiles the search indexes.

**1. Virtual Environment Setup**
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**2. Full System Rebuild**
To comprehensively build the database indexes from scratch, execute the rebuilding pipeline script. This script will sequentially load the required data, execute the fetcher, and build the MMseqs2 and FAISS databases utilizing high-thread inference.
```bash
./scripts/full_rebuild_pipeline.sh
```

**3. Launch the API Server**
Once the databases are indexed, launch the robust dispatch server to serve queries.
```bash
./scripts/start_api.sh
```

## Data Privacy and Version Control
In compliance with data protection policies and repository size limits, raw datasets, sequence FASTAs, model weights, and constructed vector indexes are strictly excluded from version control. Only the core logic and pipeline infrastructure are tracked.

## Authors & Contributions
Developed for the Plant Stress Functional Database (PSFD) to accelerate agricultural and biochemical mechanisms discovery.
