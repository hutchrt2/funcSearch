#!/bin/bash
set -e

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== [1/6] Loading new normalization data ==="
$PROJECT_DIR/scripts/load_normalization_data.sh

cd $PROJECT_DIR
PYTHON_EXE="python"
if [ -f ".venv/bin/python" ]; then
    PYTHON_EXE=".venv/bin/python"
fi

echo "=== [2/6] Running sequence fetcher pipeline ==="
$PYTHON_EXE -m psmm.fetcher.pipeline

echo "=== [3/6] Performing clean rebuild of Seq2Graph Database ==="
export PATH=/programs/mmseqs/bin:$PATH
$PYTHON_EXE -m psmm.bridges.seq2graph --init --clean

echo "=== [4/6] Performing clean rebuild of Embed2Graph Database ==="
OMP_NUM_THREADS=20 MKL_NUM_THREADS=20 $PYTHON_EXE -m psmm.bridges.embed2graph --init --clean

echo "=== [5/6] Calculating and adding enrichments ==="
$PYTHON_EXE $PROJECT_DIR/scripts/add_enrichments.py

echo "=== [6/6] Updating paper files with enrichments ==="
$PYTHON_EXE $PROJECT_DIR/scripts/update_paper_enrichments.py

echo "=== Rebuild completed successfully! ==="
