#!/bin/bash
set -e

PROJECT_DIR="/local/storage/thomas/5_PSMM"

echo "=== [1/4] Loading new normalization data ==="
$PROJECT_DIR/scripts/load_normalization_data.sh

cd $PROJECT_DIR
PYTHON_EXE="python"
if [ -f ".venv/bin/python" ]; then
    PYTHON_EXE=".venv/bin/python"
fi

echo "=== [2/4] Running sequence fetcher pipeline ==="
$PYTHON_EXE -m psmm.fetcher.pipeline

echo "=== [3/4] Performing clean rebuild of Seq2Graph Database ==="
export PATH=/programs/mmseqs/bin:$PATH
$PYTHON_EXE -m psmm.bridges.seq2graph --init --clean

echo "=== [4/4] Performing clean rebuild of Embed2Graph Database ==="
OMP_NUM_THREADS=20 MKL_NUM_THREADS=20 $PYTHON_EXE -m psmm.bridges.embed2graph --init --clean

echo "=== Rebuild completed successfully! ==="
