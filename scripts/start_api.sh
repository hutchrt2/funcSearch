#!/bin/bash
cd /local/storage/thomas/5_PSMM
PYTHON_EXE="python"
if [ -f ".venv/bin/python" ]; then
    PYTHON_EXE=".venv/bin/python"
fi

$PYTHON_EXE -m uvicorn psmm.api.server:app --host 0.0.0.0 --port 8000
