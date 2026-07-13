#!/bin/bash
PORT=8999
cd /local/storage/thomas/5_PSMM 2>/dev/null || cd "$(dirname "$0")/.."

echo "Checking for zombie processes on port $PORT..."
fuser -k ${PORT}/tcp 2>/dev/null
sleep 1

PYTHON_EXE="python"
if [ -f ".venv/bin/python" ]; then
    PYTHON_EXE=".venv/bin/python"
fi

echo "Starting Uvicorn API Dispatcher on port $PORT..."
$PYTHON_EXE -m uvicorn psmm.api.server:app --host 0.0.0.0 --port $PORT
