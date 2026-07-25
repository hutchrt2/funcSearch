#!/bin/bash

# Define paths
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SOURCE_BASE="${PROJECT_DIR}/input/gene_protein_entity_summary"
INPUT_DIR="${PROJECT_DIR}/data/input"

echo "Checking for source file..."
if [ -f "${SOURCE_BASE}.csv.zst" ]; then
    SOURCE_FILE="${SOURCE_BASE}.csv.zst"
    IS_ZST=true
elif [ -f "${SOURCE_BASE}.csv" ]; then
    SOURCE_FILE="${SOURCE_BASE}.csv"
    IS_ZST=false
else
    echo "Error: Neither ${SOURCE_BASE}.csv.zst nor ${SOURCE_BASE}.csv found!"
    exit 1
fi

echo "Cleaning input directories..."
mkdir -p "$INPUT_DIR"
rm -rf "$INPUT_DIR"/*

echo "Copying $SOURCE_FILE to target directory..."
cp "$SOURCE_FILE" "$INPUT_DIR/"

if [ "$IS_ZST" = true ]; then
    echo "Decompressing .zst file..."
    zstd -d "$INPUT_DIR/gene_protein_entity_summary.csv.zst"
    rm "$INPUT_DIR/gene_protein_entity_summary.csv.zst"
fi

echo "Done! The input database has been updated successfully."
