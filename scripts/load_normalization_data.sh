#!/bin/bash

# Define paths
SOURCE_FILE="/local/storage/alen/projects/1_PSFD/output/10_gene_protein_normalization/940_gene_protein_normalization/gene_protein_entity_summary.csv.zst"
INPUT_DIR="/local/storage/thomas/5_PSMM/data/input"

echo "Checking if source file exists..."
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file $SOURCE_FILE not found!"
    exit 1
fi

echo "Cleaning input directories..."
mkdir -p "$INPUT_DIR"
rm -rf "$INPUT_DIR"/*

echo "Copying $SOURCE_FILE to target directory..."
cp "$SOURCE_FILE" "$INPUT_DIR/"

echo "Decompressing and removing .zst files..."
if [ -f "$INPUT_DIR/gene_protein_entity_summary.csv.zst" ]; then
    zstd -d "$INPUT_DIR/gene_protein_entity_summary.csv.zst"
    rm "$INPUT_DIR/gene_protein_entity_summary.csv.zst"
fi

echo "Done! The input database has been updated successfully."
