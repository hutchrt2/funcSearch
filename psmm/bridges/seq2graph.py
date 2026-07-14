#!/usr/bin/env python3
"""
seq2graph_bridge.py

A command-line Python wrapper utility that automates sequence alignment searches
via MMseqs2 and maps the alignment hits back to a relational knowledge graph schema.
"""

import os
import sys
import shutil
import argparse
import subprocess
import glob
import pandas as pd
from typing import List

def get_mmseqs_binary() -> str:
    """
    Locates the mmseqs executable on the system path, virtualenv, or common conda paths.
    """
    env_bin = os.getenv("MMSEQS_BINARY")
    if env_bin:
        return env_bin
        
    which_bin = shutil.which("mmseqs")
    if which_bin:
        return which_bin
        
    # Check same directory as python interpreter (conda/virtualenv env)
    py_dir = os.path.dirname(sys.executable)
    conda_bin = os.path.join(py_dir, "mmseqs")
    if os.path.exists(conda_bin):
        return conda_bin
        
    for path in [
        os.path.expanduser("~/bin/mmseqs"),
        os.path.expanduser("~/.local/bin/mmseqs"),
        "/usr/local/bin/mmseqs",
        "/usr/bin/mmseqs"
    ]:
        if os.path.exists(path):
            return path
            
    return "mmseqs"

MMSEQS_BIN = get_mmseqs_binary()

def run_command(cmd: List[str], cwd: str) -> None:
    """
    Run a shell command using subprocess.run, ensuring standard output and
    standard error are captured and logged on failure.

    Args:
        cmd: List of command arguments.
        cwd: Directory where the command should be executed.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    print(f"Running command: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd
        )
        if result.stdout:
            print(result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Return code: {e.returncode}", file=sys.stderr)
        if e.stdout:
            print("--- Standard Output ---", file=sys.stderr)
            print(e.stdout, file=sys.stderr)
        if e.stderr:
            print("--- Standard Error ---", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
        raise e

def get_convert_command() -> str:
    """
    Determine whether 'convertalignment' or 'convertalis' is supported by the mmseqs binary.

    Returns:
        str: 'convertalignment' or 'convertalis'.
    """
    try:
        result = subprocess.run(
            [MMSEQS_BIN, "convertalignment", "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if result.returncode == 0:
            return "convertalignment"
    except Exception:
        pass
    return "convertalis"

def init_database(script_dir: str, clean: bool = False) -> None:
    """
    Initialize the local MMseqs2 database environment.
    Creates the database from reference sequences and builds the index.

    Args:
        script_dir: The directory of the script to resolve relative paths.
        clean: Delete old database files first.
    """
    project_dir = os.path.dirname(os.path.dirname(script_dir))
    db_dir = os.path.join(project_dir, "data", "blastdb")
    if clean:
        print("Cleaning old database files for a fresh start...")
        db_dir_abs = db_dir
        if os.path.exists(db_dir_abs):
            for f in os.listdir(db_dir_abs):
                f_path = os.path.join(db_dir_abs, f)
                if os.path.isfile(f_path):
                    try:
                        os.remove(f_path)
                    except Exception as e:
                        print(f"Warning: Failed to delete {f_path}: {e}")

    print("Initializing MMseqs2 database...")
    input_fasta = os.path.join(project_dir, "data", "build", "psfd_sequences.fasta")
    db_path = os.path.join(db_dir, "psfd_protein_db")
    tmp_dir = os.path.join(project_dir, "data", "tmp")

    # Ensure blastdb and tmp directories exist
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # 1. mmseqs createdb input_database/psfd_sequences.fasta blastdb/psfd_protein_db
    createdb_cmd = [
        MMSEQS_BIN, "createdb",
        input_fasta,
        db_path
    ]
    run_command(createdb_cmd, cwd=script_dir)

    # 2. mmseqs createindex blastdb/psfd_protein_db tmp/
    createindex_cmd = [
        MMSEQS_BIN, "createindex",
        db_path,
        tmp_dir
    ]
    run_command(createindex_cmd, cwd=script_dir)
    print("Database initialization completed successfully.")

def check_db_exists(script_dir: str) -> bool:
    """
    Check if the MMseqs2 database files exist.

    Args:
        script_dir: The directory of the script to resolve relative paths.

    Returns:
        bool: True if the database files exist, False otherwise.
    """
    project_dir = os.path.dirname(os.path.dirname(script_dir))
    db_path = os.path.join(project_dir, "data", "blastdb", "psfd_protein_db")
    return os.path.exists(db_path) or os.path.exists(db_path + ".dbtype")

def extract_uniprot_id(target_id: str) -> str:
    """
    Extract the UniProt ID from an MMseqs target identifier.
    
    If the identifier is in a pipe-delimited format (e.g., 'sp|Q9FL62|WRK30_ARATH'),
    the UniProt ID is the second element (index 1). Otherwise, returns the target_id as is.

    Args:
        target_id: The target sequence identifier string.

    Returns:
        The extracted UniProt ID.
    """
    if not isinstance(target_id, str):
        return target_id
    if "|" in target_id:
        parts = target_id.split("|")
        if len(parts) >= 2:
            return parts[1]
    return target_id

def clean_tmp_files(script_dir: str) -> None:
    """
    Clean up temporary query database and alignment files from the tmp folder.
    """
    project_dir = os.path.dirname(os.path.dirname(script_dir))
    tmp_path = os.path.join(project_dir, "data", "tmp")
    if os.path.exists(tmp_path):
        for pattern in ["query_db*", "results_db*", "output_hits.tab"]:
            for f in glob.glob(os.path.join(tmp_path, pattern)):
                try:
                    if os.path.isfile(f):
                        os.remove(f)
                except Exception as e:
                    print(f"Warning: Failed to remove temporary file {f}: {e}", file=sys.stderr)

def run_single_query(query_path: str, script_dir: str, evalue: float = None, min_seq_id: float = None) -> pd.DataFrame:
    """
    Run the MMseqs2 alignment search pipeline for a single query FASTA file.

    Args:
        query_path: Path to the query FASTA file.
        script_dir: Directory where the script runs.
        evalue: Expectation value threshold.
        min_seq_id: Minimum sequence identity threshold.

    Returns:
        DataFrame containing the raw alignment outputs.
    """
    clean_tmp_files(script_dir)

    project_dir = os.path.dirname(os.path.dirname(script_dir))
    tmp_dir = os.path.join(project_dir, "data", "tmp")
    tmp_query_db = os.path.join(tmp_dir, "query_db")
    db_path = os.path.join(project_dir, "data", "blastdb", "psfd_protein_db")
    results_db = os.path.join(tmp_dir, "results_db")
    output_hits_tab = os.path.join(tmp_dir, "output_hits.tab")

    # Ensure tmp directory exists
    os.makedirs(tmp_dir, exist_ok=True)

    # 1. Create query database
    print(f"Converting query FASTA '{os.path.basename(query_path)}' to MMseqs database...")
    createdb_cmd = [
        MMSEQS_BIN, "createdb",
        query_path,
        tmp_query_db
    ]
    run_command(createdb_cmd, cwd=script_dir)

    # 2. Run search
    search_cmd = [
        MMSEQS_BIN, "search",
        tmp_query_db,
        db_path,
        results_db,
        tmp_dir
    ]
    if evalue is not None:
        search_cmd.extend(["-e", str(evalue)])
    if min_seq_id is not None:
        search_cmd.extend(["--min-seq-id", str(min_seq_id)])
        
    print(f"Running MMseqs2 sequence search for '{os.path.basename(query_path)}'...")
    run_command(search_cmd, cwd=script_dir)

    # 3. Convert alignment
    convert_cmd_name = get_convert_command()
    convert_cmd = [
        MMSEQS_BIN, convert_cmd_name,
        tmp_query_db,
        db_path,
        results_db,
        output_hits_tab,
        "--format-output", "query,target,pident,evalue,qcov,tcov"
    ]
    print(f"Converting alignment to tab format using '{convert_cmd_name}'...")
    run_command(convert_cmd, cwd=script_dir)

    output_hits_abs = os.path.join(script_dir, output_hits_tab)
    if not os.path.exists(output_hits_abs) or os.path.getsize(output_hits_abs) == 0:
        print(f"No alignment hits found for '{os.path.basename(query_path)}'.")
        return pd.DataFrame(columns=["query", "target", "pident", "evalue", "qcov", "tcov"])

    alignment_cols = ["query", "target", "pident", "evalue", "qcov", "tcov"]
    return pd.read_csv(output_hits_abs, sep="\t", names=alignment_cols)

def process_query(query_path: str, output_path: str, script_dir: str, evalue: float = None, min_seq_id: float = None) -> None:
    """
    Process the query sequence (or directory of sequences): search against reference database,
    convert alignment, perform relational join with metadata, and serialize output.

    Args:
        query_path: Absolute path to the query FASTA file or directory containing FASTA files.
        output_path: Absolute path to the output CSV file.
        script_dir: Absolute path to the script directory.
        evalue: Expectation value threshold.
        min_seq_id: Minimum sequence identity threshold.
    """
    # 1. Collect all query files to process
    query_files = []
    if os.path.isdir(query_path):
        print(f"Scanning directory '{query_path}' for FASTA files...")
        for name in sorted(os.listdir(query_path)):
            if name.lower().endswith((".fasta", ".fa", ".faa")):
                query_files.append(os.path.join(query_path, name))
        if not query_files:
            print(f"No FASTA files (.fasta, .fa, .faa) found in directory: {query_path}", file=sys.stderr)
            # Create an empty CSV with correct columns
            empty_df = pd.DataFrame(columns=[
                "query", "uniprot_id", "pident", "evalue", "qcov", "tcov",
                "global_node_id", "selected_protein_name", "selected_gene_name", "selected_organism"
            ])
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            empty_df.to_csv(output_path, index=False)
            print(f"Empty results table saved to: {output_path}")
            return
    else:
        query_files = [query_path]

    # 2. Run sequential alignment searches
    all_alignments = []
    for f in query_files:
        print(f"Processing query file: {f}")
        df = run_single_query(f, script_dir, evalue=evalue, min_seq_id=min_seq_id)
        if not df.empty:
            all_alignments.append(df)

    # Clean up temp files after processing is completed
    clean_tmp_files(script_dir)

    # 3. Consolidate and map results
    if not all_alignments:
        print("No alignment hits found across all queries.")
        empty_df = pd.DataFrame(columns=[
            "query", "uniprot_id", "pident", "evalue", "qcov", "tcov",
            "global_node_id", "selected_protein_name", "selected_gene_name", "selected_organism"
        ])
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        empty_df.to_csv(output_path, index=False)
        print(f"Empty results table saved to: {output_path}")
        return

    alignment_df = pd.concat(all_alignments, ignore_index=True)

    # Load metadata
    project_dir = os.path.dirname(os.path.dirname(script_dir))
    metadata_csv_path = os.path.join(project_dir, "data", "build", "sequence_metadata.csv")
    if not os.path.exists(metadata_csv_path):
        print(f"Error: Metadata file not found at {metadata_csv_path}", file=sys.stderr)
        sys.exit(1)
    
    metadata_df = pd.read_csv(metadata_csv_path)
    if "target_accession" in metadata_df.columns:
        metadata_df = metadata_df.rename(columns={"target_accession": "uniprot_id"})

    # Check if alignment target corresponds to global_node_id or uniprot_id
    sample_targets = alignment_df["target"].dropna().unique()
    is_global_node_id = any(val in metadata_df["global_node_id"].values for val in sample_targets)

    if is_global_node_id:
        print("Detected global_node_id targets. Performing join on global_node_id...")
        alignment_df = alignment_df.rename(columns={"target": "global_node_id"})
        merged_df = pd.merge(alignment_df, metadata_df, on="global_node_id", how="inner")
    else:
        print("Detected standard target IDs. Performing join on uniprot_id...")
        alignment_df = alignment_df.rename(columns={"target": "uniprot_id"})
        alignment_df["uniprot_id"] = alignment_df["uniprot_id"].apply(extract_uniprot_id)
        merged_df = pd.merge(alignment_df, metadata_df, on="uniprot_id", how="inner")

    # Map pident to score (0-1 range)
    if "pident" in merged_df.columns:
        merged_df["score"] = merged_df["pident"] / 100.0
    else:
        merged_df["score"] = 1.0
    merged_df["score_type"] = "k-mer score"

    # Select and order final columns
    final_cols = [
        "query",
        "uniprot_id",
        "pident",
        "evalue",
        "qcov",
        "tcov",
        "global_node_id",
        "selected_protein_name",
        "selected_gene_name",
        "selected_organism",
        "score",
        "score_type"
    ]
    # Filter columns to only what we need, handling case where metadata columns might slightly differ
    for col in final_cols:
        if col not in merged_df.columns:
            merged_df[col] = None
    
    final_df = merged_df[final_cols]

    # Save to final output path
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_df.to_csv(output_path, index=False)
    print(f"All alignment hits successfully mapped to metadata and written to: {output_path}")

def main() -> None:
    """
    Main CLI entry point for the seq2graph bridge utility.
    """
    parser = argparse.ArgumentParser(
        description="MMseqs2 sequence search to relational knowledge graph bridge."
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Initialize/re-initialize the reference sequence database."
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete old database files before starting initialization."
    )
    parser.add_argument(
        "--query",
        type=str,
        default="input_FASTA",
        help="Path to the query protein sequence FASTA file or folder (default: input_FASTA)."
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=None,
        help="MMseqs2 search expectation value threshold (default: 0.001)."
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=None,
        help="MMseqs2 search minimum sequence identity threshold (0.0 to 1.0)."
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to the output CSV file."
    )

    args = parser.parse_args()

    # Determine script execution path
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if we should initialize database
    db_missing = not check_db_exists(script_dir)
    if args.init or db_missing:
        if db_missing and not args.init:
            print("Reference database not found. Initializing...")
        init_database(script_dir, clean=args.clean)
        # If user explicitly requested only --init (without --query), exit here
        if args.init and "--query" not in sys.argv:
            return

    # Resolve query paths (checking current dir or script dir)
    query_abs = os.path.abspath(args.query)
    if not os.path.exists(query_abs):
        query_in_script_dir = os.path.join(script_dir, args.query)
        if os.path.exists(query_in_script_dir):
            query_abs = query_in_script_dir
        else:
            print(f"Error: Query path not found at {query_abs} or {query_in_script_dir}", file=sys.stderr)
            sys.exit(1)

    if args.output:
        output_abs = os.path.abspath(args.output)
    else:
        output_abs = os.path.join(script_dir, "output", "blast_query_results.csv")

    process_query(query_abs, output_abs, script_dir, evalue=args.evalue, min_seq_id=args.min_seq_id)

if __name__ == "__main__":
    main()
