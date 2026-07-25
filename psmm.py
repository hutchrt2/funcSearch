#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

def get_python_exe():
    venv_python = os.path.join(PROJECT_DIR, ".venv", "bin", "python")
    if os.path.exists(venv_python):
        return venv_python
    return sys.executable

PYTHON_EXE = get_python_exe()

def run_command(cmd, cwd=PROJECT_DIR, env=None):
    print(f"Running: {' '.join(cmd)}")
    sys.stdout.flush()
    try:
        if env:
            current_env = os.environ.copy()
            current_env.update(env)
            subprocess.run(cmd, check=True, cwd=cwd, env=current_env)
        else:
            subprocess.run(cmd, check=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        sys.exit(e.returncode)

def main():
    ASCII_ART = """                  .::::..                                                :+.               
            :#@@@@@@@@@@@@@@@+.                  =+.                  .****=               
         +@@@+:            .-#@@%:               +***: =***+        :*******.+*:  =+       
      .%@%-                    .*@@+             -*****=*=****-    +*******=.-.-***+       
     #@%-@@@@@@@@@@@@@@@@@@@@@@@@#=@@-            ***=.:**-  ***: *********:.*******:      
   :@@- .@* -%%%%%%+   =%%%%%%= ##  #@#            :*+ :***:=**************: -******-      
  :@%.  .@* =%****%+   =+=##-*= ##   +@%.     .***=.*--.***+.=**--****-.**-+:+-*****.      
 .@@:   .@* =+    ++   =+    #= ##    *@*  :**=:***+ ++ ****+ +-.       =-****+  -*-     .:
 *@+    .@*  .+#%*:     .*#%*:  ##     %@:  :********-:+.*****--       ++*****+ =*.   :+*  
.@@.    .@@@@@@@@@@@@@@@@@@@@@@@@#     *@+    *********=*-+*+*-*=: :+--+******-.= .-***=   
.@@.    .@*   ....:-=@--:....   ##     =@*     .****-:*= +*******+.=  *******=: ******-    
.@@.    .@*   *#... :%. ...#+   ##     +@*       .-****+  +*********- ******=*******=      
 %@-    .@*   *#.::.:%.::: #+   ##     %@=  .:::--------=*+:-******+  ****=+++**+=.        
 -@%.    #%   *#-++::%.=+= #+  .@=    -@#    .************+ =**-. .+**++***:-**+-.         
  *@*    :@*  *#=##::%.### #+  ##    :@@:              .=+******++******=. =-+*********+=  
   *@#    :@#.*#--=++@+==--%+.%*    -@@:            :**********=   +-:::=+ =***********+=. 
    =@@-    *@*     :%.    .%@:    #@%            :**********+     +******=                
     .#@@=    =@%-  :%.  =%#:   .*@@-             =******+-.       :******+                
        *@@%-   .*@%+@*@%-   .+@@%:                                 -*****+                
          .*@@@@%=:.=+::-*@@@@%-                                     .****=                
              .:=*%@@@@%#*-.                                           :++                 

PSMM Unified Command-Line Tool"""

    parser = argparse.ArgumentParser(
        description=ASCII_ART,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. Serve
    serve_parser = subparsers.add_parser("serve", help="Start FastAPI Web Dispatcher")
    serve_parser.add_argument("--port", type=int, default=8999, help="Server port")
    serve_parser.add_argument("--host", default="0.0.0.0", help="Server host")
    serve_parser.add_argument("--dev", action="store_true", help="Start server with code reloading active")

    # 2. Benchmark
    bench_parser = subparsers.add_parser("benchmark", help="Run search & masking validation benchmarks")
    bench_subparsers = bench_parser.add_subparsers(dest="bench_command")
    
    # benchmark run (default)
    bench_run_parser = bench_subparsers.add_parser("run", help="Run sequence masking evaluation")
    bench_run_parser.add_argument("--query", help="FASTA file/folder containing queries")
    bench_run_parser.add_argument("--output", default="benchmark_results.csv", help="Result file path")
    bench_run_parser.add_argument("--method", choices=["seq2graph", "embed2graph", "both"], default="both", help="Search method(s) to benchmark")
    
    # benchmark evaluate
    bench_eval_parser = bench_subparsers.add_parser("evaluate", help="Evaluate benchmark results")

    # 3. Fetch
    fetch_parser = subparsers.add_parser("fetch", help="Fetch reference sequences (ETL)")
    fetch_parser.add_argument("--input", help="Custom list/dataframe of accessions to fetch")
    fetch_parser.add_argument("--no-cache", action="store_true", help="Force raw downloads")

    # 4. DB
    db_parser = subparsers.add_parser("db", help="Manage index and reference databases")
    db_parser.add_argument("--init", action="store_true", help="Build search indexes")
    db_parser.add_argument("--clean", action="store_true", help="Clear old index before build")
    db_parser.add_argument("--method", choices=["seq2graph", "embed2graph", "all"], default="all", help="Target specific search method index")

    # 5. Rebuild
    rebuild_parser = subparsers.add_parser("rebuild", help="Sequential end-to-end rebuild")
    rebuild_parser.add_argument("--sync-norm", action="store_true", help="Sync normalization zst dataset first")

    # 6. Enrich
    enrich_parser = subparsers.add_parser("enrich", help="Run pathway enrichment tests and tools")
    enrich_subparsers = enrich_parser.add_subparsers(dest="enrich_command")
    
    # enrich calculate
    enrich_calc_parser = enrich_subparsers.add_parser("calculate", help="Run Chi-Squared / Fisher pathway enrichment tests directly on knowledge graph")
    enrich_calc_parser.add_argument("--db", default=os.path.join(PROJECT_DIR, "data", "global_path_index.json"), help="Path to knowledge graph JSON")
    
    # enrich export
    enrich_export_parser = enrich_subparsers.add_parser("export", help="Export global path index enrichments to CSV")
    enrich_export_parser.add_argument("--output", default="enrichments.csv", help="Output CSV path")

    # enrich update-papers
    enrich_subparsers.add_parser("update-papers", help="Trickle down global enrichments into individual paper JSONs")

    # 7. Verify
    subparsers.add_parser("verify", help="Run integration tests")

    args = parser.parse_args()

    if args.command == "serve":
        cmd = [PYTHON_EXE, "-m", "uvicorn", "psmm.api.server:app"]
        cmd.extend(["--host", args.host])
        cmd.extend(["--port", str(args.port)])
        if args.dev:
            cmd.append("--reload")
        run_command(cmd)

    elif args.command == "benchmark":
        if getattr(args, "bench_command", "run") in ["run", None]:
            cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "benchmark_masking.py"), "benchmark"]
            if getattr(args, "query", None):
                cmd.extend(["--query", args.query])
            if getattr(args, "output", None):
                cmd.extend(["--output", args.output])
            if getattr(args, "method", None):
                if args.method == "both":
                    cmd.extend(["--methods", "seq2graph,embed2graph"])
                else:
                    cmd.extend(["--methods", args.method])
            run_command(cmd)
        elif args.bench_command == "evaluate":
            cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "evaluate_benchmarks.py")]
            run_command(cmd)

    elif args.command == "fetch":
        cmd = [PYTHON_EXE, "-m", "psmm.fetcher.pipeline"]
        if args.no_cache:
            cmd.append("--force")
        if args.input:
            cmd.extend(["--input", args.input])
        run_command(cmd)

    elif args.command == "db":
        mmseqs_env = {"PATH": f"/programs/mmseqs/bin:{os.environ.get('PATH', '')}"}
        embed_env = mmseqs_env.copy()
        embed_env.update({"OMP_NUM_THREADS": "20", "MKL_NUM_THREADS": "20"})

        if args.method in ["seq2graph", "all"]:
            cmd = [PYTHON_EXE, "-m", "psmm.bridges.seq2graph"]
            if args.init: cmd.append("--init")
            if args.clean: cmd.append("--clean")
            run_command(cmd, env=mmseqs_env)
        
        if args.method in ["embed2graph", "all"]:
            cmd = [PYTHON_EXE, "-m", "psmm.bridges.embed2graph"]
            if args.init: cmd.append("--init")
            if args.clean: cmd.append("--clean")
            run_command(cmd, env=embed_env)

    elif args.command == "rebuild":
        if args.sync_norm:
            run_command([os.path.join(PROJECT_DIR, "scripts", "load_normalization_data.sh")])
        run_command([os.path.join(PROJECT_DIR, "scripts", "full_rebuild_pipeline.sh")])

    elif args.command == "enrich":
        if getattr(args, "enrich_command", None) == "update-papers":
            db_path = os.path.join(PROJECT_DIR, "data", "global_path_index.json")
            papers_path = os.path.join(PROJECT_DIR, "data", "papers")
            cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "update_paper_enrichments.py"), "--db", db_path, "--papers", papers_path]
            run_command(cmd)
        elif getattr(args, "enrich_command", None) == "export":
            # Simple wrapper to dump global_path_index enrichments
            export_script = f"""
import json
import csv
with open('{PROJECT_DIR}/data/global_path_index.json', 'r') as f:
    data = json.load(f)
rows = []
for e in data.get('entities', []):
    if e.get('enrichments'):
        for enr in e['enrichments']:
            rows.append({{'entity_id': e['id'], 'ontology_id': e.get('ontology_id'), 'trait_concept': enr.get('trait_concept'), 'trait_label': enr.get('trait_label'), 'p_value': enr.get('p_value'), 'fdr': enr.get('fdr')}})
if rows:
    with open('{args.output}', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"Exported {{len(rows)}} enrichments to {args.output}")
else:
    print("No enrichments found to export.")
"""
            run_command([PYTHON_EXE, "-c", export_script])
        elif getattr(args, "enrich_command", None) == "calculate":
            cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "add_enrichments.py"), "--db", args.db]
            run_command(cmd)
        else:
            parser.parse_args(["enrich", "--help"])

    elif args.command == "verify":
        run_command([PYTHON_EXE, "-m", "psmm.api.verify"])

if __name__ == "__main__":
    main()
