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

def run_command(cmd, cwd=PROJECT_DIR, env=None, quiet=False):
    if not quiet:
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

funcSearch Unified Command-Line Tool

PIPELINE REBUILD COMMANDS:
  rebuild full                 Sequential end-to-end rebuild (Graph + Search)
      [--sync-norm]            Sync normalization dataset first
  
  rebuild graph                Rebuild the UI Knowledge Graph and Enrichments
      compile                  (Subcommand) Only compile raw outputs into UI database
      enrich                   (Subcommand) Only calculate and update pathway enrichments
      
  rebuild search               Rebuild Sequence Search Databases (MMseqs2/PLM)
      fetch                    (Subcommand) Only fetch physical FASTA sequences
      index                    (Subcommand) Only build the MMseqs2/PLM search indices
      [--method METHOD]        Target specific search index (seq2graph, embed2graph, all)

SERVER COMMANDS:
  serve                        Start FastAPI Web Dispatcher
      [--port PORT]            Server port (default: 8999)
      [--host HOST]            Server host (default: 0.0.0.0)
      [--dev]                  Start server with code reloading active
      [--tunnel]               Expose the server publicly via ngrok tunnel
      [--ngrok-token TOKEN]    Ngrok authentication token

EVALUATION & BENCHMARKING:
  benchmark run                Run sequence masking evaluation
      [--query FILE]           FASTA file/folder containing queries
      [--output FILE]          Result file path
      [--method METHOD]        Search method to benchmark (seq2graph, embed2graph, both)
  benchmark evaluate           Evaluate benchmark results

DATA & ENRICHMENTS:
  enrich calculate             Run Chi-Squared / Fisher Exact enrichment tests
  enrich export                Export global path index enrichments to CSV
  enrich update-papers         Trickle down global enrichments into paper JSONs
  export                       Alias for 'enrich export'
  verify                       Run integration tests on the database
  fetch                        Standalone sequence fetcher from upstream sources
  db                           Standalone DB index builder
  build-graph                  Standalone raw output compilation tool
"""

    class CustomHelpParser(argparse.ArgumentParser):
        def print_help(self, file=None):
            if file is None:
                import sys
                file = sys.stdout
            if self.description:
                file.write(self.description + "\n")
            else:
                super().print_help(file)

    parser = CustomHelpParser(
        description=ASCII_ART,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. Serve
    serve_parser = subparsers.add_parser("serve", help="Start FastAPI Web Dispatcher")
    serve_parser.add_argument("--port", type=int, default=8999, help="Server port")
    serve_parser.add_argument("--host", default="0.0.0.0", help="Server host")
    serve_parser.add_argument("--dev", action="store_true", help="Start server with code reloading active")
    serve_parser.add_argument("--tunnel", action="store_true", help="Expose the server publicly via ngrok tunnel")
    serve_parser.add_argument("--ngrok-token", help="Ngrok authentication token (required for tunnel)")

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

    # 3. Fetch (Hidden)
    fetch_parser = subparsers.add_parser("fetch", help=argparse.SUPPRESS)
    fetch_parser.add_argument("--input", help="Custom list/dataframe of accessions to fetch")
    fetch_parser.add_argument("--output", help="Output directory for FASTA and sequence metadata")
    fetch_parser.add_argument("--no-cache", action="store_true", help="Force raw downloads")

    # 4. DB (Hidden)
    db_parser = subparsers.add_parser("db", help=argparse.SUPPRESS)
    db_parser.add_argument("--init", action="store_true", help="Build search indexes")
    db_parser.add_argument("--clean", action="store_true", help="Clear old index before build")
    db_parser.add_argument("--method", choices=["seq2graph", "embed2graph", "all"], default="all", help="Target specific search method index")

    # 5. Rebuild
    rebuild_parser = subparsers.add_parser("rebuild", help="Sequential pipeline rebuilds")
    rebuild_subparsers = rebuild_parser.add_subparsers(dest="rebuild_target")
    
    # rebuild full
    rebuild_full = rebuild_subparsers.add_parser("full", help="Sequential end-to-end rebuild")
    rebuild_full.add_argument("--sync-norm", action="store_true", help="Sync normalization zst dataset first")
    
    # rebuild graph
    rebuild_graph = rebuild_subparsers.add_parser("graph", help="Rebuild the UI Knowledge Graph and Enrichments")
    rebuild_graph.add_argument("step", nargs="?", choices=["all", "compile", "enrich"], default="all", help="Which step to run (default: all)")

    # rebuild search
    rebuild_search = rebuild_subparsers.add_parser("search", help="Rebuild Sequence Search Databases (MMseqs2/PLM)")
    rebuild_search.add_argument("step", nargs="?", choices=["all", "fetch", "index"], default="all", help="Which step to run (default: all)")
    rebuild_search.add_argument("--method", choices=["all", "seq2graph", "embed2graph"], default="all", help="Which search index to build (default: all)")

    # 6. Enrich (Hidden)
    enrich_parser = subparsers.add_parser("enrich", help=argparse.SUPPRESS)
    enrich_subparsers = enrich_parser.add_subparsers(dest="enrich_command")
    
    # enrich calculate
    enrich_calc_parser = enrich_subparsers.add_parser("calculate", help="Run Chi-Squared / Fisher pathway enrichment tests directly on knowledge graph")
    enrich_calc_parser.add_argument("--db", default=os.path.join(PROJECT_DIR, "data", "global_path_index.json"), help="Path to knowledge graph JSON")
    
    # enrich export
    enrich_export_parser = enrich_subparsers.add_parser("export", help="Export global path index enrichments to CSV")
    enrich_export_parser.add_argument("--output", default=os.path.join(PROJECT_DIR, "data", "enrichments.csv"), help="Output CSV path")

    # enrich update-papers
    enrich_subparsers.add_parser("update-papers", help="Trickle down global enrichments into individual paper JSONs")

    # 7. Verify
    subparsers.add_parser("verify", help="Run integration tests")

    # 8. Build Graph (Hidden)
    build_graph_parser = subparsers.add_parser("build-graph", help=argparse.SUPPRESS)
    build_graph_parser.add_argument("--outdir", default=os.path.join(PROJECT_DIR, "data"), help="Directory to write the compiled database")

    # 9. Export
    export_parser = subparsers.add_parser("export", help="Export pathway enrichments to CSV")
    export_parser.add_argument("--output", default=os.path.join(PROJECT_DIR, "data", "enrichments.csv"), help="Output CSV path")

    args = parser.parse_args()

    if args.command == "serve":
        cmd = [PYTHON_EXE, "-m", "uvicorn", "funcSearch.api.server:app"]
        cmd.extend(["--host", args.host])
        cmd.extend(["--port", str(args.port)])
        if args.dev:
            cmd.append("--reload")

        if args.tunnel:
            try:
                from pyngrok import ngrok
                token = args.ngrok_token or os.environ.get("NGROK_AUTHTOKEN") or "3H6shVlTJKQheeO4ONFLYD22TCM_5dY3jvyE1bFygntBKx5pM"
                if token:
                    ngrok.set_auth_token(token)
                public_url = ngrok.connect(args.port).public_url
                print(f"\n=======================================================")
                print(f" [TUNNEL ACTIVE] API exposed at: {public_url}")
                print(f"=======================================================\n")
            except ImportError:
                print("Error: --tunnel requires pyngrok. Install with: pip install pyngrok")
                sys.exit(1)
            except Exception as e:
                print(f"Error starting tunnel: {e}")
                sys.exit(1)

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
        cmd = [PYTHON_EXE, "-m", "funcSearch.fetcher.pipeline"]
        if args.no_cache:
            cmd.append("--force")
        if args.input:
            cmd.extend(["--input", args.input])
        if args.output:
            cmd.extend(["--output", args.output])
        run_command(cmd)

    elif args.command == "db":
        embed_env = {"OMP_NUM_THREADS": "20", "MKL_NUM_THREADS": "20"}

        if args.method in ["seq2graph", "all"]:
            cmd = [PYTHON_EXE, "-m", "funcSearch.bridges.seq2graph"]
            if args.init: cmd.append("--init")
            if args.clean: cmd.append("--clean")
            run_command(cmd)
        
        if args.method in ["embed2graph", "all"]:
            cmd = [PYTHON_EXE, "-m", "funcSearch.bridges.embed2graph"]
            if args.init: cmd.append("--init")
            if args.clean: cmd.append("--clean")
            run_command(cmd, env=embed_env)

    elif args.command == "rebuild":
        def run_funcSearch(*subcmd):
            run_command([PYTHON_EXE, os.path.abspath(__file__)] + list(subcmd))

        target = getattr(args, "rebuild_target", None)
        
        if target in ["full", "graph"]:
            step = getattr(args, "step", "all") if target == "graph" else "all"
            if step in ["all", "compile"]:
                print("\n=== [Knowledge Graph Pipeline] ===")
                print("--- Compiling raw outputs into global_path_index.json ---")
                run_funcSearch("build-graph")
            if step in ["all", "enrich"]:
                print("--- Calculating pathway enrichments ---")
                run_funcSearch("enrich", "calculate")
                print("--- Trickling enrichments into paper bundles ---")
                run_funcSearch("enrich", "update-papers")
            
        if target in ["full", "search"]:
            step = getattr(args, "step", "all") if target == "search" else "all"
            print("\n=== [Sequence Search Pipeline] ===")
            if target == "full" and getattr(args, "sync_norm", False):
                print("--- Loading new normalization data ---")
                run_command([os.path.join(PROJECT_DIR, "scripts", "load_normalization_data.sh")])
            if step in ["all", "fetch"]:
                print("--- Running sequence fetcher ---")
                run_funcSearch("fetch")
            if step in ["all", "index"]:
                print("--- Rebuilding Seq2Graph & Embed2Graph indices ---")
                run_funcSearch("db", "--init", "--clean", "--method", getattr(args, "method", "all"))
            
        if target not in ["full", "graph", "search"]:
            rebuild_parser.print_help()

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
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs): return iterable

with open('{PROJECT_DIR}/data/global_path_index.json', 'r') as f:
    data = json.load(f)
rows = []
for e in tqdm(data.get('entities', []), desc="Extracting enrichments"):
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
            print("Preparing to extract enrichments for export...")
            run_command([PYTHON_EXE, "-c", export_script], quiet=True)
        elif getattr(args, "enrich_command", None) == "calculate":
            cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "add_enrichments.py"), "--db", args.db]
            run_command(cmd)
        else:
            parser.parse_args(["enrich", "--help"])

    elif args.command == "verify":
        run_command([PYTHON_EXE, "-m", "funcSearch.api.verify"])

    elif args.command == "build-graph":
        cmd = [PYTHON_EXE, os.path.join(PROJECT_DIR, "scripts", "build_knowledge_graph.py"), "--outdir", args.outdir]
        run_command(cmd)

    elif args.command == "export":
        run_command([PYTHON_EXE, __file__, "enrich", "export", "--output", args.output])

if __name__ == "__main__":
    main()
