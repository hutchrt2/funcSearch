import os
import re
import glob
import gzip
import time
import csv
import argparse
import hashlib
import traceback
from datetime import datetime
import pandas as pd
import requests

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs): return iterable
    tqdm.write = print

def get_sequence_hash(sequence: str) -> str:
    """Generate SHA-256 hash of the sequence string."""
    clean_seq = re.sub(r'\s+', '', sequence.upper())
    return hashlib.sha256(clean_seq.encode('utf-8')).hexdigest()

def parse_accession_from_header(header: str) -> str:
    """Extract accession ID from our custom FASTA header."""
    # Custom format: >{global_node_id} | {source_database}:{accession} | {selected_organism}
    h = header.lstrip('>').strip()
    if '|' in h:
        parts = [p.strip() for p in h.split('|')]
        if len(parts) >= 2:
            db_acc = parts[1]
            if ':' in db_acc:
                return db_acc.split(':', 1)[1].strip().upper()
            return db_acc.strip().upper()
    return h.split()[0].strip().upper()

def parse_raw_uniprot_accession(header: str) -> str:
    """Extract accession from standard UniProt/NCBI headers (e.g. >sp|Q9FL62|WRK30 or >NP_568439.1)"""
    h = header.lstrip('>').strip()
    if '|' in h:
        parts = h.split('|')
        if len(parts) >= 2:
            return parts[1].strip().upper()
    return h.split()[0].strip().upper()

def split_multi_fasta(multi_fasta: str) -> dict:
    """Helper to split a multi-FASTA string into {accession: sequence}"""
    records = {}
    current_acc = None
    current_lines = []
    
    for line in multi_fasta.strip().split("\n"):
        if line.startswith(">"):
            if current_acc:
                records[current_acc] = "\n".join(current_lines)
            # Try to get the raw accession from standard uniprot fasta header
            h = line.lstrip('>').strip()
            if '|' in h:
                parts = h.split('|')
                if len(parts) >= 2:
                    current_acc = parts[1].strip().upper()
                else:
                    current_acc = h.split()[0].strip().upper()
            else:
                current_acc = h.split()[0].strip().upper()
            current_lines = []
        else:
            if current_acc:
                current_lines.append(line.strip())
                
    if current_acc:
        records[current_acc] = "\n".join(current_lines)
        
    return records

class BaseFetcher:
    def __init__(self, output_fasta, metadata_csv, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats=None):
        self.output_fasta = output_fasta
        self.metadata_csv = metadata_csv
        self.log_filename = log_filename
        self.accession_to_hash = accession_to_hash
        self.sequence_hashes = sequence_hashes
        self.failed_accessions = failed_accessions
        self.stats = stats if stats is not None else {}
        self.session = requests.Session()
        
    def fetch_queue(self, queue_df):
        """Must be implemented by subclasses."""
        raise NotImplementedError

    def _exponential_backoff_request(self, method, url, max_attempts=5, timeout=15, **kwargs):
        """Polite exponential backoff for REST APIs."""
        for attempt in range(1, max_attempts + 1):
            try:
                if method.lower() == 'post':
                    response = self.session.post(url, timeout=timeout, **kwargs)
                else:
                    response = self.session.get(url, timeout=timeout, **kwargs)
                    
                if response.status_code == 200:
                    return response
                elif response.status_code in [429, 500, 502, 503, 504]:
                    wait_time = 2 ** attempt
                    time.sleep(wait_time)
                else:
                    return f"HTTP status {response.status_code}"
            except requests.exceptions.RequestException as error:
                if attempt == max_attempts:
                    return f"Network error: {error}"
                wait_time = 2 ** attempt
                time.sleep(wait_time)
        return "Unknown error"

    def save_sequence(self, target_db, accession, sequence, global_node_id, organism):
        """Saves physical sequence to FASTA if unique and updates the accession map."""
        if not sequence:
            self.log_failure(target_db, accession, "Empty sequence returned")
            return
            
        seq_hash = get_sequence_hash(sequence)
        self.accession_to_hash[accession] = seq_hash
        
        # Clean organism name to avoid nan/none
        clean_organism = ""
        if pd.notna(organism) and organism is not None:
            o_str = str(organism).strip()
            if o_str.lower() not in ['nan', 'none']:
                clean_organism = o_str
        
        # Check Tier 2 cache
        if seq_hash not in self.sequence_hashes:
            self.sequence_hashes.add(seq_hash)
            header = f">{global_node_id} | {target_db}:{accession} | {clean_organism}"
            fasta_record = f"{header}\n{sequence.strip()}\n"
            with open(self.output_fasta, "a", encoding="utf-8") as f:
                f.write(fasta_record)
        else:
            print(f"[{target_db}] {accession}: Sequence duplicate found via hash, skipped FASTA append.")
            if 'deduplicated_sequences' in self.stats:
                self.stats['deduplicated_sequences'] += 1
                
        if 'api_successes' in self.stats:
            self.stats['api_successes'] += 1
        if 'succeeded_list' in self.stats:
            self.stats['succeeded_list'].append((accession, target_db, global_node_id))
        if 'resolved_accessions' in self.stats:
            self.stats['resolved_accessions'][accession] = (target_db, accession, seq_hash)

    def log_failure(self, db_name, accession, error_msg):
        self.failed_accessions.add(accession)
        # Note: We do not write intermediate database attempt failures to the log file.
        # The log file only records final errors when all attempts for an entity are exhausted.
        print(f"[{db_name}] Failed to fetch {accession}: {error_msg}")
        if 'api_failures' in self.stats:
            self.stats['api_failures'] += 1
        if 'failed_list' in self.stats:
            self.stats['failed_list'].append((accession, db_name, error_msg))

class LocalCacheFetcher(BaseFetcher):
    def __init__(self, output_fasta, metadata_csv, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats=None, search_dirs=None):
        super().__init__(output_fasta, metadata_csv, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats)
        self.search_dirs = search_dirs or []
        self.local_sequences = {}
        self._load_local_fastas()

    def _load_local_fastas(self):
        print("[LocalCacheFetcher] Scanning local directories for FASTA files...")
        found_files = []
        for d in self.search_dirs:
            if not d or not os.path.exists(d):
                continue
            for root, _, files in os.walk(d):
                for f in files:
                    if f.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
                        found_files.append(os.path.join(root, f))
        
        count = 0
        for filepath in found_files:
            try:
                open_fn = gzip.open if filepath.endswith('.gz') else open
                with open_fn(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
                    current_header = None
                    seq_lines = []
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_header and seq_lines:
                                self._index_header(current_header, "".join(seq_lines))
                                count += 1
                            current_header = line[1:].strip()
                            seq_lines = []
                        else:
                            seq_lines.append(line)
                    if current_header and seq_lines:
                        self._index_header(current_header, "".join(seq_lines))
                        count += 1
            except Exception as e:
                print(f"[LocalCacheFetcher] Warning: failed to read {filepath}: {e}")
        print(f"[LocalCacheFetcher] Indexed {count} headers across {len(found_files)} FASTA files.")

    def _index_header(self, header, sequence):
        raw_id = header.split()[0].upper()
        self.local_sequences[raw_id] = sequence
        if '_' in raw_id:
            parts = raw_id.split('_', 1)
            sub_id = parts[1]
            self.local_sequences[sub_id] = sequence
            if '.' in sub_id:
                self.local_sequences[sub_id.split('.')[0]] = sequence
        if '.' in raw_id:
            self.local_sequences[raw_id.split('.')[0]] = sequence

    def fetch_queue(self, queue_df):
        for idx, row in queue_df.iterrows():
            orig_acc = str(row['target_accession']).strip()
            acc = orig_acc.upper()
            seq = self.local_sequences.get(acc)
            if not seq and '.' in acc:
                seq = self.local_sequences.get(acc.split('.')[0])
            if not seq and '_' in acc:
                seq = self.local_sequences.get(acc.split('_', 1)[1])
                
            if seq:
                self.save_sequence('LocalCache', orig_acc, seq, row.get('global_node_id', ''), row.get('selected_organism', ''))
                print(f"[LocalCacheFetcher] Success: Found {orig_acc} in local FASTA cache!")
            else:
                self.log_failure('LocalCache', orig_acc, "Not found in local FASTA files")

class UniProtFetcher(BaseFetcher):
    def fetch_queue(self, queue_df):
        batch_size = 100
        for i in tqdm(range(0, len(queue_df), batch_size), desc="Fetching UniProt batches"):
            batch_df = queue_df.iloc[i:i+batch_size]
            batch_accessions = batch_df['target_accession'].tolist()
            
            url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={','.join(batch_accessions)}&format=fasta"
            
            response = self._exponential_backoff_request('get', url, max_attempts=5, timeout=20)
            
            if isinstance(response, requests.Response) and response.status_code == 200:
                print("[UniProt] Batch successful! Parsing sequences...")
                fetched_dict = split_multi_fasta(response.text)
                
                for _, row in batch_df.iterrows():
                    acc = row['target_accession']
                    if acc in fetched_dict:
                        clean_seq = fetched_dict[acc]
                        self.save_sequence('UniProt', acc, clean_seq, row.get('global_node_id', ''), row.get('selected_organism', ''))
                    else:
                        print(f"[UniProt] Batch missed {acc}. Confirmed inactive or deleted in UniProt.")
                        self.log_failure('UniProt', acc, "Empty sequence returned (deleted/inactive)")
                time.sleep(1)
            else:
                print("[UniProt] Batch network failed! Triggering Sequential Fallback...")
                for _, row in batch_df.iterrows():
                    self._fetch_single(row)
                    
    def _fetch_single(self, row):
        acc = row["target_accession"]
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        response = self._exponential_backoff_request('get', url)
        if isinstance(response, requests.Response) and response.status_code == 200:
            seq_lines = response.text.strip().split('\n')[1:]
            clean_seq = '\n'.join(seq_lines)
            if clean_seq.strip():
                self.save_sequence('UniProt', acc, clean_seq, row.get('global_node_id', ''), row.get('selected_organism', ''))
                print(f"[UniProt] Success: Fetched {acc}")
            else:
                self.log_failure('UniProt', acc, "Empty sequence returned")
        else:
            error_info = response if isinstance(response, str) else "Failed"
            self.log_failure('UniProt', acc, error_info)
        time.sleep(1)

class NCBIFetcher(BaseFetcher):
    def fetch_queue(self, queue_df):
        ncbi_key = os.environ.get("NCBI_API_KEY")
        key_param = f"&api_key={ncbi_key}" if ncbi_key else ""
        sleep_time = 0.15 if ncbi_key else 0.4
        
        for _, row in tqdm(queue_df.iterrows(), total=len(queue_df), desc="Fetching NCBI sequences"):
            acc = row['target_accession']
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id={acc}{key_param}"
            response = self._exponential_backoff_request('get', url)
            
            if isinstance(response, requests.Response) and response.status_code == 200:
                text = response.text.strip()
                if text and text.startswith('>'):
                    seq_lines = text.split('\n')[1:]
                    clean_seq = '\n'.join(seq_lines)
                    self.save_sequence('NCBI', acc, clean_seq, row.get('global_node_id', ''), row.get('selected_organism', ''))
                    print(f"[NCBI] Success: Fetched {acc}")
                else:
                    self.log_failure('NCBI', acc, "Invalid FASTA returned")
            else:
                error_info = response if isinstance(response, str) else "Failed"
                self.log_failure('NCBI', acc, error_info)
            time.sleep(sleep_time)

class EnsemblFetcher(BaseFetcher):
    def fetch_queue(self, queue_df):
        headers = {"Content-Type": "text/x-fasta"}
        for _, row in tqdm(queue_df.iterrows(), total=len(queue_df), desc="Fetching Ensembl sequences"):
            acc = row['target_accession']
            url = f"https://rest.ensembl.org/sequence/id/{acc}?type=protein"
            response = self._exponential_backoff_request('get', url, headers=headers)
            
            if isinstance(response, requests.Response) and response.status_code == 200:
                text = response.text.strip()
                if text and text.startswith('>'):
                    seq_lines = text.split('\n')[1:]
                    clean_seq = '\n'.join(seq_lines)
                    self.save_sequence('Ensembl', acc, clean_seq, row.get('global_node_id', ''), row.get('selected_organism', ''))
                    print(f"[Ensembl] Success: Fetched {acc}")
                else:
                    self.log_failure('Ensembl', acc, "Invalid FASTA returned")
            else:
                error_info = response if isinstance(response, str) else "Failed"
                self.log_failure('Ensembl', acc, error_info)
            time.sleep(1)

class PhytozomeFetcher(BaseFetcher):
    def _normalize_query(self, query: str) -> str:
        q = query.strip()
        
        # Automated MSU (LOC_Os) to RAP-DB (Os) mapping
        loc_to_rap = {
            'LOC_Os01g01840': 'Os01g0108400',
            'LOC_Os01g06220': 'Os01g0155000',
            'LOC_Os01g16010': 'Os01g0264900',
            'LOC_Os01g22370': 'Os01g0327400',
            'LOC_Os01g24820': 'Os01g0350300',
            'LOC_Os01g50770': 'Os01g0703600',
            'LOC_Os01g56880': 'Os01g0776600',
            'LOC_Os01g63220': 'Os01g0851000',
            'LOC_Os01g64360': 'Os01g0863300',
            'LOC_Os01g74410': 'Os01g0975300',
            'LOC_Os02g05880': 'Os02g0152700',
            'LOC_Os02g15580': 'Os02g0255000',
            'LOC_Os02g28980': 'Os02g0491400',
            'LOC_Os02g30060': 'Os02g0503500',
            'LOC_Os02g33850': 'Os02g0543300',
            'LOC_Os02g40130': 'Os02g0614966',
            'LOC_Os02g40190': 'Os02g0615400',
            'LOC_Os02g45170': 'Os02g0673500',
            'LOC_Os02g46420': 'Os02g0689500',
            'LOC_Os02g53340': 'Os02g0773400',
            'LOC_Os04g01740': 'Os04g0107900',
            'LOC_Os04g12980': 'Os04g0206700',
            'LOC_Os04g30760': 'Os04g0376300',
            'LOC_Os04g32950': 'Os04g0402100',
            'LOC_Os04g33240': 'Os04g0405300',
            'LOC_Os04g42950': 'Os04g0508500',
            'LOC_Os04g58840': 'Os04g0685200',
            'LOC_Os04g59060': 'Os04g0687300',
            'LOC_Os04g59420': 'Os04g0690500',
            'LOC_Os05g01140': 'Os05g0102000',
            'LOC_Os05g05740': 'Os05g0150000',
            'LOC_Os05g06140': 'Os05g0153300',
            'LOC_Os05g06700': 'Os05g0159100',
            'LOC_Os05g09640': 'Os05g0188700',
            'LOC_Os05g34310': 'Os05g0415400',
            'LOC_Os05g41950': 'Os05g0498900',
            'LOC_Os05g43810': 'Os05g0513700',
            'LOC_Os05g43870': 'Os05g0514500',
            'LOC_Os05g51650': 'Os05g0594900',
            'LOC_Os06g08850': 'Os06g0188000',
            'LOC_Os06g10580': 'Os06g0207700',
            'LOC_Os06g11510': 'Os06g0218800',
            'LOC_Os06g38120': 'Os06g0579200',
            'LOC_Os07g02570': 'Os07g0116900',
            'LOC_Os07g03920': 'Os07g0131375',
            'LOC_Os07g09000': 'Os07g0187700',
            'LOC_Os07g46660': 'Os07g0661300',
            'LOC_Os07g48830': 'Os07g0687900',
            'LOC_Os08g04500': 'Os08g0139700',
            'LOC_Os08g30080': 'Os08g0390200',
            'LOC_Os08g30100': 'Os08g0390700',
            'LOC_Os08g31980': 'Os08g0414700',
            'LOC_Os08g35110': 'Os08g0452500',
            'LOC_Os08g44810': 'Os08g0562100',
            'LOC_Os09g18159': 'Os09g0350900',
            'LOC_Os09g26670': 'Os09g0438100',
            'LOC_Os09g27590': 'Os09g0448500',
            'LOC_Os09g27650': 'Os09g0449400',
            'LOC_Os09g28000': 'Os09g0453400',
            'LOC_Os09g28160': 'Os09g0454600',
            'LOC_Os09g28210': 'Os09g0455300',
            'LOC_Os09g34250': 'Os09g0518200',
            'LOC_Os09g38580': 'Os09g0558300',
            'LOC_Os11g27329': 'Os11g0461000'
        }
        
        # Check mapping first for LOC_Os keys (case insensitive check)
        for loc_id, rap_id in loc_to_rap.items():
            if q.upper().startswith(loc_id.upper()):
                return rap_id

        if q.lower().startswith('glyma'):
            parts = q.split('.')
            if len(parts) >= 2:
                q = f"{parts[0]}.{parts[1]}"
        elif q.lower().startswith('loc_os'):
            parts = q.split('.')
            if len(parts) >= 2:
                q = f"{parts[0]}.{parts[1]}"
        else:
            q = q.split('.')[0]
        
        # Strip transcript/protein suffixes like _P01, _T01, _P02, etc.
        q = re.sub(r'_[PT]\d+$', '', q)
        return q

    def fetch_queue(self, queue_df):
        ncbi_key = os.environ.get("NCBI_API_KEY")
        sleep_time = 0.15 if ncbi_key else 0.4
        
        for _, row in tqdm(queue_df.iterrows(), total=len(queue_df), desc="Fetching Phytozome sequences via fallback"):
            acc = row['target_accession']
            organism = row.get('selected_organism', '')
            norm_acc = self._normalize_query(acc)
            print(f"[Phytozome] Resolving and fetching sequence for {acc} (normalized: {norm_acc})...")
            
            # Step 1: Try searching UniProt
            uniprot_acc = self._search_uniprot(norm_acc)
            if uniprot_acc:
                print(f"[Phytozome] Resolved {acc} to UniProt accession {uniprot_acc}")
                url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.fasta"
                response = self._exponential_backoff_request('get', url)
                if isinstance(response, requests.Response) and response.status_code == 200:
                    seq_lines = response.text.strip().split('\n')[1:]
                    clean_seq = '\n'.join(seq_lines)
                    seq_hash = get_sequence_hash(clean_seq)
                    self.save_sequence('UniProt', uniprot_acc, clean_seq, row.get('global_node_id', ''), organism)
                    # Update accession mapping so the final metadata writer finds it
                    self.accession_to_hash[acc] = seq_hash
                    if 'resolved_accessions' in self.stats:
                        self.stats['resolved_accessions'][acc] = ('UniProt', uniprot_acc, seq_hash)
                    time.sleep(sleep_time)
                    continue
            
            # Step 2: Try searching NCBI protein
            ncbi_id = self._search_ncbi(norm_acc, organism)
            if ncbi_id:
                print(f"[Phytozome] Resolved {acc} to NCBI protein ID {ncbi_id}")
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id={ncbi_id}"
                response = self._exponential_backoff_request('get', url)
                if isinstance(response, requests.Response) and response.status_code == 200:
                    text = response.text.strip()
                    if text and text.startswith('>'):
                        seq_lines = text.split('\n')[1:]
                        clean_seq = '\n'.join(seq_lines)
                        seq_hash = get_sequence_hash(clean_seq)
                        self.save_sequence('NCBI', ncbi_id, clean_seq, row.get('global_node_id', ''), organism)
                        # Update accession mapping so the final metadata writer finds it
                        self.accession_to_hash[acc] = seq_hash
                        if 'resolved_accessions' in self.stats:
                            self.stats['resolved_accessions'][acc] = ('NCBI', ncbi_id, seq_hash)
                        time.sleep(sleep_time)
                        continue
            
            self.log_failure('Phytozome', acc, "Could not resolve or fetch sequence from UniProt or NCBI")
            time.sleep(sleep_time)

    def _search_uniprot(self, query):
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json"
        response = self._exponential_backoff_request('get', url)
        if isinstance(response, requests.Response) and response.status_code == 200:
            try:
                data = response.json()
                results = data.get('results', [])
                if results:
                    return results[0].get('primaryAccession')
            except Exception:
                pass
        return None

    def _search_ncbi(self, query, organism):
        ncbi_key = os.environ.get("NCBI_API_KEY")
        key_param = f"&api_key={ncbi_key}" if ncbi_key else ""
        sleep_time = 0.15 if ncbi_key else 0.4
        
        term = f"{query}[All Fields]"
        if organism and str(organism).lower() not in ['nan', 'none']:
            term += f" AND \"{organism}\"[Organism]"
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={term}{key_param}"
        response = self._exponential_backoff_request('get', url)
        time.sleep(sleep_time)
        
        if isinstance(response, requests.Response) and response.status_code == 200:
            import xml.etree.ElementTree as ET
            try:
                root = ET.fromstring(response.text)
                ids = [id_elem.text for id_elem in root.findall('.//IdList/Id')]
                if ids:
                    return ids[0]
            except Exception:
                pass
        # Fallback search without organism filter
        if organism and str(organism).lower() not in ['nan', 'none']:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={query}[All Fields]{key_param}"
            response = self._exponential_backoff_request('get', url)
            time.sleep(sleep_time)
            if isinstance(response, requests.Response) and response.status_code == 200:
                import xml.etree.ElementTree as ET
                try:
                    root = ET.fromstring(response.text)
                    ids = [id_elem.text for id_elem in root.findall('.//IdList/Id')]
                    if ids:
                        return ids[0]
                except Exception:
                    pass
        return None

def migrate_metadata_schema(fasta_path: str, metadata_path: str):
    """If metadata has the old schema, migrate it to the new format."""
    if not os.path.exists(metadata_path):
        return
        
    # Read the first line to verify schema
    with open(metadata_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        try:
            headers = next(reader)
        except StopIteration:
            return
            
    if "uniprot_id" in headers and "target_db" not in headers:
        print("[Migration] Old metadata schema detected. Migrating to new schema...")
        
        # 1. Parse accession sequence hashes from the existing FASTA
        acc_to_hash = {}
        if os.path.exists(fasta_path):
            current_acc = None
            current_seq = []
            with open(fasta_path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        if current_acc and current_seq:
                            acc_to_hash[current_acc] = get_sequence_hash("\n".join(current_seq))
                        # Try parsing custom header first
                        current_acc = parse_accession_from_header(line)
                        if current_acc == line.lstrip('>').strip().split()[0].upper():
                            # Fallback to raw uniprot format
                            current_acc = parse_raw_uniprot_accession(line)
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                if current_acc and current_seq:
                    acc_to_hash[current_acc] = get_sequence_hash("\n".join(current_seq))
                    
        # 2. Convert old rows to the new schema
        migrated_rows = []
        with open(metadata_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                uniprot_id = row.get("uniprot_id", "").strip()
                node_id = row.get("global_node_id", "").strip()
                
                target_db = ""
                target_accession = ""
                retrieval_status = "missing_id"
                seq_hash = ""
                
                if uniprot_id and uniprot_id.lower() != 'nan':
                    target_db = "UniProt"
                    target_accession = uniprot_id.upper()
                    retrieval_status = "success"
                    seq_hash = acc_to_hash.get(target_accession, "")
                
                migrated_row = {
                    "global_node_id": node_id,
                    "target_db": target_db,
                    "target_accession": target_accession,
                    "retrieval_status": retrieval_status,
                    "sequence_hash": seq_hash,
                    "selected_protein_name": row.get("selected_protein_name", ""),
                    "selected_gene_name": row.get("selected_gene_name", ""),
                    "selected_organism": row.get("selected_organism", "")
                }
                # Clean NaNs and 'nan'/'none' strings in migrated rows too
                for k in migrated_row:
                    val = migrated_row[k]
                    if pd.isna(val) or val is None or str(val).strip().lower() in ['nan', 'none']:
                        migrated_row[k] = ""
                migrated_rows.append(migrated_row)
                
        # 3. Rewrite metadata file with the new headers
        new_headers = [
            "global_node_id",
            "target_db",
            "target_accession",
            "retrieval_status",
            "sequence_hash",
            "selected_protein_name",
            "selected_gene_name",
            "selected_organism"
        ]
        with open(metadata_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=new_headers)
            writer.writeheader()
            for r in migrated_rows:
                writer.writerow(r)
                
        print(f"[Migration] Successfully migrated {len(migrated_rows)} metadata rows.")

def load_two_tier_cache(fasta_path: str, metadata_path: str, log_filename: str):
    """Loads processed nodes (Tier 1), accession-to-hash mappings, and unique sequence hashes (Tier 2)."""
    processed_nodes = set()
    accession_to_hash = {}
    sequence_hashes = set()
    
    # Load hashes from FASTA
    if os.path.exists(fasta_path):
        current_seq = []
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    if current_seq:
                        sequence_hashes.add(get_sequence_hash("\n".join(current_seq)))
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_seq:
                sequence_hashes.add(get_sequence_hash("\n".join(current_seq)))
                
    # Load processed node IDs and mapping from Metadata CSV
    if os.path.exists(metadata_path):
        with open(metadata_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                node_id = row.get('global_node_id')
                status = row.get('retrieval_status')
                acc = row.get('target_accession')
                seq_hash = row.get('sequence_hash')
                
                # Only mark as processed if we resolved it or explicitly skipped it (not failed)
                if node_id and status in ['success', 'family_only', 'missing_id']:
                    processed_nodes.add(node_id.strip())
                    
                if acc and status == 'success' and seq_hash:
                    accession_to_hash[acc.strip().upper()] = seq_hash
                    
    print(f"Cache: Loaded {len(processed_nodes)} processed entities, {len(accession_to_hash)} accession-to-hash matches, and {len(sequence_hashes)} unique sequence hashes.")
    return processed_nodes, accession_to_hash, sequence_hashes

def triage_row(row):
    """Applies the waterfall triage logic to determine target database and accession, plus all fallback attempts."""
    attempts = []
    # Priority 1: UniProt (direct or representative fallback)
    for col in ['selected_uniprot_accession', 'representative_uniprot_accession']:
        if pd.notna(row.get(col)) and str(row.get(col)).strip():
            vals = [v.strip().upper() for v in str(row.get(col)).split(';')]
            for val in vals:
                if val not in ['NAN', 'NONE', '']:
                    attempts.append(('UniProt', val))
        
    # Priority 2: NCBI
    for col in ['selected_refseq', 'selected_geneid']:
        if pd.notna(row.get(col)) and str(row.get(col)).strip():
            vals = [v.strip().upper() for v in str(row.get(col)).split(';')]
            for val in vals:
                if val not in ['NAN', 'NONE', '']:
                    attempts.append(('NCBI', val))
            
    # Priority 3: Ensembl
    for col in ['selected_ensembl_plants', 'selected_gramene', 'selected_tair']:
        if pd.notna(row.get(col)) and str(row.get(col)).strip():
            vals = [v.strip().upper() for v in str(row.get(col)).split(';')]
            for val in vals:
                if val not in ['NAN', 'NONE', '']:
                    attempts.append(('Ensembl', val))
            
    # Priority 4: Phytozome
    for col in ['selected_phytozome_gene_id', 'selected_phytozome_base_gene_id']:
        if pd.notna(row.get(col)) and str(row.get(col)).strip():
            vals = [v.strip() for v in str(row.get(col)).split(';')]
            for val in vals:
                if val.upper() not in ['NAN', 'NONE', ''] and val != '':
                    attempts.append(('Phytozome', val))
        
    # Deduplicate attempts while preserving order and prepending LocalCache
    seen = set()
    final_attempts = []
    for db, acc in attempts:
        if ('LocalCache', acc) not in seen:
            seen.add(('LocalCache', acc))
            final_attempts.append(('LocalCache', acc))
        if (db, acc) not in seen:
            seen.add((db, acc))
            final_attempts.append((db, acc))
            
    if final_attempts:
        target_db, target_accession = final_attempts[0]
        return pd.Series([target_db, target_accession, 'queued', final_attempts])
        
    # Guardrail: Family Only
    has_family = False
    for col in ['selected_interpro', 'selected_pfam', 'selected_family_id']:
        if pd.notna(row.get(col)) and str(row.get(col)).strip():
            has_family = True
            
    if has_family:
        return pd.Series([None, None, 'family_only', []])
        
    # Nothing found
    return pd.Series([None, None, 'missing_id', []])

def ingest_input_csvs(input_dir: str):
    """Load, combine, and apply triage to CSV files in the input directory."""
    csv_files = glob.glob(os.path.join(input_dir, "**", "*.csv"), recursive=True)
    allowed_files = {"gene_protein_entity_summary.csv", "manual_normalizations.csv"}
    csv_files = [f for f in csv_files if os.path.basename(f) in allowed_files]
    if not csv_files:
        print(f"Warning: No CSV files found in {input_dir}")
        return pd.DataFrame(), 0
    
    dfs = []
    for filepath in csv_files:
        try:
            df = pd.read_csv(filepath)
            dfs.append(df)
            print(f"Ingested {filepath} with {len(df)} rows.")
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            
    if not dfs:
        return pd.DataFrame(), 0
        
    combined_df = pd.concat(dfs, ignore_index=True)
    initial_count = len(combined_df)
    
    if combined_df.empty:
        return combined_df, 0
    
    # Run the triage waterfall
    combined_df[['target_db', 'target_accession', 'retrieval_status', 'retrieval_attempts']] = combined_df.apply(triage_row, axis=1)
    
    # Priority sorting: queued (0) > family_only (1) > missing_id (2)
    status_priority = {'queued': 0, 'family_only': 1, 'missing_id': 2}
    combined_df['_priority'] = combined_df['retrieval_status'].map(status_priority)
    
    # Sort by priority and drop duplicates on global_node_id keeping the highest priority mapping
    combined_df = combined_df.sort_values(by='_priority')
    combined_df = combined_df.drop_duplicates(subset=['global_node_id'], keep='first')
    combined_df = combined_df.drop(columns=['_priority'], errors='ignore')
    
    # Clean target_accession without turning nulls/nans into literal strings
    if 'target_accession' in combined_df.columns:
        combined_df['target_accession'] = combined_df['target_accession'].apply(
            lambda x: str(x).strip().upper() if pd.notna(x) and x is not None and str(x).strip().lower() not in ['nan', 'none'] else None
        )
    
    print(f"Ingested {initial_count} total rows. Deduplicated to {len(combined_df)} unique entities by global_node_id.")
    return combined_df, initial_count

def run_pipeline(input_dir: str, output_dir: str, log_dir: str, force: bool):
    """Execute the multi-database creator pipeline."""
    print("Starting Multi-DB FASTA Pipeline...")
    
    fasta_output = os.path.join(output_dir, "psfd_sequences.fasta")
    metadata_output = os.path.join(output_dir, "sequence_metadata.csv")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = os.path.join(log_dir, f"run_{run_timestamp}.log")

    # Initialize log file with a header
    with open(log_filename, "w", encoding="utf-8") as f:
        f.write("==================================================\n")
        f.write("PSFD SEQUENCE RETRIEVAL RUN LOG\n")
        f.write("==================================================\n")
        f.write(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input Directory: {input_dir}\n")
        f.write(f"Output Directory: {output_dir}\n")
        f.write(f"Logs Directory: {log_dir}\n")
        f.write(f"Force Re-download: {force}\n\n")
        f.write("DETAILED ERRORS:\n")

    # Initialize stats dictionary
    stats = {
        'total_ingested_rows': 0,
        'total_unique_entities': 0,
        'queued_count': 0,
        'family_only_count': 0,
        'missing_id_count': 0,
        'cache_hits': 0,
        'api_successes': 0,
        'api_failures': 0,
        'deduplicated_sequences': 0,
        'succeeded_list': [],
        'failed_list': [],
        'resolved_accessions': {}
    }

    start_time = time.time()
    try:
        # Migrate metadata schema if old format is present
        if not force:
            migrate_metadata_schema(fasta_output, metadata_output)

        df, initial_count = ingest_input_csvs(input_dir)
        if df.empty:
            print("No valid data to process. Exiting.")
            return
            
        stats['total_ingested_rows'] = initial_count
        stats['total_unique_entities'] = len(df)
        stats['queued_count'] = len(df[df['retrieval_status'] == 'queued'])
        stats['family_only_count'] = len(df[df['retrieval_status'] == 'family_only'])
        stats['missing_id_count'] = len(df[df['retrieval_status'] == 'missing_id'])

        processed_nodes, accession_to_hash, sequence_hashes = set(), {}, set()
        if not force:
            processed_nodes, accession_to_hash, sequence_hashes = load_two_tier_cache(fasta_output, metadata_output, log_filename)

        # Filter out entities that are already successfully processed
        pre_filter_len = len(df)
        df = df[~df['global_node_id'].isin(processed_nodes)]
        stats['cache_hits'] = pre_filter_len - len(df)
        
        if df.empty:
            print("All rows in the input CSV have already been successfully processed. Exiting.")
            return

        # Add missing sequence_hash column if not there
        df['sequence_hash'] = ""

        # Keep a strict set of columns for the metadata CSV
        fieldnames = [
            "global_node_id",
            "target_db",
            "target_accession",
            "retrieval_status",
            "sequence_hash",
            "selected_protein_name",
            "selected_gene_name",
            "selected_organism"
        ]
        
        # Make sure file header exists
        if not os.path.exists(metadata_output):
            with open(metadata_output, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()

        # Process "family_only" and "missing_id" rows (no fetch, omitted from success-only metadata CSV)
        no_fetch_df = df[df['retrieval_status'].isin(['family_only', 'missing_id'])]
        if not no_fetch_df.empty:
            print(f"Skipping writing {len(no_fetch_df)} family_only or missing_id rows to success-only metadata CSV.")

        # Process queued items
        queued_df = df[df['retrieval_status'] == 'queued']
        if queued_df.empty:
            print("No new sequences to fetch. Pipeline finished.")
            return

        # 1. Separate queue into already-cached accession mapping vs. uncached using retrieval_attempts
        cached_rows = []
        uncached_rows = []
        for idx, row in queued_df.iterrows():
            attempts = row['retrieval_attempts']
            cached_info = None
            for db, acc in attempts:
                if acc in accession_to_hash:
                    cached_info = (db, acc)
                    break
            if cached_info:
                row_dict = row.to_dict()
                row_dict['target_db'] = cached_info[0]
                row_dict['target_accession'] = cached_info[1]
                cached_rows.append(row_dict)
            else:
                uncached_rows.append(row.to_dict())

        cached_queued = pd.DataFrame(cached_rows)
        uncached_queued = pd.DataFrame(uncached_rows)

        # 2. Write cached queued items to metadata immediately (no network call needed)
        if not cached_queued.empty:
            print(f"Writing {len(cached_queued)} cached accession rows directly to metadata (Tier 1 hits)...")
            with open(metadata_output, "a", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                for _, row in cached_queued.iterrows():
                    row_dict = row.to_dict()
                    acc = row_dict['target_accession']
                    row_dict['sequence_hash'] = accession_to_hash[acc]
                    row_dict['retrieval_status'] = 'success'
                    write_dict = {field: row_dict.get(field, "") for field in fieldnames}
                    for k in write_dict:
                        val = write_dict[k]
                        if pd.isna(val) or val is None or str(val).strip().lower() in ['nan', 'none']:
                            write_dict[k] = ""
                    writer.writerow(write_dict)

        # 3. For uncached accessions, fetch using dynamic multi-round database fallback
        if not uncached_queued.empty:
            # Search directories for local FASTA files
            search_dirs = [
                input_dir,
                "/home/thomas/Projects/PlantStress-MechanismMap/Notes/phytozome_raw"
            ]
            
            # Initialize Fetchers
            failed_accessions = set()
            fetchers = {
                'LocalCache': LocalCacheFetcher(fasta_output, metadata_output, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats, search_dirs=search_dirs),
                'UniProt': UniProtFetcher(fasta_output, metadata_output, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats),
                'NCBI': NCBIFetcher(fasta_output, metadata_output, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats),
                'Ensembl': EnsemblFetcher(fasta_output, metadata_output, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats),
                'Phytozome': PhytozomeFetcher(fasta_output, metadata_output, log_filename, accession_to_hash, sequence_hashes, failed_accessions, stats)
            }

            active_queue = uncached_queued.to_dict('records')
            for item in active_queue:
                item['active_attempt_index'] = 0

            round_num = 1
            try:
                while active_queue:
                    print(f"\n--- Fallback Fetching Round {round_num} (Queue size: {len(active_queue)}) ---")
                    
                    # Group items by database of their current attempt
                    db_to_accession_items = {}
                    for item in active_queue:
                        attempts = item['retrieval_attempts']
                        idx = item['active_attempt_index']
                        db, acc = attempts[idx]
                        
                        if db not in db_to_accession_items:
                            db_to_accession_items[db] = {}
                        if acc not in db_to_accession_items[db]:
                            db_to_accession_items[db][acc] = []
                        db_to_accession_items[db][acc].append(item)
                    
                    # Process queries for each database in this round
                    for db_name, acc_map in db_to_accession_items.items():
                        rows_to_fetch = []
                        for acc, items in acc_map.items():
                            rep = items[0]
                            rows_to_fetch.append({
                                'target_accession': acc,
                                'target_db': db_name,
                                'global_node_id': rep['global_node_id'],
                                'selected_organism': rep.get('selected_organism', '')
                            })
                        
                        db_queue = pd.DataFrame(rows_to_fetch)
                        print(f"[{db_name}] Round {round_num}: Fetching {len(db_queue)} unique accessions...")
                        fetchers[db_name].fetch_queue(db_queue)
                    
                    # Post-round status scan: check successes and advance failures
                    next_round_queue = []
                    succeeded_items_in_round = []
                    for item in active_queue:
                        attempts = item['retrieval_attempts']
                        idx = item['active_attempt_index']
                        db, acc = attempts[idx]
                        
                        if acc in accession_to_hash or acc in stats['resolved_accessions']:
                            # Query succeeded for this accession!
                            if acc in stats['resolved_accessions']:
                                succ_db, succ_acc, succ_hash = stats['resolved_accessions'][acc]
                            else:
                                succ_db = db
                                succ_acc = acc
                                succ_hash = accession_to_hash[acc]
                            
                            # Store the resolution mapping for metadata writer
                            stats['resolved_accessions'][acc] = (succ_db, succ_acc, succ_hash)
                            
                            # Cache the successful metadata values in the item dictionary
                            item['target_db'] = succ_db
                            item['target_accession'] = succ_acc
                            item['sequence_hash'] = succ_hash
                            item['retrieval_status'] = 'success'
                            succeeded_items_in_round.append(item)
                        else:
                            # Attempt failed. Advance to next attempt
                            item['active_attempt_index'] += 1
                            if item['active_attempt_index'] < len(attempts):
                                next_round_queue.append(item)
                            else:
                                # All attempts exhausted for this entity
                                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                                err_msg = f"All retrieval attempts failed for entity {item['global_node_id']}. Attempts tried: {attempts}"
                                with open(log_filename, "a", encoding="utf-8") as f:
                                    f.write(f"[{timestamp}] ERROR - Entity: {item['global_node_id']} | Reason: {err_msg}\n")
                                print(f"[{item['global_node_id']}] Final failure: {err_msg}")
                                if 'failed_list' in stats:
                                    stats['failed_list'].append((item['global_node_id'], 'All_DBs', 'All attempts failed'))
                    
                    # Write succeeded items from this round progressively
                    if succeeded_items_in_round:
                        print(f"Saving {len(succeeded_items_in_round)} newly resolved metadata rows to CSV...")
                        with open(metadata_output, "a", newline="", encoding="utf-8") as f:
                            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                            for item in succeeded_items_in_round:
                                write_dict = {field: item.get(field, "") for field in fieldnames}
                                for k in write_dict:
                                    val = write_dict[k]
                                    if pd.isna(val) or val is None or str(val).strip().lower() in ['nan', 'none']:
                                        write_dict[k] = ""
                                writer.writerow(write_dict)

                    active_queue = next_round_queue
                    round_num += 1
                    
            except KeyboardInterrupt:
                print("\nPipeline interrupted by user. Gracefully exiting...")
                raise

        print(f"\nPipeline execution finished.")
        
    finally:
        end_time = time.time()
        run_duration = end_time - start_time
        
        # Append summary statistics and success list to the log file
        with open(log_filename, "a", encoding="utf-8") as f:
            f.write("\n==================================================\n")
            f.write("RUN STATISTICS SUMMARY\n")
            f.write("==================================================\n")
            f.write(f"End Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Ingested Rows: {stats['total_ingested_rows']}\n")
            f.write(f"Unique Entities: {stats['total_unique_entities']}\n")
            f.write(f"  - Queued for fetching: {stats['queued_count']}\n")
            f.write(f"  - Skipped (Family only): {stats['family_only_count']}\n")
            f.write(f"  - Skipped (Missing ID): {stats['missing_id_count']}\n")
            f.write(f"Cache Hits (already processed): {stats['cache_hits']}\n")
            f.write(f"API Queries Attempted: {stats['api_successes'] + stats['api_failures']}\n")
            f.write(f"  - Succeeded: {stats['api_successes']}\n")
            f.write(f"  - Failed: {stats['api_failures']}\n")
            f.write(f"  - Duplicate sequences avoided: {stats['deduplicated_sequences']}\n")
            f.write(f"Total Run Time: {run_duration:.2f} seconds ({run_duration/60:.2f} minutes)\n")
            if stats['api_successes'] + stats['api_failures'] > 0:
                avg_time = run_duration / (stats['api_successes'] + stats['api_failures'])
                f.write(f"Average Time per Query: {avg_time:.2f} seconds\n")
            
            f.write("\n==================================================\n")
            f.write(f"SUCCEEDED ACCESSIONS ({stats['api_successes']})\n")
            f.write("==================================================\n")
            if stats['succeeded_list']:
                for acc, db, node_id in stats['succeeded_list']:
                    f.write(f"Entity: {node_id} -> {db}:{acc}\n")
            else:
                f.write("No newly succeeded accessions in this run.\n")
                
        print(f"Log and run statistics written to: {log_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multi-DB Sequence Fetcher Engine")
    project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    parser.add_argument("--input", "-i", default=os.path.join(project_dir, "input"), help="Path to input directory containing CSV files (default: 'input').")
    parser.add_argument("--output", "-o", default=os.path.join(project_dir, "data", "build"), help="Path to output directory (default: 'data/build').")
    parser.add_argument("--logs", "-l", default="logs", help="Path to directory for log files (default: 'logs').")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-download of all sequences, ignoring existing cache.")
    args = parser.parse_args()

    run_pipeline(args.input, args.output, args.logs, args.force)
