#!/usr/bin/env python3
"""Build the static PSFD JSON bundle for the sequence annotation demo.

The builder reads existing PSFD pipeline outputs, removes machine-local absolute
paths, and emits the compact per-paper JSON bundle consumed by this repository's
static browser app. It intentionally writes only `data/` so the independent
frontend assets in this repository are not overwritten.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import shutil
import subprocess
import traceback
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any
from urllib.parse import quote

import scipy.stats as stats
try:
    from tqdm import tqdm
except ImportError:
    print("Warning: tqdm not installed. Install with 'pip install tqdm' for a progress bar.")
    def tqdm(iterable, *args, **kwargs):
        return iterable
    tqdm.write = print

REPO_ROOT = Path(__file__).resolve().parents[1]

import tarfile
import shutil

INPUT_DIR = Path("input")

def extract_tarball_if_needed(tar_path: Path, extract_to: Path) -> Path:
    if not extract_to.exists():
        print(f"Extracting {tar_path.name}...")
        extract_to.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path, "r:gz") as tar:
            
            import os
            
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            safe_extract(tar, extract_to)
            
        # Many tarballs wrap their contents in a single top-level directory (e.g. 720_hypergraph/...)
        # We should check if there's exactly one directory inside, and if so, return that.
        contents = list(extract_to.iterdir())
        if len(contents) == 1 and contents[0].is_dir():
            return contents[0]
    else:
        contents = list(extract_to.iterdir())
        if len(contents) == 1 and contents[0].is_dir():
            return contents[0]
    return extract_to


def resolve_input_dir(*names: str) -> Path:
    found = []
    
    for name in names:
        if not name:
            continue
            
        # Try directory match
        dir_path = INPUT_DIR / name
        if dir_path.is_dir():
            contents = list(dir_path.iterdir())
            if len(contents) == 1 and contents[0].name.endswith(".tar.gz"):
                nested_tar = contents[0]
                print(f"Extracting nested tarball {nested_tar.name} in {dir_path}...")
                import subprocess
                subprocess.run(["tar", "-xzf", nested_tar.name], cwd=dir_path, check=True)
            found.append((name, dir_path, False))
            
        # Try tarball match
        tar_path = INPUT_DIR / f"{name}.tar.gz"
        if tar_path.is_file():
            found.append((name, tar_path, True))
            
    if len(found) > 1:
        found_names = [f[0] for f in found]
        raise ValueError(f"Ambiguous input files found. Multiple valid candidates exist: {found_names}. Please remove duplicates to ensure the correct data is loaded.")
        
    if not found:
        raise FileNotFoundError(f"Could not find input directory or tarball for any of {names}. Expected one in {INPUT_DIR.resolve()}")
        
    name, path, is_tarball = found[0]
    
    if is_tarball:
        tmp_dir = INPUT_DIR / "tmp" / name
        return extract_tarball_if_needed(path, tmp_dir)
        
    # Check for single nested directory wrapper (ignoring any leftover .tar.gz files)
    dirs = [f for f in path.iterdir() if f.is_dir()]
    files = [f for f in path.iterdir() if f.is_file()]
    if len(dirs) == 1 and all(f.name.endswith(".tar.gz") for f in files):
        return dirs[0]
        
    return path

HYPERGRAPH_DIR = Path()
SENTENCE_DIR = Path()
TRIPLES_EVALUATION_DIR = Path()
NORMALIZATION_DIR = Path()
COMPOUND_CLASSIFICATION_DIR = Path()
GENE_PROTEIN_NORMALIZATION_DIR = Path()
CONTEXT_PROPAGATION_DIR = Path()

def default_resource_path(env_name: str, portable_path: Path) -> Path:
    configured = os.environ.get(env_name)
    if configured:
        return Path(configured).expanduser()
    return portable_path


# These route targets are the Phytozome browser/proteome names used by the
# current website for the short codes present in the Step 10 output.
PHYTOZOME_ROUTE_TARGETS = {
    "Osat": {"browser_name": "Osativa_v7_0", "proteome_id": "323", "label": "Oryza sativa v7.0"},
    "Slyc": {"browser_name": "Slycopersicum_ITAG4_0", "proteome_id": "691", "label": "Solanum lycopersicum ITAG4.0"},
}
COMPOUND_CLASSIFICATION_FIELDS = [
    "pmcid",
    "entity_instance_id",
    "entity_surface_id",
    "global_node_id",
    "canonical_form",
    "aliases",
    "normalization_decision",
    "normalization_status",
    "selected_ontology",
    "selected_ontology_id",
    "selected_label",
    "selected_rank",
    "compound_status",
    "classification_status",
    "chebi_id",
    "chebi_name",
    "chebi_formula",
    "chebi_inchikey",
    "pubchem_cid",
    "pubchem_link_type",
    "structure_source",
    "structure_inchikey",
    "structure_smiles",
    "classyfire_cache_hit",
    "classyfire_kingdom",
    "classyfire_superclass",
    "classyfire_class",
    "classyfire_subclass",
    "classyfire_direct_parent",
    "npclassifier_applicable",
    "npclassifier_cache_hit",
    "np_pathway",
    "np_superclass",
    "np_class",
    "np_is_glycoside",
    "npclassifier_error",
]
GENE_PROTEIN_NORMALIZATION_FIELDS = [
    "pmcid",
    "entity_instance_id",
    "entity_surface_id",
    "global_node_id",
    "canonical_form",
    "aliases",
    "decomposed_index",
    "gene_query",
    "gene_query_type",
    "lookup_query",
    "lookup_strategy",
    "decision",
    "status",
    "normalization_scope",
    "normalization_confidence",
    "selected_uniprot_accession",
    "selected_uniprot_entry",
    "selected_gene_name",
    "selected_protein_name",
    "selected_organism",
    "selected_taxon_id",
    "selected_taxon_rank",
    "reviewed_status",
    "source_database",
    "selected_refseq",
    "selected_geneid",
    "selected_ensembl_plants",
    "selected_gramene",
    "selected_tair",
    "selected_interpro",
    "selected_pfam",
    "selected_phytozome_code",
    "selected_phytozome_gene_id",
    "selected_phytozome_base_gene_id",
    "selected_phytozome_sequence_length",
    "selected_phytozome_source_file",
    "selected_family_id",
    "selected_family_database",
    "selected_family_type",
    "selected_family_name",
    "selected_family_alias",
    "selected_family_alias_type",
    "selected_family_linked_interpro_id",
    "selected_family_linked_pfam_ids",
    "representative_uniprot_accession",
    "representative_uniprot_entry",
    "representative_organism",
    "representative_source_database",
    "representative_reviewed_status",
    "representative_basis",
    "match_type",
    "taxon_context_source",
    "taxon_context_ids",
    "taxon_context_labels",
    "candidate_count",
]

GLOBAL_NORMALIZATIONS: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)

def preload_global_normalizations() -> None:
    csv_path = NORMALIZATION_DIR / "normalized_entities.csv"
    if source_exists(csv_path):
        resolved_path = source_path(csv_path)
        print(f"Preloading global normalizations from {resolved_path}...")
        rows = load_csv_rows(csv_path)
        for row in rows:
            pmcid = row.get("pmcid", "")
            node_id = row.get("node_id", "")
            if not pmcid or not node_id:
                continue
            GLOBAL_NORMALIZATIONS[pmcid][node_id] = {
                "selected_ontology": row.get("selected_ontology", ""),
                "selected_ontology_id": row.get("selected_ontology_id", ""),
                "selected_label": row.get("selected_label", ""),
                "canonical_form": row.get("canonical_form", ""),
                "status": row.get("status", ""),
            }

    manual_csv_path = INPUT_DIR / "manual_normalizations.csv"
    if source_exists(manual_csv_path):
        resolved_manual = source_path(manual_csv_path)
        print(f"Preloading manual normalizations from {resolved_manual}...")
        manual_rows = load_csv_rows(manual_csv_path)
        manual_count = 0
        for row in manual_rows:
            pmcid = row.get("pmcid", "").strip()
            node_id = row.get("node_id", "").strip()
            global_node_id = row.get("global_node_id", "").strip()
            if not pmcid or not node_id:
                if global_node_id and ":" in global_node_id:
                    pmcid, node_id = global_node_id.split(":", 1)
            if not pmcid or not node_id:
                continue

            existing = GLOBAL_NORMALIZATIONS[pmcid].get(node_id, {})
            existing.update({
                "selected_ontology": row.get("selected_ontology", existing.get("selected_ontology", "")),
                "selected_ontology_id": row.get("selected_ontology_id", existing.get("selected_ontology_id", "")),
                "selected_label": row.get("selected_label", existing.get("selected_label", "")),
                "canonical_form": row.get("canonical_form", existing.get("canonical_form", "")),
                "status": row.get("status", existing.get("status", "manual_override")),
            })
            GLOBAL_NORMALIZATIONS[pmcid][node_id] = existing
            manual_count += 1
        print(f"Injected {manual_count} manual normalization overrides.")



def source_path(path: Path) -> Path:
    """Return an existing source path, accepting storage-efficient .zst output."""
    if path.exists():
        return path
    compressed = Path(f"{path}.zst")
    if compressed.exists():
        return compressed
    return path


def source_exists(path: Path) -> bool:
    return source_path(path).exists()


def read_source_text(path: Path) -> str:
    resolved = source_path(path)
    if resolved.suffix == ".zst":
        return subprocess.check_output(["zstd", "-dc", str(resolved)], text=True)
    return resolved.read_text(encoding="utf-8")


def load_json(path: Path) -> Any:
    return json.loads(read_source_text(path))


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(read_source_text(path))))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"))


def as_list(value: Any) -> list[Any]:
    if value is None or value == "":
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, str):
        return [part.strip() for part in value.split("|") if part.strip()]
    return [value]


def truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def compact_dict(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in payload.items()
        if value not in (None, "", [], {})
    }


def relpath(path: Path) -> str:
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj.name)
    else:
        return str(path)


def unique_dicts_by_id(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[str, dict[str, Any]] = {}
    for item in items:
        item_id = str(item.get("id") or item.get("ontology_id") or "").strip()
        if item_id and item_id not in unique:
            unique[item_id] = item
    return list(unique.values())


def phytozome_route_target(code: str) -> dict[str, str]:
    return PHYTOZOME_ROUTE_TARGETS.get(str(code or "").strip(), {})


def phytozome_report_type(gene_id: str, base_gene_id: str = "") -> str:
    gene_id = str(gene_id or "").strip()
    base_gene_id = str(base_gene_id or "").strip()
    return "transcript" if gene_id and base_gene_id and gene_id != base_gene_id else "gene"


def phytozome_gene_url(code: str, gene_id: str, base_gene_id: str = "") -> str:
    if not code or not gene_id:
        return "https://phytozome-next.jgi.doe.gov/"
    target = phytozome_route_target(code).get("browser_name") or str(code)
    report_type = phytozome_report_type(gene_id, base_gene_id)
    report_id = str(gene_id if report_type == "transcript" else (base_gene_id or gene_id))
    return (
        f"https://phytozome-next.jgi.doe.gov/report/{report_type}/"
        f"{quote(str(target), safe='')}/{quote(report_id, safe='')}"
    )


def phytozome_search_url(gene_id: str, base_gene_id: str = "") -> str:
    query = str(gene_id or base_gene_id or "").strip()
    if not query:
        return "https://phytozome-next.jgi.doe.gov/"
    return f"https://phytozome-next.jgi.doe.gov/search?query={quote(query, safe='')}"


def family_ontology_id(database: str, identifier: str) -> str:
    database = str(database or "").strip()
    identifier = str(identifier or "").strip()
    if not identifier:
        return ""
    if ":" in identifier:
        return identifier
    if database.lower() == "interpro" or identifier.upper().startswith("IPR"):
        return f"InterPro:{identifier}"
    if database.lower() == "pfam" or identifier.upper().startswith("PF"):
        return f"Pfam:{identifier}"
    return f"{database}:{identifier}" if database else identifier


def family_resource_url(ontology_id: str) -> str:
    if ontology_id.startswith("InterPro:"):
        return f"https://www.ebi.ac.uk/interpro/entry/InterPro/{quote(ontology_id.split(':', 1)[1], safe='')}/"
    if ontology_id.startswith("Pfam:"):
        return f"https://www.ebi.ac.uk/interpro/entry/pfam/{quote(ontology_id.split(':', 1)[1], safe='')}/"
    return f"https://identifiers.org/{ontology_id}" if ontology_id else ""


def split_identifier_values(value: Any) -> list[str]:
    """Split compact Step 10 identifier lists while preserving gene-like tokens."""
    if not value:
        return []
    if isinstance(value, list):
        values = value
    else:
        values = re.split(r"[|;,]", str(value))
    return unique_preserve_order(
        part.strip()
        for part in values
        if str(part).strip() and str(part).strip() not in {"-", "NA", "None", "null"}
    )


def unique_preserve_order(values: Any) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        text = str(value or "").strip()
        key = text.lower()
        if not text or key in seen:
            continue
        seen.add(key)
        unique.append(text)
    return unique


def gene_database_ontology_id(database: str, identifier: str) -> str:
    database = str(database or "").strip()
    identifier = str(identifier or "").strip()
    if not identifier:
        return ""
    if ":" in identifier and identifier.split(":", 1)[0].lower() in {
        "refseq",
        "ncbigene",
        "ensemblplants",
        "gramene",
        "tair",
        "interpro",
        "pfam",
    }:
        return identifier
    prefix = {
        "RefSeq": "RefSeq",
        "NCBIGene": "NCBIGene",
        "EnsemblPlants": "EnsemblPlants",
        "Gramene": "Gramene",
        "TAIR": "TAIR",
        "InterPro": "InterPro",
        "Pfam": "Pfam",
    }.get(database, database)
    return f"{prefix}:{identifier}" if prefix else identifier


def gene_database_resource_url(database: str, identifier: str) -> str:
    database = str(database or "").strip()
    identifier = str(identifier or "").strip()
    if not identifier:
        return ""
    if database == "RefSeq":
        return f"https://www.ncbi.nlm.nih.gov/protein/{quote(identifier, safe='')}"
    if database == "NCBIGene":
        return f"https://www.ncbi.nlm.nih.gov/gene/{quote(identifier, safe='')}"
    if database == "EnsemblPlants":
        return f"https://plants.ensembl.org/Multi/Search/Results?q={quote(identifier, safe='')}"
    if database == "Gramene":
        return f"https://ensembl.gramene.org/Multi/Search/Results?q={quote(identifier, safe='')}"
    if database == "TAIR":
        return f"https://www.arabidopsis.org/servlets/Search?action=new_search&type=general&search_action=detail&method=1&show_obsolete=F&name={quote(identifier, safe='')}"
    if database in {"InterPro", "Pfam"}:
        return family_resource_url(gene_database_ontology_id(database, identifier))
    ontology_id = gene_database_ontology_id(database, identifier)
    return f"https://identifiers.org/{ontology_id}" if ontology_id else ""


def gene_database_ids(raw: dict[str, Any]) -> list[dict[str, str]]:
    database_specs = [
        ("selected_refseq", "RefSeq", r"^[NX][PMR]_\d+(?:\.\d+)?$"),
        ("selected_geneid", "NCBIGene", r"^\d+$"),
        ("selected_ensembl_plants", "EnsemblPlants", r"^[A-Za-z0-9_.-]+$"),
        ("selected_gramene", "Gramene", r"^[A-Za-z0-9_.-]+$"),
        ("selected_tair", "TAIR", r"^(?:AT[1-5CM]G\d+(?:\.\d+)?|[A-Za-z][A-Za-z0-9_.-]*)$"),
        ("selected_interpro", "InterPro", r"^IPR\d+$"),
        ("selected_pfam", "Pfam", r"^PF\d+$"),
    ]
    ids: dict[str, dict[str, str]] = {}
    for field, database, pattern in database_specs:
        for token in split_identifier_values(raw.get(field, "")):
            if not re.match(pattern, token, flags=re.IGNORECASE):
                continue
            ontology_id = gene_database_ontology_id(database, token)
            ids.setdefault(
                ontology_id,
                {
                    "ontology_id": ontology_id,
                    "database": database,
                    "identifier": token,
                    "source_field": field,
                    "resource_url": gene_database_resource_url(database, token),
                },
            )
    return list(ids.values())



def round_float(value: Any, digits: int = 4) -> float | str:
    if value in (None, ""):
        return ""
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return str(value)


def article_metadata(pmcid: str) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    path = SENTENCE_DIR / f"{pmcid}.json"
    if not source_exists(path):
        return {"pmcid": pmcid}, {}

    data = load_json(path)
    document = data[0]["documents"][0]
    passages = document.get("passages", [])
    first = passages[0] if passages else {}
    infons = first.get("infons", {})

    authors: list[str] = []
    for key, value in sorted(infons.items()):
        if not key.startswith("name_"):
            continue
        pieces = {}
        for part in str(value).split(";"):
            if ":" in part:
                k, v = part.split(":", 1)
                pieces[k] = v
        surname = pieces.get("surname", "").strip()
        given = pieces.get("given-names", "").strip()
        authors.append(", ".join([p for p in [surname, given] if p]))

    sentences: dict[str, dict[str, Any]] = {}
    for passage in passages:
        passage_infons = passage.get("infons", {})
        section = (
            passage_infons.get("section_type")
            or passage_infons.get("type")
            or passage_infons.get("name")
            or ""
        )
        for sentence in passage.get("sentences", []):
            sinfons = sentence.get("infons", {})
            sentence_id = sinfons.get("sentence_id", "")
            if not sentence_id:
                continue
            sentences[sentence_id] = {
                "id": sentence_id,
                "text": sentence.get("text", ""),
                "section": section,
                "passage_index": int(str(sinfons.get("passage_index") or 0) or 0),
                "sentence_index": int(str(sinfons.get("sentence_index") or 0) or 0),
                "offset": sentence.get("offset", ""),
            }

    article = {
        "pmcid": pmcid,
        "title": first.get("text", "").replace("\n", " ").strip(),
        "doi": infons.get("article-id_doi") or infons.get("pub-id_doi") or "",
        "pmid": infons.get("article-id_pmid", ""),
        "journal": infons.get("source", ""),
        "year": infons.get("year", ""),
        "volume": infons.get("volume", ""),
        "issue": infons.get("issue", ""),
        "pages": "-".join([x for x in [infons.get("fpage", ""), infons.get("lpage", "")] if x]),
        "license": infons.get("license", document.get("infons", {}).get("license", "")),
        "authors": authors,
    }
    return article, sentences


CONTEXT_SENTENCE_RE = re.compile(r"(PMC\d+\.p\d+\.s\d+(?:-\d+)?)\s*:")


def evidence_context_sentences(
    raw: dict[str, Any],
    sentence_map: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    text = str(raw.get("evidence_context_text") or "").strip()
    if not text:
        return []
    matches = list(CONTEXT_SENTENCE_RE.finditer(text))
    if not matches:
        return [{"role": "neighboring context", "text": text}]

    evidence_sentence_ids = set(as_list(raw.get("evidence_sentence_ids")))
    records: list[dict[str, Any]] = []
    for index, match in enumerate(matches):
        sentence_id = match.group(1)
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        extracted_text = re.sub(r"\s+", " ", text[start:end]).strip()
        sentence = sentence_map.get(sentence_id, {})
        passage_id = sentence_id.rsplit(".", 1)[0] if "." in sentence_id else raw.get("current_passage_id", "")
        records.append(
            compact_dict(
                {
                    "id": sentence_id,
                    "text": sentence.get("text") or extracted_text,
                    "section": sentence.get("section", ""),
                    "passage_id": passage_id,
                    "passage_index": sentence.get("passage_index", raw.get("current_passage_index", "")),
                    "sentence_index": sentence.get("sentence_index", ""),
                    "offset": sentence.get("offset", ""),
                    "role": "core evidence" if sentence_id in evidence_sentence_ids else "neighboring context",
                }
            )
        )
    return records


def load_relation_evaluations(pmcid: str) -> dict[str, dict[str, Any]]:
    """Load Step 4 relation evaluation metadata by source relation ID."""
    path = TRIPLES_EVALUATION_DIR / f"{pmcid}.json"
    if not source_exists(path):
        return {}
    try:
        payload = load_json(path)
    except json.JSONDecodeError:
        return {}
    evaluations: dict[str, dict[str, Any]] = {}
    for row in as_list(payload.get("evaluations")):
        if not isinstance(row, dict):
            continue
        relation_id = str(row.get("relation_id") or "").strip()
        if not relation_id:
            continue
        evaluations[relation_id] = compact_dict(
            {
                "assertion_modifier": row.get("assertion_modifier", ""),
                "verdict": row.get("verdict", ""),
                "status": row.get("status", ""),
                "confidence": round_float(row.get("confidence")),
                "passage_id": row.get("passage_id", ""),
                "passage_index": row.get("passage_index", ""),
                "evidence_sentence_ids": as_list(row.get("evidence_sentence_ids")),
            }
        )
    return evaluations


def external_links(entity: dict[str, Any]) -> list[dict[str, str]]:
    ontology = str(entity.get("selected_ontology") or entity.get("ontology") or "")
    ontology_id = str(entity.get("selected_ontology_id") or entity.get("ontology_id") or "")
    label = str(entity.get("selected_label") or entity.get("canonical_form") or "")
    links: list[dict[str, str]] = []

    def add_link(link_label: str, url: str) -> None:
        if url and not any(item["url"] == url for item in links):
            links.append({"label": link_label, "url": url})

    concepts = list(entity.get("selected_concepts") or [])
    if ontology_id and not any(concept.get("id") == ontology_id for concept in concepts):
        concepts.insert(0, {"id": ontology_id, "ontology": ontology, "label": label})

    for concept in concepts:
        concept_id = str(concept.get("id") or "").strip()
        concept_ontology = str(concept.get("ontology") or (concept_id.split(":", 1)[0] if ":" in concept_id else ontology))
        concept_label = str(concept.get("label") or label)
        if not concept_id:
            continue
        add_link(f"Identifiers.org {concept_id}", f"https://identifiers.org/{concept_id}")
        local_id = concept_id.split(":", 1)[1] if ":" in concept_id else concept_id
        if concept_ontology == "PubChem" and local_id:
            add_link("PubChem compound", f"https://pubchem.ncbi.nlm.nih.gov/compound/{local_id}")
            if concept_label:
                add_link("PubChem name search", f"https://pubchem.ncbi.nlm.nih.gov/#query={quote(concept_label, safe='')}")
        elif concept_ontology == "NCBITaxon" and local_id:
            add_link("NCBI Taxonomy", f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={local_id}")
        elif concept_ontology == "CHEBI" and concept_id:
            add_link("ChEBI", f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={concept_id}")
        else:
            add_link(f"OLS search {concept_id}", f"https://www.ebi.ac.uk/ols4/search?q={quote(concept_id, safe='')}")

    compound = entity.get("compound_classification") or {}
    pubchem_cid = str((compound.get("pubchem") or {}).get("cid") or "")
    if pubchem_cid and not any(link["url"].endswith(f"/{pubchem_cid}") for link in links):
        links.append({"label": "PubChem compound", "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}"})
    gene_profile = entity.get("gene_protein_normalization") or {}
    for accession in gene_profile.get("fasta_accessions", [])[:4]:
        acc = accession.get("accession", "")
        if not acc:
            continue
        if not any(link["url"].endswith(f"/{acc}") for link in links):
            links.append({"label": f"UniProt {acc}", "url": f"https://www.uniprot.org/uniprotkb/{acc}"})
        fasta_url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        if not any(link["url"] == fasta_url for link in links):
            links.append({"label": f"FASTA {acc}", "url": fasta_url})
    for taxon_id in as_list((gene_profile.get("best") or {}).get("selected_taxon_id"))[:2]:
        if taxon_id and not any(link["url"].endswith(f"id={taxon_id}") for link in links):
            links.append({"label": f"NCBI Taxonomy {taxon_id}", "url": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxon_id}"})
    for item in gene_profile.get("phytozome_ids", [])[:4]:
        gene_id = item.get("gene_id", "")
        report_url = item.get("gene_report_url") or phytozome_gene_url(item.get("code", ""), gene_id, item.get("base_gene_id", ""))
        search_url = item.get("search_url") or phytozome_search_url(gene_id, item.get("base_gene_id", ""))
        if gene_id and report_url and not any(link["url"] == report_url for link in links):
            links.append({"label": f"Phytozome report {gene_id}", "url": report_url})
        if gene_id and search_url and not any(link["url"] == search_url for link in links):
            links.append({"label": f"Search Phytozome {gene_id}", "url": search_url})
    for item in gene_profile.get("family_ids", [])[:6]:
        ontology_id = item.get("ontology_id", "")
        url = item.get("resource_url") or family_resource_url(ontology_id)
        if ontology_id and url and not any(link["url"] == url for link in links):
            links.append({"label": ontology_id, "url": url})
    for item in gene_profile.get("database_ids", [])[:8]:
        ontology_id = item.get("ontology_id", "")
        url = item.get("resource_url", "")
        if ontology_id and url and not any(link["url"] == url for link in links):
            links.append({"label": ontology_id, "url": url})
    return links


def slim_compound_classification(raw: dict[str, Any]) -> dict[str, Any]:
    """Keep only browser-facing compound metadata from Step 930."""
    return compact_dict(
        {
            "source": relpath(COMPOUND_CLASSIFICATION_DIR),
            "entity_instance_id": raw.get("entity_instance_id", ""),
            "entity_surface_id": raw.get("entity_surface_id", ""),
            "global_node_id": raw.get("global_node_id", ""),
            "canonical_form": raw.get("canonical_form", ""),
            "aliases": as_list(raw.get("aliases")),
            "normalization": compact_dict(
                {
                    "decision": raw.get("normalization_decision", ""),
                    "status": raw.get("normalization_status", ""),
                    "selected_ontology": raw.get("selected_ontology", ""),
                    "selected_ontology_id": raw.get("selected_ontology_id", ""),
                    "selected_label": raw.get("selected_label", ""),
                    "selected_rank": raw.get("selected_rank", ""),
                }
            ),
            "compound_status": raw.get("compound_status", ""),
            "classification_status": raw.get("classification_status", ""),
            "chebi": compact_dict(
                {
                    "id": raw.get("chebi_id", ""),
                    "name": raw.get("chebi_name", ""),
                    "formula": raw.get("chebi_formula", ""),
                    "inchikey": raw.get("chebi_inchikey", ""),
                }
            ),
            "pubchem": compact_dict(
                {
                    "cid": raw.get("pubchem_cid", ""),
                    "link_type": raw.get("pubchem_link_type", ""),
                }
            ),
            "structure": compact_dict(
                {
                    "source": raw.get("structure_source", ""),
                    "inchikey": raw.get("structure_inchikey", ""),
                    "smiles": raw.get("structure_smiles", ""),
                }
            ),
            "classyfire": compact_dict(
                {
                    "cache_hit": truthy(raw.get("classyfire_cache_hit", "")),
                    "kingdom": raw.get("classyfire_kingdom", ""),
                    "superclass": raw.get("classyfire_superclass", ""),
                    "class": raw.get("classyfire_class", ""),
                    "subclass": raw.get("classyfire_subclass", ""),
                    "direct_parent": raw.get("classyfire_direct_parent", ""),
                }
            ),
            "npclassifier": compact_dict(
                {
                    "applicable": truthy(raw.get("npclassifier_applicable", "")),
                    "cache_hit": truthy(raw.get("npclassifier_cache_hit", "")),
                    "pathway": raw.get("np_pathway", ""),
                    "superclass": raw.get("np_superclass", ""),
                    "class": raw.get("np_class", ""),
                    "is_glycoside": truthy(raw.get("np_is_glycoside", "")),
                    "error": raw.get("npclassifier_error", ""),
                }
            ),
            "raw_fields": {field: raw.get(field, "") for field in COMPOUND_CLASSIFICATION_FIELDS},
        }
    )


def load_compound_classifications(pmcid: str) -> dict[str, dict[str, Any]]:
    path = COMPOUND_CLASSIFICATION_DIR / f"{pmcid}.json"
    rows: list[dict[str, Any]] = []
    if source_exists(path):
        rows = [row for row in load_json(path).get("rows", []) if isinstance(row, dict)]
    else:
        csv_path = COMPOUND_CLASSIFICATION_DIR / f"{pmcid}.compound_classifications.csv"
        if source_exists(csv_path):
            rows = load_csv_rows(csv_path)

    by_key: dict[str, dict[str, Any]] = {}
    for raw in rows:
        slim = slim_compound_classification(raw)
        for key in (
            raw.get("entity_instance_id", ""),
            raw.get("global_node_id", ""),
            raw.get("selected_ontology_id", ""),
            raw.get("chebi_id", ""),
            f"PubChem:{raw.get('pubchem_cid', '')}" if raw.get("pubchem_cid") else "",
        ):
            if key:
                by_key.setdefault(key, slim)
    return by_key


def slim_gene_protein_row(raw: dict[str, Any]) -> dict[str, Any]:
    return compact_dict(
        {
            "entity_instance_id": raw.get("entity_instance_id", ""),
            "entity_surface_id": raw.get("entity_surface_id", ""),
            "global_node_id": raw.get("global_node_id", ""),
            "canonical_form": raw.get("canonical_form", ""),
            "aliases": as_list(raw.get("aliases")),
            "decomposed_index": raw.get("decomposed_index", ""),
            "gene_query": raw.get("gene_query", ""),
            "gene_query_type": raw.get("gene_query_type", ""),
            "lookup_query": raw.get("lookup_query", ""),
            "lookup_strategy": raw.get("lookup_strategy", ""),
            "decision": raw.get("decision", ""),
            "status": raw.get("status", ""),
            "normalization_scope": raw.get("normalization_scope", ""),
            "selected": compact_dict(
                {
                    "uniprot_accession": raw.get("selected_uniprot_accession", ""),
                    "uniprot_entry": raw.get("selected_uniprot_entry", ""),
                    "gene_name": raw.get("selected_gene_name", ""),
                    "protein_name": raw.get("selected_protein_name", ""),
                    "organism": raw.get("selected_organism", ""),
                    "taxon_id": raw.get("selected_taxon_id", ""),
                    "taxon_rank": raw.get("selected_taxon_rank", ""),
                    "reviewed_status": raw.get("reviewed_status", ""),
                    "source_database": raw.get("source_database", ""),
                }
            ),
            "phytozome": compact_dict(
                {
                    "code": raw.get("selected_phytozome_code", ""),
                    "gene_id": raw.get("selected_phytozome_gene_id", ""),
                    "base_gene_id": raw.get("selected_phytozome_base_gene_id", ""),
                    "sequence_length": raw.get("selected_phytozome_sequence_length", ""),
                    "source_file": relpath(Path(raw.get("selected_phytozome_source_file", ""))),
                    "source_database": raw.get("source_database", ""),
                }
            ),
            "family": compact_dict(
                {
                    "id": raw.get("selected_family_id", ""),
                    "ontology_id": family_ontology_id(raw.get("selected_family_database", ""), raw.get("selected_family_id", "")),
                    "database": raw.get("selected_family_database", ""),
                    "type": raw.get("selected_family_type", ""),
                    "name": raw.get("selected_family_name", ""),
                    "alias": raw.get("selected_family_alias", ""),
                    "alias_type": raw.get("selected_family_alias_type", ""),
                    "linked_interpro_id": raw.get("selected_family_linked_interpro_id", ""),
                    "linked_interpro_ontology_id": family_ontology_id("InterPro", raw.get("selected_family_linked_interpro_id", "")),
                    "linked_pfam_ids": as_list(raw.get("selected_family_linked_pfam_ids", "")),
                    "linked_pfam_ontology_ids": [
                        family_ontology_id("Pfam", item)
                        for item in as_list(raw.get("selected_family_linked_pfam_ids", ""))
                    ],
                    "source_database": raw.get("source_database", ""),
                }
            ),
            "representative": compact_dict(
                {
                    "uniprot_accession": raw.get("representative_uniprot_accession", ""),
                    "uniprot_entry": raw.get("representative_uniprot_entry", ""),
                    "organism": raw.get("representative_organism", ""),
                    "source_database": raw.get("representative_source_database", ""),
                    "reviewed_status": raw.get("representative_reviewed_status", ""),
                    "basis": raw.get("representative_basis", ""),
                }
            ),
            "database_ids": gene_database_ids(raw),
            "match_type": raw.get("match_type", ""),
            "taxon_context": compact_dict(
                {
                    "source": raw.get("taxon_context_source", ""),
                    "ids": as_list(raw.get("taxon_context_ids")),
                    "labels": as_list(raw.get("taxon_context_labels")),
                }
            ),
            "candidate_count": raw.get("candidate_count", ""),
            "ambiguity_reason": raw.get("ambiguity_reason", ""),
            "notes": raw.get("notes", ""),
            "raw_fields": {field: raw.get(field, "") for field in GENE_PROTEIN_NORMALIZATION_FIELDS},
        }
    )


def gene_protein_rank(row: dict[str, Any]) -> tuple[int, int]:
    decision = row.get("decision", "")
    has_selected = bool((row.get("selected") or {}).get("uniprot_accession"))
    has_phytozome = bool((row.get("phytozome") or {}).get("gene_id"))
    has_family = bool((row.get("family") or {}).get("id"))
    has_representative = bool((row.get("representative") or {}).get("uniprot_accession"))
    rank = {
        "match": 0,
        "sequence_match": 0,
        "representative_match": 1,
        "ambiguous": 2,
        "family_domain_match": 2,
        "family_or_class": 3,
        "protein_complex": 4,
        "generic_gene_set": 5,
        "no_match": 8,
    }.get(decision, 7)
    if has_selected:
        rank -= 1
    elif has_phytozome:
        rank -= 1
    elif has_family:
        rank -= 1
    elif has_representative:
        rank -= 0
    return (rank, int(str(row.get("decomposed_index") or 999) or 999))


def build_gene_protein_profile(rows: list[dict[str, Any]]) -> dict[str, Any]:
    sorted_rows = sorted(rows, key=gene_protein_rank)
    best = sorted_rows[0] if sorted_rows else {}
    accessions: dict[str, dict[str, str]] = {}
    phytozome_ids: dict[str, dict[str, str]] = {}
    family_ids: dict[str, dict[str, Any]] = {}
    database_ids: dict[str, dict[str, str]] = {}
    for row in sorted_rows:
        selected = row.get("selected") or {}
        phytozome = row.get("phytozome") or {}
        family = row.get("family") or {}
        representative = row.get("representative") or {}
        for item in row.get("database_ids", []):
            ontology_id = item.get("ontology_id", "")
            if ontology_id:
                database_ids.setdefault(ontology_id, item)
        if selected.get("uniprot_accession"):
            acc = selected["uniprot_accession"]
            accessions.setdefault(
                acc,
                {
                    "accession": acc,
                    "kind": "selected",
                    "entry": selected.get("uniprot_entry", ""),
                    "gene_name": selected.get("gene_name", ""),
                    "protein_name": selected.get("protein_name", ""),
                    "organism": selected.get("organism", ""),
                    "taxon_id": selected.get("taxon_id", ""),
                    "reviewed_status": selected.get("reviewed_status", ""),
                    "source_database": selected.get("source_database", ""),
                },
            )
        if phytozome.get("gene_id"):
            gene_id = phytozome["gene_id"]
            route_target = phytozome_route_target(phytozome.get("code", ""))
            phytozome_ids.setdefault(
                gene_id,
                {
                    "ontology_id": f"Phytozome:{gene_id}",
                    "gene_id": gene_id,
                    "base_ontology_id": f"PhytozomeBase:{phytozome.get('base_gene_id')}" if phytozome.get("base_gene_id") else "",
                    "base_gene_id": phytozome.get("base_gene_id", ""),
                    "code": phytozome.get("code", ""),
                    "browser_name": route_target.get("browser_name", ""),
                    "proteome_id": route_target.get("proteome_id", ""),
                    "route_label": route_target.get("label", ""),
                    "report_type": phytozome_report_type(gene_id, phytozome.get("base_gene_id", "")),
                    "sequence_length": phytozome.get("sequence_length", ""),
                    "source_file": phytozome.get("source_file", ""),
                    "source_database": phytozome.get("source_database", ""),
                    "search_url": phytozome_search_url(gene_id, phytozome.get("base_gene_id", "")),
                    "gene_report_url": phytozome_gene_url(phytozome.get("code", ""), gene_id, phytozome.get("base_gene_id", "")),
                },
            )
        if family.get("ontology_id"):
            ontology_id = family["ontology_id"]
            family_ids.setdefault(
                ontology_id,
                {
                    "ontology_id": ontology_id,
                    "id": family.get("id", ""),
                    "database": family.get("database", ""),
                    "type": family.get("type", ""),
                    "name": family.get("name", ""),
                    "alias": family.get("alias", ""),
                    "alias_type": family.get("alias_type", ""),
                    "linked_interpro_ontology_id": family.get("linked_interpro_ontology_id", ""),
                    "linked_pfam_ontology_ids": family.get("linked_pfam_ontology_ids", []),
                    "source_database": family.get("source_database", ""),
                    "resource_url": family_resource_url(ontology_id),
                },
            )
        for linked_id in [family.get("linked_interpro_ontology_id", ""), *as_list(family.get("linked_pfam_ontology_ids", []))]:
            if linked_id:
                family_ids.setdefault(
                    linked_id,
                    {
                        "ontology_id": linked_id,
                        "id": linked_id.split(":", 1)[1] if ":" in linked_id else linked_id,
                        "database": linked_id.split(":", 1)[0] if ":" in linked_id else "",
                        "type": "linked_family_domain",
                        "name": family.get("name", ""),
                        "alias": family.get("alias", ""),
                        "source_database": family.get("source_database", ""),
                        "resource_url": family_resource_url(linked_id),
                    },
                )
        if representative.get("uniprot_accession"):
            acc = representative["uniprot_accession"]
            accessions.setdefault(
                acc,
                {
                    "accession": acc,
                    "kind": "representative",
                    "entry": representative.get("uniprot_entry", ""),
                    "gene_name": "",
                    "protein_name": "",
                    "organism": representative.get("organism", ""),
                    "taxon_id": "",
                    "reviewed_status": representative.get("reviewed_status", ""),
                    "source_database": representative.get("source_database", ""),
                },
            )
    return compact_dict(
        {
            "source": relpath(GENE_PROTEIN_NORMALIZATION_DIR),
            "entity_instance_id": best.get("entity_instance_id", ""),
            "global_node_id": best.get("global_node_id", ""),
            "canonical_form": best.get("canonical_form", ""),
            "aliases": best.get("aliases", []),
            "row_count": len(rows),
            "best": best,
            "rows": sorted_rows,
            "fasta_accessions": list(accessions.values()),
            "phytozome_ids": list(phytozome_ids.values()),
            "family_ids": list(family_ids.values()),
            "database_ids": list(database_ids.values()),
        }
    )


def load_gene_protein_normalizations(pmcid: str) -> dict[str, dict[str, Any]]:
    csv_path = GENE_PROTEIN_NORMALIZATION_DIR / f"{pmcid}.gene_protein_normalizations.csv"
    rows: list[dict[str, Any]] = []
    if source_exists(csv_path):
        rows = load_csv_rows(csv_path)

    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for raw in rows:
        entity_instance_id = raw.get("entity_instance_id", "")
        if not entity_instance_id:
            continue
        groups[entity_instance_id].append(slim_gene_protein_row(raw))

    by_key: dict[str, dict[str, Any]] = {}
    for entity_instance_id, slim_rows in groups.items():
        profile = build_gene_protein_profile(slim_rows)
        best_phytozome = (profile.get("best") or {}).get("phytozome", {})
        phytozome_ids = profile.get("phytozome_ids", [])
        family_ids = profile.get("family_ids", [])
        keys = {
            entity_instance_id,
            profile.get("global_node_id", ""),
            (profile.get("best") or {}).get("selected", {}).get("uniprot_accession", ""),
            (profile.get("best") or {}).get("representative", {}).get("uniprot_accession", ""),
            best_phytozome.get("gene_id", ""),
            best_phytozome.get("base_gene_id", ""),
            ((profile.get("best") or {}).get("family") or {}).get("ontology_id", ""),
            ((profile.get("best") or {}).get("family") or {}).get("id", ""),
        }
        for item in phytozome_ids:
            keys.update(
                {
                    item.get("gene_id", ""),
                    item.get("base_gene_id", ""),
                    item.get("ontology_id", ""),
                    item.get("base_ontology_id", ""),
                }
            )
        for item in family_ids:
            keys.update(
                {
                    item.get("ontology_id", ""),
                    item.get("id", ""),
                    item.get("alias", ""),
                }
            )
        for item in profile.get("database_ids", []):
            keys.update(
                {
                    item.get("ontology_id", ""),
                    item.get("identifier", ""),
                }
            )
        for key in keys:
            if key:
                by_key.setdefault(key, profile)
    return by_key


def gene_protein_ontology_concepts(entity: dict[str, Any]) -> list[dict[str, str]]:
    profile = entity.get("gene_protein_normalization") or {}
    concepts: list[dict[str, str]] = []
    for accession in profile.get("fasta_accessions", []):
        acc = accession.get("accession", "")
        if not acc:
            continue
        concepts.append(
            {
                "id": f"UniProt:{acc}",
                "ontology": "UniProt",
                "label": accession.get("gene_name") or accession.get("entry") or acc,
                "description": "Gene/protein normalization from Step 10.",
            }
        )
    for phytozome in profile.get("phytozome_ids", []):
        gene_id = phytozome.get("gene_id", "")
        if gene_id:
            concepts.append(
                {
                    "id": phytozome.get("ontology_id") or f"Phytozome:{gene_id}",
                    "ontology": "Phytozome",
                    "label": gene_id,
                    "description": "Phytozome gene model selected by Step 10 gene/protein normalization.",
                }
            )
        base_gene_id = phytozome.get("base_gene_id", "")
        if base_gene_id:
            concepts.append(
                {
                    "id": phytozome.get("base_ontology_id") or f"PhytozomeBase:{base_gene_id}",
                    "ontology": "PhytozomeBase",
                    "label": base_gene_id,
                    "description": "Base Phytozome gene ID derived from the selected transcript/protein model.",
                }
            )
    for family in profile.get("family_ids", []):
        ontology_id = family.get("ontology_id", "")
        if not ontology_id:
            continue
        ontology = ontology_id.split(":", 1)[0] if ":" in ontology_id else family.get("database", "")
        concepts.append(
            {
                "id": ontology_id,
                "ontology": ontology,
                "label": family.get("name") or family.get("alias") or family.get("id") or ontology_id,
                "description": "Gene/protein family or domain ontology selected by Step 10 normalization.",
            }
        )
    for item in profile.get("database_ids", []):
        ontology_id = item.get("ontology_id", "")
        if not ontology_id:
            continue
        concepts.append(
            {
                "id": ontology_id,
                "ontology": item.get("database", "") or ontology_id.split(":", 1)[0],
                "label": item.get("identifier", "") or ontology_id,
                "description": "Gene/protein database identifier selected by Step 10 normalization.",
            }
        )
    unique: dict[str, dict[str, str]] = {}
    for concept in concepts:
        if concept.get("id"):
            unique.setdefault(concept["id"], concept)
    return list(unique.values())


def selected_entity_concepts(raw: dict[str, Any], selected: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    """Return every ontology assignment selected for an entity.

    Step 9 preserves composite concepts in plural fields such as
    selected_ontology_ids, and the JSON form also carries selected_candidates.
    Keep those as first-class browser metadata instead of collapsing to the
    primary selected_ontology_id.
    """
    selected = selected or {}
    concepts: list[dict[str, Any]] = []

    for candidate in as_list(raw.get("selected_candidates")):
        if not isinstance(candidate, dict):
            continue
        concepts.append(
            compact_dict(
                {
                    "id": candidate.get("ontology_id", ""),
                    "ontology": candidate.get("ontology", ""),
                    "label": candidate.get("label", ""),
                    "description": candidate.get("description", ""),
                    "rank": candidate.get("rank", ""),
                }
            )
        )

    ids = as_list(raw.get("selected_ontology_ids"))
    if not ids and (raw.get("selected_ontology_id") or selected.get("ontology_id")):
        ids = [raw.get("selected_ontology_id") or selected.get("ontology_id", "")]
    ontologies = as_list(raw.get("selected_ontologies"))
    labels = as_list(raw.get("selected_labels"))
    descriptions = as_list(raw.get("selected_descriptions"))
    ranks = as_list(raw.get("selected_ranks"))

    for index, ontology_id in enumerate(ids):
        ontology_id = str(ontology_id or "").strip()
        if not ontology_id:
            continue
        ontology = str(ontologies[index] if index < len(ontologies) else "").strip()
        label = str(labels[index] if index < len(labels) else "").strip()
        description = str(descriptions[index] if index < len(descriptions) else "").strip()
        rank = ranks[index] if index < len(ranks) else ""
        concepts.append(
            compact_dict(
                {
                    "id": ontology_id,
                    "ontology": ontology or ontology_id.split(":", 1)[0],
                    "label": label or ontology_id,
                    "description": description,
                    "rank": rank,
                }
            )
        )

    if selected.get("ontology_id"):
        concepts.append(
            compact_dict(
                {
                    "id": selected.get("ontology_id", ""),
                    "ontology": selected.get("ontology", ""),
                    "label": selected.get("label", ""),
                    "description": selected.get("description", ""),
                    "rank": selected.get("rank", ""),
                }
            )
        )

    return unique_dicts_by_id(concepts)


def slim_entity(raw: dict[str, Any], fallback: dict[str, Any] | None = None) -> dict[str, Any]:
    fallback = fallback or {}
    selected = raw.get("selected_candidate") or {}
    selected_concepts = selected_entity_concepts(raw, selected)
    primary_concept = selected_concepts[0] if selected_concepts else selected
    entity = {
        "node_id": raw.get("node_id") or fallback.get("node_id", ""),
        "global_node_id": raw.get("global_node_id") or fallback.get("global_node_id", ""),
        "canonical_form": raw.get("canonical_form") or fallback.get("label", ""),
        "normalized_label": raw.get("normalized_label") or fallback.get("normalized_label", ""),
        "entity_type": raw.get("entity_type") or fallback.get("type", ""),
        "aliases": as_list(raw.get("aliases")),
        "roles": as_list(raw.get("roles") or fallback.get("roles")),
        "mention_count": raw.get("mention_count") or fallback.get("mention_count") or 0,
        "relation_count": raw.get("relation_count") or 0,
        "decision": raw.get("decision", ""),
        "decision_source": raw.get("decision_source", ""),
        "status": raw.get("status", ""),
        "candidate_count": raw.get("candidate_count", 0),
        "selected_ontology": raw.get("selected_ontology") or primary_concept.get("ontology", ""),
        "selected_ontology_id": raw.get("selected_ontology_id") or primary_concept.get("id") or primary_concept.get("ontology_id", ""),
        "selected_label": raw.get("selected_label") or primary_concept.get("label", ""),
        "selected_description": raw.get("selected_description") or primary_concept.get("description", ""),
        "selected_rank": raw.get("selected_rank") or primary_concept.get("rank", ""),
        "selected_concepts": selected_concepts,
        "selected_ontologies": [concept.get("ontology", "") for concept in selected_concepts if concept.get("ontology")],
        "selected_ontology_ids": [concept.get("id", "") for concept in selected_concepts if concept.get("id")],
        "selected_labels": [concept.get("label", "") for concept in selected_concepts if concept.get("label")],
        "selected_descriptions": [concept.get("description", "") for concept in selected_concepts if concept.get("description")],
        "selected_count": raw.get("selected_count") or len(selected_concepts),
        "evidence_sentence_ids": as_list(raw.get("evidence_sentence_ids")),
        "evidence_preview": raw.get("evidence_preview", ""),
        "auto_accept": raw.get("auto_accept") or {},
        "candidate_generation": raw.get("candidate_generation") or {},
        "provenance_ref": raw.get("provenance_ref") or {},
    }
    entity["external_links"] = external_links(entity)
    return entity


def load_entities(pmcid: str, hypergraph: dict[str, Any]) -> dict[str, dict[str, Any]]:
    fallback_by_node = {
        row.get("node_id", ""): row
        for row in hypergraph.get("entity_nodes", [])
        if row.get("node_id")
    }
    path = NORMALIZATION_DIR / f"{pmcid}.json"
    entities: dict[str, dict[str, Any]] = {}
    compounds_by_key = load_compound_classifications(pmcid)
    gene_proteins_by_key = load_gene_protein_normalizations(pmcid)

    if source_exists(path):
        for raw in load_json(path).get("entities", []):
            entity = slim_entity(raw, fallback_by_node.get(raw.get("node_id", "")))
            attach_compound_classification(entity, compounds_by_key)
            attach_gene_protein_normalization(entity, gene_proteins_by_key)
            if entity["node_id"]:
                entities[entity["node_id"]] = entity

    for node_id, raw in fallback_by_node.items():
        if node_id in entities:
            continue
        norm = GLOBAL_NORMALIZATIONS.get(pmcid, {}).get(node_id)
        if norm:
            for k, v in norm.items():
                if v:
                    raw[k] = v
        entity = slim_entity({}, raw)
        attach_compound_classification(entity, compounds_by_key)
        attach_gene_protein_normalization(entity, gene_proteins_by_key)
        if node_id:
            entities[node_id] = entity
    return entities


def attach_compound_classification(entity: dict[str, Any], compounds_by_key: dict[str, dict[str, Any]]) -> None:
    keys = [
        entity.get("provenance_ref", {}).get("entity_instance_id", ""),
        entity.get("global_node_id", ""),
        entity.get("selected_ontology_id", ""),
        *as_list(entity.get("selected_ontology_ids")),
    ]
    compound = next((compounds_by_key[key] for key in keys if key and key in compounds_by_key), None)
    if not compound:
        return
    entity["compound_classification"] = compound
    entity["external_links"] = external_links(entity)


def attach_gene_protein_normalization(entity: dict[str, Any], gene_proteins_by_key: dict[str, dict[str, Any]]) -> None:
    keys = [
        entity.get("provenance_ref", {}).get("entity_instance_id", ""),
        entity.get("global_node_id", ""),
        entity.get("canonical_form", ""),
        entity.get("selected_ontology_id", ""),
        *as_list(entity.get("selected_ontology_ids")),
    ]
    profile = next((gene_proteins_by_key[key] for key in keys if key and key in gene_proteins_by_key), None)
    if not profile:
        return
    entity["gene_protein_normalization"] = profile
    entity["external_links"] = external_links(entity)


def slim_relation(
    raw: dict[str, Any],
    sentence_map: dict[str, dict[str, Any]] | None = None,
    evaluation_map: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    sentence_map = sentence_map or {}
    evaluation_map = evaluation_map or {}
    evaluation = evaluation_map.get(str(raw.get("source_relation_id") or "").strip(), {})
    evidence_items = [
        {
            "sentence_id": item.get("sentence_id", ""),
            "section_type": item.get("section_type", ""),
            "text": item.get("text", ""),
        }
        for item in as_list(raw.get("evidence"))
        if isinstance(item, dict)
    ]
    evidence_preview = raw.get("evidence_preview", "") or " ".join(
        item["text"] for item in evidence_items[:2] if item.get("text")
    )
    return {
        "record_id": raw.get("record_id") or raw.get("hyperedge_id", ""),
        "source_relation_id": raw.get("source_relation_id", ""),
        "record_type": raw.get("record_type", ""),
        "merge_decision": raw.get("merge_decision", ""),
        "assertion_modifier": evaluation.get("assertion_modifier", ""),
        "relation_evaluation_verdict": evaluation.get("verdict", ""),
        "relation_evaluation_status": evaluation.get("status", ""),
        "relation_evaluation_confidence": evaluation.get("confidence", ""),
        "triple": raw.get("triple", ""),
        "subject": raw.get("subject", ""),
        "subject_type": raw.get("subject_type", ""),
        "subject_node_id": raw.get("subject_node_id", ""),
        "predicate": raw.get("predicate", ""),
        "predicate_class": raw.get("predicate_class", ""),
        "object": raw.get("object", ""),
        "object_type": raw.get("object_type", ""),
        "object_node_id": raw.get("object_node_id", ""),
        "context": raw.get("context") or {},
        "context_node_ids": as_list(raw.get("context_node_ids")),
        "event_ids": as_list(raw.get("event_ids") or raw.get("event_id")),
        "evidence_sentence_ids": as_list(raw.get("evidence_sentence_ids")),
        "evidence": evidence_items,
        "evidence_preview": evidence_preview,
        "evidence_context_text": raw.get("evidence_context_text", ""),
        "evidence_context_sentences": evidence_context_sentences(raw, sentence_map),
        "current_passage_id": raw.get("current_passage_id", ""),
        "current_passage_index": raw.get("current_passage_index", ""),
        "context_enrichment_source": raw.get("context_enrichment_source", ""),
        "context_enriched": bool(raw.get("context_enriched", False)),
        "relation_event_membership_count": raw.get("relation_event_membership_count", 0),
    }


def slim_event(raw: dict[str, Any]) -> dict[str, Any]:
    return {
        "event_id": raw.get("event_id", ""),
        "event_label": raw.get("event_label", ""),
        "event_type": raw.get("event_type", ""),
        "event_scope": raw.get("event_scope", ""),
        "relation_ids": as_list(raw.get("relation_ids")),
        "relation_count": raw.get("relation_count", 0),
        "participant_node_ids": as_list(raw.get("participant_node_ids")),
        "participants": [
            {
                "node_id": item.get("node_id", ""),
                "label": item.get("label", ""),
                "entity_type": item.get("entity_type", ""),
                "participant_role": item.get("participant_role", ""),
                "source_relation_ids": as_list(item.get("source_relation_ids")),
                "evidence_sentence_ids": as_list(item.get("evidence_sentence_ids")),
            }
            for item in as_list(raw.get("participants"))
            if isinstance(item, dict)
        ],
        "context": raw.get("context") or {},
        "evidence_sentence_ids": as_list(raw.get("evidence_sentence_ids")),
        "evidence_context_text": raw.get("evidence_context_text", ""),
        "source_candidate_event_ids": as_list(raw.get("source_candidate_event_ids")),
        "candidate_strategy": as_list(raw.get("candidate_strategy")),
        "candidate_strength": raw.get("candidate_strength", ""),
        "confidence": round_float(raw.get("confidence")),
        "reason_code": raw.get("reason_code", ""),
        "quality_flags": as_list(raw.get("quality_flags")),
    }


def dependency_tier(raw: dict[str, Any], fallback: str = "") -> str:
    verdict = str(raw.get("verdict", ""))
    label = str(raw.get("label", ""))
    if verdict == "accepted" or label == "1":
        return "accepted"
    if verdict == "hypothesis":
        return "hypothesis"
    if verdict == "rejected" or label == "0":
        return "rejected"
    return fallback or "candidate"


def slim_dependency(raw: dict[str, Any], tier: str = "") -> dict[str, Any]:
    return {
        "dependency_id": raw.get("dependency_id", ""),
        "tier": dependency_tier(raw, tier),
        "verdict": raw.get("verdict", ""),
        "support_verdict": raw.get("support_verdict", ""),
        "label": raw.get("label", ""),
        "dependency_type": raw.get("dependency_type") or raw.get("inter_event_relation_type", ""),
        "dependency_scope": raw.get("dependency_scope", ""),
        "upstream_event_id": raw.get("upstream_event_id") or raw.get("source_id", ""),
        "downstream_event_id": raw.get("downstream_event_id") or raw.get("target_id", ""),
        "confidence": round_float(raw.get("confidence")),
        "reason_code": raw.get("reason_code", ""),
        "dependency_origin": raw.get("dependency_origin", ""),
        "candidate_dependency_ids": as_list(raw.get("candidate_dependency_ids")),
        "origin_candidate_dependency_ids": as_list(raw.get("origin_candidate_dependency_ids")),
        "bridge_entities": as_list(raw.get("bridge_entities")),
        "supporting_relation_pairs": as_list(raw.get("supporting_relation_pairs")),
        "evidence_sentence_ids": as_list(raw.get("evidence_sentence_ids")),
    }


def compact_sources(source_files: dict[str, Any]) -> dict[str, str]:
    return {key: relpath(Path(value)) if value else "" for key, value in source_files.items()}


def compact_context_propagation_item(item: dict[str, Any]) -> dict[str, Any]:
    return compact_dict(
        {
            "entity_id": item.get("entity_id", ""),
            "label": item.get("label", ""),
            "entity_type": item.get("entity_type", ""),
            "ontology_id": item.get("ontology_id", ""),
            "selected_label": item.get("selected_label", ""),
            "kind": item.get("kind", ""),
            "mark": item.get("mark", ""),
            "agreement_status": item.get("agreement_status", ""),
            "agreement_source": item.get("agreement_source", ""),
            "provenance_count": item.get("provenance_count", 0),
        }
    )


def compact_context_propagation_group(group: dict[str, Any]) -> dict[str, Any]:
    group = group or {}
    return {
        "display": group.get("display", ""),
        "taxa": [compact_context_propagation_item(item) for item in as_list(group.get("taxa")) if isinstance(item, dict)],
        "tissues": [compact_context_propagation_item(item) for item in as_list(group.get("tissues")) if isinstance(item, dict)],
        "assays": [compact_context_propagation_item(item) for item in as_list(group.get("assays")) if isinstance(item, dict)],
        "provenance_count": len(as_list(group.get("provenance"))),
    }


def compact_context_propagation(row: dict[str, Any] | None) -> dict[str, Any]:
    row = row or {}
    return {
        "direct": compact_context_propagation_group(row.get("direct_taxon_tissue_context") or {}),
        "event": compact_context_propagation_group(row.get("event_taxon_tissue_context") or {}),
        "entity_linked": compact_context_propagation_group(row.get("entity_linked_taxon_tissue_context") or {}),
        "agreement": compact_context_propagation_group(row.get("agreement_taxon_tissue_context") or {}),
    }


def load_context_propagation(pmcid: str) -> dict[str, dict[str, Any]]:
    path = CONTEXT_PROPAGATION_DIR / f"{pmcid}.relation_context_propagation.json"
    if not source_exists(path):
        return {}
    payload = load_json(path)
    return {
        row.get("relation_id", ""): compact_context_propagation(row)
        for row in payload.get("relations", [])
        if row.get("relation_id")
    }


def build_paper(pmcid: str) -> dict[str, Any]:
    hypergraph = load_json(HYPERGRAPH_DIR / f"{pmcid}.json")
    article, sentence_map = article_metadata(pmcid)
    entities_by_id = load_entities(pmcid, hypergraph)
    relation_evaluations = load_relation_evaluations(pmcid)
    context_propagation_by_relation = load_context_propagation(pmcid)

    relations = [
        slim_relation(row, sentence_map, relation_evaluations)
        for row in hypergraph.get("relation_hyperedges", [])
    ]
    for relation in relations:
        relation["taxon_tissue_context"] = context_propagation_by_relation.get(
            relation.get("record_id", ""),
            compact_context_propagation({}),
        )
    relation_ids = {row["record_id"] for row in relations}
    events = [slim_event(row) for row in hypergraph.get("atomic_events", [])]
    event_ids = {row["event_id"] for row in events}

    dependencies: list[dict[str, Any]] = []
    seen_dependencies: set[str] = set()
    for row in hypergraph.get("inter_event_relations", []):
        dep = slim_dependency(row, "accepted")
        if dep["upstream_event_id"] in event_ids and dep["downstream_event_id"] in event_ids:
            dependencies.append(dep)
            seen_dependencies.add(dep["dependency_id"])
    hypothesis_rows = hypergraph.get("hypothesis_inter_event_relations", [])
    for row in hypothesis_rows:
        dep = slim_dependency(row, "hypothesis")
        if dep["upstream_event_id"] in event_ids and dep["downstream_event_id"] in event_ids:
            dependencies.append(dep)
            seen_dependencies.add(dep["dependency_id"])
    for row in hypergraph.get("inter_event_relation_labels", []):
        if dependency_tier(row) != "rejected":
            continue
        dep = slim_dependency(row, "rejected")
        if dep["dependency_id"] in seen_dependencies:
            continue
        if dep["upstream_event_id"] in event_ids and dep["downstream_event_id"] in event_ids:
            dependencies.append(dep)
            seen_dependencies.add(dep["dependency_id"])

    entity_relation_counts: Counter[str] = Counter()
    entity_event_counts: Counter[str] = Counter()
    for rel in relations:
        for node_id in [rel["subject_node_id"], rel["object_node_id"], *rel["context_node_ids"]]:
            if node_id:
                entity_relation_counts[node_id] += 1
    for event in events:
        for node_id in event["participant_node_ids"]:
            if node_id:
                entity_event_counts[node_id] += 1

    for node_id, entity in entities_by_id.items():
        entity["relation_count"] = max(int(entity.get("relation_count") or 0), entity_relation_counts[node_id])
        entity["event_count"] = entity_event_counts[node_id]

    dependency_counts_by_event: dict[str, Counter[str]] = defaultdict(Counter)
    for dep in dependencies:
        for event_id in [dep["upstream_event_id"], dep["downstream_event_id"]]:
            dependency_counts_by_event[event_id][dep["tier"]] += 1

    for event in events:
        counts = dependency_counts_by_event[event["event_id"]]
        event["dependency_counts"] = {
            "accepted": int(counts.get("accepted", 0)),
            "hypothesis": int(counts.get("hypothesis", 0)),
            "rejected": int(counts.get("rejected", 0)),
        }
        event["has_accepted_dependency"] = counts.get("accepted", 0) > 0

    referenced_sentence_ids: set[str] = set()
    for collection in [events, relations, dependencies, list(entities_by_id.values())]:
        for row in collection:
            referenced_sentence_ids.update(as_list(row.get("evidence_sentence_ids")))
            for item in as_list(row.get("evidence_context_sentences")):
                if isinstance(item, dict) and item.get("id"):
                    referenced_sentence_ids.add(item["id"])
            if isinstance(row.get("participants"), list):
                for item in row["participants"]:
                    referenced_sentence_ids.update(as_list(item.get("evidence_sentence_ids")))

    sentences = [sentence_map[key] for key in sorted(referenced_sentence_ids) if key in sentence_map]

    stats = {
        "events": len(events),
        "relations": len(relations),
        "entities": len(entities_by_id),
        "dependencies": len(dependencies),
        "accepted_dependencies": sum(1 for dep in dependencies if dep["tier"] == "accepted"),
        "hypothesis_dependencies": sum(1 for dep in dependencies if dep["tier"] == "hypothesis"),
        "rejected_dependencies": sum(1 for dep in dependencies if dep["tier"] == "rejected"),
        "events_without_accepted_dependency": sum(1 for row in events if not row["has_accepted_dependency"]),
        "normalized_entities": sum(
            1
            for row in entities_by_id.values()
            if row.get("selected_ontology_id") or gene_protein_ontology_concepts(row)
        ),
        "compound_classifications": sum(1 for row in entities_by_id.values() if row.get("compound_classification")),
        "gene_protein_normalizations": sum(1 for row in entities_by_id.values() if row.get("gene_protein_normalization")),
        "gene_protein_phytozome_entities": sum(
            1
            for row in entities_by_id.values()
            if (row.get("gene_protein_normalization") or {}).get("phytozome_ids")
        ),
        "gene_protein_family_entities": sum(
            1
            for row in entities_by_id.values()
            if (row.get("gene_protein_normalization") or {}).get("family_ids")
        ),
        "gene_protein_fasta_entities": sum(
            1
            for row in entities_by_id.values()
            if (row.get("gene_protein_normalization") or {}).get("fasta_accessions")
        ),
        "relation_ids": len(relation_ids),
    }

    return {
        "pmcid": pmcid,
        "article": article,
        "stats": stats,
        "source_files": compact_sources(hypergraph.get("source_files", {})),
        "metadata": hypergraph.get("metadata", {}),
        "events": sorted(events, key=lambda row: row["event_id"]),
        "relations": sorted(relations, key=lambda row: row["record_id"]),
        "dependencies": sorted(
            dependencies,
            key=lambda row: (
                {"accepted": 0, "hypothesis": 1, "rejected": 2}.get(row["tier"], 9),
                row["dependency_id"],
            ),
        ),
        "entities": sorted(entities_by_id.values(), key=lambda row: row["canonical_form"].lower()),
        "sentences": sentences,
    }


def calculate_enrichments(entities: list[dict[str, Any]], relations: list[dict[str, Any]]) -> None:
    def is_gene(t: str) -> bool:
        t = t.lower()
        return "gene" in t or "protein" in t
    
    def is_trait(t: str) -> bool:
        t = t.lower()
        return "trait" in t or "pathway" in t or "phenotype" in t or "process" in t or "disease" in t

    gene_entities = [e for e in entities if is_gene(e.get("type", ""))]
    trait_entities = [e for e in entities if is_trait(e.get("type", ""))]

    node_to_concept = {}
    for e in entities:
        concept_id = e.get("ontology_id")
        if concept_id:
            node_to_concept[e["id"]] = concept_id

    gene_concepts = {node_to_concept[e["id"]] for e in gene_entities if e["id"] in node_to_concept}
    trait_concepts = {node_to_concept[e["id"]] for e in trait_entities if e["id"] in node_to_concept}

    if not gene_concepts or not trait_concepts:
        return

    total_relations = len(relations)

    gene_counts = Counter()
    trait_counts = Counter()
    co_counts = Counter()

    for rel in relations:
        node_ids = [rel.get("subject_entity_id"), rel.get("object_entity_id")] + rel.get("context_entity_ids", [])
        node_ids = [n for n in node_ids if n]
        
        rel_concepts = {node_to_concept[n] for n in node_ids if n in node_to_concept}
        
        g_c = rel_concepts & gene_concepts
        t_c = rel_concepts & trait_concepts
        
        for g in g_c:
            gene_counts[g] += 1
        for t in t_c:
            trait_counts[t] += 1
            
        for g in g_c:
            for t in t_c:
                co_counts[(g, t)] += 1

    p_values = []
    tests = []

    for (g, t), v1 in co_counts.items():
        v2 = gene_counts[g] - v1
        v3 = trait_counts[t] - v1
        v4 = total_relations - v1 - v2 - v3
        
        _, p = stats.fisher_exact([[v1, v2], [v3, v4]], alternative="greater")
        
        tests.append((g, t))
        p_values.append(p)
    
    if not p_values:
        return

    fdr_adjusted = stats.false_discovery_control(p_values)
    
    significant_pairs = defaultdict(list)
    concept_to_label = {e["ontology_id"]: e.get("selected_label") or e["ontology_id"] for e in trait_entities if "ontology_id" in e}

    for i, (g, t) in enumerate(tests):
        if fdr_adjusted[i] < 0.05:
            significant_pairs[g].append({
                "trait_concept": t,
                "trait_label": concept_to_label.get(t, t),
                "p_value": p_values[i],
                "fdr": fdr_adjusted[i]
            })
            
    for e in gene_entities:
        concept = e.get("ontology_id")
        if concept and concept in significant_pairs:
            sorted_enrichments = sorted(significant_pairs[concept], key=lambda x: x["fdr"])
            e["enrichments"] = sorted_enrichments


def build_global_path_index(papers: list[dict[str, Any]], generated_at: str) -> dict[str, Any]:
    """Create a compact cross-paper graph for the Pathfinder view."""
    entities: list[dict[str, Any]] = []
    events: list[dict[str, Any]] = []
    relations: list[dict[str, Any]] = []
    dependencies: list[dict[str, Any]] = []
    concepts: dict[str, dict[str, Any]] = {}

    for paper in papers:
        pmcid = paper["pmcid"]
        title = paper.get("article", {}).get("title", "")

        for entity in paper.get("entities", []):
            ontology_id = entity.get("selected_ontology_id", "")
            ontology = entity.get("selected_ontology", "")
            selected_label = entity.get("selected_label", "")
            selected_concepts = [
                concept
                for concept in as_list(entity.get("selected_concepts"))
                if isinstance(concept, dict) and concept.get("id")
            ]
            if ontology_id and not any(concept.get("id") == ontology_id for concept in selected_concepts):
                selected_concepts.insert(
                    0,
                    {
                        "id": ontology_id,
                        "ontology": ontology,
                        "label": selected_label or ontology_id,
                        "description": entity.get("selected_description", ""),
                    },
                )
            derived_concepts = gene_protein_ontology_concepts(entity)
            ontology_ids = [
                *[concept["id"] for concept in selected_concepts if concept.get("id")],
                *[concept["id"] for concept in derived_concepts if concept.get("id")],
            ]
            display_concept = selected_concepts[0] if selected_concepts else (derived_concepts[0] if derived_concepts else {})
            display_ontology_id = display_concept.get("id", "")
            display_ontology = display_concept.get("ontology", "")
            display_label = display_concept.get("label", "")
            node_id = entity.get("node_id", "")
            entities.append(
                {
                    "id": node_id,
                    "pmcid": pmcid,
                    "paper_title": title,
                    "label": entity.get("canonical_form") or entity.get("normalized_label") or node_id,
                    "type": entity.get("entity_type", ""),
                    "normalized_label": entity.get("normalized_label", ""),
                    "ontology": display_ontology,
                    "ontology_id": display_ontology_id,
                    "ontology_ids": [item for item in dict.fromkeys(ontology_ids) if item],
                    "selected_label": display_label,
                    "selected_description": display_concept.get("description", "") or entity.get("selected_description", ""),
                    "selected_concepts": selected_concepts,
                    "decision": entity.get("decision", ""),
                    "relation_count": entity.get("relation_count", 0),
                    "event_count": entity.get("event_count", 0),
                    "compound_classification": entity.get("compound_classification") or {},
                    "gene_protein_normalization": entity.get("gene_protein_normalization") or {},
                }
            )
            concept_inputs = []
            concept_inputs.extend(selected_concepts)
            concept_inputs.extend(derived_concepts)
            for concept_input in concept_inputs:
                concept_id = concept_input.get("id", "")
                if not concept_id:
                    continue
                concept = concepts.setdefault(
                    concept_id,
                    {
                        "id": concept_id,
                        "ontology": concept_input.get("ontology", ""),
                        "label": concept_input.get("label", "") or concept_id,
                        "description": concept_input.get("description", ""),
                        "entity_ids": [],
                        "papers": set(),
                        "types": Counter(),
                    },
                )
                concept["entity_ids"].append(node_id)
                concept["papers"].add(pmcid)
                concept["types"][entity.get("entity_type", "")] += 1

        for event in paper.get("events", []):
            events.append(
                {
                    "id": event.get("event_id", ""),
                    "pmcid": pmcid,
                    "label": event.get("event_label", ""),
                    "type": event.get("event_type", ""),
                    "scope": event.get("event_scope", ""),
                    "relation_count": event.get("relation_count", 0),
                    "accepted_dependency_count": event.get("dependency_counts", {}).get("accepted", 0),
                    "participant_entity_ids": event.get("participant_node_ids", []),
                    "evidence_sentence_ids": event.get("evidence_sentence_ids", []),
                }
            )

        for relation in paper.get("relations", []):
            relations.append(
                {
                    "id": relation.get("record_id", ""),
                    "pmcid": pmcid,
                    "triple": relation.get("triple", ""),
                    "predicate": relation.get("predicate", ""),
                    "predicate_class": relation.get("predicate_class", ""),
                    "assertion_modifier": relation.get("assertion_modifier", ""),
                    "relation_evaluation_verdict": relation.get("relation_evaluation_verdict", ""),
                    "subject_entity_id": relation.get("subject_node_id", ""),
                    "object_entity_id": relation.get("object_node_id", ""),
                    "context_entity_ids": relation.get("context_node_ids", []),
                    "taxon_tissue_context": relation.get("taxon_tissue_context") or compact_context_propagation({}),
                    "event_ids": relation.get("event_ids", []),
                    "evidence_sentence_ids": relation.get("evidence_sentence_ids", []),
                    "evidence_preview": " ".join(
                        item.get("text", "")
                        for item in relation.get("evidence", [])[:2]
                        if isinstance(item, dict)
                    ),
                }
            )

        for dep in paper.get("dependencies", []):
            dependencies.append(
                {
                    "id": dep.get("dependency_id", ""),
                    "pmcid": pmcid,
                    "tier": dep.get("tier", ""),
                    "type": dep.get("dependency_type", ""),
                    "source_event_id": dep.get("upstream_event_id", ""),
                    "target_event_id": dep.get("downstream_event_id", ""),
                    "confidence": dep.get("confidence", ""),
                    "reason_code": dep.get("reason_code", ""),
                    "bridge_entities": dep.get("bridge_entities", []),
                    "evidence_sentence_ids": dep.get("evidence_sentence_ids", []),
                    "supporting_relation_pairs": dep.get("supporting_relation_pairs", [])[:8],
                }
            )

    calculate_enrichments(entities, relations)

    concept_rows: list[dict[str, Any]] = []
    for concept in concepts.values():
        concept_rows.append(
            {
                "id": concept["id"],
                "ontology": concept["ontology"],
                "label": concept["label"],
                "description": concept["description"],
                "entity_ids": sorted(concept["entity_ids"]),
                "papers": sorted(concept["papers"]),
                "types": dict(concept["types"]),
            }
        )

    return {
        "generated_at": generated_at,
        "stats": {
            "entities": len(entities),
            "concepts": len(concept_rows),
            "events": len(events),
            "relations": len(relations),
            "dependencies": len(dependencies),
        },
        "entities": sorted(entities, key=lambda row: (row["label"].lower(), row["pmcid"], row["id"])),
        "concepts": sorted(concept_rows, key=lambda row: row["id"]),
        "events": sorted(events, key=lambda row: row["id"]),
        "relations": sorted(relations, key=lambda row: row["id"]),
        "dependencies": sorted(dependencies, key=lambda row: row["id"]),
    }


def build(outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    data_dir = outdir / "data"
    if data_dir.exists():
        shutil.rmtree(data_dir)
    (data_dir / "papers").mkdir(parents=True, exist_ok=True)

    papers = sorted({
        path.name.split(".", 1)[0]
        for path in HYPERGRAPH_DIR.glob("PMC*.json*")
        if path.name.endswith(".json") or path.name.endswith(".json.zst")
    })
    paper_summaries: list[dict[str, Any]] = []
    paper_payloads: list[dict[str, Any]] = []
    totals: Counter[str] = Counter()
    
    print(f"Discovered {len(papers)} papers to process.")
    for pmcid in tqdm(papers, desc="Compiling knowledge graph", unit="paper"):
        try:
            payload = build_paper(pmcid)
            write_json(data_dir / "papers" / f"{pmcid}.json", payload)
            paper_payloads.append(payload)
            paper_summaries.append(
                {
                    "pmcid": pmcid,
                    "title": payload["article"].get("title", ""),
                    "year": payload["article"].get("year", ""),
                    "doi": payload["article"].get("doi", ""),
                    "stats": payload["stats"],
                }
            )
            totals.update(payload["stats"])
        except Exception as e:
            tqdm.write(f"\n[ERROR] Failed to process paper {pmcid} due to: {type(e).__name__}: {e}")
            tqdm.write(f"Traceback for {pmcid}:\n{traceback.format_exc()}")
            continue

    global_path_index = build_global_path_index(paper_payloads, generated_at)
    write_json(data_dir / "global_path_index.json", global_path_index)

    manifest = {
        "generated_at": generated_at,
        "source": {
            "hypergraph": relpath(HYPERGRAPH_DIR),
            "sentences": relpath(SENTENCE_DIR),
            "normalization": relpath(NORMALIZATION_DIR),
            "compound_classification": relpath(COMPOUND_CLASSIFICATION_DIR),
            "gene_protein_normalization": relpath(GENE_PROTEIN_NORMALIZATION_DIR),
            "context_propagation": relpath(CONTEXT_PROPAGATION_DIR),
        },
        "papers": paper_summaries,
        "totals": {key: int(value) for key, value in totals.items()},
        "path_index": {
            "file": "data/global_path_index.json",
            "stats": global_path_index["stats"],
            "ontology_bridge_count": global_path_index["stats"]["concepts"],
        },
        "notes": [
            "Generated from existing PSFD pipeline outputs.",
            "Absolute local paths are removed from the public data bundle.",
            "Rejected dependency candidates come from inter_event_relation_labels.csv.",
            "Pathfinder uses normalized ontology IDs as cross-paper bridge nodes.",
            "Compound entities include Step 930 ClassyFire and NPClassifier metadata when available.",
            "Gene/protein entities include Step 10 UniProt normalization metadata when available.",
            "Relations include Step 11 event-propagated and entity-linked taxon/tissue context when available.",
        ],
    }
    write_json(data_dir / "manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
Build the static PSFD JSON bundle.
"""
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT,
        help="Sequence-demo repository root. The script writes outdir/data.",
    )
    return parser.parse_args()


def build_database(outdir: Path) -> None:
    global HYPERGRAPH_DIR, SENTENCE_DIR, TRIPLES_EVALUATION_DIR, NORMALIZATION_DIR
    global COMPOUND_CLASSIFICATION_DIR, GENE_PROTEIN_NORMALIZATION_DIR
    global CONTEXT_PROPAGATION_DIR

    HYPERGRAPH_DIR = resolve_input_dir("hypergraph_core_relations", "720_hypergraph", "7_relations_grouping")
    SENTENCE_DIR = resolve_input_dir("bioc_sentences_spacy", "111_bioc_sentences_spacy")
    TRIPLES_EVALUATION_DIR = resolve_input_dir("triples_evaluations", "410_LLM_triples_evaluation", "4_triples_evaluation")
    NORMALIZATION_DIR = resolve_input_dir("normalized_entities", "920_normalized_entities", "9_entity_normalization_llm")
    COMPOUND_CLASSIFICATION_DIR = resolve_input_dir("compound_classifications", "930_compound_classifications")
    GENE_PROTEIN_NORMALIZATION_DIR = resolve_input_dir("gene_protein_normalization", "940_gene_protein_normalization", "10_gene_protein_normalization")
    CONTEXT_PROPAGATION_DIR = resolve_input_dir("context_propagation", "11_context_propagation")
    
    preload_global_normalizations()
    manifest = build(outdir.resolve())
    print(f"Read outputs from {INPUT_DIR.resolve()}")
    print(f"Wrote PSFD static data to {outdir.resolve() / 'data'}")
    print(f"Papers: {len(manifest['papers'])}")
    print(f"Events: {manifest['totals'].get('events', 0)}")
    print(f"Relations: {manifest['totals'].get('relations', 0)}")
    print(f"Dependencies: {manifest['totals'].get('dependencies', 0)}")

def main() -> None:
    args = parse_args()
    build_database(args.outdir)

if __name__ == "__main__":
    main()
