import os
import subprocess
import pandas as pd
import json
import re
import uuid
import asyncio
from collections import defaultdict
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any

# ==========================================
# Caching System for Search & Extraction
# ==========================================

_search_cache = {}
_search_cache_lock = asyncio.Lock()

async def get_cached_search(sequence: str, method: str):
    normalized_seq = sequence.strip().upper()
    key = (normalized_seq, method)
    async with _search_cache_lock:
        return _search_cache.get(key)

async def set_cached_search(sequence: str, method: str, results: List[dict]):
    normalized_seq = sequence.strip().upper()
    key = (normalized_seq, method)
    async with _search_cache_lock:
        if len(_search_cache) > 2000:
            _search_cache.clear()
        _search_cache[key] = results

_extract_cache = {}
_extract_cache_lock = asyncio.Lock()

async def get_cached_extract(compounds: str, fasta: str, attributes: dict, method: str):
    frozen_attrs = tuple(sorted(attributes.items())) if attributes else ()
    key = (compounds.strip(), fasta.strip(), frozen_attrs, method)
    async with _extract_cache_lock:
        return _extract_cache.get(key)

async def set_cached_extract(compounds: str, fasta: str, attributes: dict, method: str, results: List[dict]):
    frozen_attrs = tuple(sorted(attributes.items())) if attributes else ()
    key = (compounds.strip(), fasta.strip(), frozen_attrs, method)
    async with _extract_cache_lock:
        if len(_extract_cache) > 1000:
            _extract_cache.clear()
        _extract_cache[key] = results

# ==========================================
# Helper Functions for Entity Mapping
# ==========================================

def as_list(val):
    if val is None:
        return []
    if isinstance(val, list):
        return val
    return [val]

def unique_strings_list(values) -> list:
    seen = set()
    result = []
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        low = s.lower()
        if low not in seen:
            seen.add(low)
            result.append(s)
    return result

def selected_concepts(entity) -> list:
    return [
        concept for concept in as_list(entity.get("selected_concepts"))
        if concept and isinstance(concept, dict) and (concept.get("id") or concept.get("ontology_id"))
    ]

def entity_ontology_ids(entity) -> list:
    ids = []
    if entity.get("selected_ontology_id"):
        ids.append(entity["selected_ontology_id"])
    if entity.get("ontology_id"):
        ids.append(entity["ontology_id"])
    for oid in as_list(entity.get("selected_ontology_ids")):
        if oid:
            ids.append(oid)
    for oid in as_list(entity.get("ontology_ids")):
        if oid:
            ids.append(oid)
    for concept in selected_concepts(entity):
        cid = concept.get("id") or concept.get("ontology_id")
        if cid:
            ids.append(cid)
    for oid in gene_protein_ontology_ids(entity):
        if oid:
            ids.append(oid)
    return unique_strings_list(ids)

def gene_protein_meta(entity):
    return entity.get("gene_protein_normalization") or None

def gene_protein_ontology_ids(entity) -> list:
    profile = gene_protein_meta(entity)
    if not profile:
        return []
    ids = []
    for item in as_list(profile.get("fasta_accessions")):
        acc = item.get("accession") if isinstance(item, dict) else ""
        if acc:
            ids.append(f"UniProt:{acc}")
    for item in as_list(profile.get("phytozome_ids")):
        if isinstance(item, dict):
            if item.get("ontology_id"):
                ids.append(item["ontology_id"])
            if item.get("base_ontology_id"):
                ids.append(item["base_ontology_id"])
    for item in as_list(profile.get("family_ids")):
        if isinstance(item, dict) and item.get("ontology_id"):
            ids.append(item["ontology_id"])
    for item in as_list(profile.get("database_ids")):
        if isinstance(item, dict) and item.get("ontology_id"):
            ids.append(item["ontology_id"])
    return unique_strings_list(ids)

def gene_protein_search_text(entity) -> str:
    profile = gene_protein_meta(entity)
    if not profile:
        return ""
    
    rows_text = []
    for row in as_list(profile.get("rows")):
        if isinstance(row, dict):
            row_parts = [
                str(row.get("gene_query") or ""),
                str(row.get("lookup_query") or ""),
                str(row.get("decision") or ""),
                str(row.get("status") or ""),
                str(row.get("normalization_scope") or ""),
                str(row.get("match_type") or ""),
                str(row.get("ambiguity_reason") or ""),
                json.dumps(row.get("phytozome") or {}),
                json.dumps(row.get("family") or {}),
                json.dumps(row.get("database_ids") or {}),
                json.dumps(row.get("raw_fields") or {})
            ]
            rows_text.append(" ".join(row_parts))
            
    parts = [
        str(profile.get("canonical_form") or ""),
        " ".join([str(a or "") for a in as_list(profile.get("aliases"))]),
        json.dumps(profile.get("best") or {}),
        json.dumps(profile.get("fasta_accessions") or {}),
        json.dumps(profile.get("phytozome_ids") or {}),
        json.dumps(profile.get("family_ids") or {}),
        json.dumps(profile.get("database_ids") or {}),
        " ".join(rows_text)
    ]
    return " ".join(parts)

def compound_meta(entity):
    return entity.get("compound_classification") or None

def compound_search_text(entity) -> str:
    compound = compound_meta(entity)
    if not compound:
        return ""
    cf = compound.get("classyfire") or {}
    np = compound.get("npclassifier") or {}
    chebi = compound.get("chebi") or {}
    pubchem = compound.get("pubchem") or {}
    structure = compound.get("structure") or {}
    raw = compound.get("raw_fields") or {}
    
    parts = [
        str(compound.get("compound_status") or ""),
        str(compound.get("classification_status") or ""),
        str(chebi.get("id") or ""),
        str(chebi.get("name") or ""),
        str(chebi.get("formula") or ""),
        str(chebi.get("inchikey") or ""),
        str(pubchem.get("cid") or ""),
        str(structure.get("inchikey") or ""),
        str(structure.get("smiles") or ""),
        str(cf.get("kingdom") or ""),
        str(cf.get("superclass") or ""),
        str(cf.get("class") or ""),
        str(cf.get("subclass") or ""),
        str(cf.get("direct_parent") or ""),
        str(np.get("pathway") or ""),
        str(np.get("superclass") or ""),
        str(np.get("class") or ""),
        " ".join([str(val or "") for val in raw.values()])
    ]
    return " ".join(parts)

def is_compound_like_global_entity(entity) -> bool:
    type_val = str(entity.get("type") or entity.get("entity_type") or "").lower()
    ids_str = " ".join([str(entity.get("ontology_id") or "")] + [str(x or "") for x in as_list(entity.get("ontology_ids"))])
    has_chebi_pubchem = bool(re.search(r'\b(CHEBI|ChEBI|PubChem):', ids_str, re.IGNORECASE))
    return type_val == "compound" or bool(entity.get("compound_classification")) or has_chebi_pubchem

def is_chemical_exposure_global_entity(entity) -> bool:
    type_val = str(entity.get("type") or entity.get("entity_type") or "").lower()
    if type_val != "experimental_condition":
        return False
    ontology_ids = entity_ontology_ids(entity)
    text = " ".join([
        str(entity.get("label") or ""),
        str(entity.get("selected_label") or ""),
        str(entity.get("normalized_label") or ""),
        str(entity.get("selected_description") or "")
    ] + ontology_ids).lower()
    has_peco = any(str(oid or "").upper().startswith("PECO:") for oid in ontology_ids)
    has_exposure = bool(re.search(r'\b(exposure|treatment|treated|application|hormone|acid)\b', text))
    return has_peco and has_exposure

def query_suggests_condition(query: str) -> bool:
    return bool(re.search(r'\b(exposure|treatment|treated|condition|stress|application|spray|dose|medium)\b', str(query or ""), re.IGNORECASE))

def compound_alias_search_text(entity) -> str:
    compound = entity.get("compound_classification") or {}
    chebi = compound.get("chebi") or {}
    pubchem = compound.get("pubchem") or {}
    normalization = compound.get("normalization") or {}
    raw = compound.get("raw_fields") or {}
    
    parts = [
        str(entity.get("label") or ""),
        str(entity.get("selected_label") or ""),
        str(entity.get("normalized_label") or ""),
        str(entity.get("ontology_id") or ""),
        str(compound.get("canonical_form") or ""),
        str(compound.get("aliases") or ""),
        str(normalization.get("selected_label") or ""),
        str(normalization.get("selected_ontology_id") or ""),
        str(chebi.get("id") or ""),
        str(chebi.get("name") or ""),
        f"PubChem:{pubchem.get('cid')}" if pubchem.get('cid') else "",
        str(raw.get("canonical_form") or ""),
        str(raw.get("aliases") or ""),
        str(raw.get("selected_label") or ""),
        str(raw.get("selected_ontology_id") or "")
    ]
    return " ".join(parts).lower()

def text_has_query_token(text: str, query: str) -> bool:
    if not query:
        return False
    escaped = re.escape(query)
    pattern = rf'(^|[^a-z0-9]){escaped}([^a-z0-9]|$)'
    return bool(re.search(pattern, str(text or ""), re.IGNORECASE))

def global_entity_search_text(entity) -> str:
    parts = [
        str(entity.get("id") or ""),
        str(entity.get("pmcid") or ""),
        str(entity.get("label") or ""),
        str(entity.get("normalized_label") or ""),
        str(entity.get("selected_label") or ""),
        str(entity.get("ontology_id") or ""),
        " ".join([str(oid or "") for oid in as_list(entity.get("ontology_ids"))]),
        " ".join([
            " ".join(filter(None, [concept.get("id"), concept.get("label"), concept.get("description")]))
            for concept in selected_concepts(entity)
        ]),
        str(entity.get("ontology") or ""),
        str(entity.get("type") or ""),
        str(entity.get("selected_description") or ""),
        compound_search_text(entity),
        gene_protein_search_text(entity)
    ]
    return " ".join(parts).lower()

def path_entity_name(entity) -> str:
    if not entity:
        return "Unknown entity"
    return (
        entity.get("label") or
        entity.get("canonical_form") or
        entity.get("selected_label") or
        entity.get("normalized_label") or
        entity.get("id") or
        entity.get("node_id") or
        "Unknown entity"
    )

def entity_evidence_weight(entity) -> int:
    return int(entity.get("relation_count") or 0) * 3 + int(entity.get("event_count") or 0) * 2 + (4 if len(entity_ontology_ids(entity)) else 0)

def entity_matches_annotation_category(entity, category: str = "auto") -> bool:
    if not category or category == "auto":
        return True
    categories_map = {
        "compound": ["compound"],
        "gene_protein": ["gene_protein", "gene", "protein"],
        "pathway_or_process": ["pathway_or_process"],
        "experimental_condition": ["experimental_condition"],
        "trait_function": ["plant_trait", "molecular_trait_or_function", "phenotype"],
        "anatomy": ["anatomical_structure", "anatomical_part", "cellular_component", "developmental_stage"],
        "context": ["genotype", "genetic_perturbation", "regulatory_motif"],
        "taxon": ["taxon"],
        "assay": ["assay_method", "assay_or_measurement"],
    }
    types_for_cat = categories_map.get(category, [])
    ent_type = str(entity.get("type") or entity.get("entity_type") or "").lower()
    
    if ent_type in types_for_cat:
        return True
    if category == "compound":
        return is_compound_like_global_entity(entity)
    if category == "gene_protein":
        profile = entity.get("gene_protein_normalization") or {}
        return bool(profile) or len(gene_protein_ontology_ids(entity)) > 0
    return False

def ranked_entity_matches(db, term: str, category: str = "auto") -> list:
    query = str(term or "").strip().lower()
    if not query:
        return []
    short_query = len(query) <= 3
    matches = []
    
    candidates = db.entities
    if category and category != "auto" and category in db.entities_by_category:
        candidates = db.entities_by_category[category]
        
    for entity in candidates:
        name = entity["_clean_name"]
        ontology = entity["_clean_ontology_id"]
        ontology_ids_set = entity["_clean_ontology_ids_set"]
        ontology_local_ids_set = entity["_clean_ontology_local_ids_set"]
        selected = entity["_clean_selected_label"]
        search = entity["_search_text"]
        compound_text = entity["_compound_alias_search_text"]
        score = 0
        
        if ontology == query:
            score += 120
        if query in ontology_ids_set:
            score += 120
        if query in ontology_local_ids_set:
            score += 110
        if name == query or selected == query:
            score += 100
            
        if short_query:
            if text_has_query_token(f"{name} {selected}", query):
                score += 48
            if text_has_query_token(search, query):
                score += 24
        else:
            if query in name or query in selected:
                score += 50
            if query in search:
                score += 30
                
        if entity["_is_compound"]:
            if name == query or selected == query or text_has_query_token(compound_text, query):
                score += 45
            if not short_query and query in compound_text:
                score += 30
                
        if entity["_is_chemical_exposure"] and not query_suggests_condition(query):
            score -= 40
            
        score += min(20, entity["_evidence_weight"])
        
        min_score = 44 if short_query else 20
        if score > min_score and (category == "auto" or category not in db.entities_by_category or entity_matches_annotation_category(entity, category)):
            matches.append((entity, score))
            
    matches.sort(key=lambda x: (-x[1], x[0]["_clean_name"]))
    return [m[0] for m in matches]

# ==========================================
# Database & Relationship Extraction Core
# ==========================================

class PSFDDatabase:
    def __init__(self):
        self.entities = []
        self.relations = []
        self.entity_by_id = {}
        self.relations_by_entity = defaultdict(list)
        self.entities_by_hash = defaultdict(list)
        self.entities_by_category = defaultdict(list)
        self.concepts_count = 0
        
    def load_database(self, filepath: str):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Database file not found at: {filepath}")
            
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
            
        self.entities = data.get("entities", [])
        self.relations = data.get("relations", [])
        
        self.entity_by_id = {}
        self.relations_by_entity = defaultdict(list)
        self.entities_by_hash = defaultdict(list)
        self.entities_by_category = defaultdict(list)
        
        concepts = set()
        for ent in self.entities:
            ent_id = ent.get("id")
            if not ent_id:
                continue
            self.entity_by_id[ent_id] = ent
            
            # Index Hash Suffix
            hash_suffix = ent_id.split(".")[-1]
            self.entities_by_hash[hash_suffix].append(ent)
            
            # Count concepts
            for oid in entity_ontology_ids(ent):
                concepts.add(oid)
                
            # Precompute fields for optimization
            ent["_clean_name"] = path_entity_name(ent).lower()
            ent["_clean_ontology_id"] = str(ent.get("ontology_id") or "").lower()
            
            ontology_ids = [str(oid or "").lower() for oid in as_list(ent.get("ontology_ids"))]
            ent["_clean_ontology_ids_set"] = set(ontology_ids)
            
            ontology_local_ids = []
            for oid in ontology_ids:
                if ":" in oid:
                    ontology_local_ids.append(oid.split(":", 1)[1])
                else:
                    ontology_local_ids.append(oid)
            ent["_clean_ontology_local_ids_set"] = set(ontology_local_ids)
            ent["_clean_selected_label"] = str(ent.get("selected_label") or "").lower()
            ent["_search_text"] = global_entity_search_text(ent)
            ent["_compound_alias_search_text"] = compound_alias_search_text(ent)
            
            is_comp = is_compound_like_global_entity(ent)
            ent["_is_compound"] = is_comp
            ent["_is_chemical_exposure"] = is_chemical_exposure_global_entity(ent)
            ent["_evidence_weight"] = entity_evidence_weight(ent)
            
            # Index categories to avoid scanning all entities
            if is_comp:
                self.entities_by_category["compound"].append(ent)
            
            profile = ent.get("gene_protein_normalization") or {}
            is_gene_prot = bool(profile) or len(gene_protein_ontology_ids(ent)) > 0
            if is_gene_prot:
                self.entities_by_category["gene_protein"].append(ent)
                
            ent_type = str(ent.get("type") or ent.get("entity_type") or "").lower()
            categories_map = {
                "compound": ["compound"],
                "gene_protein": ["gene_protein", "gene", "protein"],
                "pathway_or_process": ["pathway_or_process"],
                "experimental_condition": ["experimental_condition"],
                "trait_function": ["plant_trait", "molecular_trait_or_function", "phenotype"],
                "anatomy": ["anatomical_structure", "anatomical_part", "cellular_component", "developmental_stage"],
                "context": ["genotype", "genetic_perturbation", "regulatory_motif"],
                "taxon": ["taxon"],
                "assay": ["assay_method", "assay_or_measurement"],
            }
            for cat, types in categories_map.items():
                if ent_type in types:
                    self.entities_by_category[cat].append(ent)
                
        self.concepts_count = len(concepts)
            
        # Index relations
        for rel in self.relations:
            subject_ids = self._get_global_endpoint_ids(rel, "subject")
            object_ids = self._get_global_endpoint_ids(rel, "object")
            context_ids = rel.get("context_entity_ids", [])
            
            all_associated_ids = set(subject_ids + object_ids + context_ids)
            for eid in all_associated_ids:
                if eid:
                    self.relations_by_entity[eid].append(rel)

    def _get_global_endpoint_ids(self, rel: dict, role: str) -> list:
        single = rel.get(f"{role}_entity_id")
        plural = rel.get(f"{role}_entity_ids", [])
        ids = []
        if single:
            ids.append(single)
        if isinstance(plural, list):
            ids.extend(plural)
        return list(set(ids))

# Instantiate Global DB
db = PSFDDatabase()

def relation_has_endpoint(rel, entity_id: str) -> bool:
    return entity_id in relation_global_subject_ids(rel) or entity_id in relation_global_object_ids(rel)

def relation_global_subject_ids(rel) -> list:
    ids = as_list(rel.get("subject_entity_ids")) + as_list(rel.get("subject_entity_id"))
    return list(set([str(x) for x in ids if x]))

def relation_global_object_ids(rel) -> list:
    ids = as_list(rel.get("object_entity_ids")) + as_list(rel.get("object_entity_id"))
    return list(set([str(x) for x in ids if x]))

def relation_attribute_category(entity) -> dict:
    ent_type = str(entity.get("type") or entity.get("entity_type") or "")
    mapping = [
        {"key": "genes", "label": "Gene/protein", "types": ["gene_protein", "gene", "protein"]},
        {"key": "metabolites", "label": "Metabolite/compound", "types": ["compound"]},
        {"key": "pathways", "label": "Pathway/process", "types": ["pathway_or_process"]},
        {"key": "tissues", "label": "Tissue/anatomy", "types": ["anatomical_structure", "anatomical_part", "cellular_component", "developmental_stage"]},
        {"key": "species", "label": "Species/taxon", "types": ["taxon"]},
        {"key": "experimental_conditions", "label": "Experimental condition", "types": ["experimental_condition", "condition_parameter", "genotype", "genetic_perturbation"]},
        {"key": "plant_traits", "label": "Plant trait", "types": ["plant_trait", "phenotype"]},
        {"key": "molecular_traits", "label": "Molecular trait/function", "types": ["molecular_trait_or_function"]},
        {"key": "human_traits", "label": "Human trait", "types": ["human_trait"]},
    ]
    for item in mapping:
        if ent_type in item["types"]:
            return item
    return None

def clean_predicate(value: str) -> str:
    val = str(value or "unknown")
    if val.endswith("_event"):
        val = val[:-6]
    return val.replace("_", " ")

def clean_optional_display(value: str) -> str:
    val = str(value or "")
    if val.endswith("_event"):
        val = val[:-6]
    return val.replace("_", " ").strip()

def entity_stable_id(entity) -> str:
    return entity.get("node_id") or entity.get("id") or ""

def context_entity_merge_key(entity) -> str:
    if not entity:
        return ""
    ids = sorted([str(oid or "").strip().lower() for oid in entity_ontology_ids(entity) if oid])
    if not ids:
        return f"node:{str(entity_stable_id(entity)).lower()}"
    if len(ids) == 1:
        return f"ontology:{ids[0]}"
    return f"ontology-set:{'|'.join(ids)}"

def context_entity_member_ids(entity) -> list:
    ids = [entity_stable_id(entity)] + as_list(entity.get("merged_node_ids"))
    return unique_strings_list(ids)

def unique_by(items, key_fn) -> list:
    seen = set()
    result = []
    for item in items:
        k = key_fn(item)
        if k not in seen:
            seen.add(k)
            result.append(item)
    return result

def merge_context_entities_by_ontology(entities) -> list:
    groups = {}
    for entity in as_list(entities):
        if not entity:
            continue
        key = context_entity_merge_key(entity)
        if not key:
            continue
        if key not in groups:
            groups[key] = []
        groups[key].append(entity)
        
    result = []
    for key, items in groups.items():
        unique_items = unique_by(items, entity_stable_id)
        if len(unique_items) == 1:
            result.append(unique_items[0])
            continue
        
        representative = unique_items[0]
        aliases = []
        for entity in unique_items:
            aliases.extend([
                entity_name(entity),
                entity.get("canonical_form"),
                entity.get("normalized_label"),
                entity.get("selected_label")
            ] + as_list(entity.get("aliases")))
        aliases = unique_strings_list(aliases)
        
        ont_ids = []
        for entity in unique_items:
            ont_ids.extend(entity_ontology_ids(entity))
        ont_ids = unique_strings_list(ont_ids)
        
        sel_ont_ids = []
        for entity in unique_items:
            sel_ont_ids.extend(as_list(entity.get("selected_ontology_ids")))
        sel_ont_ids = unique_strings_list(sel_ont_ids)
        
        merged_ids = unique_strings_list([entity_stable_id(entity) for entity in unique_items])
        
        merged_entity = {
            **representative,
            "canonical_form": entity_name(representative),
            "aliases": aliases,
            "ontology_ids": ont_ids,
            "selected_ontology_ids": sel_ont_ids,
            "merged_node_ids": merged_ids,
            "context_merge_count": len(unique_items)
        }
        result.append(merged_entity)
    return result

def entity_name(entity) -> str:
    if not entity:
        return "Unknown entity"
    return entity.get("canonical_form") or entity.get("normalized_label") or entity.get("node_id") or ""

def is_useful_annotation_context(entity) -> bool:
    ent_type = str(entity.get("type") or entity.get("entity_type") or "")
    return ent_type in {
        "experimental_condition",
        "condition_parameter",
        "genotype",
        "genetic_perturbation",
        "taxon",
        "anatomical_structure",
        "cellular_component",
        "developmental_stage",
        "assay_method",
    }

def annotation_ontology_ids(entity) -> list:
    compound = entity.get("compound_classification") or {}
    chebi = compound.get("chebi") or {}
    pubchem = compound.get("pubchem") or {}
    normalization = compound.get("normalization") or {}
    
    ids = entity_ontology_ids(entity)
    if normalization.get("selected_ontology_id"):
        ids.append(normalization["selected_ontology_id"])
    if chebi.get("id"):
        ids.append(chebi["id"])
    if pubchem.get("cid"):
        ids.append(f"PubChem:{pubchem['cid']}")
        
    return unique_strings_list(ids)

def annotation_relation_context_entities(db, rel) -> list:
    context_entities = []
    for eid in as_list(rel.get("context_entity_ids")):
        entity = db.entity_by_id.get(eid)
        if entity and is_useful_annotation_context(entity):
            context_entities.append(entity)
            
    merged = merge_context_entities_by_ontology(context_entities)
    result = []
    for entity in merged:
        ids = entity_ontology_ids(entity)
        result.append({
            "id": entity_stable_id(entity),
            "mergeKey": context_entity_merge_key(entity),
            "memberIds": context_entity_member_ids(entity),
            "pmcid": entity.get("pmcid"),
            "label": path_entity_name(entity),
            "type": clean_predicate(entity.get("type") or entity.get("entity_type") or "context"),
            "ontologyId": ids[0] if ids else "",
            "ontologyLabel": entity.get("selected_label") or "",
            "color": "",
            "mergeCount": int(entity.get("context_merge_count") or 0),
            "aliases": as_list(entity.get("aliases")),
        })
    return result

def relation_taxon_tissue_context(rel) -> dict:
    return rel.get("taxon_tissue_context") or {}

def relation_taxon_tissue_display(rel, key: str) -> str:
    context = relation_taxon_tissue_context(rel)
    if key == "agreement":
        computed_display = clean_optional_display(format_taxon_tissue_display(
            category_aware_taxon_tissue_agreement_items(context.get("event") or {}, context.get("entity_linked") or {})
        ))
        if computed_display:
            return computed_display
        agreement_display = clean_optional_display(context.get("agreement", {}).get("display") or "")
        if agreement_display:
            return agreement_display
            
    display = format_taxon_tissue_display(relation_taxon_tissue_display_items(rel, key))
    baseline = context.get(key, {}).get("display") if isinstance(context.get(key), dict) else ""
    return clean_optional_display(display or baseline or "")

def relation_taxon_tissue_overlap_display(rel) -> str:
    agreement_display = relation_taxon_tissue_display(rel, "agreement")
    if agreement_display:
        return agreement_display
    context = relation_taxon_tissue_context(rel)
    computed_items = category_aware_taxon_tissue_agreement_items(context.get("event") or {}, context.get("entity_linked") or {})
    return clean_optional_display(format_taxon_tissue_display(computed_items))

def taxon_tissue_item_is_informative(item) -> bool:
    if not item or not isinstance(item, dict):
        return False
    label = str(item.get("label") or "").strip()
    id_str = str(item.get("ontology_id") or item.get("entity_id") or "").strip()
    norm_label = label.lower()
    norm_id = id_str.lower()
    return bool(label or id_str) and norm_label not in ("unknown", "-", "") and norm_id not in ("unknown", "-", "")

def taxon_tissue_comparable_key(item) -> str:
    if not taxon_tissue_item_is_informative(item):
        return ""
    kind = str(item.get("kind") or "").lower()
    id_str = str(item.get("ontology_id") or item.get("entity_id") or item.get("label") or "").strip().lower()
    return f"{kind}|{id_str}" if kind and id_str else ""

def format_taxon_tissue_display(items) -> str:
    informative = [item for item in items if taxon_tissue_item_is_informative(item)]
    parts = []
    for item in informative:
        label = item.get("label") or item.get("entity_id") or ""
        mark = item.get("mark") or ""
        ont = f" ({item['ontology_id']})" if item.get("ontology_id") else ""
        parts.append(f"{label}{mark}{ont}".strip())
    
    seen = set()
    result = []
    for p in parts:
        low = p.lower()
        if low not in seen:
            seen.add(low)
            result.append(p)
    return "; ".join(result)

def category_aware_taxon_tissue_agreement_items(event_group=None, entity_linked_group=None) -> list:
    if event_group is None:
        event_group = {}
    if entity_linked_group is None:
        entity_linked_group = {}
    items = []
    categories = ["taxa", "tissues", "assays"]
    for category in categories:
        event_items = [item for item in as_list(event_group.get(category)) if taxon_tissue_item_is_informative(item)]
        linked_items = [item for item in as_list(entity_linked_group.get(category)) if taxon_tissue_item_is_informative(item)]
        if not linked_items:
            continue
        event_keys = set([taxon_tissue_comparable_key(item) for item in event_items if taxon_tissue_comparable_key(item)])
        if event_keys:
            for item in linked_items:
                if taxon_tissue_comparable_key(item) in event_keys:
                    items.append({
                        **item,
                        "agreement_status": item.get("agreement_status") or "overlap",
                        "agreement_source": item.get("agreement_source") or "event_and_entity_linked"
                    })
            continue
        for item in linked_items:
            items.append({
                **item,
                "agreement_status": item.get("agreement_status") or "entity_linked_fallback",
                "agreement_source": item.get("agreement_source") or "entity_linked_category_fallback"
            })
    return unique_by(items, taxon_tissue_comparable_key)

def relation_taxon_tissue_display_items(rel, key: str) -> list:
    context = relation_taxon_tissue_context(rel)
    group = context.get(key) or {}
    baseline_keys = relation_taxon_tissue_baseline_keys(rel) if key == "entity_linked" else set()
    
    items = as_list(group.get("taxa")) + as_list(group.get("tissues")) + as_list(group.get("assays"))
    result = []
    for item in items:
        if not isinstance(item, dict):
            continue
        if key != "entity_linked":
            result.append({**item, "mark": ""})
        else:
            method = str(item.get("method") or "")
            is_context_only_propagation = method == "entity_linked_direct_context"
            is_triple_propagation = method == "entity_linked_location_relation"
            is_already_assigned = taxon_tissue_comparable_key(item) in baseline_keys
            mark = "*" if (is_context_only_propagation and not is_triple_propagation and not is_already_assigned) else ""
            result.append({**item, "mark": mark})
    return result

def relation_taxon_tissue_baseline_keys(rel) -> set:
    context = relation_taxon_tissue_context(rel)
    direct = context.get("direct") or {}
    event = context.get("event") or {}
    
    baseline_items = (
        as_list(direct.get("taxa")) + as_list(direct.get("tissues")) + as_list(direct.get("assays")) +
        as_list(event.get("taxa")) + as_list(event.get("tissues")) + as_list(event.get("assays"))
    )
    keys = set()
    for item in baseline_items:
        k = taxon_tissue_comparable_key(item)
        if k:
            keys.add(k)
    return keys

def normalized_entity_for_relation(entity) -> str:
    ids = annotation_ontology_ids(entity)
    first_id = ids[0] if ids else "unresolved"
    return f"{path_entity_name(entity)}|{first_id}"

def relation_extraction_row(db, query_name: str, query_entity: dict, rel: dict, attr_filters: dict) -> dict:
    subject = db.entity_by_id.get(rel.get("subject_entity_id"))
    obj = db.entity_by_id.get(rel.get("object_entity_id"))
    if not subject or not obj:
        return None
    
    query_entity_id = query_entity.get("id") or query_entity.get("node_id")
    query_is_subject = query_entity_id in relation_global_subject_ids(rel)
    query_is_object = query_entity_id in relation_global_object_ids(rel)
    
    other = obj if (query_is_subject and not query_is_object) else subject
    attribute = relation_attribute_category(other)
    if not attribute:
        return None
        
    attr_key = attribute["key"]
    if not attr_filters.get(attr_key, False):
        return None
        
    predicate = clean_predicate(rel.get("predicate") or rel.get("predicate_class") or "relates to")
    
    contexts = []
    for item in annotation_relation_context_entities(db, rel):
        lbl = item["label"]
        oid = item["ontologyId"]
        contexts.append(f"{lbl} ({oid})" if oid else lbl)
    context_str = "; ".join(unique_strings_list(contexts))
    
    event_taxon_tissue_context = relation_taxon_tissue_display(rel, "event")
    entity_linked_taxon_tissue_context = relation_taxon_tissue_display(rel, "entity_linked")
    overlap_taxon_tissue_context = relation_taxon_tissue_overlap_display(rel)
    
    return {
        "query_name": path_entity_name(query_entity) or query_name,
        "pmcid": rel.get("pmcid") or "",
        "relation_id": rel.get("id"),
        "attribute_entity_id": other.get("id") or other.get("node_id"),
        "attribute_key": attr_key,
        "relation": f"{path_entity_name(subject)} {predicate} {path_entity_name(obj)}",
        "context": context_str,
        "event_taxon_tissue_context": event_taxon_tissue_context,
        "entity_linked_taxon_tissue_context": entity_linked_taxon_tissue_context,
        "overlap_taxon_tissue_context": overlap_taxon_tissue_context,
        "attribute_type": attribute["label"],
        "normalized_relation": f"{normalized_entity_for_relation(subject)} {predicate} {normalized_entity_for_relation(obj)}",
    }

def relation_rows_for_query_entities(db, query_entities: list, attr_filters: dict) -> list:
    seen = set()
    rows = []
    for qe in query_entities:
        query_name = qe["query_name"]
        entity = qe["entity"]
        if not entity:
            continue
        ent_id = entity.get("id") or entity.get("node_id")
        if not ent_id:
            continue
        relations = db.relations_by_entity.get(ent_id) or []
        for rel in relations:
            if not relation_has_endpoint(rel, ent_id):
                continue
            row = relation_extraction_row(db, query_name, entity, rel, attr_filters)
            if not row:
                continue
            key = f"{row['query_name']}|{row['relation_id']}|{row['attribute_entity_id']}"
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)
            
    attribute_rank = [
        "genes", "metabolites", "pathways", "tissues", "species", 
        "experimental_conditions", "plant_traits", "molecular_traits", "human_traits"
    ]
    def sort_key(row):
        q_name = row["query_name"].lower()
        attr_key = row["attribute_key"]
        try:
            rank = attribute_rank.index(attr_key)
        except ValueError:
            rank = len(attribute_rank)
        rel_str = row["relation"].lower()
        return (q_name, rank, rel_str)
        
    rows.sort(key=sort_key)
    return rows

def clean_protein_sequence(sequence: str) -> str:
    return re.sub(r'[^A-Z*]', '', str(sequence or "").upper()).replace("*", "")

def parse_fasta_records(text: str) -> list:
    raw = str(text or "").strip()
    if not raw:
        return []
    if ">" not in raw:
        seq = clean_protein_sequence(raw)
        return [{"name": "query_1", "sequence": seq}] if seq else []
        
    records = []
    current = None
    for line in raw.splitlines():
        trimmed = line.strip()
        if not trimmed:
            continue
        if trimmed.startswith(">"):
            if current and current["sequence"]:
                current["sequence"] = clean_protein_sequence(current["sequence"])
                records.append(current)
            name = trimmed[1:].strip()
            if not name:
                name = f"query_{len(records) + 1}"
            current = {"name": name, "sequence": ""}
        elif current:
            current["sequence"] += trimmed
            
    if current and current["sequence"]:
        current["sequence"] = clean_protein_sequence(current["sequence"])
        records.append(current)
        
    return [r for r in records if len(r["sequence"]) >= 20]

app = FastAPI(
    title="PSMM API Dispatcher",
    description="Bridge API between the live frontend UI and the PlantStress-MechanismMap architecture."
)

# Production CORS Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://alenzimic.github.io",
        "http://localhost",
        "http://127.0.0.1",
        "http://localhost:3000",
        "http://localhost:8000",
        "http://127.0.0.1:8000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SearchRequest(BaseModel):
    sequence: str
    method: str  # e.g., 'embed2graph' or 'seq2graph'

class SearchResult(BaseModel):
    query: str
    uniprot_id: Optional[str] = None
    global_node_id: Optional[str] = None
    selected_protein_name: Optional[str] = None
    selected_gene_name: Optional[str] = None
    selected_organism: Optional[str] = None
    score: Optional[float] = None
    score_type: Optional[str] = None

def convert_dataframe_to_json(df: pd.DataFrame) -> List[dict]:
    """
    Translates the bridge pandas DataFrame into the JSON schema expected by the UI.
    """
    if df.empty:
        return []
        
    records = df.to_dict(orient="records")
    for r in records:
        for k, v in r.items():
            if pd.isna(v):
                r[k] = None
    return records

async def perform_search_internal(sequence: str, method: str) -> List[dict]:
    # Check cache first
    cached = await get_cached_search(sequence, method)
    if cached is not None:
        return cached

    # Generate unique filename for temporary FASTA and output CSV to avoid race conditions
    query_id = uuid.uuid4().hex
    query_file = f"query_{query_id}.fasta"
    output_csv = f"bridge_output_{query_id}.csv"
    
    try:
        with open(query_file, "w") as f:
            if not sequence.startswith(">"):
                f.write(">query_sequence\n")
            f.write(sequence)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to write query sequence: {str(e)}")

    import sys
    script_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    python_exe = sys.executable
    
    if method == "embed2graph":
        module_name = "psmm.bridges.embed2graph"
    elif method == "seq2graph":
        module_name = "psmm.bridges.seq2graph"
    else:
        raise HTTPException(status_code=400, detail="Invalid method. Must be 'embed2graph' or 'seq2graph'.")

    # 3. Safe Subprocess Execution
    try:
        # Use absolute paths because CWD changes to script_dir during execution
        abs_query_file = os.path.abspath(query_file)
        abs_output_csv = os.path.abspath(output_csv)
        
        loop = asyncio.get_event_loop()
        def run_subprocess():
            return subprocess.run(
                [python_exe, "-m", module_name, "--query", abs_query_file, "--output", abs_output_csv],
                capture_output=True,
                text=True,
                cwd=script_dir,
                timeout=120  # 2 minute timeout to prevent hanging the HPC server
            )
            
        result = await loop.run_in_executor(None, run_subprocess)
        
        if result.returncode != 0:
            raise HTTPException(status_code=500, detail=f"Bridge script failed: {result.stderr or result.stdout}")
            
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="Bridge script execution timed out.")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to execute bridge script: {str(e)}")
    finally:
        # Clean up files after processing
        if os.path.exists(query_file):
            os.remove(query_file)

    # 4. Read output and convert to JSON schema
    if not os.path.exists(output_csv):
        raise HTTPException(status_code=500, detail="Bridge script did not produce the expected output CSV.")

    try:
        df = pd.read_csv(output_csv)
        response_data = convert_dataframe_to_json(df)
        
        # Save to cache
        await set_cached_search(sequence, method, response_data)
        
        return response_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to parse bridge output: {str(e)}")
    finally:
        if os.path.exists(output_csv):
            os.remove(output_csv)

@app.post("/search", response_model=List[SearchResult])
async def search(request: SearchRequest):
    return await perform_search_internal(request.sequence, request.method)

# ==========================================
# Expose Additional API Endpoints
# ==========================================

class StatsResponse(BaseModel):
    entities: int
    concepts: int
    relations: int

class ExtractRequest(BaseModel):
    compounds: str = ""
    fasta: str = ""
    attributes: Optional[Dict[str, bool]] = None
    method: str = "embed2graph"

class ExtractResultRow(BaseModel):
    query_name: str
    pmcid: str
    relation_id: str
    attribute_entity_id: str
    attribute_key: str
    relation: str
    context: str
    event_taxon_tissue_context: str
    entity_linked_taxon_tissue_context: str
    overlap_taxon_tissue_context: str
    attribute_type: str
    normalized_relation: str

@app.get("/api/stats", response_model=StatsResponse)
async def get_stats():
    return StatsResponse(
        entities=len(db.entities),
        concepts=db.concepts_count,
        relations=len(db.relations)
    )

@app.post("/api/extract", response_model=List[ExtractResultRow])
async def extract(request: ExtractRequest):
    filters = request.attributes
    if filters is None:
        filters = {
            "genes": True,
            "metabolites": True,
            "pathways": True,
            "tissues": True,
            "species": True,
            "experimental_conditions": True,
            "plant_traits": True,
            "molecular_traits": True,
            "human_traits": True,
        }
    
    # Check cache first
    cached = await get_cached_extract(request.compounds, request.fasta, filters, request.method)
    if cached is not None:
        return cached
        
    query_entities = []
    
    # 1. Process compounds
    if request.compounds:
        terms = unique_strings_list(re.split(r'[\n;,]+', request.compounds))
        for term in terms:
            matched = ranked_entity_matches(db, term, "compound")
            if matched:
                query_entities.append({
                    "query_name": term,
                    "entity": matched[0],
                    "source": "compound"
                })
                
    # 2. Process FASTA
    if request.fasta:
        queries = parse_fasta_records(request.fasta)
        for q in queries:
            search_results = await perform_search_internal(q["sequence"], request.method)
            for item in search_results:
                global_node_id = item.get("global_node_id")
                if not global_node_id:
                    continue
                hash_suffix = global_node_id.split(".")[-1]
                matched_entities = db.entities_by_hash.get(hash_suffix, [])
                for entity in matched_entities:
                    query_entities.append({
                        "query_name": path_entity_name(entity),
                        "entity": entity,
                        "source": "fasta"
                    })
                    
    # 3. Extract relations
    rows = relation_rows_for_query_entities(db, query_entities, filters)
    
    # Cache results
    await set_cached_extract(request.compounds, request.fasta, filters, request.method, rows)
    
    return rows

@app.on_event("startup")
async def startup_event():
    db_file = os.environ.get("PSFD_DB_FILE")
    
    if not db_file:
        candidates = [
            os.path.join(os.path.dirname(__file__), "../psfd-sequence-annotation-demo/data/global_path_index.json"),
            os.path.join(os.path.dirname(__file__), "PSFD/psfd-sequence-annotation-demo-main/data/global_path_index.json"),
            os.path.join(os.path.dirname(__file__), "PSFD/psfd-sequence-annotation-demo/data/global_path_index.json"),
            "/local/storage/thomas/psfd-sequence-annotation-demo-copy/data/global_path_index.json",
            "/local/storage/thomas/psfd-sequence-annotation-demo/data/global_path_index.json",
        ]
        for path in candidates:
            if os.path.exists(path):
                db_file = path
                break
        if not db_file:
            db_file = "/local/storage/thomas/psfd-sequence-annotation-demo/data/global_path_index.json"

    print(f"Loading PSFD database from {db_file}...")
    try:
        db.load_database(db_file)
        print(f"Database loaded successfully: {len(db.entities)} entities, {db.concepts_count} concepts, {len(db.relations)} relations.")
    except Exception as e:
        print(f"Error loading database: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    # Do NOT hardcode port 8000. Accept port dynamically via env var.
    port = int(os.environ.get("PSMM_PORT", 8080))
    uvicorn.run("api_server:app", host="0.0.0.0", port=port, reload=False)
