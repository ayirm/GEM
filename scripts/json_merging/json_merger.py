#!/usr/bin/env python3

import json
import argparse
import sys

from pathlib import Path
from datetime import datetime, timezone

def load_json(path: Path):
    """
    Load a JSON file safely.
    """
    try:
        with path.open() as f:
            return json.load(f)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to load JSON file {path}: {e}")

def merge_annotations(parsed_genes: list,uniprot_to_kegg: dict,go_by_uniprot: dict,kegg_by_id: dict,) -> dict:
    """
    Merge gene records with GO annotations and KEGG data using UniProt_ID as key.
    """

    if not isinstance(parsed_genes, list):
        raise TypeError("parsed_genes must be a list of gene objects")

    merged = {}

    for gene in parsed_genes:
        uniprot_id = gene.get("UniProt_ID")

        # UniProt_ID is the only valid join key
        if not uniprot_id:
            continue

        kegg_id = uniprot_to_kegg.get(uniprot_id)

        merged[uniprot_id] = {
            **gene,
            "GO": go_by_uniprot.get(uniprot_id, []),
            "KEGG_ID": kegg_id,
            "KEGG": kegg_by_id.get(kegg_id) if kegg_id else None,
        }

    return merged

def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Merge parsed genes, GO annotations, UniProt→KEGG mapping, and KEGG data"
    )

    parser.add_argument("--parsed_json", required=True, type=Path)
    parser.add_argument("--mapping_json", required=True, type=Path)
    parser.add_argument("--go_json", required=True, type=Path)
    parser.add_argument("--kegg_json", required=True, type=Path)
    parser.add_argument("--merged_json", required=True, type=Path)

    args = parser.parse_args()

    parsed_genes = load_json(args.parsed_json)
    uniprot_to_kegg = load_json(args.mapping_json)
    go_by_uniprot = load_json(args.go_json)
    kegg_by_id = load_json(args.kegg_json)

    merged = merge_annotations(
        parsed_genes=parsed_genes,
        uniprot_to_kegg=uniprot_to_kegg,
        go_by_uniprot=go_by_uniprot,
        kegg_by_id=kegg_by_id,
    )

    if not isinstance(parsed_genes, list):
        sys.exit("[ERROR] parsed_json must be a list of gene objects")

    with args.merged_json.open("w") as f:
        json.dump(merged, f, indent=2)

    write_versions()


if __name__ == "__main__":
    main()
