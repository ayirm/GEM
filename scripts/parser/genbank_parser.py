#!/usr/bin/env python3

import argparse
import json
import sys
import tempfile
import re

from datetime import datetime, timezone
from Bio import SeqIO

def sanitize_genbank(in_path):
    """
    Fix non-standard LOCUS lines so Biopython can parse them.
    """
    fixed_lines = []
    with open(in_path, "r") as f:
        for line in f:
            if line.startswith("LOCUS"):
                # extract locus name (first non-space token after LOCUS)
                m_locus = re.match(r"LOCUS\s+(\S+)", line)
                locus = m_locus.group(1) if m_locus else "UNKNOWN"
                # try to find a bp length (digits before 'bp'), otherwise first integer in the line
                m_len = re.search(r"(\d+)\s*bp", line)
                if not m_len:
                    m_len = re.search(r"(\d+)", line)
                length = m_len.group(1) if m_len else "0"
                line = f"LOCUS       {locus:<16} {length} bp    DNA     linear   BCT 01-JAN-2000\n"
            fixed_lines.append(line)

    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
    tmp.writelines(fixed_lines)
    tmp.close()
    return tmp.name

def parse_genbank(gb_file):
    """
        Gets the following parts from the provided gff or gbk file
        - gene="thrA"
        - EC_number="2.7.1.39"
        - uni_idserence="similar to AA sequence:UniProtKB:P00561"
        - product="Bifunctional aspartokinase/homoserine dehydrogenase 1"

        Only the cds records that have UniProt id's are saved. Other areas can be empty but UniProtID will always be filled
    """
    gb_file = sanitize_genbank(gb_file)

    cds_values = []
    total_cds = 0

    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            total_cds += 1
            qualifiers = feature.qualifiers

            gene = qualifiers.get("gene", [None])[0]
            protein = qualifiers.get("product", [None])[0]
            ec_number = qualifiers.get("EC_number", [None])[0]

            uniprot_id = None
            for uni_ids in qualifiers.get("inference", []):
                if "UniProtKB:" in uni_ids:
                    uniprot_id = uni_ids.split("UniProtKB:")[-1]
                    break

            if uniprot_id:
                cds_values.append({
                    "Gene": gene,
                    "Protein": protein,
                    "EC_Number": ec_number,
                    "UniProt_ID": uniprot_id
                })

    sys.stderr.write(f"Total cds: {total_cds} \n Parsing {len(cds_values)}th value")
    return cds_values

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
        description="Parse GenBank file and extract UniProt-linked CDS entries"
    )
    parser.add_argument(
        "--gb_file",
        required=True,
        help="Input GenBank file (.gbk or .gff)"
    )
    parser.add_argument(
        "--parser_json",
        required=True,
        help="Output JSON file"
    )

    args = parser.parse_args()

    cds_values_data = parse_genbank(args.gb_file)

    with open(args.parser_json, "w") as out:
        json.dump(cds_values_data, out, indent=2)

    write_versions()

if __name__ == "__main__":
    main()
