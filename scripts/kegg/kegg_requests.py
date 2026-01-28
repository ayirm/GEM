#!/usr/bin/env python3

import time
import requests
import sys
import json
import tqdm
import argparse
import re
import os

from datetime import datetime, timezone

BASE_URL = "https://rest.kegg.jp"

def kegg_request(endpoint, rate_limit=0.35, last_time=[0]):
    elapsed = time.time() - last_time[0]
    if elapsed < rate_limit:
        time.sleep(rate_limit - elapsed)
    # Basically checks if the time since the last request is smaller than 0.35
    # if it is smaller, sleeps for that amount
    # ex: time elapsed is 0,25 then it sleeps for 0,10. 0,20 for 0,15 etc.

    url = f"{BASE_URL}{endpoint}"
    try:
        r = requests.get(url)
        r.raise_for_status()
        last_time[0] = time.time()
        return r.text
    except requests.RequestException as e:
        print(f"[KEGG ERROR] {endpoint}: {e}", file=sys.stderr)
        return None

def parse_page_entries(page, section):
    """
        Searches the entry inside the provided page, returns the part of the searched section \n
        page is the full page of the kegg request
        """
    if not page:
        return []

    lines = []
    active = False
    # All KEGG entries align section names. Content starts 12 spaces in.
    indent = " " * 12

    for line in page.splitlines():
        if line.startswith(section):
            active = True
            lines.append(line[len(section):].strip())
        elif active and line.startswith(indent):
            lines.append(line.strip())
        elif active:
            break

    return lines

def find_ko_numbers(pageResult):
    """
    ORTHOLOGY   K12524  bifunctional aspartokinase / homoserine dehydrogenase 1 [EC:2.7.2.4 1.1.1.3]
    ^ Finds this in the page and gets the K12524 part.
    """
    koLines = parse_page_entries(pageResult, "ORTHOLOGY")
    koID = []

    # Regex to find one or more K numbers
    koRegex = re.compile(r'(K\d{5,})')

    for line in koLines:
        matches = koRegex.findall(line)
        for ko in matches:
            koID.append(f"KO:{ko}")

    return koID

def find_full_brite(pageResult):
    """
    BRITE       KEGG Orthology (KO) [BR:ecj00001]
         09100 Metabolism
          09105 Amino acid metabolism
           00260 Glycine, serine and threonine metabolism
                JW0001 (thrA)
           00270 Cysteine and methionine metabolism
            JW0001 (thrA)
               00300 Lysine biosynthesis
            JW0001 (thrA)
              09110 Biosynthesis of other secondary metabolites
           00261 Monobactam biosynthesis
                JW0001 (thrA)
        Enzymes [BR:ecj01000]
         1. Oxidoreductases
          1.1  Acting on the CH-OH group of donors
           1.1.1  With NAD+ or NADP+ as acceptor
            1.1.1.3  homoserine dehydrogenase
             JW0001 (thrA)
         2. Transferases
          2.7  Transferring phosphorus-containing groups
           2.7.2  Phosphotransferases with a carboxy group as acceptor
            2.7.2.4  aspartate kinase
             JW0001 (thrA)
    ^ Finds this then just returns it
    """
    briteResults = parse_page_entries(pageResult, "BRITE")
    if not briteResults:
        return "No BRITE section"
        
    return "\n".join(briteResults)

def find_brite_numbers(pageResult):
    briteResults = parse_page_entries(pageResult, "BRITE")
    if not briteResults:
        return []
        
    briteID = []
    briteRegex = re.compile(r"\[(BR:[a-z0-9]+)\]")

    for line in briteResults:
        matches = briteRegex.findall(line)
        for brID in matches:
            briteID.append(brID)

    return briteID

def find_pathways(pageResult):
    """
    PATHWAY     ecj00260  Glycine, serine and threonine metabolism
        ecj00261  Monobactam biosynthesis
        ecj00270  Cysteine and methionine metabolism
        ecj00300  Lysine biosynthesis
        ecj01100  Metabolic pathways
        ecj01110  Biosynthesis of secondary metabolites
        ecj01120  Microbial metabolism in diverse environments
        ecj01230  Biosynthesis of amino acids
    ^ Finds this and then seperates the ecj... parts    

    Then returns the ecj part. No convertion to map happens here
    """
    pathwayResults = parse_page_entries(pageResult, "PATHWAY")
    pathwayID = []
    pathwayRegex = re.compile(r'([a-z]{3,5}\d{5})')

    for line in pathwayResults:
        # Find the first ID on the line
        match = pathwayRegex.match(line)
        if match:
            pathwayID.append(f"path:{match.group(1)}") 
    # Add 'path:' prefix for consistency, normally just map00260 works but other ones have br:, ko: etc.
    # It was purely an aesthetic choice

    return pathwayID

def find_compounds(reactionPage):
    """
    DEFINITION  ATP + L-Aspartate <=> ADP + 4-Phospho-L-aspartate
    EQUATION    C00002 + C00049 <=> C00008 + C03082
    parts

    Then sends them back as equation(DEFINITION) and cmpd_id(EQUATION)
    """

    if not reactionPage:
        return {
            "equation" : "N/A(ERROR)",
            "cmpdID" : "N/A(ERROR)"
        }
        
    readableEQ = " ".join(parse_page_entries(reactionPage, "DEFINITION"))
    cmpdEQ = " ".join(parse_page_entries(reactionPage, "EQUATION"))

    return {
        "equation" : readableEQ if readableEQ else "N/A",
        "cmpdID" : cmpdEQ if cmpdEQ else "N/A"
    }

def fetch_kegg_entry(kegg_id):
    result = {
        "KO_Terms": [],
        "pathways": [],
        "BRITE": "N/A",
        "BRITE_id": [],
        "reactions": {},
        "compounds": {}
    }

    page = kegg_request(f"/get/{kegg_id}")
    if not page:
        return result

    result["KO_Terms"] = find_ko_numbers(page)
    result["BRITE"] = find_full_brite(page)
    result["BRITE_id"] = find_brite_numbers(page)
    pathways = find_pathways(page)
    result["pathways"] = pathways

    map_cache = {}
    reactions = set()

    for p in pathways:
        m = re.search(r"(\d{5})$", p)
        if not m:
            continue

        map_id = f"map{m.group(1)}"

        if map_id not in map_cache:
            links = kegg_request(f"/link/reaction/{map_id}")
            if links:
                map_cache[map_id] = [
                    parts[1]
                    for line in links.splitlines()
                    if (parts := line.split("\t")) and len(parts) == 2
                ]
            else:
                map_cache[map_id] = []

        rxns = map_cache[map_id]
        result["reactions"][p] = rxns
        reactions.update(rxns)

    for rn in reactions:
        rpage = kegg_request(f"/get/{rn}")
        result["compounds"][rn] = find_compounds(rpage)

    return result

def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)

def run_kegg_annotation(mapping_json, out_json, cache_json=None):
    with open(mapping_json) as f:
        mapping = json.load(f)

    if cache_json and os.path.exists(cache_json):
        with open(cache_json) as f:
            results = json.load(f)
    else:
        results = {}

    unique_kegg_ids = set(mapping.values())

    for kegg_id in tqdm(unique_kegg_ids, desc="KEGG annotation"):
        if kegg_id in results:
            continue
        results[kegg_id] = fetch_kegg_entry(kegg_id)

    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping_json", required=True)
    parser.add_argument("--out_json", required=True)
    args = parser.parse_args()

    run_kegg_annotation(
        mapping_json=args.mapping_json,
        out_json=args.out_json,
    )

    write_versions()
