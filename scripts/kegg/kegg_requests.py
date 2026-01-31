#!/usr/bin/env python3

import time
import requests
import sys
import json
import argparse
import re
import sqlite3
import os

from datetime import datetime, timezone

BASE_URL = "https://rest.kegg.jp"

def init_db(db_path):
    """Initialize (or open) the SQLite DB.

    If a relative path is provided while running inside a Nextflow `work/`
    directory, this function will prefer a persistent cache location so the
    database survives Nextflow's work cleanup. It tries the following in
    order and picks the first writable candidate:
      - project root `./.cache/iumobg-model_creation`
      - `$XDG_CACHE_HOME/iumobg-model_creation`
      - `~/.cache/iumobg-model_creation`
      - `/tmp/iumobg-model_creation`
    """

    db_path = os.path.expanduser(db_path)

    if not os.path.isabs(db_path):
        cwd = os.getcwd()
        project_root = None
        sep_work = os.path.sep + 'work' + os.path.sep
        if sep_work in cwd:
            project_root = cwd.split(sep_work)[0]

        candidates = []
        if project_root:
            candidates.append(os.path.join(project_root, '.cache', 'iumobg-model_creation'))
        xdg = os.environ.get('XDG_CACHE_HOME')
        if xdg:
            candidates.append(os.path.join(xdg, 'iumobg-model_creation'))
        home = os.path.expanduser('~')
        if home and home not in ('/', ''):
            candidates.append(os.path.join(home, '.cache', 'iumobg-model_creation'))
        candidates.append(os.path.join('/tmp', 'iumobg-model_creation'))

        cache_base = None
        for cand in candidates:
            try:
                os.makedirs(cand, exist_ok=True)
                # verify writable
                testfile = os.path.join(cand, '.write_test')
                with open(testfile, 'w') as tf:
                    tf.write('ok')
                os.remove(testfile)
                cache_base = cand
                break
            except Exception:
                continue

        if cache_base:
            db_path = os.path.join(cache_base, db_path)
        else:
            db_path = os.path.abspath(db_path)
    else:
        db_path = os.path.abspath(db_path)

    db_dir = os.path.dirname(db_path)
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)

    print(f"    [DB] using database at {db_path}", file=sys.stderr)

    conn = sqlite3.connect(db_path, timeout=30)
    cur = conn.cursor()
    try:
        cur.execute("PRAGMA journal_mode=WAL;")
        cur.execute("PRAGMA synchronous=NORMAL;")
    except sqlite3.DatabaseError:
        pass

    cur.execute("""
    CREATE TABLE IF NOT EXISTS kegg_results (
        kegg_id TEXT PRIMARY KEY,
        payload TEXT NOT NULL
    )
    """)

    conn.commit()
    return conn

def kegg_is_complete(cur, kegg_id):
    # Read the payload as text and parse in Python to avoid relying on
    # SQLite's JSON1 extension (may not be available in some environments).
    cur.execute(
        "SELECT payload FROM kegg_results WHERE kegg_id = ?",
        (kegg_id,)
    )
    row = cur.fetchone()
    if not row:
        return False
    try:
        payload = json.loads(row[0])
    except Exception:
        return False
    return payload.get("status") == "complete"


def save_kegg_result(cur, conn, kegg_id, result):
    cur.execute(
        "INSERT OR REPLACE INTO kegg_results VALUES (?, ?)",
        (kegg_id, json.dumps(result))
    )
    conn.commit()

def kegg_request(endpoint, rate_limit=0.35, retries=5, backoff=1.5, last_time=[0]):
    url = f"{BASE_URL}{endpoint}"

    for attempt in range(1, retries + 1):
        elapsed = time.time() - last_time[0]
        if elapsed < rate_limit:
            time.sleep(rate_limit - elapsed)

        try:
            print(f"      → GET {endpoint} (attempt {attempt})", file=sys.stderr)
            r = requests.get(url, timeout=(5, 20))  # CONNECT, READ
            r.raise_for_status()
            last_time[0] = time.time()
            return r.text

        except requests.RequestException as e:
            print(
                f"      ! ERROR {endpoint} attempt {attempt}: {e}",
                file=sys.stderr
            )

            if attempt == retries:
                return None
            time.sleep(backoff ** attempt)

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
            "equation": "N/A(ERROR)",
            "cmpdID": [],
        }

    readableEQ = " ".join(parse_page_entries(reactionPage, "DEFINITION"))
    cmpdEQ = " ".join(parse_page_entries(reactionPage, "EQUATION"))

    # Extract compound IDs like C00002
    cmpd_ids = re.findall(r"C\d{5}", cmpdEQ)

    return {
        "equation": readableEQ if readableEQ else "N/A",
        "cmpdID": list(set(cmpd_ids)),
    }


def populate_kegg_compound(compound_id):
    """
    Fetch KEGG compound entry (Cxxxxx) and extract model-relevant fields.
    FORMULA     C4H9NO3
    EXACT_MASS  119.0582
    MOL_WEIGHT  119.12
    DBLINKS     CAS: 672-15-1
            PubChem: 3561
            ChEBI: 15699
            KNApSAcK: C00001366
            PDB-CCD: HSE
            NIKKAJI: J9.199E
    """
    page = kegg_request(f"/get/{compound_id}")
    if not page:
        return None

    formula = None
    exact_mass = None
    mol_weight = None
    dblinks = {}

    for line in page.splitlines():
        if line.startswith("FORMULA"):
            formula = line.replace("FORMULA", "").strip()

        elif line.startswith("EXACT_MASS"):
            exact_mass = line.replace("EXACT_MASS", "").strip()

        elif line.startswith("MOL_WEIGHT"):
            mol_weight = line.replace("MOL_WEIGHT", "").strip()

        elif line.startswith("DBLINKS"):
            # DBLINKS can span multiple indented lines
            current = line.replace("DBLINKS", "").strip()
            if current:
                key, val = current.split(":", 1)
                dblinks[key.strip()] = val.strip()

        elif line.startswith(" " * 12) and ":" in line and dblinks:
            key, val = line.strip().split(":", 1)
            dblinks[key.strip()] = val.strip()

    return {
        "formula": formula,
        "exact_mass": exact_mass,
        "molecular_weight": mol_weight,
        "dblinks": dblinks
    }


def fetch_kegg_entry(kegg_id):
    result = {
        "status": "incomplete",
        "KO_Terms": [],
        "pathways": [],
        "BRITE": "N/A",
        "BRITE_id": [],
        "reactions": {},
        "compounds": {},
        "compound_metadata": {}
    }

    print(f"    [ENTRY] fetching main page", file=sys.stderr)
    page = kegg_request(f"/get/{kegg_id}")
    if not page:
        print(f"    [ENTRY] FAILED main page", file=sys.stderr)
        return result

    print(f"    [KO] parsing ORTHOLOGY", file=sys.stderr)
    result["KO_Terms"] = find_ko_numbers(page)

    print(f"    [BRITE] parsing BRITE", file=sys.stderr)
    result["BRITE"] = find_full_brite(page)
    result["BRITE_id"] = find_brite_numbers(page)

    print(f"    [PATHWAY] parsing PATHWAY", file=sys.stderr)
    pathways = find_pathways(page)
    result["pathways"] = pathways
    print(f"    [PATHWAY] found {len(pathways)}", file=sys.stderr)

    map_cache = {}
    reactions = set()

    for p in pathways:
        m = re.search(r"(\d{5})$", p)
        if not m:
            continue

        map_id = f"map{m.group(1)}"

        if map_id not in map_cache:
            print(f"    [MAP] {map_id} → reactions", file=sys.stderr)
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

    print(f"    [REACTIONS] total {len(reactions)}", file=sys.stderr)

    compound_cache = {}

    for rn in reactions:
        print(f"      [RXN] {rn}", file=sys.stderr)
        rpage = kegg_request(f"/get/{rn}")
        rxn_data = find_compounds(rpage)
        result["compounds"][rn] = rxn_data

        for cid in rxn_data.get("cmpdID", []):
            if cid not in compound_cache:
                print(f"        [CMPD] {cid}", file=sys.stderr)
                compound_cache[cid] = populate_kegg_compound(cid)

    result["compound_metadata"] = compound_cache
    result["status"] = "complete"

    print(f"    [DONE] {kegg_id}", file=sys.stderr)
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

def run_kegg_annotation(mapping_json, out_json, db_path):
    with open(mapping_json) as f:
        mapping = json.load(f)

    unique_kegg_ids = sorted(set(mapping.values()))

    conn = init_db(db_path)
    cur = conn.cursor()

    total = len(unique_kegg_ids)
    for i, kegg_id in enumerate(unique_kegg_ids, 1):
        if kegg_is_complete(cur, kegg_id):
            continue

        print(f"[{i}/{total}] Annotating {kegg_id}", file=sys.stderr)

        result = fetch_kegg_entry(kegg_id)
        save_kegg_result(cur, conn, kegg_id, result)

    # EXPORT SQLITE → JSON (single pass)
    cur.execute("SELECT kegg_id, payload FROM kegg_results")
    rows = cur.fetchall()

    results = {
        kegg_id: json.loads(payload)
        for kegg_id, payload in rows
    }

    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping_json", required=True)
    parser.add_argument("--out_json", required=True)
    parser.add_argument("--db", required=True)
    args = parser.parse_args()

    run_kegg_annotation(
        mapping_json=args.mapping_json,
        out_json=args.out_json,
        db_path=args.db,
    )

    write_versions()
