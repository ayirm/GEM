import requests, time, sys, json, os

from Bio import SeqIO
from bioservices import QuickGO
from openpyxl import Workbook

def map_uniprot_to_kegg_via_rest(uniprot_ids, chunk_size=500, max_retries=3):
    """
    Maps UniProt IDs to KEGG IDs using the UniProt REST API directly.
    Includes a retry mechanism and correctly handles paginated results.
    """
    mapping_dict = {}
    
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i:i+chunk_size]
        if not chunk:
            continue

        print(f"\nSubmitting chunk {i//chunk_size + 1}/{(len(uniprot_ids) - 1)//chunk_size + 1} for mapping...")
        
        for attempt in range(max_retries):
            job_id = None
            try:
                payload = {
                    "from": "UniProtKB_AC-ID",
                    "to": "KEGG",
                    "ids": ",".join(chunk),
                }
                
                submit_response = requests.post("https://rest.uniprot.org/idmapping/run", data=payload)
                submit_response.raise_for_status()
                job_id = submit_response.json()["jobId"]
                
                # --- Polling Loop ---
                spinner_idx = 0
                job_finished = False
                while not job_finished:
                    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
                    status_response = requests.get(status_url)

                    if status_response.status_code == 404:
                        sys.stdout.write(f"\rPolling job status... (job not ready, waiting) {'/-\\|'[spinner_idx % 4]}")
                        sys.stdout.flush()
                        spinner_idx += 1
                        time.sleep(5)
                        continue

                    status_response.raise_for_status()
                    status_data = status_response.json()

                    # Check for completion signals
                    if "results" in status_data or status_data.get("jobStatus") == "FINISHED":
                        job_finished = True
                    elif status_data.get("jobStatus") in ["RUNNING", "QUEUED"]:
                        sys.stdout.write(f"\rPolling job status... {'/-\\|'[spinner_idx % 4]}")
                        sys.stdout.flush()
                        spinner_idx += 1
                        time.sleep(2)
                    else:
                        print(f"\nJob has an unexpected status. Full response below:")
                        print(json.dumps(status_data, indent=2))
                        raise requests.exceptions.RequestException("Unexpected job status")
                
                # --- Result Fetching with Pagination ---
                if job_finished:
                    sys.stdout.write("\rPolling job status... Done. Fetching all results.\n")
                    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
                    
                    page_count = 1
                    while results_url:
                        results_response = requests.get(results_url)
                        results_response.raise_for_status()

                        for item in results_response.json().get("results", []):
                            uniprot_id = item['from']
                            if uniprot_id not in mapping_dict:
                                mapping_dict[uniprot_id] = item['to']
                        
                        # Check for the 'next' page link in the response headers
                        if 'next' in results_response.links:
                            results_url = results_response.links['next']['url']
                            sys.stdout.write(f"\rFetching page {page_count+1}...")
                            sys.stdout.flush()
                            page_count += 1
                        else:
                            results_url = None # No more pages

                # If we successfully processed the job, exit the retry loop
                break 

            except requests.exceptions.RequestException as e:
                print(f"\nAttempt {attempt + 1}/{max_retries} failed for chunk starting at index {i}: {e}")
                if attempt + 1 == max_retries:
                    print(f"--- Giving up on this chunk after {max_retries} attempts. ---")
                else:
                    delay = 2 ** attempt
                    print(f"--- Retrying in {delay} seconds... ---")
                    time.sleep(delay)

    return mapping_dict

def parse_kegg_entry(entry_text):
    """
    Parses the text from a KEGG 'get' command to find pathways,
    the KO number, and the BRITE hierarchy.
    """
    pathways = []
    brite_hierarchy = []
    ko_number = None
    
    lines = entry_text.strip().split('\n')
    
    current_section = None
    
    for line in lines:
        if not line.startswith(" "):
            parts = line.split(None, 1)
            if parts:
                keyword = parts[0]
                if keyword in ["ENTRY", "NAME", "SYMBOL", "ORTHOLOGY", "ORGANISM", 
                               "PATHWAY", "MODULE", "BRITE", "POSITION", "MOTIF", 
                               "DBLINKS", "STRUCTURE", "AASEQ", "NTSEQ"]:
                    current_section = keyword
                else:
                    current_section = None 
        
        content = line[12:].strip()
        
        if current_section == "PATHWAY":
            if content:
                pathways.append(content)
        
        elif current_section == "ORTHOLOGY":
            if "[KO:" in content:
                try:
                    ko_num = content.split("[KO:")[-1].split("]")[0].strip()
                    if ko_num: 
                        ko_number = ko_num
                except Exception:
                    pass 
        
        elif current_section == "BRITE":
            if content:
                brite_hierarchy.append(content)

    return pathways, ko_number, brite_hierarchy

def parse_ko_entry(entry_text):
    """Parses the text from a KEGG KO entry to find reactions."""
    reactions = []
    lines = entry_text.strip().split('\n')
    in_reaction_section = False
    for line in lines:
        if line.startswith("REACTION"):
            in_reaction_section = True
            reaction_part = line[12:].strip()
            reactions.extend([r.strip() for r in reaction_part.split(',')])
        elif in_reaction_section and line.startswith("            "): # Continuation line (12 spaces)
            reaction_part = line[12:].strip()
            reactions.extend([r.strip() for r in reaction_part.split(',')])
        elif in_reaction_section and not line.startswith("            "):
            in_reaction_section = False
    return [r for r in reactions if r]

def uniprot_search(ann_path:str, prkName="prokkaAnnotes", cache_file="go_cache.json"):
    """
    Uses prokka results to get gene_id, protein_name, EC_number and Gene Onthology via UniProt
    """
    gbf_file = os.path.join(ann_path, f"{prkName}.gbf")
    if not os.path.exists(gbf_file):
        gbf_file = os.path.join(ann_path, f"{prkName}.gbk")

    total_CDS = 0
    try:
        total_CDS = sum(1 for record in SeqIO.parse(gbf_file, "genbank") for f in record.features if f.type == "CDS")
    except FileNotFoundError:
        print(f"Error: Could not find GBK/GBF file at {gbf_file}")
        return []

    cds_values = []

    for i , record in enumerate(SeqIO.parse(gbf_file, "genbank"), 1):
        for feature in record.features:
            if feature.type != "CDS":
                continue # Use continue instead of pass

            qualifiers = feature.qualifiers
            gene = qualifiers.get("gene", [None])[0]
            product = qualifiers.get("product", [None])[0]
            ec_number = qualifiers.get("EC_number", [None])[0]

            uni_inference = None
            for inf in qualifiers.get("inference", []):
                if "UniProtKB" in inf:
                    uni_inference = inf.split("UniProtKB:")[-1]
                    break

            if uni_inference:
                cds_values.append({
                    "gene": gene,
                    "inference": uni_inference,
                    "product": product,
                    "EC_number": ec_number
                })
            
    sys.stdout.write(f"Parsing {len(cds_values)} valid UniProt ids out of {total_CDS}")
    print(f"\nTotal with UniProt IDs: {len(cds_values)}")
    
    # GO Term Search Part
    go = QuickGO()
    uniprot_ids = [entry["inference"] for entry in cds_values if entry["inference"]]
    go_data = {}
    
    # Caching is now relative to the execution directory (the process 'work' dir)
    go_cache_path = cache_file 

    if cache_file and os.path.exists(go_cache_path):
        print(f"\nFound GO cache: '{go_cache_path}'. Using that")
        try:
            with open(go_cache_path, 'r') as f:
                go_data = json.load(f)
            print("Loaded GO data from cache.")
        except Exception as e:
            print(f"Error loading cache '{go_cache_path}': {e}. Re-fetching all GO terms.")
            go_data = {} 

    if not go_data:
        if cache_file:
            print(f"\nCache not found or empty. Fetching {len(uniprot_ids)} GO terms from QuickGO...")
        
        go = QuickGO()
        
        for i, uid in enumerate(uniprot_ids, 1):
            try:
                response = go.Annotation(
                    geneProductId= uid,
                    includeFields= "goName",
                )
                go_data[uid] = response.get("results", [])
            except Exception as e:
                print(f"\nError fetching GO for {uid}: {e}")
                go_data[uid] = [] 

            sys.stdout.write(f"\r Searching {i}th GO term out of {len(uniprot_ids)}")
            sys.stdout.flush()
        
        print("\nFinished fetching GO terms. Yey :)")

        if cache_file:
            try:
                with open(go_cache_path, 'w') as f:
                    json.dump(go_data, f, indent=2)
                print(f"Saved GO data to cache: '{go_cache_path}'")
            except Exception as e:
                print(f"Error saving cache: {e}")
    
    # Adding the values back to the dictornary
    for d in cds_values:
        uid = d["inference"]
        go_info = go_data.get(uid, [])
        d["GO"] = [f"{g['goId']}: {g['goName']}" for g in go_info]
    
    return cds_values

def kegg_search(cds_values, ann_path, gene_cache_file="kegg_gene_cache.json", ko_cache_file="kegg_ko_cache.json"):
    """
    Function that makes 3 request per second to https://https://rest.kegg.jp/conv/genes/ convert Uniprot ids into KEGG Ids
    Then uses the id to get KO id's, pathways and reactions
    """
    # Caching is now relative to the execution directory (the process 'work' dir)
    gene_cache_path = gene_cache_file
    ko_cache_path = ko_cache_file

    gene_cache = {}
    if os.path.exists(gene_cache_path):
        print(f"\nFound KEGG Gene cache: '{gene_cache_path}'. Using that.")
        with open(gene_cache_path, 'r') as f:
            gene_cache = json.load(f)

    ko_cache = {}
    if os.path.exists(ko_cache_path):
        print(f"Found KEGG KO cache: '{ko_cache_path}'. Using that.")
        with open(ko_cache_path, 'r') as f:
            ko_cache = json.load(f)
    
    uniprot_ids = [entry["inference"] for entry in cds_values]
    unique_uniprot_ids = list(set(uniprot_ids))

    print(f"Doing UniProt to KEGG transform for {len(unique_uniprot_ids)} unique IDs")
    uniprot_to_kegg = map_uniprot_to_kegg_via_rest(unique_uniprot_ids)
    print(f"\nSuccessfully mapped {len(uniprot_to_kegg)} UniProt IDs.")

    total_entries = len(cds_values)
    for i, entry in enumerate(cds_values, 1):
        uniprot_id = entry["inference"]
        kegg_id = uniprot_to_kegg.get(uniprot_id)
        
        entry["kegg_id"] = kegg_id if kegg_id else "N/A"
        entry["pathways"] = []
        entry["reactions"] = []
        entry["brite_hierarchy"] = []

        if not kegg_id:
            continue

        try:
            if kegg_id not in gene_cache:
                time.sleep(0.35) # (max 3 requests/sec)
                response = requests.get(f"https://rest.kegg.jp/get/{kegg_id}")
                response.raise_for_status()
                gene_cache[kegg_id] = response.text
            
            gene_entry_text = gene_cache[kegg_id]
            pathways, ko_number, brite_data = parse_kegg_entry(gene_entry_text)
            
            entry["pathways"] = pathways if pathways else ["N/A"]
            entry["brite_hierarchy"] = brite_data if brite_data else ["N/A"]

            if ko_number:
                if ko_number not in ko_cache:
                    time.sleep(0.35) # (max 3 requests/sec)
                    response = requests.get(f"https://rest.kegg.jp/get/{ko_number}")
                    response.raise_for_status()
                    ko_cache[ko_number] = response.text
                
                ko_entry_text = ko_cache[ko_number]
                reactions = parse_ko_entry(ko_entry_text)
                entry["reactions"] = reactions if reactions else ["N/A"]
            else:
                entry["reactions"] = ["N/A (No KO number found)"]

        except requests.exceptions.HTTPError as e:
            print(f"\nError fetching data for {kegg_id} (from UniProt {uniprot_id}): {e}")
            entry["pathways"] = [f"Error: {e}"]
            entry["reactions"] = [f"Error: {e}"]
            entry["brite_hierarchy"] = [f"Error: {e}"]
        
        sys.stdout.write(f"\rProcessing {i}/{total_entries} entries...")
        sys.stdout.flush()
    
    print("\n\nFinished KEGG search.")
    try:
        with open(gene_cache_path, 'w') as f:
            json.dump(gene_cache, f, indent=2)
        print(f"Saved KEGG Gene data to cache: '{gene_cache_path}'")
    except Exception as e:
        print(f"Error saving gene cache: {e}")

    try:
        with open(ko_cache_path, 'w') as f:
            json.dump(ko_cache, f, indent=2)
        print(f"Saved KEGG KO data to cache: '{ko_cache_path}'")
    except Exception as e:
        print(f"Error saving KO cache: {e}")

    return cds_values

def excel_creation(enriched_cds_values, output_file="annotation_results.xlsx"):
    """
    Writes the final, enriched CDS annotations into a single Excel sheet.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Annotations"

    if not enriched_cds_values:
        print("Warning: No data to write to Excel file.")
        ws.append(["No data found"])
        wb.save(output_file)
        return

    headers = list(enriched_cds_values[0].keys())
    ws.append(headers)

    for entry in enriched_cds_values:
        row = [
            "\n".join(map(str, entry.get(h, []))) if isinstance(entry.get(h), list) else entry.get(h, "")
            for h in headers
        ]
        ws.append(row)

    wb.save(output_file)
    print(f"✅ Results written to '{output_file}'")