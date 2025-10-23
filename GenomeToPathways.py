import subprocess, os, sys, json, requests, time

from Bio import SeqIO
from openpyxl import Workbook
from bioservices import *

def run_cmd(cmd:str, cmd_name="No Name"):
    """Helper for subprocess command, helps to remove repetitive try-except blocks"""
    print(f"---Running {cmd_name}---\n")
    print(f"Full command: \n {''.join(cmd)}")

    try:
        process = subprocess.run(
            cmd,
            capture_output=True,
            check=True,
            shell = True, # For allowing pipes(|) and redirection arrows(>)
        )
        print(f"---Run of {cmd_name} was successfull--- \n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"---Error happened with {cmd_name}--- \n")
        print(f"Error return: \n {e.returncode}\n")
        print(f"Stdout: \n {e.stdout} \n")
        print(f"Stderr: \n {e.stderr} \n")

def check_existing_file(file_path: str, step_name: str) -> bool:
    """
    Checks if a file already exists and asks the user if it should be reused.
    Returns True if the step should be skipped.
    """
    if os.path.exists(file_path):
        resp = input(f"Found existing {step_name} at {file_path}. Skip this step and reuse it? [Y/n]: ")
        if resp.strip().lower() in ["", "y", "yes"]:
            print(f"Skipping {step_name}, using existing file.\n")
            return True
    return False

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
        # Check for a new section start
        if not line.startswith(" "):
            # Find the section keyword (first word)
            # Use split(None, 1) to split only on the first block of whitespace
            parts = line.split(None, 1)
            if parts:
                keyword = parts[0]
                # Check against a list of known KEGG keywords
                if keyword in ["ENTRY", "NAME", "SYMBOL", "ORTHOLOGY", "ORGANISM", 
                               "PATHWAY", "MODULE", "BRITE", "POSITION", "MOTIF", 
                               "DBLINKS", "STRUCTURE", "AASEQ", "NTSEQ"]:
                    current_section = keyword
                else:
                    # If it's not a recognized keyword, it's not a new section
                    current_section = None 
        
        # Process based on the current section
        # The content of each section starts at column 13 (index 12)
        content = line[12:].strip()
        
        if current_section == "PATHWAY":
            if content:
                # content is e.g., "ecj00260  Glycine, serine and threonine metabolism"
                pathways.append(content)
        
        elif current_section == "ORTHOLOGY":
            if "[KO:" in content:
                try:
                    # Extract the KO number
                    ko_num = content.split("[KO:")[-1].split("]")[0].strip()
                    if ko_num: # Only assign if we found one
                        ko_number = ko_num
                except Exception:
                    pass # Ignore parsing errors, keep ko_number as None or as last found
        
        elif current_section == "BRITE":
            if content:
                # content is e.g., "KEGG Orthology (KO) [BR:ecj00001]"
                # or "09100 Metabolism"
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
            # Extract reactions from the first line
            reaction_part = line[12:].strip()
            reactions.extend([r.strip() for r in reaction_part.split(',')])
        elif in_reaction_section and line.startswith("            "): # Continuation line (12 spaces)
             # Extract reactions from continuation lines
            reaction_part = line[12:].strip()
            reactions.extend([r.strip() for r in reaction_part.split(',')])
        elif in_reaction_section and not line.startswith("            "):
            # No longer in the reaction section if indentation stops
            in_reaction_section = False
    # Filter out empty strings that might result from splitting
    return [r for r in reactions if r]

# --- Main Codes for Alignment--->Annotation--->Protein search in UniProt ---> KEGG pathway search
# !TODO: Add quality control steps beofre alignment

def run_alignement(ref_path:str, fastq1_loc:str, fastq2_loc:str, out_loc:str, threads=4):
    """
    Indexes the reference genome then aligns the fastq reads to the reference as a sam file

    Code in terminal:
        bowtie2-build --threads 4 ref_path/file.fna ref_path/bowtie2 <-Name of the files that it will create
        bowite2 -x ref_path/bowtie -1 fastq1_loc/file.fastq -2 fastq2_loc/file.fasta -S sam_loc/file.sam
    """

    os.makedirs(out_loc, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(ref_path))[0]
    sam_file = os.path.join(out_loc, f"{base_name}.sam")

    ref_index_dir = os.path.join(out_loc, "ref_index")
    os.makedirs(ref_index_dir, exist_ok=True)
    ref_index_prefix = os.path.join(ref_index_dir, base_name)

    run_cmd(f"bowtie2-build --threads {str(threads)} {ref_path} {ref_index_prefix}","Indexing Reference File") 

    run_cmd(f"bowtie2 -x {ref_index_prefix} -1 {fastq1_loc} -2 {fastq2_loc} -S {sam_file}","Aligning FastQ reads")

    return sam_file # This is the full path used for run_samtools()'s sam_path

def run_samtools(sam_path:str, threads:int, out_loc:str):
    """
    Runs the SAM to BAM conversion and indexing

    Code in terminal:
        samtools view -uS -o bam_loc/file.bam sam_loc/file.sam
        samtools sort -@ 4 -T temp_file.tmp.sort -o bam_loc/file_sorted.bam
        samtools index bam_loc/file_sorted.bam
    """
    os.makedirs(out_loc, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(sam_path))[0] #Basically gets the file name, strips the directory name if it's a directory and removes the . thingie
    bam_loc = os.path.join(out_loc, f"{base_name}.bam")
    sorted_loc = os.path.join(out_loc, f"{base_name}_sorted.bam")
    tmp_file = os.path.join(out_loc, f"{base_name}.tmp.sort")

    run_cmd(f"samtools view -uS -o {bam_loc} {sam_path}","SAM to BAM conversion")
    run_cmd(f"samtools sort -@ {str(threads)} -T {tmp_file} -o {sorted_loc} {bam_loc}", "BAM file sorting")
    run_cmd(f"samtools index {sorted_loc}","İndexing the sorted BAM file")

    return sorted_loc # Used in run_bcftools() as sorted_path

def run_bcftools(ref_path:str, sorted_path:str ,out_loc:str):
    """
    Collapses the reads according to reference and writes variant information to a VCF file. 
    Also compresses the VCF file with bgzip and indexes it with tabix

    Code in terminal:
        bcftools mpileup -f ref_loc/file.fa bam_loc/file_sorted.bam | bcftools call -mv -Ob -o vcf_loc/file.bcf
        bcftools convert -O v -o vcf_loc/file.vcf vcf_loc/file.bcf
        bgzip -c vcf_loc/file.vcf > gz_loc/file.vcf.gz
        tabix -p vcf gz_loc/file.vcf.gz
    """

    os.makedirs(out_loc, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(sorted_path))[0]
    bcf_loc = os.path.join(out_loc, f"{base_name}.bcf")
    vcf_loc = os.path.join(out_loc, f"{base_name}.vcf")
    gz_loc = os.path.join(out_loc, f"{base_name}.vcf.gz")

    run_cmd(f"bcftools mpileup -f {ref_path} {sorted_path} | bcftools call -mv -Ob -o {bcf_loc}", "Collapsing reads to a BCF file")
    run_cmd(f"bcftools convert -O v -o {vcf_loc} {bcf_loc}","Converting BCF to VCF")
    run_cmd(f"bgzip -c {vcf_loc} > {gz_loc}", "Compressing VCF")
    run_cmd(f"tabix -p vcf {gz_loc}", "Indexing compressed VCF")

    return gz_loc # Used in run_reference_assembly as gz_loc

def run_reference_assembly(ref_path:str, gz_loc:str, out_loc:str):
    """
    Creates a concensus sequence with bcftools. 

    Code in terminal:
        bcftools concensus -f ref_loc/file.fa gz_loc/file.vcf.gz > file.fasta
    """

    os.makedirs(out_loc, exist_ok=True)

    # base_name = os.path.splitext(os.path.basename(gz_loc))[0] #!IMPORTANT: base name here has the .vcf. I dont think it would cause an error though
    # Aaaand, it caused an error. Hurray!
    base_name = os.path.basename(gz_loc).replace(".vcf.gz","")
    cons_loc = os.path.join(out_loc, f"{base_name}.fasta")

    run_cmd(f"bcftools consensus -f {ref_path} {gz_loc} > {cons_loc}","Creating consensus sequence")

    return cons_loc # Used in prokka to annote

def run_annotation(cons_path:str, out_loc:str,prkName="prokkaAnnotes"):
    """
    Runs prokka to annote the reference assembled genome

    Code in terminal:
        prokka --outdir prk_loc --prefix genomeName file.fasta
    """

    os.makedirs(out_loc, exist_ok=True)

    run_cmd(f"prokka --outdir {out_loc} --prefix {prkName} --force {cons_path}","Annotation with Prokka")

    return prkName, out_loc

def uniprot_search(ann_path:str, prkName="prokkaAnnotes", cache_file="go_cache.json"):
    """
    Uses prokka results to get gene_id, protein_name, EC_number and Gene Onthology via UniProt
    Basically reads the prokka.gbf file and selects the CDS parts and then separetes hypothetical proteins that doesnt have any EC number or UniProt id
    After getting UniProt id, searches it to find Gene Onthologies \n
    Result is writted into a list \n

    Also GO data is cached because it takes toooo long!
    """

    gbf_file = os.path.join(ann_path, f"{prkName}.gbf")
    if not os.path.exists(gbf_file):
        gbf_file = os.path.join(ann_path, f"{prkName}.gbk")

    total_CDS = sum(1 for record in SeqIO.parse(gbf_file, "genbank") for f in record.features if f.type == "CDS")

    cds_values = []

    for i , record in enumerate(SeqIO.parse(gbf_file, "genbank"), 1):
        for feature in record.features:
            if feature.type != "CDS":
                pass

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

    print(f"Total with UniProt IDs: {len(cds_values)}")
    
    # GO Term Search Part
    go = QuickGO()

    uniprot_ids = [entry["inference"] for entry in cds_values if entry["inference"]]

    go_data = {}
    go_cache_path = os.path.join(ann_path, cache_file)

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
    Essentially the same thing as mapping() but doesn't crash the whole VSCode because it forgets to notify the server and makes VSCode think that it crashed
    Then uses the id to get KO id's, pathways and reactions
    """
    gene_cache_path = os.path.join(ann_path, gene_cache_file)
    ko_cache_path = os.path.join(ann_path, ko_cache_file)

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
        
        # Initialize with default empty values
        entry["kegg_id"] = kegg_id if kegg_id else "N/A"
        entry["pathways"] = []
        entry["reactions"] = []
        entry["brite_hierarchy"] = [] # <-- ADDED THIS

        if not kegg_id:
            continue

        try:
            if kegg_id not in gene_cache:
                time.sleep(0.35) # (max 3 requests/sec)
                response = requests.get(f"https://rest.kegg.jp/get/{kegg_id}")
                response.raise_for_status()
                gene_cache[kegg_id] = response.text
            
            gene_entry_text = gene_cache[kegg_id]
            # <-- MODIFIED CALL TO NEW PARSER -->
            pathways, ko_number, brite_data = parse_kegg_entry(gene_entry_text)
            
            entry["pathways"] = pathways if pathways else ["N/A"]
            entry["brite_hierarchy"] = brite_data if brite_data else ["N/A"] # <-- ADDED THIS

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
            entry["brite_hierarchy"] = [f"Error: {e}"] # <-- ADDED THIS
        
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

    # Write headers dynamically from the keys of the first entry
    headers = list(enriched_cds_values[0].keys())
    ws.append(headers)

    # Write the data rows
    for entry in enriched_cds_values:
        # Convert list values to comma-separated strings for cleaner Excel output
        row = [
            # Use newline (\n) as a separator for lists, which Excel can handle
            "\n".join(map(str, entry.get(h, []))) if isinstance(entry.get(h), list) else entry.get(h, "")
            for h in headers
        ]
        ws.append(row)

    wb.save(output_file)
    print(f"✅ Results written to '{output_file}'")

# We need to input two Fastq file name and one reference file
def run_pipeline(ref_path, fastq1, fastq2, out_dir="results", threads=4, prk_prefix="prokkaAnnotes", excel_file="annotation_results.xlsx"):
    os.makedirs(out_dir, exist_ok=True)

    # Folders in folders part
    align_and_sam = os.path.join(out_dir, "sam")
    bcf_and_gz = os.path.join(out_dir, "bcf")
    prokka_out = os.path.join(out_dir, "prk")

    # Names for stopping the repeptition
    base_name = os.path.splitext(os.path.basename(ref_path))[0]
    sam_file = os.path.join(align_and_sam, f"{base_name}.sam")
    sorted_bam = os.path.join(align_and_sam, f"{base_name}_sorted.bam")
    gz_vcf = os.path.join(bcf_and_gz, f"{base_name}.vcf.gz")
    cons_fasta = os.path.join(bcf_and_gz, f"{base_name}.fasta")
    prk_gbk = os.path.join(prokka_out, f"{prk_prefix}.gbf")

    # Checks for file
    if not check_existing_file(sam_file, "Alignemnt and Sam file creation"):
        sam_file = run_alignement(ref_path, fastq1, fastq2, out_loc=align_and_sam, threads=threads)

    if not check_existing_file(sorted_bam, "Sorting the Bam file"):
        sorted_bam = run_samtools(sam_file, threads=threads, out_loc=align_and_sam)

    if not check_existing_file(gz_vcf, "Collapsing Reads"):
        gz_vcf = run_bcftools(ref_path, sorted_bam, out_loc=bcf_and_gz)

    if not check_existing_file(cons_fasta, "Consensus File Creation"):
        cons_fasta = run_reference_assembly(ref_path, gz_vcf, out_loc=bcf_and_gz)

    if not os.path.exists(prk_gbk):
        prk_gbk = os.path.join(prokka_out, f"{prk_prefix}.gbk")
    if not check_existing_file(prk_gbk, "annotation (Prokka)"):
        prkName, prk_path = run_annotation(cons_fasta, out_loc=prokka_out, prkName=prk_prefix)
    else:
        prkName, prk_path = prk_prefix, prokka_out

    cds_values = uniprot_search(prkName=prk_prefix, ann_path=prk_path)
    
    if cds_values:
        print("\n--- Running KEGG Pathway and Reaction Search ---")
        enriched_results = kegg_search(cds_values, ann_path=prk_path)

        print("\n--- Creating Excel Report ---")
        excel_creation(enriched_results, output_file=os.path.join(out_dir, excel_file))

    print("✅ Pipeline finished successfully!")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Full genome-to-pathways pipeline")
    parser.add_argument("--ref", required=True, help="Path to reference genome FASTA")
    parser.add_argument("--fastq1", required=True, help="Path to FastQ R1")
    parser.add_argument("--fastq2", required=True, help="Path to FastQ R2")
    parser.add_argument("--out_dir", default="results", help="Output directory (default: results)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for alignment (default: 4)")
    parser.add_argument("--prk_prefix", default="prokkaAnnotes", help="Prefix for Prokka annotation (default: prokkaAnnotes)")
    parser.add_argument("--excel_file", default="annotation_results.xlsx", help="Excel output file name (default: annotation_results.xlsx)")

    args = parser.parse_args()

    run_pipeline(
        ref_path=args.ref,
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        out_dir=args.out_dir,
        threads=args.threads,
        prk_prefix=args.prk_prefix,
        excel_file=args.excel_file
    )
