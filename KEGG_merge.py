#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import json
import requests
import time
import argparse

# Gerekli kütüphaneler (Biopython, Openpyxl, Bioservices)
try:
    from Bio import SeqIO
    from openpyxl import Workbook
    from bioservices import QuickGO
except ImportError as e:
    print(f"Hata: Gerekli bir kütüphane eksik: {e.name}", file=sys.stderr)
    print(f"Lütfen 'pip install {e.name}' komutu ile kurun.", file=sys.stderr)
    sys.exit(1)

# -----------------------------------------------------------------
# YARDIMCI FONKSİYONLAR (PIPELINE'DAN)
# -----------------------------------------------------------------

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
        print(f"---Error happened with {cmd_name}--- \n", file=sys.stderr)
        print(f"Error return: \n {e.returncode}\n", file=sys.stderr)
        print(f"Stdout: \n {e.stdout} \n", file=sys.stderr)
        print(f"Stderr: \n {e.stderr} \n", file=sys.stderr)
        return False # Hata durumunda False dön

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
                print(f"\nAttempt {attempt + 1}/{max_retries} failed for chunk starting at index {i}: {e}", file=sys.stderr)
                if attempt + 1 == max_retries:
                    print(f"--- Giving up on this chunk after {max_retries} attempts. ---", file=sys.stderr)
                else:
                    delay = 2 ** attempt
                    print(f"--- Retrying in {delay} seconds... ---")
                    time.sleep(delay)

    return mapping_dict

# -----------------------------------------------------------------
# YENİ KEGGCLIENT SINIFI (BİRLEŞTİRİLDİ)
# -----------------------------------------------------------------

class KEGGClient:
    """
    KEGG veritabanıyla web servisi (API) üzerinden iletişim kurmak
    için kullanılan istemci sınıfı.
    
    Bu sınıf, 'requests' kütüphanesini ve saniyede 3 istek
    kuralına uymak için 'time.sleep' kullanır.
    'Bio.KEGG.REST' KULLANMAZ.
    """
    BASE_URL = "https://rest.kegg.jp"
    
    def __init__(self, rate_limit_sec=0.35):
        """Sınıf başlatıldığında bir mesaj yazdırır. (3 req/sn için ~0.35sn)"""
        print("[KEGGClient] İstemci başlatıldı. KEGG API'ye bağlanmaya hazır.")
        self.rate_limit = rate_limit_sec
        self.last_request_time = 0

    def _make_request(self, endpoint):
        """
        Rate-limiting (hız sınırı) uygulayan özel istek fonksiyonu.
        """
        # Hız sınırını kontrol et
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        
        url = f"{self.BASE_URL}{endpoint}"
        try:
            response = requests.get(url)
            self.last_request_time = time.time()
            response.raise_for_status() # Hata varsa (404, 500) exception fırlat
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"\n[KEGGClient Hata] URL '{url}' sorgulanırken hata: {e}", file=sys.stderr)
            return None # Hata durumunda None dön

    def get_linked_entries(self, target_db, source_id):
        """
        KEGG 'link' operasyonu.
        örn: /link/pathway/hsa:7157
        """
        endpoint = f"/link/{target_db}/{source_id}"
        # print(f"[KEGGClient Link Sorgusu] '{source_id}' -> '{target_db}' bağlantıları aranıyor...")
        raw_result = self._make_request(endpoint)
        
        if not raw_result:
            # print("[KEGGClient Sonuç] Bağlantı bulunamadı.")
            return []
        
        links = []
        for line in raw_result.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) == 2:
                links.append(parts[1]) # Hedef ID (ikinci sütun)
        
        # print(f"[KEGGClient Sonuç] {len(links)} adet bağlantı bulundu.")
        return links

# -----------------------------------------------------------------
# PIPELINE FONKSİYONLARI (PIPELINE'DAN)
# -----------------------------------------------------------------

def run_alignement(ref_path:str, fastq1_loc:str, fastq2_loc:str, out_loc:str, threads=4):
    """
    Indexes the reference genome then aligns the fastq reads to the reference as a sam file
    """
    os.makedirs(out_loc, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(ref_path))[0]
    sam_file = os.path.join(out_loc, f"{base_name}.sam")
    ref_index_dir = os.path.join(out_loc, "ref_index")
    os.makedirs(ref_index_dir, exist_ok=True)
    ref_index_prefix = os.path.join(ref_index_dir, base_name)

    run_cmd(f"bowtie2-build --threads {str(threads)} {ref_path} {ref_index_prefix}","Indexing Reference File") 
    run_cmd(f"bowtie2 -x {ref_index_prefix} -1 {fastq1_loc} -2 {fastq2_loc} -S {sam_file}","Aligning FastQ reads")
    return sam_file

def run_samtools(sam_path:str, threads:int, out_loc:str):
    """
    Runs the SAM to BAM conversion and indexing
    """
    os.makedirs(out_loc, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(sam_path))[0]
    bam_loc = os.path.join(out_loc, f"{base_name}.bam")
    sorted_loc = os.path.join(out_loc, f"{base_name}_sorted.bam")
    tmp_file = os.path.join(out_loc, f"{base_name}.tmp.sort")

    run_cmd(f"samtools view -uS -o {bam_loc} {sam_path}","SAM to BAM conversion")
    run_cmd(f"samtools sort -@ {str(threads)} -T {tmp_file} -o {sorted_loc} {bam_loc}", "BAM file sorting")
    run_cmd(f"samtools index {sorted_loc}","İndexing the sorted BAM file")
    return sorted_loc

def run_bcftools(ref_path:str, sorted_path:str ,out_loc:str):
    """
    Collapses the reads according to reference and writes variant information to a VCF file. 
    Also compresses the VCF file with bgzip and indexes it with tabix
    """
    os.makedirs(out_loc, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(sorted_path))[0].replace("_sorted", "")
    bcf_loc = os.path.join(out_loc, f"{base_name}.bcf")
    vcf_loc = os.path.join(out_loc, f"{base_name}.vcf")
    gz_loc = os.path.join(out_loc, f"{base_name}.vcf.gz")

    run_cmd(f"bcftools mpileup -f {ref_path} {sorted_path} | bcftools call -mv -Ob -o {bcf_loc}", "Collapsing reads to a BCF file")
    run_cmd(f"bcftools convert -O v -o {vcf_loc} {bcf_loc}","Converting BCF to VCF")
    run_cmd(f"bgzip -c {vcf_loc} > {gz_loc}", "Compressing VCF")
    run_cmd(f"tabix -p vcf {gz_loc}", "Indexing compressed VCF")
    return gz_loc

def run_reference_assembly(ref_path:str, gz_loc:str, out_loc:str):
    """
    Creates a concensus sequence with bcftools. 
    """
    os.makedirs(out_loc, exist_ok=True)
    base_name = os.path.basename(gz_loc).replace(".vcf.gz","")
    cons_loc = os.path.join(out_loc, f"{base_name}.fasta")
    run_cmd(f"bcftools consensus -f {ref_path} {gz_loc} > {cons_loc}","Creating consensus sequence")
    return cons_loc

def run_annotation(cons_path:str, out_loc:str,prkName="prokkaAnnotes"):
    """
    Runs prokka to annote the reference assembled genome
    """
    os.makedirs(out_loc, exist_ok=True)
    run_cmd(f"prokka --outdir {out_loc} --prefix {prkName} --force {cons_path}","Annotation with Prokka")
    return prkName, out_loc

def uniprot_search(ann_path:str, prkName="prokkaAnnotes", cache_file="go_cache.json"):
    """
    Uses prokka results to get gene_id, protein_name, EC_number and Gene Onthology via UniProt
    """
    gbf_file = os.path.join(ann_path, f"{prkName}.gbf")
    if not os.path.exists(gbf_file):
        gbf_file = os.path.join(ann_path, f"{prkName}.gbk")

    total_CDS = 0
    cds_values = []

    try:
        for record in SeqIO.parse(gbf_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    total_CDS += 1
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
    except FileNotFoundError:
        print(f"Hata: Prokka GBK/GBF dosyası bulunamadı: {gbf_file}", file=sys.stderr)
        return []

    sys.stdout.write(f"Parsing {len(cds_values)} valid UniProt ids out of {total_CDS}\n")
    print(f"Total with UniProt IDs: {len(cds_values)}")
    
    # GO Term Search Part
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
            print(f"Error loading cache '{go_cache_path}': {e}. Re-fetching all GO terms.", file=sys.stderr)
            go_data = {} 

    new_ids_to_fetch = [uid for uid in uniprot_ids if uid not in go_data]
    if new_ids_to_fetch:
        print(f"\nCache found, but {len(new_ids_to_fetch)} new GO terms need to be fetched from QuickGO...")
        go = QuickGO()
        
        for i, uid in enumerate(new_ids_to_fetch, 1):
            try:
                response = go.Annotation(
                    geneProductId= uid,
                    includeFields= "goName",
                )
                go_data[uid] = response.get("results", [])
            except Exception as e:
                print(f"\nError fetching GO for {uid}: {e}", file=sys.stderr)
                go_data[uid] = [] 

            sys.stdout.write(f"\r Searching {i}/{len(new_ids_to_fetch)} new GO terms...")
            sys.stdout.flush()
        
        print("\nFinished fetching new GO terms.")

        if cache_file:
            try:
                with open(go_cache_path, 'w') as f:
                    json.dump(go_data, f, indent=2)
                print(f"Saved updated GO data to cache: '{go_cache_path}'")
            except Exception as e:
                print(f"Error saving cache: {e}", file=sys.stderr)
    else:
        print("All GO terms were loaded from cache.")

    # Adding the values back to the dictornary
    for d in cds_values:
        uid = d["inference"]
        go_info = go_data.get(uid, [])
        d["GO_Terms"] = "\n".join([f"{g['goId']}: {g['goName']}" for g in go_info]) if go_info else "N/A"
    
    return cds_values

# -----------------------------------------------------------------
# YENİ KEGG_SEARCH FONKSİYONU (BİRLEŞTİRİLDİ)
# -----------------------------------------------------------------

def kegg_search(cds_values):
    """
    KEGGClient sınıfını ve "Gene -> Pathway -> Reaction" iş akışını kullanarak
    UniProt ID'lerini zenginleştirir.
    
    Pipeline'dan gelen `map_uniprot_to_kegg_via_rest` fonksiyonunu kullanır.
    Pipeline'dan gelen `KEGGClient` sınıfını kullanır.
    """
    
    # --- 1. UniProt -> KEGG ID Eşleştirmesi ---
    uniprot_ids = [entry["inference"] for entry in cds_values]
    unique_uniprot_ids = list(set(uniprot_ids))

    print(f"\n--- Running KEGG Pathway and Reaction Search ---")
    print(f"Doing UniProt to KEGG transform for {len(unique_uniprot_ids)} unique IDs")
    uniprot_to_kegg = map_uniprot_to_kegg_via_rest(unique_uniprot_ids)
    print(f"\nSuccessfully mapped {len(uniprot_to_kegg)} UniProt IDs.")

    # --- 2. KEGG API İstemcisini Başlat ---
    kegg_api_client = KEGGClient() # Yeni, requests-tabanlı istemcimiz
    
    # --- 3. Her Geni İşle (İş Akışı) ---
    total_entries = len(cds_values)
    print(f"Processing {total_entries} entries for KEGG pathways and reactions...")
    
    for i, entry in enumerate(cds_values, 1):
        sys.stdout.write(f"\rProcessing entry {i}/{total_entries}...")
        sys.stdout.flush()
        
        uniprot_id = entry["inference"]
        kegg_id = uniprot_to_kegg.get(uniprot_id)
        
        # Başlangıç değerlerini ayarla
        entry["kegg_id"] = kegg_id if kegg_id else "N/A"
        entry["pathways"] = "N/A" # String olarak başlat
        entry["reactions_by_pathway"] = "N/A" # String olarak başlat
        
        if not kegg_id:
            continue # Eşleşme yoksa atla

        try:
            # --- YENİ İŞ AKIŞI (GENE -> PATHWAY -> REACTION) ---
            
            # 1. Gene bağlı pathway'leri bul
            pathway_idler = kegg_api_client.get_linked_entries("pathway", kegg_id)
            if not pathway_idler:
                entry["pathways"] = "No pathway found"
                continue
            
            # Pathway'leri \n ile ayrılmış metin yap
            entry["pathways"] = "\n".join(sorted(list(set(pathway_idler))))
            
            # 2. Bu pathway'lerin reaksiyonlarını bul
            pathway_reaksiyon_iliskisi = {}
            
            for p_id in pathway_idler:
                # 'path:ecj00260' -> 'map00260'
                map_number = p_id[-5:]
                reference_map_id = f"map{map_number}"
                
                # Referans haritadaki reaksiyonları çek
                reaksiyonlar = kegg_api_client.get_linked_entries("reaction", reference_map_id)
                
                if reaksiyonlar:
                    pathway_reaksiyon_iliskisi[p_id] = sorted(list(set(reaksiyonlar)))
            
            # 3. Sonuçları JSON formatında kaydet (Excel için en iyisi)
            if pathway_reaksiyon_iliskisi:
                # JSON'ı Excel hücresine sığdırmak için indent=None
                entry["reactions_by_pathway"] = json.dumps(pathway_reaksiyon_iliskisi)
            else:
                entry["reactions_by_pathway"] = "No reactions found for these pathways"
        
        except Exception as e:
            print(f"\nError processing KEGG workflow for {kegg_id}: {e}", file=sys.stderr)
            entry["pathways"] = f"Error: {e}"
            entry["reactions_by_pathway"] = f"Error: {e}"
    
    print("\n\nFinished KEGG search.")
    return cds_values

# -----------------------------------------------------------------
# EXCEL VE ANA PİPELİNE FONKSİYONLARI (PIPELINE'DAN)
# -----------------------------------------------------------------

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
        row = []
        for h in headers:
            value = entry.get(h, "")
            # Listeleri (artık sadece GO_Terms) \n ile ayır
            if isinstance(value, list):
                row.append("\n".join(map(str, value)))
            # Diğer her şeyi (string, None vb.) olduğu gibi ekle
            else:
                row.append(value)
        ws.append(row)

    wb.save(output_file)
    print(f"✅ Results written to '{output_file}'")

def run_pipeline(ref_path, fastq1, fastq2, out_dir="results", threads=4, prk_prefix="prokkaAnnotes", excel_file="annotation_results.xlsx"):
    os.makedirs(out_dir, exist_ok=True)

    # Folders in folders part
    align_and_sam = os.path.join(out_dir, "sam_bam") # Samtools de burayı kullanıyor
    bcf_and_gz = os.path.join(out_dir, "vcf_consensus") # BCF/VCF ve Consensus
    prokka_out = os.path.join(out_dir, "prokka_annotation")

    # Names for stopping the repeptition
    base_name = os.path.splitext(os.path.basename(ref_path))[0]
    sam_file = os.path.join(align_and_sam, f"{base_name}.sam")
    sorted_bam = os.path.join(align_and_sam, f"{base_name}_sorted.bam")
    gz_vcf = os.path.join(bcf_and_gz, f"{base_name}.vcf.gz")
    cons_fasta = os.path.join(bcf_and_gz, f"{base_name}.fasta")
    
    # Prokka'nın gbf veya gbk oluşturma durumunu kontrol et
    prk_gbf = os.path.join(prokka_out, f"{prk_prefix}.gbf")
    prk_gbk = os.path.join(prokka_out, f"{prk_prefix}.gbk")

    # Checks for file
    if not check_existing_file(sam_file, "Alignemnt and Sam file creation"):
        sam_file = run_alignement(ref_path, fastq1, fastq2, out_loc=align_and_sam, threads=threads)

    if not check_existing_file(sorted_bam, "Sorting the Bam file"):
        sorted_bam = run_samtools(sam_file, threads=threads, out_loc=align_and_sam)

    if not check_existing_file(gz_vcf, "Collapsing Reads"):
        gz_vcf = run_bcftools(ref_path, sorted_bam, out_loc=bcf_and_gz)

    if not check_existing_file(cons_fasta, "Consensus File Creation"):
        cons_fasta = run_reference_assembly(ref_path, gz_vcf, out_loc=bcf_and_gz)

    # İki potansiyel Prokka çıktı dosyasından birini kontrol et
    if not (check_existing_file(prk_gbf, "annotation (Prokka)") or check_existing_file(prk_gbk, "annotation (Prokka)")):
        prkName, prk_path = run_annotation(cons_fasta, out_loc=prokka_out, prkName=prk_prefix)
    else:
        prkName, prk_path = prk_prefix, prokka_out

    # --- ANOTASYON VE ZENGİNLEŞTİRME ---
    cds_values = uniprot_search(prkName=prkName, ann_path=prk_path)
    
    if cds_values:
        # Kodu basitleştirmek için kegg_search çağrısından gereksiz parametreleri kaldırdık
        enriched_results = kegg_search(cds_values)

        print("\n--- Creating Excel Report ---")
        excel_creation(enriched_results, output_file=os.path.join(out_dir, excel_file))
    else:
        print("Prokka'dan UniProt ID'si olan CDS bulunamadı. KEGG ve Excel adımları atlandı.", file=sys.stderr)

    print("✅ Pipeline finished successfully!")

# -----------------------------------------------------------------
# ANA ÇALIŞTIRMA BLOĞU (PIPELINE'DAN)
# -----------------------------------------------------------------

if __name__ == "__main__":
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