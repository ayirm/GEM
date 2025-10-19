import subprocess, os, sys

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

    run_cmd(f"bowtie2-build --threads {str(threads)} {ref_path} {base_name}","Indexing Reference File") 
    # Just as a nitpick, this creates the indexes in the current diretory. I don't like the cluster it creates, possibly modify the {base_name} to group with others

    run_cmd(f"bowtie2 -x {base_name} -1 {fastq1_loc} -2 {fastq2_loc} -S {sam_file}","Aligning FastQ reads")

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

prokka_bin = os.environ.get("PROKKA_BIN", "prokka")
prokka_db = os.environ.get("PROKKA_DB", "./prokka_db")

def run_annotation(cons_path:str, out_loc:str,prkName="prokkaAnnotes"):
    """
    Runs prokka to annote the reference assembled genome

    Code in terminal:
        prokka --outdir prk_loc --prefix genomeName file.fasta
    """

    os.makedirs(out_loc, exist_ok=True)

    run_cmd(f"{prokka_bin} --outdir {out_loc} --prefix {prkName} --kingdom Bacteria --force {cons_path}","Annotation with Prokka")

    return prkName, out_loc

def uniprot_search(ann_path:str, prkName="prokkaAnnotes"):
    """
    Uses prokka results to get gene_id, protein_name, EC_number and Gene Onthology via UniProt
    Basically reads the prokka.gbf file and selects the CDS parts and then separetes hypothetical proteins that doesnt have any EC number or UniProt id
    After getting UniProt id, searches it to find Gene Onthologies \n
    Result is writted into a list \n

    NOTE: genes can also be used for UniProt but would need to specify the organism, and i don't think .gbf has it. \n
    NOTE: Also .gff file can be read too since it has the same information(at least the things we need), but i find .gbf easier on the eyes
    """

    gbf_file = os.path.join(ann_path, f"{prkName}.gbf")

    total_CDS = sum(1 for record in SeqIO.parse(gbf_file, "genbank") for f in record.features if f.type == "CDS")

    cds_values = []

    for i , record in enumerate(SeqIO.parse(gbf_file, "genbank"), 1):
        for feature in record.features:
            if feature != "CDS":
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
    # !TODO: Make another function to check the empty slots for each value and then hit UniProt to fill them
    # Only values with the UniProt id is written and this includes the hypothetical proteins
    # Also some proteins have id's but no gene or EC number
    
    # GO Term Search Part
    go = QuickGO()

    uniprot_ids = [entry["inference"] for entry in cds_values if entry["inference"]]

    go_data = {}

    for i, uid in enumerate(uniprot_ids, 1):
        response = go.Annotation(
            geneProductId= uid,
            includeFields= "goName",
        )
   
        go_data[uid] = response.get("results", [])

        sys.stdout.write(f"\r Searching {i}th GO term out of {len(uniprot_ids)}")
        sys.stdout.flush()

    # UniProt ids to KEGG ids part
    u = UniProt(verbose=False)

    try:
        mapping_result = u.mapping("UniProtKB_AC-ID", "KEGG", uniprot_ids)
    except Exception as e:
        print(f"Warning: failed KEGG mapping: {e}")
        mapping_result = {}

    # Adding the values to dictonary in order to put it in excel file later
    for d in cds_values:
        uid = d["inference"]
        go_info = go_data.get(uid, [])
        d["KEGG"] = mapping_result.get(uid)
        d["GO"] = [f"{g['goId']}: {g['goName']}" for g in go_info]
    
    return mapping_result, cds_values

# !FIXME : It downloads something it think? I haven't let it run fully because it showed up as 7 hours long download
# It's disables as of now
def kegg_search(KeggIDs):
    """
    Finds the pathways that has the input protein involved, also results are filtered on organism codes
    KeggID value should be similiar to {'results': [{'from': 'P43403', 'to': 'hsa:7535'}]}
    Only to keys' values will be used for pathway creation
    """
    k = KEGG(verbose=False)

    to_list = [r.get('to') for r in KeggIDs.get('results', []) if r.get('to')]
    
    orgs = []
    genes = []
    for entry in to_list:
        if ':' in entry:
            org, gene_id = entry.split(':', 1)
            orgs.append(org)
            genes.append(gene_id)
        else:
            # If malformed, skip or fill with None
            orgs.append(None)
            genes.append(entry)
    
    kegg_results = k.get_pathway_by_gene(gene=genes, organism=orgs)

    return kegg_results

def excel_creation(cds_values, kegg_results, output_file="annotation_results.xlsx"):
    """
    Writes CDS annotations and KEGG pathway results into an Excel file.
    """
    wb = Workbook()

    # --- Sheet 1: CDS values ---
    ws1 = wb.active
    ws1.title = "CDS_Annotations"

    # Write headers dynamically
    headers = list(cds_values[0].keys()) if cds_values else []
    ws1.append(headers)

    for entry in cds_values:
        # Convert lists (like GO) to comma-separated strings
        row = [
            ", ".join(v) if isinstance(v, list) else v
            for v in entry.values()
        ]
        ws1.append(row)

    # --- Sheet 2: KEGG pathways ---
    ws2 = wb.create_sheet("KEGG_Pathways")
    ws2.append(["KEGG_ID", "Pathway_ID", "Pathway_Name"])

    for gene, data in kegg_results.items():
        if not data:
            continue
        pathways = data.get('pathways', {})
        for pid, pname in pathways.items():
            ws2.append([gene, pid, pname])

    wb.save(output_file)
    print(f"✅ Results written to '{output_file}'")

# We need to input two Fastq file name and one reference file
def run_pipeline(ref_path, fastq1, fastq2, out_dir="results", threads=4, prk_prefix="prokkaAnnotes", excel_file="annotation_results.xlsx"):
    os.makedirs(out_dir, exist_ok=True)

    # Folders in folders part
    align_and_sam = os.path.join(out_dir, "sam")
    bcf_and_gz = os.path.join(out_dir, "bcf")
    prokka_out = os.path.join(out_dir, "prk")

    # Alignment
    sam_file = run_alignement(ref_path, fastq1, fastq2, out_loc=align_and_sam, threads=threads)
    sorted_bam = run_samtools(sam_file, threads=threads, out_loc=align_and_sam)
    gz_vcf = run_bcftools(ref_path, sorted_bam, out_loc=bcf_and_gz)
    cons_fasta = run_reference_assembly(ref_path, gz_vcf, out_loc=bcf_and_gz)
    prkName, prk_path = run_annotation(cons_fasta, out_loc=prokka_out, prkName=prk_prefix)

    # Protein search & KEGG mapping
    mapping_result, cds_values = uniprot_search(prkName=prkName, ann_path=prk_path)
    
    #kegg_results = kegg_search(mapping_result)
    kegg_results = {}

    # Excel output
    excel_creation(cds_values, kegg_results, output_file=os.path.join(out_dir, excel_file))

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
