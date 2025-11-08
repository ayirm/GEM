import subprocess, os, sys, argparse
import pandas as pd

from uniprotGO import UniprotGO
from uniprotMapping import UniprotMapping
from keggAPI import KeggAPI

def run_cmd(cmd:str, cmd_step="Unknown Step"):
    """
    Mainly a helper for removing try-except blocks and making my brain burn less
    """
    print(f"---Running {cmd_step}---\n")
    print(f"Full Command: \n {''.join(cmd)}\n")

    try:
        process = subprocess.run(
            cmd,
            capture_output=True,
            check=True,
            shell=True # For allowing | and > 
        )
        print(f"---Command {cmd_step} was successfull---\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"---Error happened with {cmd_step}--- \n", file=sys.stderr)
        print(f"Error return: \n {e.returncode}\n", file=sys.stderr)
        print(f"Stdout: \n {e.stdout} \n", file=sys.stderr)
        print(f"Stderr: \n {e.stderr} \n", file=sys.stderr)
        return False
    
def check_output_files(out_file:str, step_name:str) -> bool:
    """
    Doesn't stop the subsequent steps in multiple shell command functions.
    Mainly used for prokka
    """
    if os.path.exists(out_file):
        resp = input(f"Found existing {step_name} at {out_file}. Skip this step and reuse it? [Y/n]: ")
        if resp.strip().lower() in ["", "y", "yes"]:
            print(f"Skipping {step_name}, using existing file.\n")
            return True
    return False

def flatten_entry(entry: dict) -> dict:
    """
    Flattens nested lists and dicts in a KEGG annotation entry for Excel export.
    """
    flat = entry.copy()

    # Convert lists to comma-separated strings
    flat["KO_Terms"] = ", ".join(entry.get("KO_Terms", []))
    flat["Pathways"] = ", ".join(entry.get("Pathways", []))
    flat["BRITE_id"] = ", ".join(entry.get("BRITE_id", []))

    # Flatten reactions: pathID -> [RXXXX]
    reactions = entry.get("Reactions", {})
    flat["Reactions"] = "; ".join(
        f"{path}: {', '.join(rxs)}" for path, rxs in reactions.items()
    )

    # Flatten compounds: RXXXX -> {equation, cmpdID}
    compounds = entry.get("Compounds", {})
    flat["Compounds"] = "; ".join(
        f"{rid}: {vals.get('equation','N/A')} ({', '.join(vals.get('cmpdID', []))})"
        for rid, vals in compounds.items()
    )

    return flat


def export_to_excel(cds_values, output_path="cds_kegg_results.xlsx"):
    """
    Converts cds_values (list of dicts) into a nicely flattened Excel file.
    """
    flattened = [flatten_entry(entry) for entry in cds_values]
    df = pd.DataFrame(flattened)
    df.to_excel(output_path, index=False)
    print(f"✅ KEGG results saved to {output_path}")

# ---- Terminal Pipeline Functions ---

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

# --- Other file calling and excel writing

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
    if not check_output_files(sam_file, "Alignemnt and Sam file creation"):
        sam_file = run_alignement(ref_path, fastq1, fastq2, out_loc=align_and_sam, threads=threads)

    if not check_output_files(sorted_bam, "Sorting the Bam file"):
        sorted_bam = run_samtools(sam_file, threads=threads, out_loc=align_and_sam)

    if not check_output_files(gz_vcf, "Collapsing Reads"):
        gz_vcf = run_bcftools(ref_path, sorted_bam, out_loc=bcf_and_gz)

    if not check_output_files(cons_fasta, "Consensus File Creation"):
        cons_fasta = run_reference_assembly(ref_path, gz_vcf, out_loc=bcf_and_gz)

    if not os.path.exists(prk_gbk):
        prk_gbk = os.path.join(prokka_out, f"{prk_prefix}.gbk")
    if not check_output_files(prk_gbk, "annotation (Prokka)"):
        prkName, prk_path = run_annotation(cons_fasta, out_loc=prokka_out, prkName=prk_prefix)
    else:
        prkName, prk_path = prk_prefix, prokka_out

    print(f"---Control: Starting parsing and GO Term finding")
    goTerms = UniprotGO(annPath=prk_path, prkName=prk_prefix)
    cds_values = goTerms.organize_GO()
    print(f"---Control: Parsing and GO term finding ended")

    print(f"---Control: Turning UniProtIDs to KEGG IDs")
    uniprotID = [entry["UniProt_ID"] for entry in cds_values if entry["UniProt_ID"]]
    mapping = UniprotMapping()
    mapping_dict = mapping.organize_Mapping(uniProt_ids=uniprotID)
    
    for entry in cds_values:
        uid = entry.get("UniProt_ID")
        entry["KEGG_ID"] = mapping_dict.get(uid)

    print(f"---Control: Get the pathway and the linked info from KEGG")
    # We need to run keggAPI for every seperated value, for loop can be run either in here or in KeggAPI
    keggRequests = KeggAPI(cdsValues=cds_values)

    for i, keggID in enumerate(cds_values, 1):
        kID = entry.get("KEGG_ID")
        if not kID:
            print(f"No id for {kID}, skipping")
            continue

        results_dict = keggRequests.organizeKeggResults(keggID=kID)

        if results_dict:
            keggID.update({
                "KO_Terms": results_dict.get("KO_Terms", []),
                "Pathways": results_dict.get("pathways", []),
                "BRITE": results_dict.get("BRITE", "N/A"),
                "BRITE_id": results_dict.get("BRITE_id", []),
                "Reactions": results_dict.get("reactions", {}),
                "Compounds": results_dict.get("compunds", {})
            })
        else:
            print(f"[{i}] No KEGG data found for {kID}")

    print(f"---Control: Create the excel file")
    export_to_excel(cds_values, output_path="annotation_results.xlsx")


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