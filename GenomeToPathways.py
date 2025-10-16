import subprocess, os, shutil, sys

from Bio import Blast, SeqIO
from openpyxl import workbook

def run_cmd(cmd:list, cmd_name="No Name"):
    """Helper for subprocess command, helps to remove repetitive try-except blocks"""
    print(f"---Running {cmd_name}---\n")
    print(f"Full command: \n {''.join(cmd)}")

    try:
        process = subprocess.run(
            cmd,
            capture_output=True,
            check=True,
        )
        print(f"---Run of {cmd_name} was successfull--- \n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"---Error happened with {cmd_name}--- \n")
        print(f"Error return: \n {e.returncode}\n")
        print(f"Stdout: \n {e.stdout} \n")
        print(f"Stderr: \n {e.stderr} \n")

# --- Main Codes for Alignment--->Annotation--->Protein search in UniProt ---> KEGG pathway search
def run_alignement(sam_path:str, threads=4, out_loc="sam"):
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

    run_cmd(["samtools","view","-uS","-o",bam_loc, sam_path],"SAM to BAM conversion")
    run_cmd(["samtools","sort","-@",str(threads),"-T",tmp_file, "-o", sorted_loc, bam_loc], "BAM file sorting")
    run_cmd(["samtools","index",sorted_loc])

def run_bcftools(ref_path:str, sorted_path:str ,out_loc="vcf"):
    """
    Collapses the reads according to reference and writes variant information to a VCF file. 
    Also compresses the VCF file with bgzip and indexes it with tabix

    Code in terminal:
        bcftools mpileup -f ref_loc/file.fa bam_loc/file_sorted.bam | bcftools call -mv -Ob -o vcf_loc/file.bcf
        bcftools convert -O v -o vcf_loc/file.vcf vcf_loc/file.bcf
    """

    os.makedirs(out_loc, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(sorted_path))[0]
    bcf_loc = os.path.join(out_loc, f"{base_name}.bcf")
    vcf_loc = os.path.join(out_loc, f"{base_name}.vcf")
    gz_loc = os.path.join(out_loc, f"{base_name}.vcf.gz")

    run_cmd(["bcftools","mpileup","-f", ref_path, sorted_path, "|", "bcftools","call","-mv","-Ob","-o",bcf_loc], "Collapsing reads to a BCF file")
    run_cmd(["bcftools","convert","-O","v",vcf_loc,bcf_loc],"Converting BCF to VCF")
    run_cmd(["bgzip","-c",vcf_loc, ">",gz_loc], "Compressing VCF")
    run_cmd(["tabix","-p","vcf", gz_loc], "Indexing compressed VCF")

def run_reference_assembly():
    """
    Creates a concensus sequence with bcftools. 

    Code in terminal:
        bcftools concensus -f ref_loc/file.fa gz_loc/file.vcf.gz > file.fasta
    """

def run_annotation():
    """
    Runs prokka to annote the reference assembled genome

    Code in terminal:
        prokka --outdir prk_loc/ --prefix genomeName file.fasta
    """

def uniprot_search():
    """
    Uses the concensus genome to find the gene onthology, also translates protein ids into KEGG ids to send to KEGG function
    """

def kegg_search():
    """
    Finds the pathways that has the input protein involved, also results are filtered on organism codes
    """

def excel_creation():
    pass