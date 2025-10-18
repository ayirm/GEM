import subprocess, os, shutil, sys

from Bio import Blast, SeqIO
from openpyxl import workbook

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
def run_alignement(ref_path:str, fastq1_loc:str, fastq2_loc:str, out_loc="sam", threads=4):
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
    run_cmd(f"bowtie2 -x {base_name} -1 {fastq1_loc} -2 {fastq2_loc} -S {sam_file}","Aligning FastQ reads")

    return sam_file # This is the full path used for run_samtools()'s sam_path

def run_samtools(sam_path:str, threads=4, out_loc="sam"):
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

def run_bcftools(ref_path:str, sorted_path:str ,out_loc="vcf"):
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
    run_cmd(f"bcftools convert -O v {vcf_loc} {bcf_loc}","Converting BCF to VCF")
    run_cmd(f"bgzip -c {vcf_loc} > {gz_loc}", "Compressing VCF")
    run_cmd(f"tabix -p vcf {gz_loc}", "Indexing compressed VCF")

    return gz_loc # Used in run_reference_assembly as gz_loc

def run_reference_assembly(ref_path:str, gz_loc:str, out_loc="cons"):
    """
    Creates a concensus sequence with bcftools. 

    Code in terminal:
        bcftools concensus -f ref_loc/file.fa gz_loc/file.vcf.gz > file.fasta
    """

    os.makedirs(out_loc, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(gz_loc))[0]
    cons_loc = os.path.join(out_loc, f"{base_name}.fasta")

    run_cmd(f"bcftools consensus -f {ref_path} {gz_loc} > {cons_loc}","Creating consensus sequence")

    return cons_loc # Used in prokka to annote

def run_annotation(cons_path:str, out_loc="prk",prkName="prokkaAnnotes"):
    """
    Runs prokka to annote the reference assembled genome

    Code in terminal:
        prokka --outdir prk_loc --prefix genomeName file.fasta
    """

    os.makedirs(out_loc, exist_ok=True)

    run_cmd(f"prokka --outdir {out_loc} --prefix {prkName} {cons_path}","Annotation with Prokka")

    return prkName, out_loc

def uniprot_search(prkName="prokkaAnnotes", ann_path="prk"):
    """
    Uses prokka results to get gene_id, protein_name, EC_number and Gene Onthology via UniProt \n
    Basically reads the prokka.gbf file and selects the CDS parts and then separetes hypothetical proteins that doesn't have any EC number or UniProt id
    After getting UniProt id, searches it to find Gene Onthologies \n
    Result is writted into a list \n

    NOTE: genes can also be used for UniProt but would need to specify the organism, and i don't think .gbf has it. \n
    NOTE: Also .gff file can be read too since it has the same information(at least the things we need), but i find .gbf easier on the eyes
    """

    gbf_file = os.path.join(ann_path, f"{prkName}.gbf")

def kegg_search():
    """
    Finds the pathways that has the input protein involved, also results are filtered on organism codes
    """

def excel_creation():
    pass