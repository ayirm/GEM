

class Controller(object):
    def __init__(self):
        pass

    def Quality_Control(self, fastq1, fastq2):
        """
        What it needs to do: Check the quality then send the files to Alignment
        Terminal Code: 
            fastqc -o out_dir/ *.fastq <-- * makes it so that any fastq file is checked
        """

    def Alignment_Refs(self, ref_fasta, fastq1, fastq2, out_dir):
        """
        What it needs to do: Index the reference and then create a .sam file which then will be send to Samtools
        Terminal Code:
            bowtie2-build --threads 4 ref_path/file.fna ref_path/bowtie2 <-Name of the files that it will create
            bowite2 -x ref_path/bowtie -1 fastq1_loc/file.fastq -2 fastq2_loc/file.fasta -S sam_loc/file.sam
        """

    def Samtools_Sorting(self, threads, sam_file, out_dir):
        """
        What it needs to do: Sort the .bam file and then send the file to the bcftools
        Terminal Code:
            samtools view -uS -o bam_loc/file.bam sam_loc/file.sam
            samtools sort -@ 4 -T temp_file.tmp.sort -o bam_loc/file_sorted.bam
            samtools index bam_loc/file_sorted.bam
        """

    def Bcftools_Piling(self, ref_fasta, sorted_bam, out_dir):
        """
        What it needs to do: Collapse the reads, zip the resulting .vcf file and send the zipped file to assembly
        Terminal Code:
            bcftools mpileup -f ref_loc/file.fa bam_loc/file_sorted.bam | bcftools call -mv -Ob -o vcf_loc/file.bcf
            bcftools convert -O v -o vcf_loc/file.vcf vcf_loc/file.bcf
            bgzip -c vcf_loc/file.vcf > gz_loc/file.vcf.gz
            tabix -p vcf gz_loc/file.vcf.gz
        """

    def Reference_Assembly(self, ref_fasta, zip_file, out_dir):
        """
        What it needs to do: Creates the consensus sequence then sends the fasta file to Annotation
        Terminal Code:
            bcftools concensus -f ref_loc/file.fa gz_loc/file.vcf.gz > file.fasta
        """

    def Annotation(self, consensus_path, prokkaName, out_dir):
        """
        What it needs to do: Annote the file and send the .gbf or .gbk file to Reading
        Terminal Code:
            prokka --outdir prk_loc --prefix genomeName file.fasta
        """

    def Reading_Annoted(self, ann_file, prokkaName):
        """
        What it needs to do: Read the file and get the
            - gene_id
            - protein_name
            - EC_number
            - UniProt ID
        then send the UniProt ID to UniProt Search and EC_number to KEGG Search.
        """

    def UniProt_Search(self, cds_values):
        """
        What it needs to do: Uses QuickGO from bioservices to send the UniProt IDs. Results are then added to the cds_values dictionary
        """

    def Kegg_Pathway_Search(self, cds_values):
        """
        What it needs to do: Get the UniProt ID's from the cds_values dictionary then use mapping() to convert them to Kegg IDs which look like this: \n
        {'results': [{'from': 'P43403', 'to': 'hsa:7535'}]} \n
        Then it needs to strip the organism id(hsa) and gene_id(7535) part and get the pathways with get_pathway_by_gene() \n
        Pathways then added to the dictionary.
        """

    def Kegg_Reaction_Search(self, pathway_info):
        """
        What it needs to do: Get ALL of the reactions inside a pathway(which is in the dictionary)
        """

    