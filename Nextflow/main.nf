#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
READS - TO - PATHWAYS - PIPLINE
=====================================================
Reference:    ${params.ref}
FastQ 1:      ${params.fastq1}
FastQ 2:      ${params.fastq2}
Outdir:       ${params.outdir}
Threads:      ${params.threads}
"""

// Channel for the reference FASTA, used by multiple processes
ch_ref_fasta = Channel.value( file(params.ref) )

// Get a sample ID from the FastQ name
// e.g., "my_sample_R1.fastq.gz" -> "my_sample"
def fastq1_file = file(params.fastq1)
def sample_id = fastq1_file.simpleName.replaceAll('(?i)_R1.*', '')
log.info "Using sample ID: ${sample_id}"

// Channel for a single sample, matching your script's structure
ch_sample = Channel.of( [ sample_id, file(params.fastq1), file(params.fastq2) ] )

workflow {
    // 1. Index the reference genome
    ch_bowtie_index = Bowtie_Index( ch_ref_fasta )

    // 2. Align reads to the indexed reference
    ch_sam = Bowtie_Align( ch_bowtie_index.combine(ch_sample) )

    // 3. Convert SAM to sorted, indexed BAM
    ch_bam = Run_Samtools( ch_sam )

    // 4. Call variants with BCFtools and collapse them into vcf
    ch_vcf = Run_Bcftools( ch_ref_fasta.combine(ch_bam) )

    // 5. Create consensus sequence
    ch_consensus = Reference_Assembly( ch_ref_fasta.combine(ch_vcf) )

    // 6. Annotate consensus with Prokka
    ch_prokka = Annotation( ch_consensus )

    // 7. Parse Prokka output, search UniProt to get gene,EC,Protein name,UniProt id
    ch_cds_json = Uniprot_Search( ch_prokka )

    // 8. Search KEGG, get pathways/reactions/compunds
    ch_kegg_json = Kegg_Search( ch_cds_json )

    Create_Excel( ch_kegg_json )
}

process Bowtie_Index {
    publishDir "${params.outdir}/bowtie2_index", mode: 'copy', overwrite: false

    input:
    path ref_file

    output:
    path "bowtie2_index"

    script:
    def ref_prefix = ref_fasta.simpleName
    """
    mkdir bowtie2_index
    bowtie2-build --threads ${task.cpus} ${ref_fasta} bowtie2_index/${ref_prefix}
    """
}

process Bowtie_Align {
    publishDir "${params.outdir}/sam", mode: 'copy', overwrite: true

    input:
    tuple path(index_dir), val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    def ref_prefix = file(params.ref).simpleName
    """
    bowtie2 -x ${index_dir}/${ref_prefix} \
        -1 ${fq1} \
        -2 ${fq2} \
        -S ${sample_id}.sam \
        --threads ${task.cpus}
    """
}

process Run_Samtools {
    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    samtools view -uS -o ${sample_id}.bam ${sam_file}
    samtools sort -@ ${task.cpus} -T ${sample_id}.tmp.sort -o ${sample_id}_sorted.bam ${sample_id}.bam
    samtools index ${sample_id}_sorted.bam
    """
}

process Run_Bcftools {
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true

    input:
    tuple path(ref_fasta), val(sample_id), path(bam_file), path(bai_file)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${bam_file} | \
        bcftools call -mv -Oz -o ${sample_id}.vcf.gz
    
    tabix -p vcf ${sample_id}.vcf.gz
    """
}

process Reference_Assembly {
    publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

    input:
    tuple path(ref_fasta), val(sample_id), path(vcf_file), path(tbi_file)

    output:
    tuple val(sample_id), path("${sample_id}.consensus.fasta")

    script:
    """
    bcftools consensus -f ${ref_fasta} ${vcf_file} > ${sample_id}.consensus.fasta
    """
}

process Annotation {
    publishDir "${params.outdir}/prokka", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(consensus_fasta)

    output:
    tuple val(sample_id), path(prokka_dir)

    script:
    def prokka_dir = "${sample_id}_prokka"
    """
    prokka --outdir ${prokka_dir} \
        --prefix ${sample_id} \
        --force \
        ${consensus_fasta}
    """
}

process Uniprot_Search {
    publishDir "${params.outdir}/annotation_data", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(prokka_dir)

    output:
    tuple val(sample_id), path(prokka_dir), path("${sample_id}_cds_values.json.gz")

    script:
    """
    #!/usr/bin/env python3
    import sys, json, os, gzip
    from annotation_helpers import uniprot_search

    sample_id = "${sample_id}"
    prokka_dir = "${prokka_dir}"
    out_file = f"{sample_id}_cds_values.json"

    # Caching will happen inside the process work directory
    cds_values = uniprot_search(ann_path=prokka_dir, prkName=sample_id, cache_file="go_cache.json")
    
    with open(out_file, 'w') as f:
        json.dump(cds_values, f, indent=2)
    
    gzip ${out_file}
    """
}

process Kegg_Search {
    """
    Takes the UniProt IDs from the previous step, maps them to KEGG IDs,
    and fetches pathway, reaction, and BRITE data. Caches KEGG results
    and outputs a final gzipped JSON.
    """
    publishDir "${params.outdir}/annotation_data", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(prokka_dir), path(cds_json_gz)

    output:
    tuple val(sample_id), path("${sample_id}_enriched_cds.json.gz")

    script:
    """
    #!/usr/bin/env python3
    import sys, json, os, gzip
    from annotation_helpers import kegg_search

    sample_id = "${sample_id}"
    prokka_dir = "${prokka_dir}"
    cds_json_file_gz = "${cds_json_gz}"
    out_file = f"{sample_id}_enriched_cds.json"

    # Read the gzipped JSON file
    with gzip.open(cds_json_file_gz, 'rt', encoding='utf-8') as f:
        cds_values = json.load(f)

    # Caching will happen inside the process work directory
    enriched_cds = kegg_search(cds_values, ann_path=prokka_dir, gene_cache_file="kegg_gene_cache.json", ko_cache_file="kegg_ko_cache.json")
    
    with open(out_file, 'w') as f:
        json.dump(enriched_cds, f, indent=2)
    
    gzip ${out_file}
    """
}

process Create_Excel {
    """
    Reads the final gzipped JSON with all enriched annotation data
    and creates a human-readable Excel report.
    """
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(enriched_json_gz)

    output:
    path "${sample_id}_annotation_results.xlsx"

    script:
    """
    #!/usr/bin/env python3
    import sys, json, os, gzip
    from annotation_helpers import excel_creation

    sample_id = "${sample_id}"
    enriched_json_file_gz = "${enriched_json_gz}"
    output_file = f"${sample_id}_annotation_results.xlsx"

    # Read the gzipped JSON file
    with gzip.open(enriched_json_file_gz, 'rt', encoding='utf-8') as f:
        enriched_data = json.load(f)
    
    excel_creation(enriched_data, output_file=output_file)
    """
}