#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    iumobg/model_creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/iumobg/model_creation
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MODEL_CREATION  } from './workflows/model_creation'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_model_creation_pipeline'
include { PIPELINE_COMPLETION      } from './subworkflows/local/utils_nfcore_model_creation_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_model_creation_pipeline'

nextflow.enable.dsl=2

// Modülleri projeye dahil ediyoruz
include { SPADES          } from './modules/nf-core/spades/main'
include { SAMTOOLS_SORT   } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX  } from './modules/nf-core/samtools/index/main'
// Kendi yazdığın veya nf-core'dan indirdiğin Prokka:
include { PROKKA          } from './modules/nf-core/prokka/main'
// --- SENİN MATRİS SUBWORKFLOW'UN BURADA ---
include { MATRIX_GENERATION_WF } from './subworkflows/nf-core/matrix/matrix_wf.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow IUMOBG_MODEL_CREATION {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    MODEL_CREATION (
        samplesheet
    )
    emit:
    multiqc_report = MODEL_CREATION.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
    // 1. HAZIRLIK KISMINI KAPATIYORUZ (Hata veren yer burası)
    /*
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        ...
    )
    */

    // 2. ANA ANALİZİ KAPATIYORUZ (Çünkü samplesheet bekliyor)
    /*
    IUMOBG_MODEL_CREATION (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    */

    // 3. SENİN MATRİS KISMIN (Burası açık kalsın)
    // params.input yerine params.matrix_input kullanalım ki nf-core karışmasın
    if (params.matrix_input) {
        ch_matrix_json = Channel.fromPath(params.matrix_input, checkIfExists: true)
        MATRIX_GENERATION_WF ( ch_matrix_json )
    }

    // 4. TAMAMLAMA KISMINI DA KAPATALIM (Rapor verecek veri olmayacak)
    /*
    PIPELINE_COMPLETION ( ... )
    */
}
 
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
