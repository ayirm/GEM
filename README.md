# iumobg/model_creation


[![GitHub Actions CI Status](https://github.com/iumobg/model_creation/actions/workflows/nf-test.yml/badge.svg)](https://github.com/iumobg/model_creation/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/iumobg/model_creation/actions/workflows/linting.yml/badge.svg)](https://github.com/iumobg/model_creation/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/iumobg/model_creation)

## Introduction

**iumobg/model_creation** is a bioinformatics pipeline that converts raw FASTQ sequencing reads into a genome-scale metabolic model (GEM) and performs constraint-based simulations. The full pipeline conducts quality control, genome assembly(de novo or reference-guided) and annotation then integrates UniProt, GO and KEGG information to construct a CobraPy model.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/))
3. De Novo Assembly ([`SPAdes`](https://github.com/ablab/spades))
4. Reference-guided assembly ([`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/)), ([`SAMtools`](https://www.htslib.org/)) and ([`BCFtools`](https://samtools.github.io/bcftools/))
5. Annotation ([`Bakta`](https://bakta.readthedocs.io/)) and ([`Prodigal`](https://github.com/hyattpd/Prodigal))
6. Custom python scripts using Biopython, openpyxl, bioservices

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. 

First, prepare a samplesheet with your input data that looks as follows:
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,RXN_ID,max_min
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,R02738,min
```
Each row represents one biological sample.
The columns are defined as follows:
- `sample`: Unique sample identifier
- `fastq_1`: FASTQ file for single-end reads or read 1 for paired-end data
- `fastq_2`: FASTQ file for read 2 (leave empty for single-end data)
- `RXN_ID`: KEGG reaction ID to be optimized in the metabolic model
- `max_min`: Specify whether the reaction should be maximized or minimized
> [!NOTE]
> An optional reference genome (ref_genome) can be added as an additional column when performing reference-guided assembly.

Now, you can run the pipeline using:

```bash
nextflow run iumobg/model_creation \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

iumobg/model_creation was originally written by Doga Yasemen Testere, İrem Ay.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use iumobg/model_creation for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
