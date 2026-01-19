---

# Whole Genome Sequencing Analysis Pipelines (Nextflow)

**Modular short-read and long-read genome assembly, annotation, and taxonomic profiling**

---

## Project Overview

This repository contains a collection of **reproducible Whole Genome Sequencing (WGS) workflows** implemented using **Nextflow DSL2**. The pipelines support:

* **Illumina short-read data** (single-end and paired-end)
* **Long-read data** from **Oxford Nanopore** and **PacBio**
* **Taxonomic profiling (Kraken2)**, **assembly**, **annotation**, and **QC**

The workflows were developed through hands-on bioinformatics research and pipeline development work at **Brigham and Women’s Hospital (BWH)**, with a focus on **real-world usability, modular design, and biological correctness** rather than single-dataset optimization.

---
## Why This Project Exists

WGS projects often require:

* Different handling for **short-read vs long-read** technologies
* Reproducible execution across **local machines, Conda, and Docker**
* Clear separation between **workflow logic** and **compute configuration**
* Transparency about **biological and data-driven limitations**

This project addresses those needs using:

* **Nextflow DSL2** for workflow orchestration
* **Technology-specific pipelines** instead of one “one-size-fits-all” workflow
* **Selective containerization** (e.g., Docker for Prokka only)
* Explicit documentation of **coverage-dependent behavior**

---

## What Makes This Project Different

* **Technology-aware design**
  Separate workflows are provided for Illumina, Nanopore, and PacBio data, reflecting the biological and computational differences between sequencing platforms.

* **Modular configuration strategy**
  Compute resources and execution logic are separated into:

  * `base.config` — shared defaults
  * `resources.config` — process-level resource tuning
  * `nextflow.config` — execution profiles (local / conda / docker)

* **Selective container usage**
  Only tools that benefit from containerization (e.g., **Prokka**) run in Docker; all other steps can run via Conda.

* **Real-world constraints acknowledged**
  Long-read assemblies are explicitly **coverage-dependent**. Low-coverage public datasets may fail to assemble complete genomes, and this behavior is documented rather than hidden.

* **Grounded in research experience**
  The pipeline design reflects hands-on experience with microbial genome assembly, annotation, and comparative genomics, including a **published WGS study**.

---

## Repository Structure

```
Whole_Genome_Sequencing_Nextflow/
├── environment_setup.md
├── README.md
├── conf/
│   ├── base.config                  # Shared defaults (CPUs, memory, time, outdir)
│   └── resources.config             # Process-specific resource tuning
├── data/                            # Example input structure (no real data)
│   ├── illumina_PE/
│   │   ├── sample_1.fastq.gz
│   │   └── sample_2.fastq.gz
│   ├── illumina_SE/
│   │   └── sample.fastq.gz
│   ├── nanopore/
│   │   └── sample.fastq.gz
│   ├── pacbio/
│   │   └── sample.fastq.gz
│   └── kraken_library_fasta/
│       └── example_genome.fna
├── nextflow.config                  # Root config with execution profiles
└── workflows/
    ├── Kraken2.nf                   # Taxonomic profiling
    ├── WGS_illumina_PE.nf            # Illumina paired-end pipeline
    ├── WGS_illumina_SE.nf            # Illumina single-end pipeline
    ├── WGS_nanopore.nf               # Oxford Nanopore pipeline
    └── WGS_pacBio.nf                 # PacBio pipeline

```

---

## Workflow Architecture

### Short-read (Illumina) workflows

1. Read QC (FastQC)
2. Read trimming (Trimmomatic)
3. Genome assembly (SPAdes)
4. Assembly QC (QUAST)
5. Genome annotation (Prokka)
6. Optional taxonomic classification (Kraken2)

### Long-read (Nanopore / PacBio) workflows

1. Taxonomic profiling (Kraken2)
2. User-defined genome size
3. Long-read assembly
4. Genome annotation (Prokka via Docker)
5. Assembly QC (QUAST)

Each workflow can be run **independently**, depending on sequencing technology and experimental design.

---

## How to Run

> **Note:** `nextflow.config` automatically loads `conf/base.config` and `conf/resources.config`.
> You only need to specify **inputs and high-level options** at runtime.

## Environment Setup

Detailed installation and dependency setup steps for macOS (Apple Silicon)
are provided in [`environment_setup.md`](environment_setup.md).

### Illumina Paired-End

```bash
nextflow run workflows/WGS_illumina_PE.nf \
  --reads "data/illumina_PE_data/*_{1,2}.fastq.gz" \
  --kraken_library_fasta "kraken_library_fasta" \
  --buildKrakenDB true \
  --kingdom "Bacteria" \
  --genus "Escherichia" \
  --outdir "Results" \
  -profile docker
```

### Illumina Single-End

```bash
nextflow run workflows/WGS_illumina_SE.nf \
  --reads "data/illumina_SE_data/*.fastq.gz" \
  --kraken_library_fasta "kraken_library_fasta" \
  --buildKrakenDB true \
  --kingdom "Bacteria" \
  --genus "Escherichia" \
  --outdir "Results" \
  -profile docker
```

### Kraken2 (Quick Species Identification Only)

```bash
nextflow run workflows/Kraken2.nf \
  --reads "data/*.fastq.gz" \
  --kraken_library_fasta "kraken_library_fasta" \
  --buildKrakenDB true
  --outdir "Results"
```

### PacBio Long-Read Assembly

```bash
nextflow run workflows/WGS_pacBio.nf \
  --reads "data/pacbio_data/*.fastq.gz" \
  --genome_size 4.6m \
  --kingdom "Bacteria" \
  --genus "Escherichia" \
  --species "coli" \
  --outdir "Results" \
  -profile docker
```

### Oxford Nanopore Assembly

```bash
nextflow run workflows/WGS_nanopore.nf \
  --reads "data/nanopore_data/*.fastq.gz" \
  --genome_size 4.6m \
  --kingdom "Bacteria" \
  --genus "Escherichia" \
  --species "coli" \
  --outdir "Results" \
  -profile docker
```

## Validation and Known Limitations

The workflows in this repository were validated by running them end-to-end on representative test datasets. The following behaviors were observed during validation and are documented to ensure transparent interpretation of results.

### Kraken2 Database Handling

- The pipeline executes successfully when `--buildKrakenDB true` is used.
  - On the first run, the Kraken2 database is built as expected.
  - On subsequent runs, Nextflow correctly detects the existing database and skips rebuilding, proceeding directly to classification.
- When running with `--buildKrakenDB false`, Kraken2 classification fails with:
despite the required `.k2d` files being present on disk.
- This behavior was consistently observed during validation for both short-read and long-read workflows.
- **Current validated usage:** run the workflows with `--buildKrakenDB true`.

### Assembly Quality and QUAST Behavior

- Assembly outcomes depend on sequencing depth and input data quality.
- For some test datasets, assemblies produced only short contigs.
- In these cases, QUAST exited with status 4, consistent with its default minimum contig length requirement (≥500 bp).
- This reflects the biological characteristics of the input data rather than a pipeline execution failure.

### Long-read Test Datasets

- Long-read workflows were validated using publicly available test datasets.
- For datasets with insufficient coverage, assemblies did not produce usable contigs, which is expected behavior for long-read assemblers under low-coverage conditions.

### Summary

- ✅ Workflows execute end-to-end under validated configurations  
- ⚠️ Kraken2 database reuse with `--buildKrakenDB false` requires refinement  
- ⚠️ Assembly and QC outcomes are dependent on sequencing depth and data quality  

These observations reflect real-world data constraints and are documented to guide correct usage and future improvements.

---

## Relation to Prior Work

This project is informed by:

* Bioinformatics research experience at **Brigham and Women’s Hospital**
* A published WGS study:

  * *Assembly, annotation, and comparative whole genome sequence of* **Fusarium verticillioides**
* Practical experience with Illumina, Nanopore, and PacBio microbial datasets

---

## Author

**Varun Gowda**
Bioinformatics | Genomics | Workflow Development

---