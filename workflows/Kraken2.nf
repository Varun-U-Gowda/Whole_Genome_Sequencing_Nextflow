#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------- Params ----------
params.reads         = "data/*.fastq.gz"  // Nanopore reads
params.outdir = params.outdir ?: "Results"
kraken_outdir = "${params.outdir}/QuickID"

params.kraken_library_fasta  = "genome_fasta"
params.kraken_dbname = "kraken_db"
params.buildKrakenDB = true  // Set false after db is built

// ---------- Processes ----------

process BuildKrakenDB {
  storeDir params.kraken_dbname

  input:
    val fasta_dir 

  output:
    path "db", type: 'dir', emit: db

  script:
  """
  set -euo pipefail

  DB_DIR="db"
  mkdir -p "\$DB_DIR"

  echo ">>> Downloading taxonomy"
  kraken2-build --download-taxonomy --db "\$DB_DIR"

  echo ">>> Adding Bacteroides genomes"
  for genome in "${fasta_dir}"/*.fna "${fasta_dir}"/*.fa "${fasta_dir}"/*.fasta; do
    [ -e "\$genome" ] || continue
    echo "  - Adding \$(basename "\$genome")"
    kraken2-build --add-to-library "\$genome" --db "\$DB_DIR"
  done

  echo ">>> Building Kraken2 DB"
  kraken2-build --build --db "\$DB_DIR" --threads ${task.cpus}

  echo ">>> Cleaning intermediate files"
  kraken2-build --clean --db "\$DB_DIR"
  
  echo ">>> Custom Bacteroides DB build complete"
  """
}

process Kraken2QuickID {
  tag "$sample"
  publishDir kraken_outdir, mode: 'copy'

  input:
    tuple val(sample), val(dbdir), path(reads)

  output:
    path "${sample}_quickID_report.txt"
    path "${sample}_quickID_output.txt"

  script:
  """
  set -euo pipefail

  echo ">>> Quick species identification: ${sample}"

  kraken2 \\
    --db "${dbdir}" \\
    --threads ${task.cpus} \\
    --report ${sample}_quickID_report.txt \\
    --output ${sample}_quickID_output.txt \\
    ${reads}
  
  echo ">>> Classification complete for ${sample}"
  echo ">>> Check ${sample}_quickID_report.txt to identify species"
  """
}

// ---------- Workflow ----------
workflow {
  // Read samples
  reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
    .map { file -> 
      def sample = file.baseName.replaceAll(/\.fastq\.gz$/, '')
      tuple(sample, file)
    }

  // Handle Kraken2 DB
  def kraken_db_ch
  
  if (params.buildKrakenDB as boolean) {
    kraken_db_ch = BuildKrakenDB(Channel.value(params.kraken_library_fasta)).db
  } else {
    if( !file("${params.kraken_dbname}/db/hash.k2d").exists()
     || !file("${params.kraken_dbname}/db/taxo.k2d").exists()
     || !file("${params.kraken_dbname}/db/opts.k2d").exists() ) {
      error "Custom Kraken2 DB missing at ${params.kraken_dbname}/db"
    }
    kraken_db_ch = Channel.value("${params.kraken_dbname}/db")
  }

  // Combine DB with reads for classification
  kraken_input = reads_ch
    .map { sample, reads -> tuple(sample, reads) }
    .combine(kraken_db_ch)
    .map { sample, reads, dbdir -> tuple(sample, dbdir, reads) }
  
  // Run Quick ID
  Kraken2QuickID(kraken_input)
}
