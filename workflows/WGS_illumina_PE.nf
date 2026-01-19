#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------- Params ----------
params.reads          = "data/*_{1,2}.fastq.gz"
params.outdir = params.outdir ?: "Results"
qc_outdir     = "${params.outdir}/QC"
trim_outdir   = "${params.outdir}/Trim_results"
spades_outdir = "${params.outdir}/Spades"
quast_outdir  = "${params.outdir}/Quast"
prokka_outdir = "${params.outdir}/Prokka"
kraken_outdir = "${params.outdir}/Kraken2"

params.kingdom        = "Bacteria"
params.genus          = "Escherichia"
params.prefix         = "my_genome"

params.kraken_library_fasta   = "genome_fasta"				
params.kraken_dbname  = "kraken_db"
params.buildKrakenDB  = true


// ---------- Processes ----------

process Fastqc {
  publishDir qc_outdir, mode: 'copy'

  input:
    path read

  output:
    path '*_fastqc.zip'
    path '*_fastqc.html'

  script:
  """
  fastqc -o ./ "$read"
  """
}

process Trim {
  tag "$sample"
  publishDir trim_outdir, mode: 'copy'

  input:
    tuple val(sample), path(read1), path(read2)

  output:
    tuple val(sample),
          path("${sample}/${sample}Trimmed_1P.fastq.gz"),
          path("${sample}/${sample}Trimmed_2P.fastq.gz")

  script:
  """
  set -euo pipefail
  mkdir -p ${sample}

  jar_path=\$(ls -1 "\${CONDA_PREFIX}"/share/trimmomatic*/trimmomatic*.jar 2>/dev/null | head -n1 || true)
  adapters_path=\$(ls -1 "\${CONDA_PREFIX}"/share/trimmomatic*/adapters/TruSeq3-PE.fa 2>/dev/null | head -n1 || true)

  if [[ -z "\$jar_path" ]]; then
    echo "ERROR: Trimmomatic JAR not found" >&2
    exit 1
  fi
  if [[ -z "\$adapters_path" ]]; then
    echo "ERROR: TruSeq3-PE.fa not found" >&2
    exit 1
  fi

  java -jar "\$jar_path" PE \\
    -threads ${task.cpus} \\
    "$read1" "$read2" \\
    ${sample}/${sample}Trimmed_1P.fastq.gz ${sample}/${sample}Trimmed_1U.fastq.gz \\
    ${sample}/${sample}Trimmed_2P.fastq.gz ${sample}/${sample}Trimmed_2U.fastq.gz \\
    ILLUMINACLIP:"\$adapters_path":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  """
}

process Spades {
  tag "$sample"
  publishDir spades_outdir, mode: 'copy'

  input:
    tuple val(sample), path(r1p), path(r2p)

  output:
    tuple val(sample),
          path("spades_${sample}", type: 'dir'),
          path("spades_${sample}/contigs.fasta")

  script:
  """
  set -euo pipefail
  outdir="spades_${sample}"
  spades.py -1 "$r1p" -2 "$r2p" \\
            -t ${task.cpus} -m ${task.memory.toGiga()} \\
            -o "\$outdir"
  """
}

process Quast {
  tag "$sample"
  publishDir quast_outdir, mode: 'copy'

  input:
    tuple val(sample), path(spades_dir), path(contigs)

  output:
    tuple val(sample), path("quast_${sample}", type: 'dir')

  script:
  """
  set -euo pipefail
  quast.py "$contigs" -o quast_${sample}
  """
}

process Prokka {
  tag "$sample"
  container 'staphb/prokka:1.14.6'  // Keep container for Prokka
  publishDir prokka_outdir, mode: 'copy'

  input:
    tuple val(sample), path(spades_dir), path(contigs)

  output:
    tuple val(sample), path("prokka_${sample}", type: 'dir')

  script:
  """
  set -euo pipefail
  prokka "$contigs" \\
    --outdir prokka_${sample} \\
    --prefix ${params.prefix}_${sample} \\
    --kingdom ${params.kingdom} \\
    --genus ${params.genus} \\
    --usegenus
  """
}

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
  kraken2-build --download-taxonomy --db "\$DB_DIR" --use-ftp

  echo ">>> Adding genomes"
  for genome in "${fasta_dir}"/*.fna "${fasta_dir}"/*.fa "${fasta_dir}"/*.fasta; do
    [ -e "\$genome" ] || continue
    echo "  - Adding \$(basename "\$genome")"
    kraken2-build --add-to-library "\$genome" --db "\$DB_DIR"
  done

  echo ">>> Building Kraken2 DB"
  kraken2-build --build --db "\$DB_DIR" --threads ${task.cpus}

  echo ">>> Cleaning intermediate files"
  kraken2-build --clean --db "\$DB_DIR"
  """
}

process Kraken2Classify {
  tag "$sample"
  publishDir kraken_outdir, mode: 'copy'

  input:
    tuple val(sample), val(dbdir), path(contigs)

  output:
    tuple val(sample), 
          path("${sample}_kraken_report.txt"),
          path("${sample}_kraken_output.txt")

  script:
  """
  set -euo pipefail

  echo ">>> Classifying ${contigs} with DB: ${dbdir}"

  kraken2 \\
    --db "${dbdir}" \\
    --threads ${task.cpus} \\
    --report ${sample}_kraken_report.txt \\
    --output ${sample}_kraken_output.txt \\
    ${contigs}
  """
}

workflow {
  // 1) Pair raw reads by sample - SIMPLIFIED
  paired = Channel
    .fromFilePairs(params.reads, flat: true, checkIfExists: true)

  // 2) FastQC on raw reads (per file)
  raw_files = paired.flatMap { sid, r1, r2 -> [r1, r2] }
  Fastqc(raw_files)

  // 3) Trim → emits paired-trimmed FASTQs (sample, r1p, r2p)
  trimmed = Trim(paired)

  // 4) SPAdes on paired-trimmed
  assembled = Spades(trimmed)

  // 5) QUAST on contigs
  Quast(assembled)

  // 6) Prokka on contigs
  Prokka(assembled)

  // 7) Kraken2 classification
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

  contigs_ch = assembled.map { sid, spdir, contigs -> tuple(sid, contigs) }

  db_contigs_ch = kraken_db_ch.combine(contigs_ch)
    .map { dbdir, sid, contigs -> tuple(sid, dbdir, contigs) }
  
  Kraken2Classify(db_contigs_ch)
}