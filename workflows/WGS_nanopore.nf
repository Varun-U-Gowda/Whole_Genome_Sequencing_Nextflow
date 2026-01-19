#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------- Params ----------
params.reads 	= "data/*.{fastq,fastq.gz,fq,fq.gz}"
params.outdir 	= params.outdir ?: "Results"
assembly_outdir = "${params.outdir}/Assembly"
quast_outdir  	= "${params.outdir}/Quast"
prokka_outdir 	= "${params.outdir}/Prokka"

// Genome size - set based on Quick ID results
params.genome_size     = "5.3m"  // Adjust per species if known

params.kingdom = "Bacteria"
params.genus   = "Escherichia"
params.species = null  // Optional: set if known from Quick ID
params.prefix  = "my_genome"

//Download reference genome
process DownloadReference {
  tag "$accession"
  storeDir "References"  // Cache downloaded references
  
  input:
    val accession
  
  output:
    path "*.fna", emit: reference
  
  when:
    accession != null
  
  script:
  """
  echo ">>> Downloading reference: ${accession}"
  
  datasets download genome accession ${accession} \
    --filename ${accession}.zip
  
  unzip ${accession}.zip
  
  # Find and rename the genomic.fna file
  GENOMIC=\$(find ncbi_dataset/data -name "*_genomic.fna")
  cp \$GENOMIC ${accession}_genomic.fna
  
  echo ">>> Reference downloaded: ${accession}_genomic.fna"
  """
}

// Assembly with Flye
process Flye {
  tag "$sample"
  publishDir assembly_outdir, mode: 'copy'

  input:
    tuple val(sample), path(reads)

  output:
    tuple val(sample),
          path("flye_${sample}", type: 'dir'),
          path("flye_${sample}/assembly.fasta")

  script:
  """
  echo ">>> Assembling ${sample}"
  echo ">>> Genome size: ${params.genome_size}"
  
  flye --nano-raw $reads \\
       --out-dir flye_${sample} \\
       --genome-size ${params.genome_size} \\
       --asm-coverage 100 \\
       --threads ${task.cpus}
  
  echo ">>> Assembly complete"
  """
}

// Quality assessment with optional reference
process Quast {
  tag "$sample"
  publishDir quast_outdir, mode: 'copy', pattern: "quast_*"

  input:
    tuple val(sample), path(flye_dir), path(assembly), path(reference)
  
  output:
    tuple val(sample), path("quast_${sample}", type: 'dir')
  
  script:
  def ref_arg = reference.name != 'NO_REF' ? "-r ${reference} --fragmented" : ""
  """
  set -euo pipefail
  
  echo ">>> Running Quast on ${assembly}"
  
  quast.py "$assembly" \\
    ${ref_arg} \\
    -o quast_${sample} \\
    --threads ${task.cpus}
  
  echo ">>> Quast complete"
  """
}

// Annotation
process Prokka {
  tag "$sample"
  container 'staphb/prokka:1.14.6'
  publishDir prokka_outdir, mode: 'copy', pattern: "prokka_*"

  input:
    tuple val(sample), path(flye_dir), path(assembly)

  output:
    path "prokka_${sample}/*"

  script:
  def species_arg = params.species ? "--species ${params.species}" : ""
  """
  set -euo pipefail
  prokka "$assembly" \\
    --outdir prokka_${sample} \\
    --prefix ${params.prefix}_${sample} \\
    --kingdom ${params.kingdom} \\
    --genus ${params.genus} \\
    ${species_arg} \\
    --usegenus
  """
}

workflow {
  reads_ch = Channel
    .fromPath(params.reads, checkIfExists: true)
    .map { file -> 
      def sample = file.simpleName
      tuple(sample, file)
    }
  
  // Download reference if accession provided
  if (params.reference_accession) {
    reference_ch = DownloadReference(Channel.value(params.reference_accession))
  } else {
    // Create dummy reference channel
    reference_ch = Channel.value(file('NO_REF'))
  }
  
  // Assembly
  assembled = Flye(reads_ch)
  
  // Add reference to assembled channel
  assembled_with_ref = assembled
    .map { sample, flye_dir, assembly -> tuple(sample, flye_dir, assembly) }
    .combine(reference_ch)
    .map { sample, flye_dir, assembly, ref -> tuple(sample, flye_dir, assembly, ref) }
  
  // Quast 
  Quast(assembled_with_ref)
  
  // Prokka
  Prokka(assembled)
}