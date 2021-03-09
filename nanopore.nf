#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
    nanopore-nf
    ================================================================
    DESCRIPTION
    Usage:
    nextflow run metagenlab/nanopore
    Options:
        --fast5       	    Input directory of fastq or fast5 files.
        --fastq             Input directory of fastq or fastq files.
        --outdir        	  Output directory for pipeline analysis results.
    Profiles:
        local               local execution
        slurm               SLURM execution with both singularity and Docker
    Author:
    Farid Chaabane (faridchaabn@gmail.com)
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

project_dir = projectDir

process fast5_batching {
  input: path f5path from params.fast5
  
  output: file '*.batch' into batches_channel
  
  script:
  """
  $project_dir/scripts/batcher.py -p $f5path -o . -f fast5 -w True  
  """
}

/*
fastqpath=Channel.fromPath(params.fastq)
process fastq_batching {
  input: fastqpath
  
  output:
  
  script:
  """
  
  """
}

process guppy_basecalling  {
  container 'genomicpariscentre/guppy:4.4.2'
  
  input: file batch from batches_channel
  
  output: fastq
  
  script:
  """
  
  """
}
*/

workflow.onComplete {
	  RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'

    log.info "nanopore-nf has finished."
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Time:     ${workflow.complete}"
    log.info "Duration: ${workflow.duration}\n"
}
