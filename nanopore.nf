#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:
    nextflow run metagenlab/nanopore [Profiles] [Options]
    
    Profiles:
        local               local execution
        slurm               SLURM execution with either Singularity or Docker
    
    Options:
        --outdir        	  Output directory for pipeline analysis results
        --input             Input directory of files to process

    ================================================================
    File batching
    ================================================================
    Options:
        --ftype       	    file type [fastq,fast5] (default=fast5)
        --watch             Set to true to watch a directory for files
        --batch             Number of files to process in batch (default=5000)
        --wait              Time in seconds to wait for fast5/q files (default=1)
    ================================================================
    Basecalling
    ================================================================
    Options:
        --basecall       	  Use guppy cpu/gpu basecalling (default=cpu)
        --model             Basecalling model for guppy (default=dna_r9.4.1_450bps_hac.cfg)
        --runid             fastq.gz name after basecalling            

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

/*
================================================================
    File batching
================================================================
*/

params.ftype='fast5'
params.watch=false
params.batch=5000
params.wait=1
fileName='done.fast5'

if (params.watch){
  process stop_watch { // porcess that creates a text file within a time limit. 
                      // The text file presence stops nextflow's watch for new files
    script:
    """
    sleep $params.wait
    echo "spent $params.wait seconds watching $params.ftype dir" > $params.input/$fileName
    """
    }
  
  Channel.watchPath("$params.input/*.$params.ftype").until{file->file.name == fileName}.set{files_ch}
}

else{
  files_ch = Channel.fromPath("$params.input/*.$params.ftype")
  }

/*
================================================================
    Basecalling
================================================================
*/

params.model='dna_r9.4.1_450bps_hac.cfg'
params.basecall=false
params.runid='runid_0'

process guppy_basecalling  {
  container 'genomicpariscentre/guppy:4.4.2'
  
  cpus = 2

  input: file fast5List from files_ch.buffer(size: params.batch, remainder: true)
  
  output: file '*.fastq.gz' into basecalled_fastq_ch
          file '*.txt' into basecall_summary_ch
  
  when: params.basecall

  script:
  if (params.basecall == 'cpu')
    """
    guppy_basecaller -i . --save_path . \
    --config $params.model --cpu_threads_per_caller ${task.cpus} \
    --num_callers 1 --compress_fastq
    """
  else if (params.basecall == 'gpu')
    """
    guppy_basecaller -i . --save_path . \
    --config $params.model --num_callers 1 -x "cuda:0" --compress_fastq 
    """
}

basecall_summary_ch.collectFile(name: "summary.txt", keepHeader:true, skip:1, storeDir:"$params.outdir/basecall")
basecalled_fastq_ch.collectFile(name: "$params.runid\.fastq.gz", storeDir:"$params.outdir/basecall")


workflow.onComplete {
	  RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'

    log.info "nanopore-nf has finished."
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Time:     ${workflow.complete}"
    log.info "Duration: ${workflow.duration}\n"
}
