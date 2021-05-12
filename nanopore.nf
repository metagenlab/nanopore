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

    ========================================================================================
    File batching
    ========================================================================================
    Options:
        --ftype             input files fromat [fastq,fast5] (default=fast5)
        --watch             Set to true to watch a directory for files
        --batch             Number of files to process in batch (default=5000)
        --wait              Time in seconds to wait for fast5/q files (default=1)

    ========================================================================================
    Basecalling
    ========================================================================================
    Options:
        --basecall       	  Use guppy cpu/gpu basecalling (default=cpu)
        --model             Basecalling model for guppy (default=dna_r9.4.1_450bps_fast.cfg)
        --runid             fastq.gz name after basecalling

    ========================================================================================
    QC
    ========================================================================================
    Options:
        --qc       	  Path to summary table for quality control 

    ========================================================================================
    Taxonomy classification
    ========================================================================================
    Options:
        --tax       	  Taxonomic classifier of choice (default=centrifuge) 
    
    ========================================================================================
    Assembly
    ========================================================================================
    Options:
        --assembly      Assembler of choice (default=flye)
    
    ========================================================================================
    Map to reference genome, plot coverage
    ========================================================================================
    Options:
        --map          mapper of choice (default=minimap2)
        --reference    path to reference genome (default=$params.input/*.fna)



    Author:
    Farid Chaabane (faridchaabn@gmail.com)
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

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
if (params.watch) {
  input_files_ch=Channel.watchPath("$params.input/*.$params.ftype").
    until{file->file.name == fileName}
    }
else {
  input_files_ch=Channel.fromPath("$params.input/*.$params.ftype")
}

input_files_ch.into{fast5_ch; fastq_ch}

// porcess that creates a text file within a time limit. 
// The text file presence stops nextflow's watch for new files

if (params.watch) and (params.basecall){

process stop_watch {
script:
"""
sleep $params.wait
echo "spent $params.wait seconds watching fast5 dir" > $params.input/$fileName
"""
}
}



/*
================================================================
    Basecalling, adapter and barcode trimming
================================================================
*/

params.model='dna_r9.4.1_450bps_fast.cfg'
params.basecall=false
params.runid='runid_0'
params.trim=false

process guppy_basecalling  {
  container 'genomicpariscentre/guppy:4.5.3'

  cpus 10

  input: file fast5List from fast5_ch.buffer(size: params.batch, remainder: true)
  
  when: params.basecall
  
  output: file '*.fastq' into basecalled_fastq_ch, collect_fastq_ch
          file '*.txt' into basecall_summary_ch
  
  script:
  if (params.basecall == 'cpu')
      if (params.trim)
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model --cpu_threads_per_caller 1 \
        --num_callers ${task.cpus} --trim_strategy "dna" --trim_barcodes --compress_fastq
        """
      else
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model --cpu_threads_per_caller 1 \
        --num_callers ${task.cpus} 
        """

  else if (params.basecall == 'gpu')
     if (params.trim)
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model -x "cuda:0" \
        --num_callers ${task.cpus} --gpu_runners_per_device ${task.cpus} --trim_strategy "dna" --trim_barcodes 
        """
    else
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model -x "cuda:0" \
        --num_callers ${task.cpus} --gpu_runners_per_device ${task.cpus}
        """
}
basecall_summary_ch.collectFile(name: "${params.runid}_summary.txt", keepHeader:true, 
  skip:1, storeDir:"$params.outdir/basecall").set{merged_summary_ch}

//basecalled_fastq_ch.into{fastq_collect_ch,bfastq_ch}

collect_fastq_ch.collectFile(name: "${params.runid}.fastq", storeDir:"$params.outdir/basecall").set{merged_fastq_ch}


/*
================================================================
    Quality Control
================================================================
*/
params.qc=false

if (params.basecall) {
  summaries_ch=Channel.fromPath("$params.outdir/basecall/*summary.txt")
}
else{
  summaries_ch=Channel.fromPath("$params.input/*summary.txt")
}

summaries_ch.into{summary_qc_ch; summary_aln_ch}

process pyco_quality_control {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'

  publishDir "$params.outdir/QC", mode: 'symlink', pattern: "*.html"
  
  input: file summaryTxt from summary_qc_ch
  
  when: params.qc
  
  output: file "${params.runid}.html" into html_report_ch
  
  script:
  """
  pycoQC -f $summaryTxt -o ${params.runid}.html
  """
}

/*
================================================================
    Taxonomy classification
================================================================
*/
all_fastqs_ch=basecalled_fastq_ch.mix(fastq_ch)

all_fastqs_ch.into{fastq_tax_ch; fastq_assembly_ch; fastq_mapping_ch}

params.tax=false

params.meta=false

process centrifuge_fastqs {
  container 'quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5'

  cpus 20
  
  input: file fastqs from fastq_tax_ch
  
  when: params.tax == 'centrifuge'
  
  output: file "report.txt" into centrifuge_reports_ch
  
  script:
  """
  centrifuge -p ${task.cpus} -x $params.db -U $fastqs -S centrifuge_out.txt --report-file report.txt
  """
}
centrifuge_reports_ch.collectFile(name:"${params.runid}_centrifuge_report.txt", keepHeader:true, 
  skip:1, storeDir:"$params.outdir/taxonomy")


/*
================================================================
    Assembly
================================================================
*/
params.assembly=false

process assembly_with_flye {
  container 'quay.io/biocontainers/flye:2.8.3--py27h6a42192_1'
  
  cpus 30
  
  publishDir "$params.outdir/assembly", mode: 'symlink', pattern: "*"

  input: file fastqs from fastq_assembly_ch
  
  when: params.assembly == 'flye'  
  
  output: file "*.fasta" into fasta_ch
          file "*.gfa" into gfa_ch
          file "*.txt" into assembly_info_ch
  
  script:
  if (params.meta)
    """
    flye --nano-raw $fastqs --out-dir . --threads ${task.cpus} --meta
    """
  else
    """
    flye --nano-raw $fastqs --out-dir . --threads ${task.cpus}
    """
}

/*
================================================================
    Map to reference genome, plot coverage
================================================================
*/

params.map=false
params.reference="$params.input/*.fna"

process minimap2_reads_to_reference {
  container 'quay.io/biocontainers/minimap2:2.18--h5bf99c6_0'
  
  cpus 20
  
  when: params.map=='minimap2'

  publishDir "$params.outdir/minimap2", mode: 'symlink', pattern: "*"

  input: 
  file fasta from Channel.fromPath("$params.reference")
  file fastq from fastq_mapping_ch
  
  output: file "*.sam" into sam_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont $fasta $fastq > ${params.runid}.sam
  """
}

process sam_to_sorted_bam {
  container 'quay.io/biocontainers/samtools:0.1.19--h270b39a_9'
  
  cpus 10

  publishDir "$params.outdir/coverage", mode: 'copy', pattern: "*"

  input: file sam from sam_ch
  
  output: 
  file "*.bam" into bam_ch
  file "*.bai" into bai_ch

  script:
  """
  samtools view -@ ${task.cpus} -bS $sam | samtools sort -@ ${task.cpus} - sorted_aln 
  samtools index sorted_aln.bam
  """
}

process coverage_plot {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  publishDir "$params.outdir/coverage", mode: 'copy', pattern: "*.html"
  
  input: 
  file bam from bam_ch
  file bai from bai_ch
  file summary from summary_aln_ch

  output: file "*.html" into cov_ch
  
  script:
  """
  pycoQC --summary_file $summary -a $bam -o ${params.runid}_cov.html
  """
}

workflow.onComplete {
	  RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'

    log.info "nanopore-nf has finished."
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Time:     ${workflow.complete}"
    log.info "Duration: ${workflow.duration}\n"
}
