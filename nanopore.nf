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
        --reference    path to reference genome
    
    ========================================================================================
    Resistance gene identifier
    ========================================================================================
    Options:
        --res       resistance gene tool (default=rgi)
        --card      path to CARD database
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

if (params.watch) {
  input_files_ch=Channel.watchPath("$params.fast5").
    until{file->file.name == fileName}
    }
else {
  input_files_ch=Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.fastq))}
}

input_files_ch.into{fast5_ch; fastq_ch}

// porcess that creates a text file within a time limit. 
// The text file presence stops nextflow's watch for new files

if (params.watch && params.basecall){

process stop_watch {
script:
"""
sleep $params.wait
echo "spent $params.wait seconds watching fast5 dir" > $params.fast5/$fileName
"""
}
}



/*
================================================================
    Basecalling, adapter and barcode trimming
================================================================
*/

process guppy_basecalling  {
  container 'genomicpariscentre/guppy:4.5.3'

  cpus 10

  when: params.basecall

  input: file fast5List from fast5_ch.buffer(size: params.batch, remainder: true)

  output: 
      file "fail/*.fastq" into failed_fastq_ch
      file "pass/*.fastq" into passed_fastq_ch
      file '*.txt' into basecall_summary_ch
      
  script:
  if (params.basecall == 'cpu')
      if (params.trim)
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model --cpu_threads_per_caller ${task.cpus} \
        --num_callers 1  --trim_strategy "dna" --trim_barcodes --compress_fastq
        """
      else
        """
        guppy_basecaller -i . --save_path . \
        --config $params.model --cpu_threads_per_caller ${task.cpus} \
        --num_callers 1
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

guppy_fastqs_ch=Channel.empty()
if (params.pass){
  guppy_fastqs_ch << passed_fastq_ch
} else {
  guppy_fastqs_ch << passed_fastq_ch.mix(failed_fastq_ch)
}

 guppy_fastqs_ch.into{basecalled_fastq_ch; collect_fastq_ch} 

if (params.collect) {
  collect_fastq_ch.collectFile(name: "${params.runid}.fastq", storeDir:"$params.outdir/basecall").set{merged_fastq_ch}
  }

all_fastqs_ch=basecalled_fastq_ch.mix(fastq_ch)

all_fastqs_ch.into{fastq_quality_ch; fastq_tax_ch; fastq_assembly_ch; fastq_mapping_ch; fastq_resistance_ch}

/*
================================================================
    Quality Control
================================================================
*/

summaries_ch = Channel.empty()
if (params.basecall) {
  summaries_ch=Channel.fromPath("$params.outdir/basecall/*summary.txt")
}
if (!params.basecall && params.qc=='pycoqc') {
  summaries_ch=Channel.fromPath("$params.summary")
}

summaries_ch.into{summary_qc_ch; summary_aln_ch}

process pyco_quality_control {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  when: params.qc=='pycoqc'
  
  publishDir "$params.outdir/qc", mode: 'symlink', pattern: "*"
  
  input: file summaryTxt from summary_qc_ch
  
  output: file "${params.runid}.html" into html_report_ch
  
  script:
  """
  pycoQC -f $summaryTxt -o ${params.runid}.html
  """
}

process nanoplot_qc {
  container 'quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0'
  
  tag "${sampleId}"

  cpus 10
  
  publishDir "$params.outdir/qc/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.qc=='nanoplot'

  input: input: set sampleId, file(fastq) from fastq_quality_ch

  output: tuple(sampleId, file("${sampleId}*")) into nanoplot_ch
  
  script:
  """
  NanoPlot -t ${task.cpus} -p $sampleId -o . --fastq $fastq
  """
}

/*
================================================================
    Taxonomy classification
================================================================
*/

process centrifuge_fastqs {
  container 'quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6'

  tag "${sampleId}"
  
  cpus 20
  
  publishDir "$params.outdir/centrifuge/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.tax == 'centrifuge'

  input: set sampleId, file(fastq) from fastq_tax_ch
  
  output: tuple(sampleId, file("${sampleId}*")) into centrifuge_reports_ch
  
  script:
  """
  centrifuge -p ${task.cpus} -x $params.cendb/$params.cenind -U $fastq -S ${sampleId}-out.txt --report-file ${sampleId}-report.txt
  """
}
if (params.collect){
  centrifuge_reports_ch.collectFile(name:"${params.runid}_centrifuge_report.txt", keepHeader:true, 
  skip:1, storeDir:"$params.outdir/taxonomy")
}



/*
================================================================
    Assembly
================================================================
*/

process assembly_with_flye {
  container 'quay.io/biocontainers/flye:2.8.3--py27h6a42192_1' 

  tag "${sampleId}"

  cpus 30

  publishDir "$params.outdir/flye/$sampleId", mode: 'copy', pattern: "*"

  when: params.assembly == 'flye' 

  input: set sampleId, file(fastq) from fastq_assembly_ch
  
  output: tuple(sampleId, file("assembly.fasta")) into fasta_ch
          tuple(sampleId, file("assembly_graph.gfa")) into gfa_ch
          tuple(sampleId, file("assembly_info.txt")) into assembly_info_ch
  
  script:
  if (params.meta)
    """
    flye --nano-raw $fastq --out-dir . --threads ${task.cpus} --meta
    """
  else
    """
    flye --nano-raw $fastq --out-dir . --threads ${task.cpus}
    """
}

assembled_fasta_ch = fasta_ch


/*
================================================================
    Map to reference genome, plot coverage
================================================================
*/

process minimap2_reads_to_reference {
  container 'quay.io/biocontainers/minimap2:2.18--h5bf99c6_0'
  
  tag "${sampleId}"

  cpus 20

  publishDir "$params.outdir/minimap2/$sampleId", mode: 'symlink', pattern: "*"

  when: params.map=='minimap2'

  input: 
  set sampleId, file(fastq) from fastq_mapping_ch
  
  output: tuple(sampleId, file("${sampleId}.sam")) into sam_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont ${params.reference} $fastq > ${sampleId}.sam
  """
}

process sam_to_sorted_bam {
  container 'quay.io/biocontainers/samtools:0.1.19--h270b39a_9'
  
  tag "${sampleId}"

  cpus 10

  publishDir "$params.outdir/coverage/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.map
  
  input: set sampleId, file(sam) from sam_ch
  
  output: 
  tuple(sampleId, file("${sampleId}.bam")) into bam_ch
  tuple(sampleId, file("${sampleId}.txt")) into map_info_ch
  tuple(sampleId, file("${sampleId}.bam.bai")) into bai_ch

  script:
  """
  samtools flagstat $sam > ${sampleId}.txt
  samtools view -@ ${task.cpus} -bS $sam | samtools sort -@ ${task.cpus} - ${sampleId} 
  samtools index ${sampleId}.bam
  """
}

all_bams_ch = bam_ch
all_bams_ch.into{pyco_bam_ch; bed_bam_ch}

process bam_to_bed{
  container 'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1'

  tag "${sampleId}"

  publishDir "$params.outdir/coverage/$sampleId", mode: 'copy', pattern: "*"

  when: params.map
  
  input: set sampleId, file(bam) from bed_bam_ch

  output: 
  tuple(sampleId, file("${sampleId}.tsv")) into cov_tsv_ch
  
  script:
  """
  bedtools genomecov -d -ibam $bam > ${sampleId}.tsv
  """
  }


process pycoQC_coverage_plot {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  publishDir "$params.outdir/coverage", mode: 'copy', pattern: "*.html"
  
  when: params.coverage=='pycoqc'
  
  input: 
  file bam from pyco_bam_ch
  file bai from bai_ch
  file summary from summary_aln_ch

  output: file "*.html" into cov_ch
  
  script:
  """
  pycoQC --summary_file $summary -a $bam -o ${params.runid}_cov.html
  """
}

/*
================================================================
    Resistance gene identifier
================================================================
*/

process rgi {
  container 'quay.io/biocontainers/rgi:5.2.0--pyhdfd78af_0'

  tag "${sampleId}"

  publishDir "$params.outdir/resistance/$sampleId", mode: 'copy', pattern: "*"

  cpus 20
  
  when params.res == 'rgi'

  input: set sampleId, file(fasta) from assembled_fasta_ch
  
  output: file "*.txt" into rgi_txt_ch
          file "*.json" into rgi_json_ch
          file "*.png" into rgi_heatmap_ch
  
  script:
  """
  rgi load -i $params.card/card.json --card_annotation $params.card/card_database_*.fasta \
  --wildcard_annotation $params.card/wildcard_database_*.fasta \
  --wildcard_index $params.card/wildcard/index-for-model-sequences.txt
  rgi main -i $fasta -o $params.runid -t contig -a BLAST -n ${task.cpus} --split_prodigal_jobs --clean
  rgi heatmap --input . --output ${sampleId}_heatmap.png
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
