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
} else {
  fast5_ch = Channel.empty()
  fastq_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.fastq))}
}

if (params.reference) {
  references_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.reference))}
} else {
    references_ch = Channel.empty()
} 

if (params.summary) {
  summaries_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.summary))}
} else {
  summaries_ch=Channel.empty()
}

if (params.hybrid) {
  illumina_fastq_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.r1), file(row.r2))}
} else {
  illumina_fastq_ch=Channel.empty()
}

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

  cpus params.cores

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

raw_fastqs_ch = basecalled_fastq_ch.mix(fastq_ch)
raw_fastqs_ch.into{filter_fastqs_ch; fastqs_ch}

/*
================================================================
    Quality filtering
================================================================
*/

process quality_reads_filtering {
  container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'

  conda 'bioconda::nanofilt=2.8.0'
  
  tag "${sampleId}"
  
  publishDir "$params.outdir/filtered/$sampleId", mode: 'copy', pattern: "*"

  when: params.filter

  input: set sampleId, file(fastq) from filter_fastqs_ch
  
  output: tuple(sampleId, file("${sampleId}.fastq")) into quality_fastqs_ch
  
  script:
  """
  NanoFilt $fastq -q $params.phred > ${sampleId}.fastq
  """
}

all_fastqs_ch = quality_fastqs_ch.mix(fastqs_ch)
all_fastqs_ch.into{fastq_qc_ch; fastq_tax_ch; fastq_assembly_ch; fastq_hybrid_assembly_ch; fastq_mapping_reads_ch; fastq_mapping_assembly_ch ; fastq_resistance_ch}


/*
================================================================
    Quality Control
================================================================
*/


summaries_ch.into{summary_qc_ch; summary_aln_ch}

process pyco_qc {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  conda 'bioconda::pycoqc=2.5.2'
  
  tag "${sampleId}"

  publishDir "$params.outdir/qc/$sampleId", mode: 'copy', pattern: "*"

  when: (params.qc=='pycoqc' && params.summary)
  
  input: set sampleId, file(summary) from summary_qc_ch
  
  output: tuple(sampleId, file("${sampleId}.html")) into html_report_ch
  
  script:
  """
  pycoQC -f $summary -o ${sampleId}.html
  """
}

process nanoplot_qc {
  container 'quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0'
  
  conda 'bioconda::nanoplot=1.38.0'
  
  tag "${sampleId}"

  cpus params.cores
  
  publishDir "$params.outdir/qc/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.qc=='nanoplot'

  input: set sampleId, file(fastq) from fastq_qc_ch

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

  conda 'bioconda::centrifuge=1.0.4_beta'

  tag "${sampleId}"
  
  cpus params.cores
  
  publishDir "$params.outdir/centrifuge/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.tax

  input: set sampleId, file(fastq) from fastq_tax_ch
  
  output: tuple(sampleId, file("${sampleId}*")) into centrifuge_reports_ch
  
  script:
  """
  centrifuge -p ${task.cpus} -x $params.cendb/$params.cenin -U $fastq -S ${sampleId}-out.txt
  centrifuge-kreport -x $params.cendb/$params.cenin ${sampleId}-out.txt > ${sampleId}-kreport.txt
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

  conda 'bioconda::flye=2.8.3'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/assembly/$sampleId", mode: 'copy', pattern: "*"

  when: params.assembly

  input: set sampleId, file(fastq) from fastq_assembly_ch
  
  output: 
  tuple(sampleId, file("assembly.fasta")) into fasta_ch
  tuple(sampleId, file("assembly_graph.gfa")) into gfa_ch
  tuple(sampleId, file("assembly_info.txt")) into assembly_info_ch
  
  script:
  if (params.meta)
    """
    flye --nano-raw $fastq -t ${task.cpus} --meta -o . 
    """
  else
    """
    flye --nano-raw $fastq -t ${task.cpus} -o .
    """
}




/*
================================================================
    Hybrid assembly
================================================================
*/

process hybrid_assembly_unicycler {
  container 'quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3'

  conda 'bioconda::unicycler=0.4.8'

  tag "${sampleId}"   

  cpus params.cores

  publishDir "$params.outdir/hybrid-assembly/$sampleId", mode: 'copy', pattern: "*"

  when: params.hybrid 

  input: 
  set sampleId, file(fastq) from fastq_hybrid_assembly_ch
  set sampleId, file(r1), file(r2) from illumina_fastq_ch
  
  output: 
  tuple(sampleId, file("assembly.fasta")) into hybrid_fasta_ch
  tuple(sampleId, file("assembly.gfa")) into hybrid_gfa_ch

  script:
  """
  unicycler -1 $r1 -2 $r2 -l $fastq -t ${task.cpus} -o .
  """
}


hybrid_fasta_ch.mix(fasta_ch).set{assembled_fasta_ch} 
assembled_fasta_ch.into{assemblies_resistance_ch; assemblies_map_ch; assemblies_prokka_ch}

/*
================================================================
    Mapping reads against assembly
================================================================
*/
process mapping_reads_against_assembly {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'
  
  conda 'bioconda::minimap2=2.20'
  
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"
  
  when: (params.map && params.assembly)
  
  input: 
  set sampleId, file(fastq) from fastq_mapping_assembly_ch
  set sampleId, file(fasta) from assemblies_map_ch
  
  output: 
  tuple(sampleId, file("${sampleId}-assembly.bam")) into bam_assembly_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont $fasta $fastq | samtools sort -@ ${task.cpus} -o ${sampleId}-assembly.bam
  """
}

/*
================================================================
    Map to reference genome, plot coverage
================================================================
*/

process minimap2_reads_to_reference {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'
  
  conda 'bioconda::minimap2=2.20'
  
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"

  when: (params.map && params.reference)

  input: 
  set sampleId, file(fastq) from fastq_mapping_reads_ch
  set sampleId, file(reference) from references_ch 

  output: 
  tuple(sampleId, file("${sampleId}-ref.bam")) into bam_ref_ch
  tuple(sampleId, file("${sampleId}-ref.bam.bai")) into bai_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont $reference $fastq | samtools sort -@ ${task.cpus} -o ${sampleId}-ref.bam
  samtools index ${sampleId}-ref.bam
  samtools view -c -F 260 ${sampleId}-ref.bam > ${sampleId}-ref.bam
  """
}

cov_bam_ch = bam_ref_ch
cov_bam_ch.into{pyco_bam_ch; bed_bam_ch}

process bam_to_bed{
  container 'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1'

  conda 'bioconda::bedtools=2.30.0'
  
  tag "${sampleId}"

  publishDir "$params.outdir/coverage/$sampleId", mode: 'copy', pattern: "*"

  when: (params.map && params.reference)
  
  input: set sampleId, file(bam) from bed_bam_ch

  output: tuple(sampleId, file("${bam.baseName}.tsv")) into cov_tsv_ch
  
  script:
  """
  bedtools genomecov -d -ibam $bam > ${bam.baseName}.tsv
  """
}


process pycoQC_coverage_plot {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  conda 'bioconda::pycoqc=2.5.2'
  
  tag "${sampleId}"

  publishDir "$params.outdir/coverage/$sampleId", mode: 'copy', pattern: "*.html"
  
  when: (params.coverage && params.reference && params.summary)
  
  input: 
  set sampleId, file(bam) from pyco_bam_ch
  set sampleId, file(bai) from bai_ch
  set sampleId, file(summary) from summary_aln_ch

  output: tuple(sampleId, file("*_cov.html")) into cov_ch
  
  script:
  """
  pycoQC --summary_file $summary -a $bam -o ${bam.baseName}_cov.html
  """
}

/*
================================================================
    Resistance gene identifier
================================================================
*/

process rgi {
  container 'docker-daemon:metagenlab/rgi:5.2.0-3.1.2'

  conda 'bioconda::rgi=5.2.0'
  
  tag "${sampleId}"

  publishDir "$params.outdir/resistance", mode: 'copy', pattern: "*"

  cpus params.cores
  
  when params.res

  input: set sampleId, file(fasta) from assemblies_resistance_ch
  
  output: 
  tuple(sampleId, file("${sampleId}.txt")) into rgi_txt_ch
  tuple(sampleId, file("${sampleId}.json")) into rgi_json_ch
  
  script:
  """
  rgi load -i $params.card/card.json --card_annotation $params.card/card_database_*.fasta \
  --wildcard_annotation $params.card/wildcard_database_*.fasta \
  --wildcard_index $params.card/wildcard/index-for-model-sequences.txt
  rgi main -i $fasta -o ${sampleId} -t contig -a BLAST -n ${task.cpus} --split_prodigal_jobs --clean
  """
}



process rgi_heatmap {
  container 'docker-daemon:metagenlab/rgi:5.2.0-3.1.2'

  conda 'bioconda::rgi=5.2.0'

  publishDir "$params.outdir/resistance", mode: 'copy', pattern: "*"

  input: file json from rgi_json_ch.collect()
  
  output: file "heatmap-all.png" into rgi_heatmap_ch
  
  script:
  """
  rgi heatmap -i . -o heatmap-all.png -cat drug_class -clus samples
  """
}

/*
================================================================
    Genome annotation
================================================================
*/

process prokka_annotation {
  container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'
  
  conda 'bioconda::prokka=1.14.6'
  
  tag "${sampleId}"

  publishDir "$params.outdir/annotation/$sampleId", mode: 'copy', pattern: "*"

  cpus params.cores
  
  when params.res

  input: set sampleId, file(fasta) from assemblies_prokka_ch
  
  output: tuple(sampleId, file("PROKKA/${sampleId}")) into prokka_out_ch
  
  script:
  """
  prokka --outdir PROKKA --prefix ${sampleId} $fasta
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
