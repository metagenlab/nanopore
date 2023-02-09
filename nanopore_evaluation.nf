#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:
    nextflow run metagenlab/nanopore [Options]
    
    Options:
        --outdir        	  Output directory for pipeline analysis results
        --input             Input directory of files to process


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
	--givenAssembly	notify if the assembly is already provided 
			(temporary: to save time during testing)
    
    ========================================================================================
    Mapping to Human
    ========================================================================================
    Options:
	--mapping	remove reads mapping to the human genome

    ========================================================================================
    Polishing
    ========================================================================================
    Options:
	--polish	Initiate all polishing (medaka, pepper, 
			homopolish, racon, pep_med, med_hom, pep_hom)
			(temporary: final version will only have the chosen polishing method)
	--qc_assembly	Quality check and comparison of assembly and all polished genomes 

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

    ========================================================================================
    Annotation 
    ========================================================================================
    Options:
        --annotation	Allow prokka annotation on every assembly/polished genome
			(currently used for comparing gene numbers)


Unused options that might be useful if someone implements it (originally from Farid)			

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
	
    Author:
    Farid Chaabane (faridchaabn@gmail.com)
    Alexandre Jann (alexandre.jann@unil.ch) 
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

/*
The general design used in this pipeline usually follows a structure where
a process outputs a single channel that is then copied into multiple subchannels :

input_channel_A

process X {
input A
output B}

output_channel_B.into{input_channel_B_for_process_Y;
                      input_channel_B_for_process_Z}
*/


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
  ref_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.reference))}
  anot_tmp_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, "reference",file(row.reference))}
  anot_tmp_ch.into{anot_reference_ch;res_reference_ch}
} else {
  ref_ch = Channel.empty()
  res_reference_ch=Channel.empty()
  anot_reference_ch=Channel.empty()
} 
ref_ch.into{references_ch;qc_reference_ch}

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

// DEBUG to avoid assembly during testing
if (params.givenAssembly) {
  given_assembly_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, file(row.assembly))}
  anot_tmp_assembly_ch=Channel.fromPath("$params.samples")
  .splitCsv(header:true)
  .map{ row-> tuple(row.sampleId, "flye",file(row.assembly))}
  anot_tmp_assembly_ch.into{anot_given_assembly_ch;res_assembly_ch;res_flye_ch}
} else {
  res_flye_ch=Channel.empty()
  res_assembly_ch=Channel.empty()
  given_assembly_ch=Channel.empty()
  anot_given_assembly_ch=Channel.empty()
}
given_assembly_ch.into{qc_flye_assembly;
	    medaka_assembly_ch;
	    fasta_ch;
            racon_mapping_assembly_ch;
            racon_assembly_ch;
            homo_assembly_ch;
            pepper_index_assembly_ch;
            pepper_assembly_ch;
            assembly_polish_ch;
            assembly_tax_ch;
	    checkm_fasta_ch;
            fasta_unclassified_ch}


Channel.fromPath(params.homopolish_db).into{bacteria_db_ch;bacteria_db_med_hom_ch;bacteria_db_pep_hom_ch}

// process that creates a text file within a time limit. 
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

  //conda 'bioconda::nanofilt=2.8.0'
  
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
all_fastqs_ch.into{medaka_fastq_ch; 
	racon_fastq_ch; 
	racon_mapping_fastq_ch; 
	pepper_fastq_ch; 
	pepper_index_fastq_ch; 
	fastq_filter_ch; 
	fastq_map_ch; 
	fastq_qc_ch; 
	fastq_tax_ch; 
	fastq_hybrid_assembly_ch; 
	fastq_mapping_reads_ch; 
	fastq_mapping_assembly_ch ; 
	fastq_resistance_ch;
	pep_med_fastq_ch}


/*
================================================================
    Quality Control
================================================================
*/

summaries_ch.into{summary_qc_ch; summary_aln_ch; summary_pycoQC_ch}

/*
This process enables sequencing data quality check of an entier 
sequencing run
*/

process pyco_qc {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  tag "${sampleId}"

  publishDir "$params.outdir/qc/$sampleId", mode: 'copy', pattern: "*"

  when: (params.qc=='pycoqc' && params.summary)
  
  input: set sampleId, file(summary) from summary_qc_ch
  
  output: tuple(sampleId, file("${sampleId}.html")) into html_report_ch
  output: tuple(sampleId, file("${sampleId}.json")) into json_report_ch
  
  script:
  """
  pycoQC -f $summary -o ${sampleId}.html -j ${sampleId}.json
  """
}

process nanoplot_qc {
  container 'quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0'
  
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
    Mapping to Human
================================================================
*/

/*
This process maps reads to the human genome (reference needs to
be set in the nextflow.config file)
*/
process minimap_to_human {
    
    container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'
  
    tag "${sampleId}"

    cpus params.cores
    publishDir "$params.outdir/blast/$sampleId", mode: 'copy', pattern: "*"
    input: set sampleId, file(fastq) from fastq_map_ch

    when: params.mapping

    output: tuple(sampleId, file("${sampleId}*")) into minimap2_aln_ch

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont '/media/IMU/GEN/PROJECTS/MINBC_19_015/Analysis/Kp_seq/GCA_000001405.15_GRCh38_full_analysis_set.fna/GCA_000001405.15_GRCh38_full_analysis_set.fna' $fastq > ${sampleId}-mapping.sam
    """
}

/*
This process creates a list of IDs of reads not mapped to the
 human genome
*/
process unmapped_to_human {

  container 'quay.io/staphb/samtools:1.14'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/samtools/$sampleId", mode: 'copy', pattern: "*"

  input: set sampleId, file(sam) from minimap2_aln_ch

  output: 
  tuple(sampleId, file("${sampleId}-unmapped.lst")) into unmapped_lst_ch

  script:
  """
  samtools view -f 4 $sam > ${sampleId}-unmapped.sam 
  cut -f 1 ${sampleId}-unmapped.sam | sort  | uniq > ${sampleId}-unmapped.lst
  """
}

/*
This process creates a new fastq file with all reads not mapped
to the human genome
*/
process filter_human_fastq {

  container 'quay.io/biocontainers/seqtk:1.3--h7132678_4'

  tag "${sampleId}"

  publishDir "$params.outdir/samtools/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(lst), file(fastq) from unmapped_lst_ch.join(fastq_filter_ch)

  output:
  tuple(sampleId, file("${sampleId}-unmapped.fastq")) into fastq_assembly_ch

  script:
  """
  seqtk subseq $fastq $lst > ${sampleId}-unmapped.fastq
  """
}


/*
================================================================
    Taxonomy classification
================================================================
*/

/*
This process performs species identification through Centrifuge
*/
process centrifuge_fastqs {
  container 'quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6'

  tag "${sampleId}"
  
  cpus params.cores
  
  publishDir "$params.outdir/centrifuge/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.tax

  input: set sampleId, file(fastq) from fastq_tax_ch
  
  output: tuple(sampleId, file("raw-reads*")) into centrifuge_reports_ch
  
  script:
  """
  centrifuge -p ${task.cpus} -x $params.cendb/$params.cenin -U $fastq -S raw-reads-out.txt
  centrifuge-kreport -x $params.cendb/$params.cenin raw-reads-out.txt > raw-reads.txt
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

/*
This process performs genome assembly through Flye

The --meta argument (for assembly from metagenomic reads) can be
used if necessary. It might be useful for polymicrobial samples
*/
process assembly_with_flye {
  container 'quay.io/biocontainers/flye:2.8.3--py27h6a42192_1' 

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/assembly/$sampleId", mode: 'copy', pattern: "*"

  when: params.assembly

  input: set sampleId, file(fastq) from fastq_assembly_ch
  
  output: 
  tuple(sampleId, file("assembly.fasta")) into flye_fasta_ch
  tuple(sampleId, name, file("assembly.fasta")) into out_anot_flye_ch 
  tuple(sampleId, file("assembly_graph.gfa")) into gfa_ch
  tuple(sampleId, file("assembly_info.txt")) into assembly_info_ch
  
  script:
  name="flye"
  if (params.meta)
    """
    flye --nano-raw $fastq -t ${task.cpus} --meta -o . 
    """
  else
    """
    flye --nano-raw $fastq -t ${task.cpus} -o .
    """
}
if(params.givenAssembly==false){
out_anot_flye_ch.into{anot_flye_assembly_ch;res_flye_ch}
flye_fasta_ch.into{fasta_ch;
            assembly_tax_ch;
            medaka_assembly_ch;
            racon_mapping_assembly_ch;
            racon_assembly_ch;
            homo_assembly_ch;
            pepper_index_assembly_ch;
            pepper_assembly_ch;
	    qc_flye_assembly;
            assembly_polish_ch;
	    checkm_fasta_ch;
            assemblies_map_ch}
} else {
anot_flye_assembly_ch=Channel.empty()
res_flye_ch=Channel.empty()
}

/*
This process performs species identification through Centrifuge
based on the genome assembly
*/
process centrifuge_assembly {
  container 'quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6'

  tag "${sampleId}"
  
  cpus params.cores
  
  publishDir "$params.outdir/centrifuge/$sampleId", mode: 'copy', pattern: "*"
  
  when: params.tax

  input: set sampleId, file(fasta) from assembly_tax_ch
  
  output: tuple(sampleId, file("flye-assembly.txt")) into centrifuge_assembly_reports_ch
  
  script:
  """
  centrifuge -p ${task.cpus} -x $params.cendb/$params.cenin -f $fasta -S flye-out.txt
  centrifuge-kreport -x $params.cendb/$params.cenin flye-out.txt > flye-assembly.txt
  """
}
if (params.collect){
  centrifuge_assembly_reports_ch.collectFile(name:"${params.runid}_centrifuge_assembly_report.txt", keepHeader:true, 
  skip:1, storeDir:"$params.outdir/taxonomy")
}

/*
================================================================
       Retrieve unclassified reads (multiplexed samples  only)
================================================================
*/

if (params.keep_unclassified){
fastq_unclassified_ch=Channel.fromPath("/media/IMU/GEN/PROJECTS/MINBC_19_015/SeqData/multiplex_SA_KP_EC_EF_AN_N_20220718_2nd/all_reads_unclassified.fastq")
fastq_unclassified_ch.into{fastq_unclass_mapping_ch;fastq_unclass_attributed_ch}

process mapping_unclassified_to_assembly {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'
  tag "${sampleId}"
  cpus params.cores
  publishDir "$params.outdir/unassigned/$sampleId", mode: 'copy', pattern: "*"
  when: params.keep_unclassified
  input: set sampleId,file(assembly),file(unclassified_fastq) from fasta_unclassified_ch.combine(fastq_unclass_mapping_ch)
  output: tuple(sampleId, file("*.sam")) into unclass_to_assembly_ch
  
  script:
  """
  minimap2 -t  ${task.cpus} -ax map-ont $assembly $unclassified_fastq > ${sampleId}_mapped.sam
  """
}

process create_read_list {
  container 'quay.io/staphb/samtools:1.14'
  tag "${sampleId}"
  cpus params.cores
  publishDir "$params.outdir/unassigned/$sampleId", mode: 'copy', pattern: "*"
  input: set sampleId, file(sam) from unclass_to_assembly_ch
  output: 
  tuple(sampleId, file("*.lst")) into mapped_unclass_list_ch
  file("*.lst") into mapped_list_ch
 
  script:
  """
  samtools view -F 4 $sam | cut -f 1 | sort  | uniq > ${sam.baseName}_mapped.lst
  """
}

// Be sure to handle this for any number of multiplexed samples (input section)
process merge_read_lists {
  cpus params.cores
  publishDir "$params.outdir/unassigned/", mode: 'copy', pattern: "duplicated.lst"
  input: file ('*') from mapped_list_ch.collect()
  output: file("duplicated.lst") into unclass_dup_ch
  
  script:
  """
  cat *.lst | sort | uniq -d > duplicated.lst
  #cat *.fastq > merged.fastq
  #awk -F' ' '{if (NR % 4 == 1) {print}}' merged.fastq | cut -c 2-37 > merged.lst
  #sort merged.lst | uniq -d > duplicated.lst
  """
}

process get_attributed_reads {
  container 'quay.io/biocontainers/seqtk:1.3--h7132678_4'
  tag "${sampleId}"
  cpus params.cores
  publishDir "$params.outdir/unassigned/$sampleId", mode: 'copy', pattern: "*.fastq"
  input: 
  set sampleId, file(mapped_lst),file(duplicated_lst), file(unclassified_fastq) from mapped_unclass_list_ch.combine(unclass_dup_ch).combine(fastq_unclass_attributed_ch)
//  file(duplicated_lst) from unclass_dup_ch 
//  file(unclassified_fastq) from fastq_unclass_attributed_ch
  output: tuple(sampleId, file("*.fastq")) into unclass_saved_fastq_ch
  
  script:
  """
  grep -v $mapped_lst -f $duplicated_lst > notDup.lst
  seqtk subseq $unclassified_fastq notDup.lst> attributed_no_dup.fastq
  """
}
}
/*
================================================================
    Polishing
================================================================
*/
if(params.polish){
/*
This process performs genome polishing on Flye's assembly with
Medaka
*/
process medaka_polishing {
  container 'quay.io/biocontainers/medaka:1.6.0--py38h84d2cc8_0'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/medaka/$sampleId", mode: 'copy', pattern: "*"

  input: 
  set sampleId, file(assembly), file(fastq) from medaka_assembly_ch.join(medaka_fastq_ch)

  output: 
  tuple(sampleId, file("medaka.fasta")) into out_medaka_ch
  tuple(sampleId, name, file("medaka.fasta")) into out_anot_medaka_ch


  script:
  name="medaka"
  """
  medaka_consensus -i $fastq -d $assembly -o . -t ${task.cpus} -m r941_min_high_g303
  mv consensus.fasta medaka.fasta
  """
}

out_medaka_ch.into{qc_medaka_ch;med_homo_assembly_ch}
out_anot_medaka_ch.into{;res_medaka_ch;anot_medaka_ch}

/*
This process performs a mapping of reads on Flye's assembly, 
a mandatory step to prepare for Racon polishing
*/
process racon_minimap {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/polishing/minimap2/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly), file(fastq) from racon_mapping_assembly_ch.join(racon_mapping_fastq_ch) 

  output:
  tuple(sampleId, file("${sampleId}-ONT*")) into racon_mapping_ch

  script:
  """
  minimap2 -t ${task.cpus} -d ${sampleId}-assembly.mmi $assembly
  minimap2 -t ${task.cpus} -a ${sampleId}-assembly.mmi $fastq > ${sampleId}-ONT-to-assembly-sorted.sam
  """
}

/*
This process performs Racon polishing
*/
process racon_polishing {
  container 'quay.io/biocontainers/racon:1.5.0--h7ff8a90_0'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/racon/$sampleId", mode: 'copy', pattern: "*"
  
  input:
  set sampleId, file(assembly),file(fastq), file(index) from racon_assembly_ch.join(racon_fastq_ch.join(racon_mapping_ch))

  output:
  tuple(sampleId, file("racon.fasta")) into out_racon_ch
  tuple(sampleId, name, file("racon.fasta")) into out_anot_racon_ch

  script:
  name="racon"
  """
  racon -t ${task.cpus} $fastq $index $assembly > ${sampleId}-racon-assembly.fasta
  mv ${sampleId}-racon-assembly.fasta racon.fasta
  """
}

out_racon_ch.set{qc_racon_ch}
out_anot_racon_ch.into{anot_racon_ch;res_racon_ch}

/*
This process performs genome polishing by Homopolish on Medaka 
polished assembly

It doesn't allow to use the --genus argument
*/
process homopolish_polishing {
  //container 'quay.io/biocontainers/homopolish:0.3.3--pyh5e36f6f_0'
  conda 'bioconda::homopolish=0.3.3' 

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/homopolish/$sampleId", mode: 'copy', pattern: "*"

  input:
    set sampleId, file(assembly),file(db) from homo_assembly_ch.combine(bacteria_db_ch)

  output:
  tuple(sampleId, file("homopolish.fasta")) into out_hom_ch
  tuple(sampleId, name, file("homopolish.fasta")) into out_anot_homo_ch

  script:
  name="homopolish"
  """
  homopolish polish -t ${task.cpus} -a $assembly -s $db -m R9.4.pkl -o .
  mv assembly_homopolished.fasta homopolish.fasta
  """
}

out_hom_ch.set{qc_homo_ch}
out_anot_homo_ch.into{anot_homo_ch;res_hom_ch}

/*
This process performs a mapping of reads on Flye's assembly,
a mandatory step to prepare for Pepper polishing
*/
process PEPPER_mapping {

  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/polish/pepper/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly), file(fastq) from pepper_index_assembly_ch.join(pepper_fastq_ch)

  output:
  tuple(sampleId, file("${sampleId}-pepper-ref.bam")) into samtools_pepper_ch

  """
  minimap2 -t ${task.cpus} -ax map-ont $assembly $fastq > ${sampleId}-pepper-ref.bam
  """
}

/*
This process performs indexing on from the mapping of
the process PEPPER_mapping
*/
process PEPPER_indexing {

  container 'quay.io/staphb/samtools:1.14'

  tag "${sampleId}"

  cpus params.cores

  //publishDir "$params.outdir/samtools/$sampleId", mode: 'copy', pattern: "*"

  input: set sampleId, file(bam) from samtools_pepper_ch

  output: 
  tuple(sampleId, file("${sampleId}-pepper-ref.bam")) into bam_pepper_ch
  tuple(sampleId, file("${sampleId}-pepper-ref.bam.bai")) into bai_pepper_ch

  script:
  """
  samtools sort -@ ${task.cpus} -o ${sampleId}-pepper-ref.bam $bam
  samtools index ${sampleId}-pepper-ref.bam
  """
}

pepper_db_ch=Channel.fromPath(params.pepper_db)
//pepper_db_ch=Channel.fromPath("database/pepper/pepper_r941_guppy305_microbial.pkl")

process PEPPER_polishing {

  container 'docker://kishwars/pepper_deepvariant:r0.8'  

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/PEPPER/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly), file(bam), file(bai),file(db) from pepper_assembly_ch.join(bam_pepper_ch.join(bai_pepper_ch)).combine(pepper_db_ch)

  output:
  tuple(sampleId, file("pepper.fasta")) into out_pepper_ch
  tuple(sampleId, name, file("pepper.fasta")) into out_anot_pepper_ch

  script:
  name="pepper"
  """
  pepper polish -t ${task.cpus} -b $bam -f $assembly -m $db -o .
  mv _pepper_polished.fa pepper.fasta
  """
}

out_pepper_ch.into{qc_pepper_ch;pep_homo_assembly_ch ;pep_med_assembly_ch}
out_anot_pepper_ch.into{anot_pepper_ch;res_pepper_ch}
/*
This process performs genome polishing by Homopolish  on Medaka 
polished assembly

It doesn't allow to use the --genus argument
*/
process homopolish_polishing_from_medaka {
  //container 'quay.io/biocontainers/homopolish:0.3.3--pyh5e36f6f_0'
  conda 'bioconda::homopolish=0.3.3' 
  
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/med_hom/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly), file(db) from med_homo_assembly_ch.combine(bacteria_db_med_hom_ch)

  output:
  tuple(sampleId, file("med_hom.fasta")) into out_med_hom_ch
  tuple(sampleId, name, file("med_hom.fasta")) into out_anot_med_homo_ch

  script:
  name="med_hom"
  """
  homopolish polish  -t ${task.cpus} -a $assembly -s $db -m R9.4.pkl -o .
  mv medaka_homopolished.fasta med_hom.fasta
  """
}

out_med_hom_ch.into{qc_med_homo_ch;rac_med_hom_cov}
out_anot_med_homo_ch.into{anot_med_homo_ch;res_med_hom_ch}

/*
This process performs genome polishing by Homopolish on Pepper
polished assembly

It doesn't allow to use the --genus argument
*/
process homopolish_polishing_from_pepper {
  //container 'quay.io/biocontainers/homopolish:0.3.3--pyh5e36f6f_0'

  conda 'bioconda::homopolish=0.3.3' 
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/pep_hom/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly),file(db) from pep_homo_assembly_ch.combine(bacteria_db_pep_hom_ch)
  //file(db) from bacteria_db_pep_hom_ch

  output:
  //tuple(sampleId, file("*homopolished.fasta")) into out_pep_hom_ch
  tuple(sampleId, file("pep_hom.fasta")) into out_pep_hom_ch
  tuple(sampleId, name, file("pep_hom.fasta")) into out_anot_pep_homo_ch

  script:
  name="pep_hom"
  """
  homopolish polish -t ${task.cpus} -a $assembly -s $db -m R9.4.pkl -o .
  mv pepper_homopolished.fasta pep_hom.fasta
  """
}

out_pep_hom_ch.set{qc_pep_homo_ch}
out_anot_pep_homo_ch.into{anot_pep_homo_ch;res_pep_hom_ch}

/*
This process allows for genome polishing by Medaka on Pepper
polished assembly
*/
process medaka_polishing_from_pepper {
  container 'quay.io/biocontainers/medaka:1.6.0--py38h84d2cc8_0'

  tag "${sampleId}"

  cpus params.cores
  publishDir "$params.outdir/pep_med/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(assembly),file(fastq) from pep_med_assembly_ch.join(pep_med_fastq_ch)

  output:
  tuple(sampleId,file("pep_med.fasta")) into out_pep_med_ch
  tuple(sampleId,name,file("pep_med.fasta")) into out_anot_pep_med_ch

  script:
  name="pep_med"
  """
  medaka_consensus -t ${task.cpus} -i $fastq -d $assembly -o . -m r941_min_high_g303
  mv consensus.fasta pep_med.fasta
  """
}

out_pep_med_ch.set{qc_pep_med_ch}
out_anot_pep_med_ch.into{anot_pep_med_ch;res_pep_med_ch}
}

/*
================================================================
    Quality Check Assembly
================================================================
*/

if(params.qc_assembly){

/*
This process performs a comparison between the different genome
assemblies and polishing

It can be used to compare the different polishing tools
*/
process quality_check_assembly {
  container 'quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7'

  tag "${sampleId}"
  
  publishDir "$params.outdir/QC_Assembly/$sampleId", mode: 'copy', pattern: "*"

  when: params.qc_assembly

  input:
  set sampleId, file(flye),file(racon) , file(homo), file(medaka), file(pepper) , file(pep_med), file(med_hom), file(pep_hom), file(ref) from qc_flye_assembly.join(qc_racon_ch).join(qc_homo_ch).join(qc_medaka_ch).join(qc_pepper_ch).join(qc_pep_med_ch).join(qc_med_homo_ch).join(qc_pep_homo_ch).join(qc_reference_ch)

  output: 
  tuple(sampleId, file("${sampleId}_quast*")) into quast_results
 
  when: (params.polish && params.reference)
 
  script:
  """
  quast $flye $racon $homo $medaka $pepper $pep_med $med_hom $pep_hom \
        -l "flye,racon,homopolish,medaka,pepper,pep_med,med_hom,pep_hom" \
        -r $ref \
        -o ${sampleId}_quast
  """
}
}

/*
================================================================
    Mapping reads against assembly
================================================================
*/

/*
This process maps reads to the genome assembled by Flye
*/
process mapping_reads_against_assembly {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'
  
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"
  
  when: (params.map && params.assembly)
  
  input: 
  set sampleId, file(fastq),file(fasta) from fastq_mapping_assembly_ch.join(assemblies_map_ch)
  
  output: 
  tuple(sampleId, file("${sampleId}-assembly.paf")) into mapping_assembly_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont $fasta $fastq > ${sampleId}-assembly.paf
  """
}

/*
This process creates a .bam file from the previous mapping
*/
process mapping_reads_against_assembly_samtools {
  container 'quay.io-staphb-samtools-1.14'
  
  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"
  
  when: (params.map && params.assembly)
  
  input: 
  set sampleId, file(paf) from mapping_assembly_ch
  
  output: 
  tuple(sampleId, file("${sampleId}-assembly.bam")) into bam_assembly_ch
  
  script:
  """
  samtools sort -@ ${task.cpus} -o ${sampleId}-assembly.bam $paf
  """
}
/*
================================================================
    Map to reference genome, plot coverage
================================================================
*/

/*
This process maps reads to a reference genome given as argument
*/
process minimap2_reads_to_reference {
  container 'quay.io/biocontainers/minimap2:2.20--h5bf99c6_0'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"

  when: (params.map && params.reference)

  input: 
  set sampleId, file(fastq), file(reference) from fastq_mapping_reads_ch.join(rac_med_hom_cov)
  
  output: 
  tuple(sampleId, file("${sampleId}-paf.sam")) into paf_ref_ch
  
  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont $reference $fastq > ${sampleId}-paf.sam
  """
}

/*
This process creates a .bam file from the previous mapping
*/
process samtools_reads_to_reference {
  container 'quay.io-staphb-samtools-1.14'

  tag "${sampleId}"

  cpus params.cores

  publishDir "$params.outdir/mapping/$sampleId", mode: 'copy', pattern: "*"

  when: (params.map && params.reference)

  input: 
  set sampleId, file(paf) from paf_ref_ch

  output: 
  tuple(sampleId, file("${sampleId}-ref.bam")) into bam_ref_ch
  tuple(sampleId, file("${sampleId}-ref.bam.bai")) into bai_ch
  
  script:
  """
  head $paf
  samtools sort -@ ${task.cpus} -o ${sampleId}-ref.bam $paf
  samtools index ${sampleId}-ref.bam
  #samtools view -c -F 260 ${sampleId}-ref.bam > ${sampleId}-ref.bam
  """
}

cov_bam_ch = bam_ref_ch
cov_bam_ch.into{pyco_bam_ch; bed_bam_ch}

// This process converts the .bam file to the .bed format
process bam_to_bed{
  container 'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1'

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

/*
This process uses pycoQC to allow for a visual representation
of the coverage of the genome
*/
process pycoQC_coverage_plot {
  container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
  
  //conda 'bioconda::pycoqc=2.5.2'
  
  tag "${sampleId}"

  publishDir "$params.outdir/coverage/$sampleId", mode: 'copy', pattern: "*.html"
  
  when: (params.coverage && params.reference && params.summary)
  
  input: 
  set sampleId, file(bam), file(bai), file(summary) from pyco_bam_ch.join(bai_ch).join(summary_aln_ch)

  output: 
  tuple(sampleId, file("*_cov.html")) into cov_ch
  tuple(sampleId, file("*_cov.json")) into cov_json_ch
  
  script:
  """
  pycoQC --summary_file $summary -a $bam -o ${sampleId}_cov.html -j ${bam.baseName}_cov.json
  """
}

/*
================================================================
    Resistance gene identifier
================================================================
*/
if(params.res){
res_rgi_ch=res_reference_ch.mix(res_assembly_ch,
	res_flye_ch,
	res_hom_ch,
	res_medaka_ch,
	res_pepper_ch,
	res_med_hom_ch,
	res_pep_hom_ch,
	res_pep_med_ch,
	res_racon_ch)

/*
This process perfoms resistance gene identification with RGI
it doesn't need to have a direct reference to the CARD database
*/
process rgi {
  container 'quay.io/biocontainers/rgi:5.2.1--pyha8f3691_2'
  containerOptions '--bind /data/databases'

  tag "${fasta.baseName}"

  publishDir "$params.outdir/resistance/$sampleId", mode: 'copy', pattern: "*"

  cpus params.cores
  
  when params.res

  input: set sampleId,name,file(fasta) from res_rgi_ch
  
  output: 
  tuple(sampleId, file("${fasta.baseName}.txt")) into rgi_txt_ch
  tuple(sampleId, file("${fasta.baseName}.json")) into rgi_json_ch
  
  script:
  """
  #rgi load -i $params.card/card.json --card_annotation $params.card/card_database_*.fasta \
  #--wildcard_annotation $params.card/wildcard_database_*.fasta \
  #--wildcard_index $params.card/wildcard/index-for-model-sequences.txt
  rgi main -i $fasta -o ${fasta.baseName} -t contig -a BLAST -n ${task.cpus} --split_prodigal_jobs --clean
  """
}


/*
This process creates a heatmap of resistance genes for a 
visualisation between polishing tools
*/
process rgi_heatmap {
  container 'quay.io/biocontainers/rgi:5.2.1--pyha8f3691_2'

  publishDir "$params.outdir/resistance/$sampleId", mode: 'copy', pattern: "*"

  input: 
  set sampleId,file(json) from rgi_json_ch.groupTuple()
  
  output: 
  tuple(sampleId,file("heatmap*")) into rgi_heatmap_ch
  
  script:
  """
  rgi heatmap -i . -o heatmap
  """
}

}
/*
================================================================
    Genome annotation
================================================================
*/
if(params.annotation){
anot_prokka_ch=anot_reference_ch.mix(anot_flye_assembly_ch,
	anot_given_assembly_ch,
	anot_medaka_ch,
	anot_homo_ch,
	anot_racon_ch,
	anot_pepper_ch,
	anot_med_homo_ch,
	anot_pep_homo_ch,
	anot_pep_med_ch)

/*
This process predicts CDS with prokka, it can be used to assess
the quality of the assembly
*/
process prokka_annotation {
  container 'quay.io-biocontainers-prokka-1.14.6--pl5321hdfd78af_4'
  
  tag "${sampleId}"

  publishDir "$params.outdir/annotation/$sampleId", mode: 'copy', pattern: "*"

  cpus params.cores
  
  when params.annotation

  input: set sampleId, name,file(fasta) from anot_prokka_ch
  
  //output: tuple(sampleId, file("${fasta.baseName}*.txt")) into prokka_out_ch
  output: tuple(sampleId, file("${fasta.baseName}*.txt")) into prokka_out_ch
  
  script:
  """
  prokka --cpus ${task.cpus} --outdir . --prefix ${fasta.baseName} $fasta --force
  """
}
}
/*
================================================================
    CheckM
================================================================
*/

/*
This process performs a contamination and completeness check using
CheckM
*/
process checkm {

  container 'metagenlab/checkm:1.0.20'
  containerOptions '--bind /data/databases'

  tag "${sampleId}"

  publishDir "$params.outdir/checkM/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file(fasta) from checkm_fasta_ch

  output:
  file "*checkm_qa.tsv" into checkM_output 

  script:
  """
  checkm lineage_wf --genes -t ${task.cpus} -x fasta . .
  checkm qa -q lineage.ms . -o 2 --tab_table > ${sampleId}_checkm_qa.tsv
  """
}

/*
================================================================
    Report
================================================================
*/
if(params.res && params.annotation &&  params.tax &&  params.coverage){

/*
This process creates a multiqc report composed of outputs from :
- Centrifgue from raw reads
- Centrifuge from Flye's assembly
- Quast
- Prokka
*/
process multiqc {

  conda 'bioconda::multiqc=1.12'

  tag "${sampleId}"

  publishDir "$params.outdir/multiQC/$sampleId", mode: 'copy', pattern: "*"

  input:
  set sampleId, file ('*') ,file ('*') ,file ('*') ,file ('*') ,file ('*') ,file ('*') from json_report_ch.join(centrifuge_reports_ch).join(centrifuge_assembly_reports_ch).join(quast_results).join(cov_json_ch).join(prokka_out_ch.groupTuple())

  output:
  file "multiqc_report.html" into multiqc_report
  file "multiqc_data"

  script:
  """
  multiqc . --cl_config prokka_fn_snames:True
  """
}
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
