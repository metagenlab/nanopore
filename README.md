# Nanopore nextflow pipeline

Clinical pipeline for real time analysis of ONT sequencing data.

This pipeline includes species identification (using Centrifuge) and resistance gene identification (using RGI-CARD).

This pipeline assumes the use of a multiplexing sequencing kit as it will gather fasta files based on the barcodes provided :

The following inputs are expected :
- fastq     : the folder in which sequencing data were gathered (absolute path, it should end by "/barcodeXX", e.g. XX = 01)
- sampleID  : a unique name for your sample
- summary   : a basecalling summary file from guppy (absolute path)  

an optionnal input can be provided to try to use reads for which demultiplexing failed with `--keep_unclassified`.
- unclassified : a *single* ONT fastq file containing reads for which demultipexing failed 

These information have to be gather in a csv formated file, with the following header (in any order):
`sampleId,fastq,summary,reference(,unclassified)`

A complete run looks like the following :
`nextflow run nanopore_main.nf --samples <your_input_file.csv> --outdir <output_folder_name> --summary --tax`

Currently, the whole pipeline needs to be run to create a multiQC report (found in the folder multiqc).
But you can also search for tool specific outputs in folders : centrifuge, resistance, coverage or checkM.



## Polishing tools evaluation

To reproduce the polishing evaluation results, you can use the nanopore_evaluation.nf pipeline.

The inputs expected are a little different from the above pipeline :
- fastq     : a *single* ONT fastq file. (absolute path, ´cat´ all passed reads is a way to create such a file)
- sampleID  : a unique name for your sample
- summary   : a basecalling summary file from guppy (absolute path)
- reference : a reference genome (absolute path)

These information have to be gather in a csv formated file, with the following header (in any order):
`sampleId,fastq,summary,reference`

A complete run looks like the following :
`nextflow run nanopore_evaluation.nf --samples <your_input_file.csv> --outdir <output_folder_name> --summary --map --coverage --tax --reference --res --annotation`

Currently, the whole pipeline needs to be run to create a multiQC report (found in the folder multiqc).
But you can also search for tool specific outputs in folders : centrifuge, resistance, coverage or checkM.




