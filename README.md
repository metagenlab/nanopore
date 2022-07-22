# Nanopore nextflow pipeline
Clinical workflow for real time analysis of ONT sequencing data. (still work in progress)

This includes species identification (using centrifuge) and resistance gene identification (using RGI-CARD).

The pipeline expects the folowing inputs :
- fastq     : a *single* ONT fastq file. (absolute path, cat all passed reads is a way to create such a file)
- sampleID  : a unique name for your sample
- summary   : a basecalling summary file from guppy (optional, absolute path)
- reference : a reference genome (optional)

These information have to be gather in a csv formated file, with the following header (in any order):
`sampleId,fastq,summary,reference`

A complete run looks like the following :
`nextflow run nanopore.nf --samples <your_input_file.csv> --outdir <unique_output_folder_name> --summary --map --coverage --tax --reference --res --annotation`

Currently, the whole pipeline needs to be run to create a multiQC report (found in the folder multiqc).
But you can also search for information in folders : centrifuge, resistance, coverage or checkM.


Because some databases are only found in my working directory I would suggest to run it from there.
(I think the pipeline currently only works with a reference genome)
