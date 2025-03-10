[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)

Kinnex 16S deconcatenation and demultiplexing 
==========

## Revio run folder structure

A typical run folder contains the following subfolders:

* fail_reads
* hifi_reads
* metadata
* pb_formats
* statistics

The pb_formats folder contains the files required for importing the run in SMRTLink.
The hifi_reads folder contains bam files (one per kinnex barcode) required for the code below.

## Method

Each HiFi bam file is first deconcatenated into single Kinnex units using Skera (https://skera.how/). This process uses the list of 12 barcode pairs (6 for cDNA or 16 for single-cell)  corresponding to the assembled subunits.
Deconcatenated reads are demultiplexed using the Lima package (https://lima.how/) and the sample-sheet provided by the lab to produce separate bam files for each barcoded 16S sample (max 384 per barcode).
The resulting BAM files are converted to fastq using bam2fastq.
Lima counts are used to produce a QC mosaic plot, provided to the customer when the project corresponds to the full barcode subset.

## Command

A bash workflow was prepared to streamline the different steps explained above. It relies on a config file that provides all specific arguments and naming necessary during the process.

An example yaml file is shown next:

```
runfolder: "/data/pacbio_data/250204.Revio1.FCC"
movie: "m84247_250204_101305_s4"
hififolder: "hifi_reads"
bcM0001:
  bam: "m84247_250204_101305_s4.hifi_reads.bcM0001.bam"
  samplesheet: "Exp4886_4898_4906_bc01_SMRTLink_Barcodefile_LDS20250128_edits.csv"
bcM0002:
  bam: "m84247_250204_101305_s4.hifi_reads.bcM0002.bam"
  samplesheet: "Exp4921_1_bc02_SMRTLink_Barcodefile_LDS20250128.csv"
bcM0003:
  bam: "m84247_250204_101305_s4.hifi_reads.bcM0003.bam"
  samplesheet: "Exp4921_2_bc03_SMRTLink_Barcodefile_LDS20250128.csv"
bcM0004:
  bam: "m84247_250204_101305_s4.hifi_reads.bcM0004.bam"
  samplesheet:  "Exp4921_3_bc04_SMRTLink_Barcodefile_LDS20250128.csv"
outfolder: "/data/pacbio_data/250204.Revio1.FCC/demux_analysis"
qc_format: "html"
inputs: "inputs"
skera_results: "skera_results"
lima_results: "lima_results"
fastq_results: "fastq_results"
log_skera: "INFO"
log_lima: "INFO"
nthr_skera: 48
nthr_lima: 48
par_bam2fastq: 12
nthr_bam2fastq: 4
mincnt: 12000
final_results: "final_results"
```

The workflow command **Kinnex_16S_decat_demux.sh -c config.yaml** is run from /opt/biotools/Kinnex_16S_decat_demux_bash where the code is stored.

Results are stored in the following folder structure:

* inputs
* skera_results
* lima_results
* fastq_results
* runlog.txt

The fastq_results folder contains the final fastq files, ready for archiving and transfer when the full barcode set corresponds to a single project. In case several projects are mixed within a single barcode sublibrary, additinal work is required to isolate each project data based on sample-sheet information (not automatized in this code version)

## Conclusion

This aggregated analysis replaces a several-step SMRT-Link process which required presence of the operator at several points in time and was not simple to run. Resource optimisation for each subtask allows rapid analysis of a typical run and delivers data ready for delivery.

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
