# ubiquitous-rotary-phone

As of the switch from DSL1 -> DSL2 this pipeline will no longer work with newer versions of Nextflow

## Analysis of FB5P-seq data using Nextflow:
### Project Goals:
* Analyze FB5P-seq data from T cells SRR10099498 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP221379) 
* Analyze FB5P-seq data from B cells SRR10099491 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10099491)
* Use Nextflow
* Align using STARsolo
* Estimate library complexity
* Output to fill in a template HTML
* Optionally perform immune repertoire reconstruction on human or mouse data
* Alternatively, calculate library complexity information and generate HTML file from data previously analyzed with STARsolo - DataT1_fromSTARouts is an example of the same dataset as DataT1, but run on a STARsolo output folder instead of fastqs

### Usage instructions:
To Run on FBP5 data:
  1. Download the SRA data or begin with fastq paths to your own data
  2. If using the SRA data or data in which the barcode read is longer than the BC+UMI, a full example of a slurm script for this can be found in trimFB5P_R2fastqs.sh
  ```
  # To trim the barcode fastq (R2) from the T-cell Data for the nextflow pipeline:
  cutadapt -u -3 -o output.fastq.gz  /path/to/ubiquitous-rotary-phone/fastqs/SRR10099498_2.fastq.gz
  mv output.fastq.gz /path/to/ubiquitous-rotary-phone/fastqs/SRR10099498_2.fastq.gz
  ```
  3. Run the nextflow main script. An example script for this is nextflow.sh. This script is written for a slurm workload manager running in local mode.
  ```
  #Run this nextflow command to run the full script after setting your inputs in the NFinput.config, with or without the -resume option, 
  nextflow run Nextflow_script_main.nf -c NFinput.config # -resume
  ```
  4. In order to run TRUST4 and produce VDJ output for your human or mouse samples, install TRUST4 from github: https://github.com/liulab-dfci/TRUST4 and use the path to the directory for params.trust4_dir.

### Sample Output:

Data from https://www.frontiersin.org/articles/10.3389/fimmu.2020.00216/full#h7

Examples of the Barcode Rank plot generated for the T cell datasets:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/Barcoderank_plot.png?raw=trues" >

Examples of the UMI based sequencing saturation plots for T cell datasets:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/UMIsat_plot.png?raw=trues" >

Example of UMI counts across the plate:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/Platelayout_umis.png?raw=trues" >


More plots for both datasets can be viewed in the summary.html file for both datasets

### Required:
* nextflow/20.11.0-edge
* star/intel/2.7.6a 
* fastqc/0.11.9
* samtools/intel/1.12
* picard/2.23.8
* multiqc/1.9
* umi-tools
* sra-tools/2.10.9 (optional if downloading SRA data)
* cutadapt/3.1 (optional if barcode read is longer than BC+UMI)
* python 3.9: 
  * datetime
  * pandas
  * base64
  * sys
  * numpy
  * re
  * operator
  * scipy.optimize
  * matplotlib.pyplot
  * altair
