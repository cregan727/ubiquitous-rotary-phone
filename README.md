# ubiquitous-rotary-phone
## Analysis of FB5P-seq data using Nextflow:
Project Goals:
* Analyze FB5P-seq data from T cells SRR10099498 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP221379) 
* Analyze FB5P-seq data from B cells SRR10099491 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10099491)
* Use Nextflow
* Align using STARsolo
* Estimate library complexity
* Output to fill in a template HTML
* Alternatively, calculate library complexity information and generate HTML file from data previously analyzed with STARsolo - DataT2 is an example of the same dataset as DataT1, but run on a STARsolo output folder instead of fastqs
* TODO: Add MultiQC stats/plots to HTML if it starts with fastqs
* TODO: Allow switching of which read is the barcode read in the STARsolo input command - in current form it's set up for read 2 as the barcode read and read 1 as the cDNA read and read 2 as the barcode read


Data from https://www.frontiersin.org/articles/10.3389/fimmu.2020.00216/full#h7

Examples of the Barcode Rank plot generated for the T cell datasets:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/Barcoderank_plot.png?raw=trues" >

Examples of the UMI based sequencing saturation plots for T cell datasets:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/UMIsat_plot.png?raw=trues" >

Example of UMI counts across the plate:
<img src="https://github.com/cregan727/ubiquitous-rotary-phone/blob/main/DataT1/Platelayout_umis.png?raw=trues" >


More plots for both datasets can be viewed in the summary.html file for both datasets

