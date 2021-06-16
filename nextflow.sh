#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --time=12:00:00
#SBATCH --mem=70GB
#SBATCH --job-name=nextflowss


## run nextflow

module load nextflow/20.11.0-edge
module load star/intel/2.7.6a 
module load fastqc/0.11.9
module load samtools/intel/1.12
module load picard/2.23.8
module load multiqc/1.9

nextflow run Nextflow_script_main.nf -c NFinput.config # -resume
