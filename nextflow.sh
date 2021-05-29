#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mem=160GB
#SBATCH --job-name=nextflowss
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cmr736@nyu.edu

## run nextflow

module load nextflow/20.11.0-edge
module load star/intel/2.7.6a 
module load fastqc/0.11.9
module load samtools/intel/1.12
module load picard/2.23.8

nextflow run Nextflow_script2.nf -resume
