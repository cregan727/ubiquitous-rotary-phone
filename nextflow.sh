#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=nextflowss
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cmr736@nyu.edu

## run nextflow

module load nextflow/20.11.0-edge
module load star/intel/2.7.6a 

nextflow run Nextflow_script.nf
