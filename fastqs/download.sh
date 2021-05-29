#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=download
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cmr736@nyu.edu

module load sra-tools/2.10.9 

fasterq-dump SRR10099491

gzip SRR10099491*.fastq
