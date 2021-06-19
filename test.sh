#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=nextflowss

module load cutadapt/3.1
# run cutadapt and then replace the original fastq with the trimmed one
cutadapt -u -3 -o output.fastq.gz  /path/to/ubiquitous-rotary-phone/fastqs/SRR10099498_2.fastq.gz
mv output.fastq.gz /path/to/ubiquitous-rotary-phone/fastqs/SRR10099498_2.fastq.gz
