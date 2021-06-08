#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=trust4

TRUST4_DIR=/scratch/cmr736/TRUST4/


${TRUST4_DIR}run-trust4 -b Aligned.sortedByCoord.out.bam  -f ${TRUST4_DIR}hg38_bcrtcr.fa --ref ${TRUST4_DIR}human_IMGT+C.fa --barcode CB
