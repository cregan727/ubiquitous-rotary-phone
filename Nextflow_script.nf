#!/usr/bin/env nextflow

params.reads = '/scratch/cmr736/ubiquitous-rotary-phone/fastqs/*_R{1,2}.fastq.gz'
params.path1 = '/scratch/cmr736/ubiquitous-rotary-phone/'
params.reference = '/scratch/cmr736/references/refdata-gex-GRCh38-2020-A/star'

Channel
	.fromList(['1','.75', '.5', '.25', '.125', '.0625'])
	.into {percents; percents2}


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch } 


process fastqc {
    tag "FASTQC on $sample_id"

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    module load fastqc/0.11.9
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
} 

process starsolo {
tag "STARsolo"


publishDir './data/'

input:
    set pair_id, path(reads) from read_pairs_ch
output:
    
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody
    file "Aligned.sortedByCoord.out.bam" into bamfile_ch

script:

    gene = reads[0]
    barcode = reads[1]

    """
STAR --genomeDir /scratch/cmr736/references/star_human/ \
--runThreadN 6 \
--readFilesIn $gene $barcode  \
--soloCBwhitelist /scratch/cmr736/ubiquitous-rotary-phone/brbseq.wlist.txt \
--limitBAMsortRAM 20000000000 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--soloType Droplet \
--soloCBlen 6 \
--soloUMIstart 7 \
--soloUMIlen 14 

module load samtools/intel/1.12

samtools index Aligned.sortedByCoord.out.bam
    """



}

process downsample {
tag "Downsample on $percent"

input:
file "Aligned.sortedByCoord.out.bam" from bamfile_ch
each(percent) from percents

output:
file "ds_${percent}_output.bam" into ds_bam_ch
file "ds_${percent}_output.bam.bai" into ds_bamind_ch
file "ds_${percent}_counts.tsv.gz" into ds_count_ch
file "ds_${percent}_Reads_per_CB.txt" into CBs_ch

script:

"""

module load picard/2.23.8

if  [ $percent == 1 ] ; then

cp Aligned.sortedByCoord.out.bam ds_${percent}_output.bam

else 
java -jar /share/apps/picard/2.23.8/picard.jar DownsampleSam \
I= Aligned.sortedByCoord.out.bam \
O=ds_${percent}_output.bam \
P=$percent

fi

module load samtools/intel/1.12
samtools index ds_${percent}_output.bam

##  Use UMI tools to count
umi_tools count --umi-tag=UB \
--cell-tag=CB \
--gene-tag=GX \
--extract-umi-method=tag \
--per-gene \
--wide-format-cell-counts \
--per-cell \
-I ds_${percent}_output.bam \
-S ds_${percent}_counts.tsv.gz


samtools view ds_${percent}_output.bam | grep -o  CB:Z:...... > CBUMI.txt
sort -o CBUMI.txt CBUMI.txt
uniq -c CBUMI.txt > ds_${percent}_Reads_per_CB.txt
rm CBUMI.txt

echo "End of Script"

"""

}


process outsatstats {

publishDir './data/'

input: 
val CB from CBs_ch.collect()
val count from ds_count_ch.collect()
val percent from percents2.collect()

output:
file "Barcoderank_plot.png" into brp_ch
file "Genesat_plot.png" into genesat_ch
file "UMIsat_plot.png" into umisat_ch


script:
"""
python /scratch/cmr736/ubiquitous-rotary-phone/write_outstats.py $CB $count $percent

"""
}
