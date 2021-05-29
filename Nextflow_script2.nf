#!/usr/bin/env nextflow

params.reads = '/scratch/cmr736/ubiquitous-rotary-phone/fastqs/SRR10099498*_{1,2}.fastq.gz'
params.bclist = '/scratch/cmr736/ubiquitous-rotary-phone/barcode_seq_2ndSet.txt'
params.reference = '/scratch/cmr736/references/star_human/'
params.author = 'Claire Regan'
params.stranded = 'Reverse'
params.barcode_length = 8 
params.UMI_length = 8


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
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
} 

process starsolo {
tag "STARsolo"


publishDir './data2/', mode: 'copy', overwrite: false

input:
    set pair_id, path(reads) from read_pairs_ch
    file reference from params.reference
output:
    
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody
    file "Aligned.sortedByCoord.out.bam" into bamfile_ch

script:

    gene = reads[1]
    barcode = reads[0]

    """
STAR --genomeDir /scratch/cmr736/references/star_human/ \
--runThreadN 8 \
--readFilesIn $barcode $gene  \
--soloCBwhitelist /scratch/cmr736/ubiquitous-rotary-phone/barcode_seq_2ndSet.txt \
--limitBAMsortRAM 20000000000 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--soloStrand Reverse \
--soloType Droplet \
--soloCBlen 8 \
--soloUMIstart 9 \
--soloUMIlen 8 

samtools index Aligned.sortedByCoord.out.bam
    """



}

process downsample {
tag "Downsample on $percent"

input:
file "Aligned.sortedByCoord.out.bam" from bamfile_ch
each(percent) from percents
val bclength from params.barcode_length

output:
file "ds_${percent}_output.bam" into ds_bam_ch
file "ds_${percent}_output.bam.bai" into ds_bamind_ch
file "ds_${percent}_counts.tsv.gz" into ds_count_ch
file "ds_${percent}_Reads_per_CB.txt" into CBs_ch

script:

"""


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


samtools view ds_${percent}_output.bam | grep -oE  "CB:Z:.{${bclength}}" > CBUMI.txt
sort -o CBUMI.txt CBUMI.txt
uniq -c CBUMI.txt > ds_${percent}_Reads_per_CB.txt
rm CBUMI.txt

echo "End of Script"

"""

}


process outsatstats {

publishDir './data2/', mode: 'copy', overwrite: true

input: 
val CB from CBs_ch.collect()
val count from ds_count_ch.collect()
val bclist from params.bclist

output:
file "*.png" into plots_ch


script:
"""
python /scratch/cmr736/ubiquitous-rotary-phone/write_outstats.py $CB $count $bclist

"""
}

process html {

publishDir './data2/', mode: 'copy', overwrite: true

input: 
val images from plots_ch.collect()
val logs from alignment_logs.collect()

output:
file "summary.html" into html_outs_ch

script:
"""
python /scratch/cmr736/ubiquitous-rotary-phone/write_html.py $logs $images

"""
}
