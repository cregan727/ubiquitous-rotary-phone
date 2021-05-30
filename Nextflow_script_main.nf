#!/usr/bin/env nextflow

/* 
Nextflow script to perform primary analysis of plate based single cell RNAseq or plate based low input RNAseq.
In its current form is is specifically designed for FB5P-seq data (https://doi.org/10.3389/fimmu.2020.00216) 
But if can be easily modified to other plate based scRNAseq methods with modifications to the input parameters.
Since FB5P-seq is a 5' chemistry based method one specific change is that the params.stranded = 'Reverse' which
is used to set the strandedness in the STARsolo command may need to be changed to Forward or Unstranded for 3' or 
Smartseq applications. 

Notes: Params can be specified with a config file, required modules must either be installed and or loaded before beginning.

Params:

params.reads = path to the fastqs
params.bclist = path to barcodes file for your method
params.reference = path to STAR reference genome
params.pubdir = directory to which the outputs will go
params.author = your name/ sample name/ anything you want to appear alongside the date in the header of the html file - note that commas will be removed
params.stranded = 'Forward' (most 3' methods), 'Reverse' (FB5P-seq), 'Unstranded' (smartseq methods)
params.barcode_length = barcode length 
params.UMI_length = umi length
params.picard_path = path to the picard.jar. In some cases this is just picard.jar

CLI tools:
module load nextflow/20.11.0-edge
module load star/intel/2.7.6a 
module load fastqc/0.11.9
module load samtools/intel/1.12
module load picard/2.23.8
umi_tools

Python libraries:
datetime
pandas
base64
sys
numpy
re
operator
scipy.optimize
matplotlib.pyplot

*/



/* Making channel for the downsampling - if you'd like to change the percentages to which the output bamfile is downsampled
the percentages list can be changed. If the percentage is set to 1, the file will be copied instead of downsampled. */

Channel
	.fromList(['1','.75', '.5', '.25', '.125', '.0625'])
	.set {percents}



// Modification below this line should not be necessary for typical pipeline use



//Channel for the fastqs - 1 for the fastqc and the other for the STAR input

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch } 


// Run FASTQC and send output to specified pubdir

process fastqc {
    tag "FASTQC on $sample_id"
    
    publishDir ${pubdir}, mode: 'copy', overwrite: True

    input:
    set sample_id, file(reads) from read_pairs2_ch
    val pubdir from params.pubdir

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
} 


// Map using STARsolo

process starsolo {
	tag "STARsolo"

	publishDir ${pubdir}, mode: 'copy', overwrite: True

	input:
    set pair_id, path(reads) from read_pairs_ch
    file reference from params.reference
	file bclist from params.bclist
	val pubdir from params.pubdir
	val barcode_length from params.barcode_length
	val UMI_length from params.UMI_length
	val stranded from params.stranded
    
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
	STAR --genomeDir ${reference} \
	--runThreadN 8 \
	--readFilesIn $barcode $gene  \
	--soloCBwhitelist ${bclist} \
	--limitBAMsortRAM 20000000000 \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	--soloStrand ${stranded} \
	--soloType Droplet \
	--soloCBlen ${barcode_length} \
	--soloUMIstart ${barcode_length}+1 \
	--soloUMIlen ${UMI_length} 

	samtools index Aligned.sortedByCoord.out.bam
	
	gzip ./Solo.out/Gene/*/*
	
    """



}


/* In order to assess the quality of the library beyond metrics like seqeuncing saturation alone,
downsampling of the bamfile output from STAR can be done and then metrics calc'd. This step produces
a downsampled counts file, and a file which has a count of the number of reads in the bam file for each
barcode. */

process downsample {
	tag "Downsample on $percent"

	input:
	file "Aligned.sortedByCoord.out.bam" from bamfile_ch
	each(percent) from percents
	val bclength from params.barcode_length
	val picard_path from params.picard_path

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
	java -jar ${picard_path} DownsampleSam \
	I= Aligned.sortedByCoord.out.bam \
	O=ds_${percent}_output.bam \
	P=$percent

	fi

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

/* Merge the results of the downsampling into one file of the counts of Reads, UMIs and Genes for each barcode at 
each downsampling percentage and produce output plots that half calculated the half saturation point for the library
based on number of genes and number of UMIs */

process outsatstats {

	publishDir ${pubdir}, mode: 'copy', overwrite: True

	input: 
	val CB from CBs_ch.collect()
	val count from ds_count_ch.collect()
	val bclist from params.bclist
	val pubdir from params.pubdir

	output:
	file "*.png" into plots_ch
	file "outstats.csv"

	script:
	"""
	python /scratch/cmr736/ubiquitous-rotary-phone/write_outstats.py $CB $count $bclist
	"""
	}


/* Produce a Summary HTML for the run. This includes a header with the author and the date the file was generated
Information about the mapping rates to the genome and transcriptome, and stats about the cells based on the STARsolo output.
It also includes a barcode rank plot, and the two saturation plots generated earlier. */

process html {

	publishDir ${pubdir}, mode: 'copy', overwrite: True

	input: 
	val images from plots_ch.collect()
	val logs from alignment_logs.collect()
	val pubdir from params.pubdir
	val author from params.author

	output:
	file "summary.html" into html_outs_ch

	script:
	"""
	python /scratch/cmr736/ubiquitous-rotary-phone/write_html.py $logs $images $author
	"""
}
