#Python script to combine the outstats from the downsampling

import pandas as pd
import sys

CBFiles = sys.argv[1]
Countfiles = sys.argv[2]
percents = sys.argv[3]

inputbcs = '/scratch/cmr736/ubiquitous-rotary-phone/brbseq.wlist.txt'


def outsatstats(percent, Reads_per_CB, counts, inputbcs):
    """
    Inputs:
    
    percent -  the percent of reads targeted
    
    Reads_per_CB -  the file containing the number of reads in the bam 
file per cell barcode after filtereing 
    
    couts - the counts.tsv file generated from the UMI tools count 
funtion for the downsampled bam
    
    inputbcs
    
    Outputs:
    
    outstats -  a table of the mean and median reads, UMIs and genes per 
input barcode
    """
    UMIs = pd.read_csv(counts, delimiter='\t', index_col='gene')
    #fill in an empty column for any barcodes that have no UMIs at this 
read depth
    for i in inputbcs['Barcode'].tolist():
        if i not in UMIs.columns:
            UMIs[i] = UMIs.shape[0] * [0]
    #remove barcodes not in the list
    UMIs = UMIs[inputbcs['Barcode']]
    #make reads df
    reads = pd.read_csv(Reads_per_CB, delimiter= '\s', header=None)
    reads[1] = reads[1].str.split(':', expand=True)[2]
    reads.columns = ["Reads", "Barcodes"]
    reads.index = reads['Barcodes']
    reads = reads[['Reads']]
    reads = reads[reads.index.isin(inputbcs['Barcode'])]
    reads.reindex(inputbcs['Barcode'])
    #count number of genes for each barcode and UMIs per barcode
    reads['Genes'] = np.count_nonzero(UMIs, axis=0)
    reads['UMI'] = UMIs.sum(axis=0)
    #calc mean and median read/UMIs/Genes per BC
    outstats= pd.DataFrame(reads.mean(axis=0)).T
    outstats.columns = ['Mean '+ i for i in outstats.columns]
    outstats2 = pd.DataFrame(reads.median(axis=0)).T
    outstats2.columns = ['Median '+ i for i in outstats2.columns]
    outstats[outstats2.columns] = outstats
    outstats['percent'] = percent
    return(outstats)
    
def outsatstats_all(percent, Reads_per_CB, counts, inputbcs):
    UMIs = pd.read_csv(counts, delimiter='\t', index_col='gene')
    #fill in an empty column for any barcodes that have no UMIs at this 
read depth
    for i in inputbcs['Barcode'].tolist():
        if i not in UMIs.columns:
            UMIs[i] = UMIs.shape[0] * [0]
    #remove barcodes not in the list
    UMIs = UMIs[inputbcs['Barcode']]
    #make reads df
    reads = pd.read_csv(Reads_per_CB, delimiter= '\s', header=None)
    reads[1] = reads[1].str.split(':', expand=True)[2]
    reads.columns = ["Reads", "Barcodes"]
    reads.index = reads['Barcodes']
    reads = reads[['Reads']]
    reads = reads[reads.index.isin(inputbcs['Barcode'])]
    reads.reindex(inputbcs['Barcode'])
    #count number of genes for each barcode and UMIs per barcode
    reads['Genes'] = np.count_nonzero(UMIs, axis=0)
    reads['UMI'] = UMIs.sum(axis=0)
    return(reads)

df = pd.DataFrame(index=inputbcs['Barcode'])
for i in range(0,len(percents)):
    df_temp = outsatstats_all(percents[i],CBFiles[i], Countfiles[i]  inputbcs)
    df_temp.columns = [percents[i]+" "+z for z in df_temp.columns]
    df = df.join(df_temp)
    print(i)
df.to_csv("fullstats.csv")
print('script done')


