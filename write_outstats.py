#Python script to combine the outstats from the downsampling

import pandas as pd
import numpy as np
import sys

system_input = str(sys.argv[1:]).replace("]", "']").replace("', '[", "'").replace(",'", "'").replace("'","").replace(" ", "").split(']')

CBFiles = system_input[0].split(",")
CBFiles = [x.replace("[","") for x in CBFiles]
Countfiles = system_input[1].split(",")
Countfiles = [x.replace("[","") for x in Countfiles]
percents = system_input[2].split(",")

inputbcs = '/scratch/cmr736/ubiquitous-rotary-phone/brbseq.wlist.txt'
inputbcs = pd.read_csv(inputbcs, header=None)
inputbcs.columns = ['Barcode']
print(inputbcs[0:10])

def outsatstats_all(percent, Reads_per_CB, counts, inputbcs):
    UMIs = pd.read_csv(counts, delimiter='\t', index_col='gene', compression='gzip')
    #fill in an empty column for any barcodes that have no UMIs at this read depth
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
    print(reads.iloc[:,0:5]
    #count number of genes for each barcode and UMIs per barcode
    #reads['Genes'] = np.count_nonzero(UMIs, axis=0)
    #reads['UMI'] = UMIs.sum(axis=0)
    #return(reads)
print(CBFiles)
print(Countfiles)

df = pd.DataFrame(index=inputbcs['Barcode'])
for i in range(0,len(percents)):
    df_temp = outsatstats_all(percents[i],CBFiles[i], Countfiles[i],  inputbcs)
    df_temp.columns = [percents[i]+" "+z for z in df_temp.columns]
    df = df.join(df_temp)
    print(i)
df.to_csv("fullstats.csv")
print('script done')


