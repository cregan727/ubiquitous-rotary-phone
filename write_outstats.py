#Python script to combine the outstats from the downsampling

import pandas as pd
import numpy as np
import sys
import re
import operator
import scipy.optimize
import matplotlib.pyplot as plt

system_input = str(sys.argv[1:]).replace("]", "']").replace("', '[", "'").replace(",'", "'").replace("'","").replace(" ", "").split(']')

CBFiles = system_input[0].split(",")
CBFiles = [x.replace("[","") for x in CBFiles]
Countfiles = system_input[1].split(",")
Countfiles = [x.replace("[","") for x in Countfiles]
percents = [re.search("ds_.*_", x).group(0).replace("ds_","").replace("_","") for x in  Countfiles]
inputbcs = system_input[2].replace(",","").replace(" ","")
print(inputbcs)

inputbcs = pd.read_csv(inputbcs, header=None)
inputbcs.columns = ['Barcode']

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
    for i in inputbcs['Barcode'].tolist():
        if i not in reads.index:
            reads.loc[i] = 0
    reads = reads.reindex(inputbcs['Barcode'], copy=False)
    #count number of genes for each barcode and UMIs per barcode
    reads['Genes'] = np.count_nonzero(UMIs, axis=0)
    reads['UMI'] = UMIs.sum(axis=0)
    return(reads)

df = pd.DataFrame(index=inputbcs['Barcode'])

#force the df columns to be in order:
sorted_pairs = sorted(enumerate(percents), key=operator.itemgetter(1))
indices = [index for index, element in sorted_pairs]

#Combine into one df
for i in indices:
    df_temp = outsatstats_all(percents[i],CBFiles[i], Countfiles[i],  inputbcs)
    df_temp.columns = [percents[i]+" "+z for z in df_temp.columns]
    df = df.join(df_temp)
    print(percents[i])

df_means = df.mean(axis = 0)
df_median = df.median(axis = 0)


#write out
df.to_csv("outstats.csv")
print('done writing stats')


#make barcode rank plot
df_bcrankplot = df[["1 UMI"]]
plt.scatter(range(1,97), df_bcrankplot.sort_values('1 UMI', ascending=False)['1 UMI'])
plt.xscale('log')
plt.yscale('log')
plt.ylabel("UMI counts")
plt.xlabel("Barcodes")
plt.ylim(1,10000000)
plt.title("Barcode Rank Plot")
plt.savefig("Barcoderank_plot.png")
plt.close()

#UMI halfsat

def model(x, vm, km):
    return (vm * x / (x + km))
paramaters, covar = scipy.optimize.curve_fit(model, 
                                             df_means.filter(regex='Reads').tolist(), 
                                             df_median.filter(regex='UMI').tolist(),
					     p0=[1000,1000])

plt.scatter(df_means.filter(regex='Reads'),df_median.filter(regex='UMI'))
x_coords = range(0, int(max(paramaters[1]*1.5,max(df_means.filter(regex='Reads')))) , int(paramaters[1]*.1))
y_coords = [model(x, paramaters[0], paramaters[1]) for x in x_coords]
plt.plot(x_coords, y_coords)
plt.axhline(paramaters[0], color='black', linestyle='-', label='Max Median UMI Count: ' + str(round(paramaters[0])))
plt.axvline(paramaters[1],color='black', linestyle='--', label='Half Saturation Point: ' + str(round(paramaters[1])))
plt.legend()
plt.ylabel("Median UMIs per Cell/Well")
plt.xlabel("Median Reads per Cell/Well")
plt.title("Mean UMI's per Barcode")
plt.savefig("UMIsat_plot.png")
plt.close()

#Gene halfsat

paramaters, covar = scipy.optimize.curve_fit(model, 
                                             df_means.filter(regex='Reads').tolist(), 
                                             df_median.filter(regex='Gene').tolist(),
                                             p0=[1000,1000])
print(paramaters)
plt.scatter(df_means.filter(regex='Reads'),df_median.filter(regex='Gene'))
x_coords = range(0, int(max(paramaters[1]*1.5,max(df_means.filter(regex='Reads')))) , int(paramaters[1]*.1))
y_coords = [model(x, paramaters[0], paramaters[1]) for x in x_coords]
plt.plot(x_coords, y_coords)
plt.axhline(paramaters[0], color='black', linestyle='-', label='Max Median Gene Count: ' + str(round(paramaters[0])))
plt.axvline(paramaters[1],color='black', linestyle='--', label='Half Saturation Point: ' + str(round(paramaters[1])))
plt.legend()
plt.ylabel("Median Genes per Cell/Well")
plt.xlabel("Mean Reads per Cell/Well")
plt.title("Median Genes per Barcode")
plt.savefig("Genesat_plot.png")
plt.close()

print("Done with script")
