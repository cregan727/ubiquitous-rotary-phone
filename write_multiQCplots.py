# Python script to make plots for the html from the multiqc output

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

# open the file, extract the json, and grab the general stats
f = open("./multiqc_data/multiqc_data.json", 'r')
multiqc = json.loads(f.read())
f.close()
pd.set_option('display.float_format', '{:.2f}'.format)
df2 = pd.DataFrame(multiqc['report_general_stats_data'][0])
df2 = df2.T


# Make the bar graphs figure
fig, ax = plt.subplots(1, len(df2.columns), figsize=(15, 4))
for z in range(0, len(df2.columns)):
    ax[z].barh(df2.index, df2[df2.columns[z]], color='cornflowerblue')
    ax[z].set_title(df2.columns[z])
    for i in range(0, len(df2.index)):
        if 'percent' in df2.columns[z]:
            ax[z].text(.45*max(df2[df2.columns[z]]),
                       i,
                       str(round(df2[df2.columns[z]][i], 2)) + "%",
                       multialignment='left',
                       color='black',
                       size='x-large')
        # if the number is really large, shift the middle point of the number back
        elif max(df2[df2.columns[z]]) < 100:
            ax[z].text(.45*max(df2[df2.columns[z]]),
                       i,
                       f"{int(df2[df2.columns[z]][i]):,d}",
                       multialignment='left',
                       color='black',
                       size='x-large')

        else:
            ax[z].text(.3*max(df2[df2.columns[z]]),
                       i,
                       f"{int(df2[df2.columns[z]][i]):,d}",
                       multialignment='left',
                       color='black',
                       size='x-large')
    if z != 0:
        # Remove y-axis label
        ax[z].set_ylabel('')
        ax[z].set_yticks(())
fig.tight_layout()
fig.suptitle('MultiQC Stats', y=1.05)
plt.savefig("MultiQCstats.png", facecolor='white')
plt.close()
print("done with Stats plot")

# Create the per base sequence quality plot
# starting with pulling the data from the json

plot_data = multiqc['report_plot_data'][
    'fastqc_per_base_sequence_quality_plot'][
    'datasets'][0]
list_of_fastqs_to_plot = []

for i in range(0, len(plot_data)):
    temp = pd.DataFrame(plot_data[i]['data'])
    list_of_fastqs_to_plot.append(temp)

plt.figure(figsize=(7.5, 5))
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20(np.linspace(0, 1, 20)))
for i in range(0, len(list_of_fastqs_to_plot)):
    plt.plot(list_of_fastqs_to_plot[i][0],
             list_of_fastqs_to_plot[i][1],
             label=plot_data[i]['name'],
             linewidth=3)

plt.ylim(3, 61)
plt.legend(fontsize='large')
plt.ylabel("Phred score", fontsize='x-large')
plt.xlabel("Base", fontsize='x-large')
plt.title("Per Base Quality Score", fontsize='xx-large')
plt.savefig("Multiqc_pbq.png", facecolor='white')
plt.close()
print("Done with per base qual plot")


# Add Status Check plot

tempdf = pd.DataFrame(multiqc['report_plot_data'][
    'fastqc-status-check-heatmap']['data'])
tempdf.columns = ['item', 'sample', 'data']
stat_check = tempdf.pivot_table(index='item', columns='sample')
stat_check.columns = multiqc['report_plot_data'][
    'fastqc-status-check-heatmap']['ycats']
stat_check.index = multiqc['report_plot_data'][
    'fastqc-status-check-heatmap']['xcats']

plt.figure(figsize=(7.5, 5))
custom_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='steelblue', lw=4)]

sns.heatmap(stat_check, cmap="RdYlBu", vmax=1.1, vmin=.15, cbar=False)
plt.legend(custom_lines, ['Fail', "Pass"])
plt.title("MultiQC Status", fontsize='xx-large')
plt.savefig("Multiqc_status.png", facecolor='white')
plt.close()
print("Done with Status Plot")

print("Done with script")

