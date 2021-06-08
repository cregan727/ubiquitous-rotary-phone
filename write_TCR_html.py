import altair as alt
import pandas as pd
import numpy as np
import sys

TRUSTTSV = sys.argv[1]
CELLNUM = sys.argv[2]


TCRTEMP = """
<div style="width: 100%; height: 700px; position: relative; clear:both;overflow: hidden; white-space:nowrap; padding-top: 10px;clear:both;"> 
<div style="width: 110%, position: absolute; clear:both;"> 
      <p style="text-align: left; color: #094D92; font-size: 30px"> TCR STATS: </p>
      </div>
FREQHIST
    <div style="width: 50%; float: right; padding-top: 50px;"> 
      TABLE
      </div>
DotPLOT

</div>
"""

# colnames for chain 1 and 2 - for B cells, the heavy chain is chain 1, and light chain is chain 2
TRUST_columns = ["V_gene",
                 "D_gene",
                 "J_gene",
                 "C_gene",
                 "cdr3_nt",
                 "cdr3_aa",
                 "read_cnt",
                 "consensus_id",
                 "CDR3_germline_similarity",
                 "consensus_full_length"]

TRUSTaligned = pd.read_csv(TRUSTTSV,
    delimiter='\t',
    index_col='#barcode')
Tcells = TRUSTaligned[(TRUSTaligned['cell_type'] == 'abT') | (TRUSTaligned['cell_type'] == 'gdT')]

print(CELLNUM)
CELLNUM = int(CELLNUM)
print(type(CELLNUM))

# Calculate the percentages of alpha, beta, and paired chains found in the T cells
No_TCR = CELLNUM - len(Tcells.index)
A_only = len(Tcells[(Tcells['chain1'] == "*") & (Tcells['chain2'] != "*")])
B_only = len(Tcells[(Tcells['chain2'] == "*") & (Tcells['chain1'] != "*")])
paired = len(Tcells[(Tcells['chain1'] != "*") & (Tcells['chain2'] != "*")])
TCR_stats = pd.DataFrame([No_TCR, A_only, B_only, paired, CELLNUM],
                         index=['No TCR',
                                "Alpha/Gamma Only",
                                "Beta/Delta only",
                                "Paired TCR",
                                "Total"],
                         columns=['Number of Cells'])
TCR_stats['Percent of Cells'] = (TCR_stats['Number of Cells']*100/CELLNUM)
TCR_stats['Percent of Cells'] = round(TCR_stats['Percent of Cells'],
                                      2).astype("str")+"%"
TCRSTATSTABLE = TCR_stats.to_html()

TCRSTATSTABLE = TCRSTATSTABLE.replace(
    """ border="1" """, " ").replace(
    "text-align: right", "text-align: left")

# split the heavy and light chain info out of its csv form
Tcells = Tcells.join(pd.DataFrame(Tcells.chain1.str.split(",").tolist(),
                                  columns=['VDJ_'+x for x in TRUST_columns],
                                  index=Tcells.index))
Tcells = Tcells.join(pd.DataFrame(Tcells.chain2.str.split(",").tolist(),
                                  columns=['VJ_'+x for x in TRUST_columns],
                                  index=Tcells.index))

Tcells = Tcells.drop(columns=['chain1',
                              'chain2',
                              'secondary_chain1',
                              'secondary_chain2'])

# calculate frequencies for freq histogram
vjchainaa = pd.DataFrame(Tcells.groupby(
    'VJ_cdr3_aa').size(), columns=['freq'])
vjchainaa['chain'] = 'VJ'
vdjchainaa = pd.DataFrame(Tcells.groupby(
    'VDJ_cdr3_aa').size(), columns=['freq'])
vdjchainaa['chain'] = 'VDJ'
aa_freq = pd.concat([vjchainaa, vdjchainaa])

freqhist = alt.Chart(aa_freq).mark_bar(
    color='reds').encode(
    alt.X('freq',
          bin=alt.Bin(step=.9999),
          title='Number of Cells the CDR3 Amino Acid Sequence is Observed in'),
    y=alt.Y('count()',
            title="Number of CDR3 AA sequences"),
    color=alt.Color('chain',
                    scale=alt.Scale(scheme='viridis'))).properties(
    width='container',
    height=200,
    title='Frequency Histogram')

# create the dataframe with the dotplot data
# create the dataframe with the dotplot data
pair_T = Tcells[(Tcells['VDJ_V_gene'] != "*") & (Tcells['VJ_V_gene'] != "*")]
dotplot_data = pair_T[['VDJ_C_gene', 'VJ_C_gene']]
dotplot_data = pd.DataFrame(dotplot_data.groupby(
    ["VDJ_C_gene", "VJ_C_gene"]).size())
dotplot_data = dotplot_data.unstack()
NaN = np.nan
dotplot_data['TRGC1'] = [NaN] * len(dotplot_data.index)
dotplot_data['TRGC2'] = [NaN] * len(dotplot_data.index)
dotplot_data.loc['TRDC1'] = [NaN] * len(dotplot_data.columns)
dotplot_data.loc['TRDC2'] = [NaN] * len(dotplot_data.columns)

# Alpha/Beta/Gamma/Delta dotplot
source = dotplot_data.melt(ignore_index=False, col_level='VJ_C_gene')
source["VDJ_C_gene"] = source.index
source['T-Cell Type'] = np.where(["TRA" in x for x in source['VJ_C_gene']],
                                 "AB T-Cell", "GD T-Cell/other")

dotplot_TCR = alt.Chart(source).mark_circle().encode(
    x=alt.X('VDJ_C_gene',
            scale=alt.Scale(
                domain=np.sort(source.index.tolist())),
            title='Beta/Delta chain'),
    y=alt.Y('VJ_C_gene',
            scale=alt.Scale(
                domain=np.sort(['*', "TRAC", "TRGC1", "TRGC2"])),
            title='Alpha/Gamma Chain'),
    size='value',
    color=alt.Color('T-Cell Type', scale=alt.Scale(scheme='viridis')),
    tooltip=['value', 'T-Cell Type']
).properties(
    width=1000,
    height=200,
    title="Chain pairing for cells with Paired Data",
).interactive()

tcr_dp_html = dotplot_TCR.to_html()
tcr_dp_html = tcr_dp_html[tcr_dp_html.find("<div id"):tcr_dp_html.find("</body>")]
tcr_dp_html = tcr_dp_html.replace("vis", "vis_dp_t")
tcr_dp_html = tcr_dp_html.replace(
    """<div id="vis_dp_t"></div>""",
    """<div id="vis_dp_t" style="width: 100%; position: absolute; clear:both; overflow: hidden; white-space:nowrap; padding-top: 300px; padding-bottom: 20px"></div>""")

freq_hist_html = freqhist.to_html()
freq_hist_html = freq_hist_html[freq_hist_html.find("<div id"):freq_hist_html.find("</body>")]
freq_hist_html = freq_hist_html.replace("vis", "vis_fq_t")
freq_hist_html = freq_hist_html.replace(
    """<div id="vis_fq_t"></div>""",
    """<div id="vis_fq_t" style="width: 50%; height: 50%; float: left; position: absolute; clear:both; overflow: hidden; white-space:nowrap"></div>""")

TCRTEMP = TCRTEMP.replace("FREQHIST", freq_hist_html)
TCRTEMP = TCRTEMP.replace("TABLE", TCRSTATSTABLE)
TCRTEMP = TCRTEMP.replace("DotPLOT", tcr_dp_html)

# final output
f = open("summary.html", "r")
htmlfile = f.read()
f.close()
htmlfile = htmlfile.replace('ADD_TCR_INFO', TCRTEMP)
f = open("summary.html", 'w')
f.write(htmlfile)
f.close()
print(np.sort(source.index.tolist()))
print("Done adding to html")
print(htmlfile)

