import altair as alt
import pandas as pd
import numpy as np
import sys

TRUSTTSV = sys.argv[1]
CELLNUM = sys.argv[2]


BCRTEMP = """
<div style="width: 100%; height: 700px; position: relative; clear:both;overflow: hidden; white-space:nowrap; padding-top: 10px;clear:both;"> 
<div style="width: 110%, position: absolute; clear:both;"> 
      <p style="text-align: left; color: #094D92; font-size: 30px"> BCR STATS: </p>
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
Bcells = TRUSTaligned[TRUSTaligned['cell_type'] == 'B']

print(CELLNUM)
CELLNUM = int(CELLNUM)
print(type(CELLNUM))

# Calculate the percentages of heavy, light, and paired chains found in the B cells
No_BCR = CELLNUM - len(Bcells.index)
L_only = len(Bcells[(Bcells['chain1'] == "*") & (Bcells['chain2'] != "*")])
H_only = len(Bcells[(Bcells['chain2'] == "*") & (Bcells['chain1'] != "*")])
paired = len(Bcells[(Bcells['chain1'] != "*") & (Bcells['chain2'] != "*")])
BCR_stats = pd.DataFrame([No_BCR, L_only, H_only, paired, CELLNUM],
                         index=['No BCR',
                                "Light Chain Only",
                                "Heavy Chain only",
                                "Paired Chains",
                                "Total"],
                         columns=['Number of Cells'])
BCR_stats['Percent of Cells'] = (BCR_stats['Number of Cells']*100/CELLNUM)
BCR_stats['Percent of Cells'] = round(BCR_stats['Percent of Cells'],
                                      2).astype("str")+"%"
BCRSTATSTABLE = BCR_stats.to_html()

BCRSTATSTABLE = BCRSTATSTABLE.replace(
    """ border="1" """, " ").replace(
    "text-align: right", "text-align: left")

# split the heavy and light chain info out of its csv form
Bcells = Bcells.join(pd.DataFrame(Bcells.chain1.str.split(",").tolist(),
                                  columns=['H_'+x for x in TRUST_columns],
                                  index=Bcells.index))
Bcells = Bcells.join(pd.DataFrame(Bcells.chain2.str.split(",").tolist(),
                                  columns=['L_'+x for x in TRUST_columns],
                                  index=Bcells.index))

Bcells = Bcells.drop(columns=['chain1',
                              'chain2',
                              'secondary_chain1',
                              'secondary_chain2'])

# calculate frequencies for freq histogram
lightchainaa = pd.DataFrame(Bcells.groupby(
    'L_cdr3_aa').size(), columns=['freq'])
lightchainaa['chain'] = 'light'
heavychainaa = pd.DataFrame(Bcells.groupby(
    'H_cdr3_aa').size(), columns=['freq'])
heavychainaa['chain'] = 'heavy'
aa_freq = pd.concat([lightchainaa, heavychainaa])

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
pair_B = Bcells[(Bcells['H_V_gene'] != "*") & (Bcells['L_V_gene'] != "*")]
dotplot_data = pair_B[['H_C_gene', 'L_C_gene']]
dotplot_data = pd.DataFrame(dotplot_data.groupby(
    ["H_C_gene", "L_C_gene"]).size())
dotplot_data = dotplot_data.unstack()
NaN = np.nan
dotplot_data.loc['IGHG3'] = [NaN] * len(dotplot_data.columns)
dotplot_data.loc['IGHG4'] = [NaN] * len(dotplot_data.columns)
dotplot_data.loc['IGHA'] = [NaN] * len(dotplot_data.columns)
dotplot_data.loc['IGHD'] = [NaN] * len(dotplot_data.columns)
dotplot_data.loc['IGHE'] = [NaN] * len(dotplot_data.columns)

# Heavy/Light dotplot
source = dotplot_data.melt(col_level='L_C_gene')
source["H_C_gene"] = source.index
source['Light Chain'] = np.where(["IGK" in x for x in source['L_C_gene']],
                                 "Kappa", "Lambda")

dotplot_BCR = alt.Chart(source).mark_circle().encode(
    x=alt.X('H_C_gene',
            scale=alt.Scale(
                domain=np.sort(source.index.tolist())),
            title='Heavy Chain'),
    y=alt.Y('L_C_gene',
            title='Light Chain'),
    size='value',
    color=alt.Color('Light Chain', scale=alt.Scale(scheme='viridis')),
    tooltip=['value', 'Light Chain']
).properties(
    width='container',
    height=200,
    title="Heavy and Light Chain usage for cells with Paired Data",
).interactive()

bcr_dp_html = dotplot_BCR.to_html()
bcr_dp_html = bcr_dp_html[bcr_dp_html.find("<div id"):bcr_dp_html.find("</body>")]
bcr_dp_html = bcr_dp_html.replace("vis", "vis_dp_b")
bcr_dp_html = bcr_dp_html.replace(
    """<div id="vis_dp_b"></div>""",
    """<div id="vis_dp_b" style="width: 100%; position: absolute; clear:both; overflow: hidden; white-space:nowrap; padding-top: 300px; padding-bottom: 20px"></div>""")

freq_hist_html = freqhist.to_html()
freq_hist_html = freq_hist_html[freq_hist_html.find("<div id"):freq_hist_html.find("</body>")]
freq_hist_html = freq_hist_html.replace("vis", "vis_fq_b")
freq_hist_html = freq_hist_html.replace(
    """<div id="vis_fq_b"></div>""",
    """<div id="vis_fq_b" style="width: 50%; height: 50%; float: left; position: absolute; clear:both; overflow: hidden; white-space:nowrap"></div>""")

BCRTEMP = BCRTEMP.replace("FREQHIST", freq_hist_html)
BCRTEMP = BCRTEMP.replace("TABLE", BCRSTATSTABLE)
BCRTEMP = BCRTEMP.replace("DotPLOT", bcr_dp_html)

# final output
f = open("summary.html", "r")
htmlfile = f.read()
f.close()
htmlfile = htmlfile.replace('ADD_BCR_INFO', BCRTEMP)
f = open("summary.html", 'w')
f.write(htmlfile)
f.close()

print("Done adding to html")
