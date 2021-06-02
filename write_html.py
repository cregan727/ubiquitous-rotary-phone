from datetime import date
import pandas as pd
import base64
import sys

# INPUTS:
# import nextflow input
system_input = str(
    sys.argv[1:]).replace(
    "]", "']").replace(
    "', '[", "'").replace(
    ",'", "'").replace(
    " /", "/").replace(
    "'", "").split(']')

logs_list = system_input[0].replace(" ", "").split(",")
logs_list = [x.replace("[", "") for x in logs_list]

Log_final_out_path = [x for x in logs_list if x.endswith('Log.final.out')][0]
summary_csv_path = [x for x in logs_list if x.endswith('Solo.out')][0]
summary_csv_path = str(summary_csv_path) + "/Gene/Summary.csv"

# Add images input in order
images_input = system_input[1].split(",")
images_input = [x.replace("[", "") for x in images_input]
# force order of images for downstream process:
images = []
for i in ['Barcoderank_plot.png', 'Genesat_plot.png', 'UMIsat_plot.png', 'Platelayout_cells.png', 'Platelayout_umis.png']:
    images.append([x.replace(" ", "") for x in images_input if x.endswith(i)][0])

# load inputs as dfs
Log_final_out = pd.read_csv(Log_final_out_path,
                            delimiter="\t", header=None, index_col=0)
Log_final_out.index = Log_final_out.index.str.replace('|', '').str. strip()
Summary_csv = pd.read_csv(summary_csv_path, header=None, index_col=0)
Summary_csv = Summary_csv.astype(str)

# Write HTML
Author = system_input[2].replace(",", "")
Sample = system_input[3].replace(",", "")

# Get date
x = date.today()
DATE_M_Y = str(
    x.strftime("%B") + " " +
    x.strftime("%d") + ", " +
    x.strftime("%Y"))

# helper functions


def fixtableforhtml(tablecolumn):
    """
    fix the formatting to make the decimals into percents and
    add comma seps to long numbers
    """
    new_table = []
    for x in tablecolumn:
        if float(x) < 1:
            x = str(round(float(x) * 100, 2))+"%"
        else:
            x = int(float(x))
            x = f"{x:,d}"
        new_table.append(x)
    return(new_table)


def outformathtml(pandasdf):
    """
    change a few formating things to prettify and make it match
    """
    pandas_table = pandasdf.to_html()
    pandas_table = pandas_table.replace(""" border="1" """, " ")
    pandas_table = pandas_table.replace("""<tr style="text-align: right;">\n      <th></th>\n      <th></th>\n    </tr>\n""",
                                        "")
    pandas_table = pandas_table.replace("""<thead>""",
                                        """<thead style="text-align: left; color: #094D92; font-size: 30px"> """)
    return(pandas_table)


# ##WRITE SEQUENCING TABLE from other tables one row at a time, out to HTML
# Then fix the HTML to fill in the html file

# Sequencing:
sequencing_pd = pd.DataFrame()
# Number of Reads:
sequencing_pd = sequencing_pd.append(Summary_csv.loc[["Number of Reads"]])
# Number of Short Reads Skipped
sequencing_pd = sequencing_pd.append(
    Log_final_out.loc[['Number of reads unmapped: too short']])
# Reads with Valid Barcodes
sequencing_pd = sequencing_pd.append(
    Summary_csv.loc[["Reads With Valid Barcodes"]])
# Sequencing Saturation
sequencing_pd = sequencing_pd.append(
    Summary_csv.loc[["Sequencing Saturation"]])
# Q30 Bases in CB+UMI
sequencing_pd = sequencing_pd.append(Summary_csv.loc[["Q30 Bases in CB+UMI"]])
# Q30 Bases in RNA Read
sequencing_pd = sequencing_pd.append(
    Summary_csv.loc[["Q30 Bases in RNA read"]])
# change column name, and fix formatting
seqpd_i = sequencing_pd.index.tolist()
seqpd_i[1] = "Number of Short Reads Skipped"
sequencing_pd.index = seqpd_i
sequencing_pd[1] = fixtableforhtml(sequencing_pd[1])
sequencing_pd.index.name = "Sequencing"
sequencing_pd.columns = [""]
Sequencing_table = outformathtml(sequencing_pd)
print("made sequencing table")


# ##WRITE MAPPING TABLE from other tables one row at a time, out to HTML
# Then fix the HTML to fill in the html file
mapping_pd = pd.DataFrame()
# Reads Mapped to Genome
mapping_pd = mapping_pd.append(
    Summary_csv.loc[["Reads Mapped to Genome: Unique+Multiple"]])
# Reads Mapped Confidently to Genome
mapping_pd = mapping_pd.append(
    Summary_csv.loc[["Reads Mapped to Genome: Unique"]])
# Reads Mapped Confidently to Transcriptome
mapping_pd = mapping_pd.append(
    Summary_csv.loc[["Reads Mapped to Transcriptome: Unique Genes"]])
mapping_pd.index = ['Reads Mapped to Genome',
                    'Reads Mapped Confidently to Genome',
                    'Reads Mapped Confidently to Transcriptome']
mapping_pd[1] = fixtableforhtml(mapping_pd[1])
mapping_pd.index.name = "Mapping"
mapping_pd.columns = [""]

Mapping_table = outformathtml(mapping_pd)
print("made mapping table")


# ##WRITE CELLS TABLE from other tables one row at a time, out to HTML
# Then fix the HTML to fill in the html file
cell_pd = pd.DataFrame()
# Estimated Number of Cells
cell_pd = cell_pd.append(
    Summary_csv.loc[['Estimated Number of Cells']])
# Fraction of Reads in Cells
cell_pd = cell_pd.append(
    Summary_csv.loc[['Fraction of Reads in Cells']])
# Mean Reads per Cell
cell_pd = cell_pd.append(
    Summary_csv.loc[['Mean Reads per Cell']])
# Median Genes per Cell
cell_pd = cell_pd.append(
    Summary_csv.loc[['Median Genes per Cell']])
# Total Genes Detected
cell_pd = cell_pd.append(
    Summary_csv.loc[['Total Genes Detected']])
# Median UMI per Cell
cell_pd = cell_pd.append(
    Summary_csv.loc[['Median UMI per Cell']])
cell_pd[1] = fixtableforhtml(cell_pd[1])
cell_pd.index.name = "Cells"
cell_pd.columns = [""]

Cell_table = outformathtml(cell_pd)
print("made cells table")

# encode the images from the previous plots
str_files = []
for image in images:
    encoded_string = ""
    with open(image, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        str_file = str(encoded_string)
        str_files.append(str_file[2:len(str_file)-1])
print("Base64 enocoded plots")


Template = """
<!DOCTYPE html>
<html>
<head>
<style>
    .tablemainhead {
        color: #094D92;
        font-size: 35px;
        }
    .tablemain {
        text-align:center;
        }
    .dataframe{
      width: 100%
    }
</style>
</head>
<body>
<div style="background-color: Black; text-align:left; padding: 20px">
    <b style="color: White; font-size: 50px"> Plate Based Sequencing Results  </b>
    <p style="color: White;">AUTHOR_NAME -  DATE_M_Y - SAMPLE_NAME</p>
    <HR WIDTH="100%" COLOR="#17BECF" SIZE="4">
</div>

<div style="align:center; width: 100%">
<table class="tablemain" style="align:center; width: 100%;">
    <thead class="tablemainhead">
        <th>CELLNUMBER</th>
        <th>MEANREADSPC</th>
        <th>MEDIANGENESPC</th>
      </thead>
      <tbody class="tablemain">
          <tr>
              <td>Estimated Number of Cells/Wells</td>
              <td>Mean Reads per Cell/Well</td>
              <td>Median Genes per Cell/Well</td>
          </tr>
      </tbody>
</table>
    </div>
<div style="float: left; width: 38%">
    SEQUENCINGTABLE
    MAPPINGTABLE
    CELLSTABLE
</div>
<div style="width: 60%; float: right; overflow:hidden; padding: .2%;">
    <img src="data:image/png;base64, BARCODERANKPLOT" alt="BarcordeRankPlot" style="width: 80%; float: right;">
      </div>
<div style="width: 110%, position: absolute; clear:both;"> 
<p style="text-align: left; color: #094D92; font-size: 30px"> Seqeuncing Saturation Curves: </p>
</div>
    <img src="data:image/png;base64, UMISSATPLT" alt="UMIs plot" style="width: 45%; float: left;">
    <img src="data:image/png;base64, GENESATPLT" alt="BarcordeRankPlot" style="width: 45%; float: right;">
<div style="width: 110%, position: absolute; clear:both;"> 
      <p style="text-align: left; color: #094D92; font-size: 30px"> Plate Layouts: </p>
      </div>
    <img src="data:image/png;base64, PLATELAYCELL" alt="UMIs plot" style="width: 44%; float: left;">
    <img src="data:image/png;base64, PLATELAYUMI" alt="BarcordeRankPlot" style="width: 55%; float: right;">
    
</body>
</html>


"""

# File in the HTML

Template = Template.replace("AUTHOR_NAME", Author)
Template = Template.replace("DATE_M_Y", DATE_M_Y)
Template = Template.replace("SAMPLE_NAME", Sample)
Template = Template.replace(
    "CELLNUMBER",
    Summary_csv.loc["Estimated Number of Cells"][1].replace(".0", ""))
Template = Template.replace(
    "MEANREADSPC",
    f"""{int(float(Summary_csv.loc["Mean Reads per Cell"][1])):,d}""")
Template = Template.replace(
    "MEDIANGENESPC",
    f"""{int(float(Summary_csv.loc["Median Genes per Cell"][1])):,d}""")
Template = Template.replace("SEQUENCINGTABLE", Sequencing_table)
Template = Template.replace("MAPPINGTABLE", Mapping_table)
Template = Template.replace("CELLSTABLE", Cell_table)
Template = Template.replace("BARCODERANKPLOT", str_files[0])
Template = Template.replace("UMISSATPLT", str_files[2])
Template = Template.replace("GENESATPLT", str_files[1])
Template = Template.replace("PLATELAYCELL", str_files[3])
Template = Template.replace("PLATELAYUMI", str_files[4])

# write to a file
f = open("summary.html", "a")
f.write(Template)
f.close()

print("End of script")
