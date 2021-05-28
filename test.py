import sys
import pandas as pd
the_files = str(sys.argv).replace("]", "']").replace("', '[", "'").replace("[", "").replace(",'", "'").replace("'","").replace(" ", "").split(']')
CBs = the_files[0].split(",")
print(CBs[1])
thing = pd.read_csv(CBs[1])
print(thing.iloc[0,0:10])
