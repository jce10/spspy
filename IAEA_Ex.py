import pandas as pd

# the service URL
livechart = "https://nds.iaea.org/relnsd/v1/data?"

# There have been cases in which the service returns an HTTP Error 403: Forbidden use this workaround
import urllib.request
def lc_pd_dataframe(url):
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    return pd.read_csv(urllib.request.urlopen(req))

# load data into a dataframe 
df = lc_pd_dataframe(livechart + "fields=levels&nuclides=13C")

# df = df[pd.to_numeric(df['energy'],errors='coerce').notna()] # remove blanks (unknown intensities)
# df.energy = df['energy'].astype(float) # convert to numeric. Note how one can specify the field by a

print(df['energy'])