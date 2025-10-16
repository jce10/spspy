import numpy as np
import requests as req
import lxml.html as xhtml
from dataclasses import dataclass
import time

PATH_TO_MASSFILE = "./etc/amdc2016_mass.txt"

@dataclass
class NucleusData:
    mass: float = 0.0
    elementSymbol: str = ""
    isotopicSymbol: str = ""
    prettyIsotopicSymbol: str = ""
    Z: int = 0
    A: int = 0

    def __str__(self):
        return self.isotopicSymbol

    def get_latex_rep(self):
        return "$^{" + str(self.A) + "}$" + self.elementSymbol

def generate_nucleus_id(z: np.uint32, a: np.uint32) -> np.uint32 :
    return z*z + z + a if z > a else a*a + z
# def generate_nucleus_id(z: int, a: int) -> int:
#     return z*z + z + a if z > a else a*a + z

class NuclearDataMap:
    U2MEV: float = 931.493614838475
    ELECTRON_MASS: float = 0.000548579909

    def __init__(self):
        self.map = {}

        with open(PATH_TO_MASSFILE) as massfile:
            massfile.readline()
            massfile.readline()
            for line in massfile:
                entries = line.split()
                data = NucleusData()
                data.Z = int(entries[1])
                data.A = int(entries[2])
                data.mass = (float(entries[4])  + 1.0e-6 * float(entries[5]) - float(data.Z) * self.ELECTRON_MASS) * self.U2MEV
                data.elementSymbol = entries[3]
                data.isotopicSymbol = f"{data.A}{entries[3]}"
                data.prettyIsotopicSymbol = f"<sup>{data.A}</sup>{entries[3]}"
                self.map[generate_nucleus_id(data.Z, data.A)] = data

    def get_data(self, z: np.uint32, a: np.uint32) -> NucleusData:
        return self.map[generate_nucleus_id(z, a)]

global_nuclear_data = NuclearDataMap()

# def get_excitations(Z: int, A: int) -> list[float]:
#     levels = []
#     text = ''
#     symbol = global_nuclear_data.get_data(Z, A).isotopicSymbol
#     site = req.get(f"https://www.nndc.bnl.gov/nudat3/getdatasetClassic.jsp?nucleus={symbol}&unc=nds")
#     contents = xhtml.fromstring(site.content)
#     tables = contents.xpath("//table")
#     rows = tables[2].xpath("./tr")
#     for row in rows[1:-2]:
#         entries = row.xpath("./td")
#         if len(entries) != 0:
#             entry = entries[0]
#             data = entry.xpath("./a")
#             if len(data) == 0:
#                 text = entry.text
#             else:
#                 text = data[0].text
#             text = text.replace('?', '')
#             text = text.replace('\xa0\xa0≈','')
#             levels.append(float(text)/1000.0) #convert to MeV
#     return levels

# def get_excitations(Z: int, A: int) -> list[float]:
#     levels = []
#     symbol = global_nuclear_data.get_data(Z, A).isotopicSymbol
#     url = f"https://www.nndc.bnl.gov/nudat3/getdatasetClassic.jsp?nucleus={symbol}&unc=nds"
#     site = req.get(url)

#     if site.status_code != 200:
#         raise RuntimeError(f"Failed to retrieve data for {symbol}. Status code: {site.status_code}")

#     if not site.content.strip():
#         raise RuntimeError(f"Empty document received for {symbol}")

#     try:
#         contents = xhtml.fromstring(site.content)
#         tables = contents.xpath("//table")
#         rows = tables[2].xpath("./tr")
#         for row in rows[1:-2]:
#             entries = row.xpath("./td")
#             if len(entries) != 0:
#                 entry = entries[0]
#                 data = entry.xpath("./a")
#                 text = data[0].text if data else entry.text
#                 if text:
#                     text = text.replace('?', '').replace('\xa0\xa0≈', '')
#                     levels.append(float(text) / 1000.0)  # convert to MeV
#     except Exception as e:
#         raise RuntimeError(f"Error parsing data for {symbol}: {e}")

#     return levels

def get_excitations(Z: int, A: int) -> list[float]:
    levels = []
    symbol = global_nuclear_data.get_data(Z, A).isotopicSymbol
    url = f"https://www.nndc.bnl.gov/nudat3/getdatasetClassic.jsp?nucleus={symbol}&unc=nds"

    for attempt in range(3):  # Try up to 3 times
        site = req.get(url)
        if site.status_code == 200:
            break
        elif site.status_code == 429:
            print(f"Rate limited (429). Waiting before retrying... (attempt {attempt+1})")
            time.sleep(5)  # back off before retrying
        else:
            raise RuntimeError(f"Request failed with status code {site.status_code}")
    else:
        raise RuntimeError("Too many requests – giving up after 3 tries.")

    contents = xhtml.fromstring(site.content)
    tables = contents.xpath("//table")
    rows = tables[2].xpath("./tr")
    for row in rows[1:-2]:
        entries = row.xpath("./td")
        if len(entries) != 0:
            entry = entries[0]
            data = entry.xpath("./a")
            if len(data) == 0:
                text = entry.text
            else:
                text = data[0].text
            text = text.replace('?', '').replace('\xa0\xa0≈','')
            levels.append(float(text)/1000.0)
    return levels


                

        