import numpy as np
import pandas as pd
import requests as req
import lxml.html as xhtml
from dataclasses import dataclass

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


def get_excitations(Z: int, A: int) -> list[float]:
    """
    Get nuclear excitation levels (in MeV) using the IAEA LiveChart API.
    Returns only excitation energies (no spin/parity), sorted in MeV.
    """
    symbol = global_nuclear_data.get_data(Z, A).isotopicSymbol.lower()

    # the service URL
    livechart = "https://nds.iaea.org/relnsd/v1/data?"

    # There have been cases in which the service returns an HTTP Error 403: Forbidden
    # use this workaround
    import urllib.request
    def lc_pd_dataframe(url):
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
        return pd.read_csv(urllib.request.urlopen(req))

    # load data into a dataframe 
    df = lc_pd_dataframe(livechart + f"fields=levels&nuclides={symbol}")

    levels = (df['energy'].dropna().astype(float)/1000).tolist() #convert keV to MeV

    return sorted(levels)

