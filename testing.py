import requests
from bs4 import BeautifulSoup
import pandas as pd
import functools
from io import StringIO
from urllib.parse import urljoin
import numpy as np 


@functools.cache
def atomicformfactor_nist(Z):
    urltmplate = 'https://physics.nist.gov/cgi-bin/ffast/ffast.pl?gtype=4&Formula={Z}'
    url = urltmplate.format(Z=Z)
    session = requests.Session()

    r = session.get(url).content

    soup = BeautifulSoup(r)
    tabledata = soup.select('body > pre')[0].text.splitlines()[3:]
    tabledata = '\n'.join(tabledata)

    data = np.genfromtxt(StringIO(tabledata))
    eV = data[:,0] * 1e3
    f = data[:,1] + 1j*data[:,2]

    return eV, f


eV, f = atomicformfactor_nist("Ga")

print(eV)
print(f)