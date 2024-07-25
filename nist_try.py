import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import scipy.constants as const

import functions

materials_data = functions.read_materials_data('data_for_different_materials.txt')

table_1 = functions.atomicformfactor_nist('Ga')

print(table_1)

energy = const.h * const.c / (materials_data['GaP']['wavelength'] * 1e-10 * const.electron_volt) 

atoms = functions.get_atoms(materials_data, 'GaP')

count_1 = atoms.count(atoms[0])
count_2 = atoms.count(atoms[-1])

# find energy value closest to the wavelength we consider
differences_1 = np.abs(table_1[0,:] - energy)
print(differences_1)

#differences_2 = np.abs(table_2[:,0] - energy)

min_index_1 = np.argmin(differences_1)
print(min_index_1)
#min_index_2 = np.argmin(differences_2)

form_factor_1 = np.full(count_1, table_1[1, min_index_1])
#form_factor_2 = np.full(count_2, table_2[min_index_2, 1])

print(form_factor_1)
print(energy)
print(table_1[0, min_index_1])