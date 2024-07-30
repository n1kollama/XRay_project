import numpy as np 
import matplotlib.pyplot as plt 
import requests
from bs4 import BeautifulSoup

import functions 

""" 
To test the functions defined in functions.py and,
most importantly, to compare with Yariv paper to ensure
that the code gives correct (seemingly) results.
Example in the paper that we can compare with:

GaP, (111) direction
N_0 = 8.1 x 10**23 cm-3
lambda = 6.154 A

deltan = 1.72 x 10**-4 cm-1
kappa = 8.8 x 10**-3 cm-1

calculated with this method by hand an on nanometer
~12.x nm-1
"""

# parameters and changing units
N_0 = 8.1e23 * 1e6# me-1
lambd = 6.154e-10 # m
L_1 = np.array([1e-2, 1e-2, 1e-2]) # m




# all parameters and definitions for the example Gallium Phosphide (as considered in the Yariv paper)
miller_GaP = np.array([1,1,1])
N = 8
unitcell_l_GaP = np.array([5.45, 5.45, 5.45]) # angstrom
unitcell_a_GaP = np.array([90, 90, 90]) # degrees
at_pos_GaP = np.array([[0, 0, 0],[1/2, 1/2, 0],[1/2, 0, 1/2],[0, 1/2, 1/2],[1/4, 1/4, 1/4],[3/4, 3/4, 1/4],[3/4, 1/4, 3/4],[1/4, 3/4, 3/4]])
at_el_GaP = ['Ga', 'Ga', 'Ga', 'Ga', 'P', 'P', 'P', 'P']
params = {
    'Ga': ([15.2354, 6.7006, 4.3591, 2.9623],[3.0669, 0.2412, 10.7805, 61.4135],1.7189),
    'P': ([6.4345, 4.1791, 1.78, 1.4908], [1.9067, 27.157, 0.526, 68.1645], 1.1149)
}


# unit cell parameter and ksi examplary value
a_0 = 5.45e-10
ksi = a_0 / np.sqrt(3)


# try out to read the data from my lovely little file
materials_data = functions.read_materials_data('data_for_different_materials.txt')

unitcell_l_GaP = materials_data['GaP']['dimensions']
unitcell_a_GaP = materials_data['GaP']['angles']
miller_GaP = materials_data['GaP']['miller_indices']
at_el_GaP = functions.get_atoms(materials_data, 'GaP')
#at_pos_GaP = functions.get_atom_positions_as_numpy_array(materials_data, 'GaP')



# plot equations for n with respect to ksi (distance along 111 direction), both according to your and Yarivs formula
x = np.linspace(0, a_0*np.sqrt(3), 100)

n_array = np.array([functions.refractive_index_ksi(xi, N_0, lambd, miller_GaP, unitcell_l_GaP, unitcell_a_GaP, at_pos_GaP, at_el_GaP, params) for xi in x])
deln_array = np.array([functions.deltan_paper_eq(N_0, lambd, xi, a_0) for xi in x])


plt.plot(x, n_array)
plt.plot(x, deln_array)
plt.show()

print(f'a_0 * sqrt3 : {a_0 * np.sqrt(3)}')
print(f'a_0 / sqrt3 : {2 * a_0 / np.sqrt(3)}')


# calcualte threshold gain and other parameters, for checking and further consideration
g_GaP = functions.threshold_gain(L_1, N_0, lambd, miller_GaP, unitcell_l_GaP, unitcell_a_GaP, at_pos_GaP, at_el_GaP, params)
n_GaP_yariv = functions.deltan_paper_eq(N_0, lambd, ksi, a_0)

print(f'----------------------------------------------------')
print(f'Threshold Gain: {g_GaP}')

print(f'----------------------------------------------------')
print(f'n as in Yariv paper: {n_GaP_yariv}')

print(f'----------------------------------------------------')

max_n = np.max(n_array)
max_n_yariv = np.max(deln_array)
delta_max = np.abs(max_n - max_n_yariv)

print(f'Max. n as calculated by me: {max_n}')
print(f'Max. n as calculated from Yariv paper: {max_n_yariv}')
print(f'Difference of Maxima: {delta_max}')
