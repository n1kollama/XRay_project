import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import os

import functions
"""
Create plots of the gain tolerance with respect to deviations of the wavelength and the lattice parameter from the 
ideal values.
Uses the gain estimate as proposed in Yarivs paper as peak gain which is then modulated by an assumed gaussian 
decay away from the ideal value. The full width at half maximum is assumed to be the spectral width of the system, 
analogous to a Bragg mirror.
Output is given as plots in "gain_tolerance_wavelength.pdf" and "gain_tolerance_a0.pdf", respectively, as well as in
text files of [gain, wavelength] in "gain_tolerance_wavelength_table.txt", and [gain, lattice parameter] in 
"gain_tolerance_a0_table.txt".
"""

# assumed gaussian curve
def exp_gain(g_0, lam, lam_peak, sigma):
    return g_0 * np.exp(-(lam - lam_peak)**2/(2 * sigma**2))


# read materials data, clean up, and extract corresponding form factors for each material
materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors_local(materials_data, "formfactor_data")

print(f'\n goodsoup \n Data all done and cleaned up! Starting with the rest of the shit... \n \n ')

output_dir = 'results'

# material for which to determine these plots, and extract this materials ideal values of lattice parameter and 
# resonance wavelength
material_name = 'CuO'
a_0 = materials_data[material_name]['dimensions'][0]
res_lam = materials_data[material_name]['wavelength']
n_bulk = 1   # bulk refractive index of the material, is assumed to be approx. 1 for xrays. can be specified if necessary.
L = 1e-7
N = 1000 # number of data points for the plots. reduce for faster evaluation.

# modification parameters
l = np.linspace(4.85, 4.9, N) # array of values for the wavelength to iterate through, format (min, max, number of values)

# arrays to save results i
gain = np.empty(N)
del_lam = np.empty(N)


# iterate over the array of wavelengths
for i, li in enumerate(l):
    n = functions.find_max_n_lam(materials_data, material_name, li)
    kappa = functions.coupling_constant_db(n, li)
    gain[i] = functions.threshold_gain_db(L, kappa)

    # determine spectral width for each wavelength value
    del_lam[i] = res_lam * np.arcsin(n/n_bulk)

sigma = del_lam/(2 * np.sqrt(2 * np.log(2)))
gain_f = exp_gain(gain, l, res_lam, del_lam)


# for outputs: convert wavelength to nanometer
l *= 1e-1

# write to text file:

# Output file name
output_file_wavelength = os.path.join(output_dir, 'gain_tolerance_wavelength_table.txt')

# Write to the file
with open(output_file_wavelength, 'w') as f:    

    for li, gain in zip(l, gain_f):
        # write the header
        f.write(f'# {material_name}: Wavelength Gain \n')

        # write the data
        f.write(f'{li} {gain}\n')



plt.plot(l, gain_f)
plt.xlabel(f'Wavelength [nm]')
plt.ylabel(f'Gain [1/nm]')
plt.savefig('results/gain_tolerance_wavelength.pdf')
plt.show()

# for further calculations transform back to nanometers
l *= 1e1


# array for the lattice parameter to iterate over, centered at ideal a_0
a_ar = np.linspace(a_0 - 0.01, a_0 + 0.01, N)
gain_ar = np.empty(N)
n_ar = np.empty(N)
orig = materials_data[material_name]['dimensions'][0] #original lattice parameter

for i, ai in enumerate(a_ar):
    # set lattice parameter to detuned value, determine respective new resonant wavelength
    materials_data[material_name]['dimensions'][0] = ai 
    materials_data[material_name]['wavelength'] = 2 * ai / np.sqrt(3)

    # determine results
    n = functions.find_max_n(materials_data, material_name)
    k = functions.coupling_constant_db(n, res_lam)
    g = functions.threshold_gain_db(L, k)
    del_lam = 2 * ai / np.sqrt(3) * np.arcsin(n/n_bulk)

    gain_ar[i] = exp_gain(g, res_lam, 2 * ai / np.sqrt(3), del_lam)
    n_ar[i] = n



# write to text file:
# Output file name
output_file_a = os.path.join(output_dir, 'gain_tolerance_lattice_table.txt')

# Write to the file
with open(output_file_a, 'w') as f:    
    # Write the data
    for ai, gain in zip(a_ar, gain_ar):
        # write the header
        f.write(f'# {material_name}: Lattice parameter Gain \n')

        # write the data
        f.write(f'{ai} {gain}\n')


plt.figure(figsize=(10, 8))
plt.plot(a_ar, gain_ar)
plt.xlabel(f'Lattice parameter [Angstrom]')
plt.ylabel(f'Gain estimate [1/nm]')
plt.savefig('results/gain_tolerance_a0.pdf')
plt.show()
