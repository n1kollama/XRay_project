import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import os

import functions

"""
Create a plot of all materials with valid atom positions from the database (data_for_different_materials.txt). Uses
functions from functions.py. 
Output as a plot saved in "gain_vs_wavelength_scatter.pdf", as a table of [material name, wavelength, gain] in 
"gain_vs_wavelength_table.txt", and as a list of [material, wavelength, refractive index modulation, coupling 
constant, gain] in "gain_vs_wavelength_list.txt".
"""

# probe dimensions in x, y, z direction
L = 1e-7

# get material data from text file, remove invalid entries (where the atom positions are not specified), and 
# get form factors from the folder that saves the form factor data for all materials.
materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors_local(materials_data, "formfactor_data")

print(f'goodsoup \n Data all done and cleaned up! Starting with gain estimates now... ')

# prepare output arrays
output_dir = 'results'
gain_estimates = np.empty(len(materials_data))
kappas = np.empty(len(materials_data))
ref_indices = np.empty(len(materials_data))
material_wavelength = np.empty(len(materials_data))
i = 0

results = {}

# loop through all materials
for mat, properties in materials_data.items():

    # for output text file
    results[mat] = {
        'ref_index' : None,
        'kappa': None,
        'gain_est': None,
        'wavelength': None
    }

    # obtain results for the material
    n_max = functions.find_max_n(materials_data, mat)
    kappa = functions.coupling_constant_db(n_max, properties['wavelength'])
    g = functions.threshold_gain_db(L, kappa)

    # save for output later on
    results[mat]['ref_index'] = n_max
    results[mat]['kappa'] = kappa 
    results[mat]['gain_est'] = g 
    results[mat]['wavelength'] = properties['wavelength'] * 1e-1 # to save data with wavelength in nanometer

    gain_estimates[i] = g
    material_wavelength[i] = properties['wavelength']
    kappas[i] = kappa 
    ref_indices[i] = n_max 

    i += 1



# write to text files:

# Output file name
output_file_list = os.path.join(output_dir, 'gain_vs_wavelength_list.txt')

# Write to the file
with open(output_file_list, 'w') as f:    
    # Write the data
    for mat, result in results.items():
        lam = result['wavelength']
        gain = result['gain_est']
        kappa = result['kappa']
        n = result['ref_index']

        f.write(f'Data for {mat}:\n')
        f.write(f'Refractive index modulation maximum: {n}\n')
        f.write(f'Coupling Constant: {kappa} \n')
        f.write(f'Gain estimate: {gain} \n')
        f.write(f'Resonance wavelength: {lam} \n \n')


# Output file name
output_file_table = os.path.join(output_dir, 'gain_vs_wavelength_table.txt')

# Write to the file
with open(output_file_table, 'w') as f:
    # Write the header
    f.write('# Material Wavelength Gain\n')
    
    # Write the data
    for mat, result in results.items():
        lam = result['wavelength']
        gain = result['gain_est']
        f.write(f'{mat} {lam} {gain} \n')


# for outputs: convert wavelength to nanometer
material_wavelength *= 1e-1

# create scatter plot for each material with the corresponding material name
plt.figure(figsize=(10,6))
plt.scatter(material_wavelength, gain_estimates)
plt.yscale('log')

for i, material in enumerate(materials_data):
    plt.annotate(material, (material_wavelength[i], gain_estimates[i]), textcoords="offset points", xytext=(5,5), ha='center')

plt.title('Gain estimates')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Gain [1/nm]')

plt.savefig("results/gain_vs_wavelength.pdf")
plt.show()