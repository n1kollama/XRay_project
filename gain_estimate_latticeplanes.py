import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import os

import functions
"""
Obtain plots for several materials for all their atom planes to see which ones give the highest gain. 
Output are scatter plots of gain versus wavelength for the respective materials. they are saved as 
"<material>_planes_gain.pdf", and the corresponding tabular text files of [plane, wavelength, gain] are 
saved as "<material>_planes_gain.txt".
"""

# read data from the database, clean up, obtain form factors.
materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors_local(materials_data, "formfactor_data")

output_dir = 'results'

# define volume of the probe
L_x = 1e-7
L_y = 1e-7
L_z = 1e-7

L = np.array([L_x, L_y, L_z])

# obtain data for copper oxide. if material is changed, the array with miller indices and d spacing have to be 
# updated as well.
material_name = 'CuO'
lattice_planes = np.array([[1, 1, 1], [2, 0, 0], [2, 2, 0], [3, 1, 1], [2, 2, 2], [4, 0, 0], [3, 3, 1], [4, 2, 0], [4, 2, 2], [5, 1, 1], [4, 4, 0]])
d_spacing = np.array([2.44, 2.117, 1.497, 1.276, 1.222, 1.058, 0.971, 0.947, 0.864, 0.815, 0.75])

wavelengths = 2 * d_spacing

gains_planes = np.empty(len(d_spacing))

# iterate through all planes to obtain respective gains.
for i, gi in enumerate(gains_planes):
    materials_data[material_name]['miller_indices'] = lattice_planes[i, :]
    materials_data[material_name]['wavelength'] = wavelengths[i]

    n = functions.find_max_n(materials_data, material_name)
    k = functions.coupling_constant_db(n, wavelengths[i])
    gains_planes[i] = functions.threshold_gain_db(L, k)


# for outputs: convert wavelength to nanometer
wavelengths *= 1e-1

# write to text file:
# Output file name
output_file_cuo = os.path.join(output_dir, 'cuo_planes_gain.txt')

# Write to the file
with open(output_file_cuo, 'w') as f:    
    # Write the data
    for plane, lam, gain in zip(lattice_planes, wavelengths, gains_planes):
        # write the header
        f.write(f'# {material_name}: Lattice plane Wavelength Gain \n')

        # write the data
        f.write(f'{plane} {lam} {gain}\n')

# create plot
plt.plot(wavelengths, gains_planes, 'r.')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Gain estimate [1/nm]')
plt.savefig('results/cuo_planes_gain.pdf')
plt.show()


# same procedure as before but for a different material.
material_name = 'MgO'

lattice_planes = np.array([[1, 1, 1], [2, 0, 0], [2, 2, 0], [3, 1, 1], [2, 2, 2], [4, 0, 0], [3, 3, 1], [4, 2, 0], [4, 2, 2], [5, 1, 1]])
d_spacing = np.array([2.421, 2.097, 1.483, 1.265, 1.211, 1.049, 0.962, 0.938, 0.856, 0.807])
wavelengths = 2 * d_spacing

gains_planes = np.empty(len(d_spacing))

for i, gi in enumerate(gains_planes):
    materials_data[material_name]['miller_indices'] = lattice_planes[i, :]
    materials_data[material_name]['wavelength'] = wavelengths[i]

    n = functions.find_max_n(materials_data, material_name)
    k = functions.coupling_constant_db(n, wavelengths[i])
    gains_planes[i] = functions.threshold_gain_db(L, k)



# for outputs: convert wavelength to nanometer
wavelengths *= 1e-1

# write to text file:
# Output file name
output_file_mgo = os.path.join(output_dir, 'mgo_planes_gain.txt')

# Write to the file
with open(output_file_mgo, 'w') as f:    
    # Write the data
    for plane, lam, gain in zip(lattice_planes, wavelengths, gains_planes):
        # write the header
        f.write(f'# {material_name}: Lattice plane Wavelength Gain \n')

        # write the data
        f.write(f'{plane} {lam} {gain}\n')

plt.plot(wavelengths, gains_planes, 'r.')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Gain estimate [1/nm]')
plt.savefig('results/mgo_planes_gain.pdf')
plt.show()


# same procedure as before but for a different material.
material_name = 'CdTe'

lattice_planes = np.array([[1, 1, 1], [2, 2, 0], [3, 1, 1], [4, 0, 0], [3, 3, 1], [4, 2, 2], [5, 1, 1], [4, 4, 0], [5, 3, 1], [6, 2, 0], [5, 3, 3], [4, 4, 4], [7, 1, 1], [6, 4, 2], [5, 5, 3], [8, 0, 0], [7, 3, 3], [8, 2, 2]])
d_spacing = np.array([3.79, 2.321, 1.979, 1.641, 1.506, 1.34, 1.263, 1.16, 1.11, 1.038, 1.001, 0.947, 0.919, 0.877, 0.855, 0.821, 0.802, 0.774])
wavelengths = 2 * d_spacing

gains_planes = np.empty(len(d_spacing))

for i, gi in enumerate(gains_planes):
    materials_data[material_name]['miller_indices'] = lattice_planes[i, :]
    materials_data[material_name]['wavelength'] = wavelengths[i]

    n = functions.find_max_n(materials_data, material_name)
    k = functions.coupling_constant_db(n, wavelengths[i])
    gains_planes[i] = functions.threshold_gain_db(L, k)



# for outputs: convert wavelength to nanometer
wavelengths *= 1e-1

# write to text file:
# Output file name
output_file_cdte = os.path.join(output_dir, 'cdte_planes_gain.txt')

# Write to the file
with open(output_file_cdte, 'w') as f:    
    # Write the data
    for plane, lam, gain in zip(lattice_planes, wavelengths, gains_planes):
        # write the header
        f.write(f'# {material_name}: Lattice plane Wavelength Gain \n')

        # write the data
        f.write(f'{plane} {lam} {gain}\n')

plt.plot(wavelengths, gains_planes, 'r.')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Gain estimate [1/nm]')
plt.savefig('results/cdte_planes_gain.pdf')
plt.show()