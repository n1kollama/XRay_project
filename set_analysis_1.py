import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests

import functions

"""
Try out the database approach to the refractive index & lasing gain calculation.

"""

L = 1e-7

# get material data and format it in the desired way
materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors(materials_data)

print(f'\n goodsoup \n Data all done and cleaned up! Starting with gain estimates now... \n \n ')

gain_estimates = np.empty(len(materials_data))
kappas = np.empty(len(materials_data))
ref_indices = np.empty(len(materials_data))
material_wavelength = np.empty(len(materials_data))
i = 0

results = {}

#x = np.linspace(0, 5.45 * np.sqrt(3), 100)
#xi = 1.75
#
#n = functions.refractive_index_db(N_0, materials_data, 'GaP', xi)
#
#n_array = np.array([functions.refractive_index_db(materials_data, 'GaP', xi) for xi in x])
#n_max = functions.find_max_n(materials_data, 'GaP')
#print(n_max)
#
#plt.plot(x, n_array)
#plt.show()
#
#kappa = functions.coupling_constant_db(n_max, materials_data['GaP']['wavelength'])
#gain = functions.threshold_gain_db(L, kappa)
#
#print(f'kappa: {kappa}')
#print(f'gain estimate: {gain}')

for mat, properties in materials_data.items():
    results[mat] = {
        'ref_index' : None,
        'kappa': None,
        'gain_est': None,
        'wavelength': None
    }

    n_max = functions.find_max_n(materials_data, mat)
    kappa = functions.coupling_constant_db(n_max, properties['wavelength'])
    g = functions.threshold_gain_db(L, kappa)

    results[mat]['ref_index'] = n_max
    results[mat]['kappa'] = kappa 
    results[mat]['gain_est'] = g 
    results[mat]['wavelength'] = properties['wavelength']

    gain_estimates[i] = g
    material_wavelength[i] = properties['wavelength']
    kappas[i] = kappa 
    ref_indices[i] = n_max 

    i += 1

for mat, result in results.items():
    print(f'Data for {mat} :')
    print('Refractive index modulation maximum:', result['ref_index'])
    print(f'Coupling Constant:', result['kappa'])
    print(f'Gain estimate:', result['gain_est'])
    print(f'Resonance wavelength:', result['wavelength'], '\n')


#plt.plot(material_wavelength, gain_estimates, '*')
#plt.yscale('log')
#plt.xlabel(f'Wavelength [Angstrom]')
#plt.ylabel(f'Lasing gain [m**-3]')
#plt.show()

# lets make a nice plot. i think we can

plt.figure(figsize=(10,6))
plt.scatter(material_wavelength, gain_estimates)
plt.yscale('log')

for i, material in enumerate(materials_data):
    plt.annotate(material, (material_wavelength[i], gain_estimates[i]), textcoords="offset points", xytext=(5,5), ha='center')

plt.title('Gain estimates')
plt.xlabel('Wavelength [Angstrom]')
plt.ylabel('Gain [nm**(-1)]')

plt.show()