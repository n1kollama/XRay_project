import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests

import functions

def exp_gain(g_0, lam, lam_peak, sigma):
    return g_0 * np.exp(-(lam - lam_peak)**2/(2 * sigma**2))


materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors_local(materials_data, "formfactor_data")

print(f'\n goodsoup \n Data all done and cleaned up! Starting with the rest of the shit... \n \n ')

material_name = 'CuO'
a_0 = materials_data[material_name]['dimensions'][0]
res_lam = materials_data[material_name]['wavelength']
n_GaP = 3
n_CuO = 1
L = 1e-7
N = 1000

l = np.linspace(4.85, 4.9, N)
del_a = 0.1 # angstrom
gain = np.empty(N)
del_lam = np.empty(N)


for i, li in enumerate(l):
    n = functions.find_max_n_lam(materials_data, material_name, li)
    kappa = functions.coupling_constant_db(n, li)
    gain[i] = functions.threshold_gain_db(L, kappa)

    del_lam[i] = res_lam * np.arcsin((n)/(n_CuO))

sigma = del_lam/(2 * np.sqrt(2 * np.log(2)))
gain_f = exp_gain(gain, l, res_lam, del_lam)

print(sigma)
print(gain_f)
print(l - res_lam)


plt.plot(l, gain_f)
plt.xlabel(f'Wavelength [Angstrom]')
plt.ylabel(f'Gain [nm**(-1)]')
plt.show()


a_ar = np.linspace(a_0 - 4*del_a, a_0 + 4*del_a, N)
gain_ar = np.empty(N)
orig = materials_data[material_name]['dimensions'][0] #original lattice parameter

for i, ai in enumerate(a_ar):
    materials_data[material_name]['dimensions'][0] = ai 
    print(materials_data[material_name]['dimensions'][0])

    n = functions.find_max_n_lam(materials_data, material_name, res_lam)
    k = functions.coupling_constant_db(n, res_lam)
    gain_ar[i] = functions.threshold_gain_db(L, k)

plt.figure(figsize=(10, 8))
plt.plot(a_ar, gain_ar)
plt.xlabel(f'Lattice parameter detuning')
plt.ylabel(f'Gain estimate')
plt.show()




materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors_local(materials_data, "formfactor_data")

print(materials_data['CuO440'])
print(materials_data['CuO'])

n_4 = functions.find_max_n(materials_data, 'CuO440')
k_4 = functions.coupling_constant_db(n, materials_data['CuO440']['wavelength'])
g_4 = functions.threshold_gain_db(L, k_4)

print(n_4)
print(k_4)
print(g_4)