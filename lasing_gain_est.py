import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests

import functions

def exp_gain(g_0, lam, lam_peak, sigma):
    return g_0 * np.exp(-(lam - lam_peak)**2/(2 * sigma**2))


materials_data = functions.read_materials_data('data_for_different_materials.txt')
functions.remove_invalid_entries(materials_data)
functions.get_form_factors(materials_data)

print(f'\n goodsoup \n Data all done and cleaned up! Starting with the rest of the shit... \n \n ')

material_name = 'CuO'
a_0 = materials_data[material_name]['dimensions'][0]
res_lam = materials_data[material_name]['wavelength']
n_GaP = 3
n_CuO = 1
L = 1e-7

l = np.linspace(res_lam - 0.02, res_lam + 0.02, 1000)
del_a = 0.01 # angstrom
gain = np.empty(1000)
del_lam = np.empty(1000)

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

a = np.array([a_0 - del_a, a_0, a_0 + del_a])
gain_ar = np.empty((1000, 3))
del_lam_ar = np.empty((1000, 3))
gain_f_ar = np.empty((1000, 3))
sigma_ar = np.empty((1000, 3))

for j, aj in enumerate(a):
    materials_data[material_name]['dimensions'][0] = aj
    materials_data[material_name]['wavelength'] = 2 * materials_data[material_name]['dimensions'][0] / np.sqrt(3)

    print(materials_data[material_name]['wavelength'])

    for i, li in enumerate(l):
        n = functions.find_max_n_lam(materials_data, material_name, li)
        kappa = functions.coupling_constant_db(n, li)
        gain_ar[i, j] = functions.threshold_gain_db(L, kappa)

        del_lam_ar[i, j] = materials_data[material_name]['wavelength'] * np.arcsin((n)/(n_CuO))

    sigma_ar[:, j] = del_lam_ar[:, j]/(2 * np.sqrt(2 * np.log(2)))
    gain_f_ar[:, j] = exp_gain(gain_ar[:, j], l, materials_data[material_name]['wavelength'], del_lam_ar[:,j])

plt.plot(l, gain_f_ar[:,0], 'r-', label='a = a_0 - del_a')
plt.plot(l, gain_f_ar[:,1], 'b.', label='a = a_0')
plt.plot(l, gain_f_ar[:,2], 'g-.', label='a = a_0 + del_a')
plt.xlabel(f'Wavelength [Angstrom]')
plt.ylabel(f'Gain [nm**(-1)]')
plt.legend()
plt.show()