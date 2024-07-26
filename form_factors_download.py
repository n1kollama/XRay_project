import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
import scipy.constants as const
import os

import functions

materials_data = functions.read_materials_data('data_for_different_materials.txt')

#functions.download_form_factors(materials_data, "formfactor_data")

functions.get_form_factors_local(materials_data, "formfactor_data")

print(materials_data['CuO'])