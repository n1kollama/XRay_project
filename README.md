Contains:

formfactor_data : directory with several text files; each text file has the name of one atomic element and contains the data downloaded from Henke for the atomic form factors.

results : directory with all results summed up in text files and pdf versions of plots.

data_for_different_materials.txt : text file with a list off all materials that are to be studied. new materials can be added as long as they have the exact same format as the already existing ones. the python files all read from this text file to perform their calculations.

form_factors_download.py : python script that accesses the henke database, downloads all atomic form factors for the atoms we are considering, and saves them in the formfactor_data folder. does not need to be executed anymore, since the folder with the data is already there.

functions.py : python script that contains all functions that are used in the other scripts. is accessed from all other python scripts.

gain_estimate_latticeplanes.py : python script that perform the gain estimate calculation for all reflective lattice planes of three different materials: CuO, MgO, CdTe. 

gain_estimate_materials.py : python script that iterates through the whole list of materials and estimates the gain for each material at the given resonance wavelength.

gain_estimate_tolerance.py : python script that gives an estimate for the tolerance of the received gain from a material for both a detuning in wavelength as well as in the lattice parameter.

example_analysis_GaP.py : script that compares this method of calculating the gain to the Yariv paper in order to verify as far as possible.