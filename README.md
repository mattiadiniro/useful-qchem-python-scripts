# useful-qchem-python-scripts
This repository contains several python scripts and executable files created with the intent of facilitating the workflow of submitting quantum mechanical computations and analyzing their results in a quick and systematic way. Most of these scripts can be particularly useful for carbene-metal-amide (CMA) complexes.

XYZ_to_INP
This script takes a folder containing several XYZ files and converts them to Q-Chem input files (INP) and performs a geometry optimization using the B3LYP/LACVP functional and basis set and then TD-DFT using CAM-B3LYP/LACVP. If a solvent is selected, it will use the PCM solvent model. 

OUT_to_CSV
This script takes a folder containing several OUT files and analyzes the main properties of the chemical structure. This script only works completely on CMA complexes. It will generate a CSV file containing important geometrical properties, as well as MO energies, and excited state energies.

OUT_to_ABS
This script takes a folder containing several OUT files and generates absorbance spectra based on the excited state energies present in the output files. It gives an option to plot and save the figure as a png as well. 

Spectra_to_CIE
This script takes a folder containing several TXT files containing emission spectra and calculates the CIE coordinates of each spectra and plots them on a CIE diagram.
