import numpy as np
import os
import matplotlib.pyplot as plt
import math
import PySimpleGUI as sg

def GaussSpectrum(x_eV, band, strength, sigma):
    "Return a normalized Gaussian, input of energy is eV"
    GaussBand = strength * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((x_eV - band) ** 2) / (2 * sigma ** 2))
    return GaussBand

# change the layout
layout = [
    [sg.Text("Convert your Q-Chem output files (.out) into absorbance spectra.")],
    [sg.Text("Enter folder path")],
    [sg.Input(key="-FOLDER-"), sg.Button("Convert")],
    [sg.Text("", key="-OUTPUT-")],
    [sg.Text("Wavelength range (nm)", size=20), sg.Input(key="-START-", default_text="200", size=5),
     sg.Input(key="-FINISH-", default_text="900", size=5)],
    [sg.Text("FWHM (eV)", size=20), sg.Input(key="-FWHM-", size=5, default_text="0.2")],
    [sg.Checkbox("Write CSV file", key="-WRITECSV-", default=True)],
    [sg.Checkbox("Save PNG image", key="-SAVEPNG-", default=False)],
    [sg.Checkbox("Show plots", key="-PLOT-", default=False)],

]

window = sg.Window("Absorbance Spectrum Calculator", layout)

# Initialize arrays to store all spectra and their legends
all_x_nm = []
all_spectra = []
all_legends = []

while True:
    event, values = window.read()

    if event == sg.WIN_CLOSED:
        break

    if event == "Convert":
        folder_path = values["-FOLDER-"]
        # Validate and convert input values
        try:
            start = int(values["-START-"])
            finish = int(values["-FINISH-"])
            width = float(values["-FWHM-"])
        except ValueError:
            sg.popup_ok("Invalid input values. Please enter valid numbers.", title="Error")
            continue  # Skip the conversion if the inputs are invalid
        if start < 100 or start > 1200 or finish < 100 or finish > 1200 or start > finish:
            sg.popup_ok("Invalid Wavelength range. Please enter valid numbers (100-1200 nm).", title="Error")
            continue  # Skip the conversion if the inputs are invalid  
        if width < 0.01 or width > 1:
            sg.popup_ok("Invalid FWHM. Please enter a valid number (0.01-1.00 eV).", title="Error")
            continue  # Skip the conversion if the inputs are invalid
        plotSpectrum = values["-PLOT-"]
        writeCSV = values["-WRITECSV-"]
        savePNG = values["-SAVEPNG-"]        



        # Check if the specified folder path exists
        if not os.path.exists(folder_path):
            sg.popup_ok("Invalid folder path. Please enter a valid path.", title="Error")
            continue  # Skip the conversion if the path is invalid
        
        

        for file_name in os.listdir(folder_path):
            if file_name.endswith('.out'):
                file_path = os.path.join(folder_path, file_name)

                # Setting the files names
                plottitle = "Calculated Absorbance Spectrum"
                compound = file_name.split('.')
                legend1 = compound[0]
                csvfilename = folder_path+"/spectra_CSV/"+legend1 + "_Calc_Abs.csv"
                PNGfilename = compound[0] + "_Calc_Abs.png"


                x_nm = [*range(start, finish + 1, 1)]
                # Convert in eV
                x_eV = np.asarray([1240 / z for z in x_nm])

                # Define the FWHM (width) in eV
                sigma = width / (2 * np.sqrt(math.log(2)))
                div = 10

                # Excitation energies in EeV and oscillator strengths
                ecc075 = []
                osc075 = []
                spectrum075 = 0

                # Read Q-Chem output and store RPA excitations with strengths
                o = open(file_path, "r")
                rl = o.readlines()
                count = 0
                for line in rl:
                    count += 1
                    if "TDDFT Excitation " in line:
                        # print("here")
                        for i in range(count, len(rl)):
                            if "excitation energy" in rl[i]:
                                ecc075.append(float(rl[i].strip().split()[-1]))
                                # print(i)
                            if "Strength   : " in rl[i]:
                                osc075.append(float(rl[i].strip().split()[-1]))
                o.close()

                # Initialize the spectrum for this file
                spectrum075 = np.zeros(len(x_nm))

                for count, peak in enumerate(ecc075):
                    GaussPeak = GaussSpectrum(x_eV, peak, osc075[count], sigma) / len(ecc075) / div
                    spectrum075 += GaussPeak

                # Append the current spectrum and legend to the respective arrays
                all_spectra.append(spectrum075)
                all_legends.append(legend1)

                # Export data in CSV for each file
                if writeCSV:
                    # Define the output directory for CSV files
                    output_csv_dir = os.path.join(folder_path, 'spectra_CSV')
                    # Check if the directory exists, and create it if it doesn't
                    if not os.path.exists(output_csv_dir):
                        os.makedirs(output_csv_dir)

                    w = open(csvfilename, "w")
                    w.write("Wavelength(nm),Absorbance(arb.u.)\n")
                    for i in range(len(x_nm)):
                        w.write("{:.1f},{:.10f}\n".format(x_nm[i], spectrum075[i]))
                    w.close()


        # Plot all spectra together
        if all_spectra and plotSpectrum:
            plt.figure(figsize=(8, 6))
            plt.xlabel('$\lambda$ / nm')
            plt.ylabel('Intensity')
            for i, (spectrum, legend) in enumerate(zip(all_spectra, all_legends)):
                plt.plot(x_nm, spectrum, label=legend)
            plt.legend()
            plt.title("Calculated Absorbance Spectra")
            if savePNG:
                plt.savefig(folder_path+"/Combined_Absorbance_Spectra.png", format='png', dpi=300, bbox_inches='tight')
            plt.show()
        
        output_text = "Conversion complete!"
        window['-OUTPUT-'].update(output_text)

window.close()