import numpy as np
import matplotlib.pyplot as plt
import os

from colour.plotting import plot_chromaticity_diagram_CIE1931
from colour.plotting import render
from PIL import Image



def interpolate_spectrum(x, y):
    # Find the nearest integer wavelengths
    x_int = np.round(x).astype(int)

    # Remove duplicate integer wavelengths and corresponding intensities
    unique_x, unique_indices = np.unique(x_int, return_index=True)
    unique_y = y[unique_indices]

    # Perform linear interpolation
    interpolated_x = np.arange(unique_x[0], unique_x[-1] + 1)
    interpolated_y = np.interp(interpolated_x, unique_x, unique_y)

    return interpolated_x, interpolated_y

def calculate_cie_coordinates(wavelength, intensity, x_bar, y_bar, z_bar):
    # Normalize intensity
    normalized_intensity = intensity / np.max(intensity)

    # Calculate tristimulus values
    x = np.sum(normalized_intensity * x_bar)
    y = np.sum(normalized_intensity * y_bar)
    z = np.sum(normalized_intensity * z_bar)

    # Calculate CIE coordinates
    cie_x = x / (x + y + z)
    cie_y = y / (x + y + z)
    cie_z = z / (x + y + z)

    return cie_x, cie_y, cie_z

def read_emission_spectrum(file_path, cmf_file_path, triangles):
    try:
        # Read the data from the emission spectrum file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find the starting line of the spectrum data
        try:
            start_index = lines.index("X\tY\n") + 1
        except ValueError:
            raise ValueError("Invalid emission spectrum file format. Missing header line 'X\tY'.")

        # Find the ending line of the spectrum data
        try:
            end_index = lines.index("\n", start_index)
        except ValueError:
            end_index = len(lines)

        # Extract the wavelength and intensity data
        try:
            data = [line.split() for line in lines[start_index:end_index] if len(line.split()) == 2]
            wavelength, intensity = zip(*data)
            wavelength = np.array(wavelength, dtype=float)
            intensity = np.array(intensity, dtype=float)
        except (ValueError, IndexError):
            raise ValueError("Invalid emission spectrum file format. Unable to parse data.")

        # Interpolate the spectrum
        interpolated_wavelength, interpolated_intensity = interpolate_spectrum(wavelength, intensity)

        # Read the color matching functions from the file
        try:
            cmf_data = np.loadtxt(cmf_file_path, delimiter=',', skiprows=1)
            cmf_wavelength = cmf_data[:, 0]
            cmf_red = cmf_data[:, 1]
            cmf_green = cmf_data[:, 2]
            cmf_blue = cmf_data[:, 3]
        except FileNotFoundError:
            raise FileNotFoundError(f"Color matching functions file '{cmf_file_path}' not found.")
        except ValueError:
            raise ValueError("Invalid color matching functions file format.")

        # Interpolate the color matching functions to match the interpolated spectrum wavelengths
        interpolated_cmf_red = np.interp(interpolated_wavelength, cmf_wavelength, cmf_red)
        interpolated_cmf_green = np.interp(interpolated_wavelength, cmf_wavelength, cmf_green)
        interpolated_cmf_blue = np.interp(interpolated_wavelength, cmf_wavelength, cmf_blue)

        # Calculate CIE coordinates
        cie_x, cie_y, cie_z = calculate_cie_coordinates(interpolated_wavelength, interpolated_intensity,
                                                        interpolated_cmf_red, interpolated_cmf_green,
                                                        interpolated_cmf_blue)


        # Plot the original and interpolated spectra
        # plt.plot(wavelength, intensity, label='Original Spectrum')
        # plt.plot(interpolated_wavelength, interpolated_intensity, label='Spectrum')
        # plt.xlabel('Wavelength (nm)')
        # plt.ylabel('Intensity')
        # plt.title('Spectrum')
        # plt.legend()
        # plt.show()

        # Print the CIE coordinates
        print("CIE Coordinates:")
        print(f"x: {cie_x}")
        print(f"y: {cie_y}")
        print(f"z: {cie_z}")

        # Plotting the *CIE 1931 Chromaticity Diagram*.
        # The argument *standalone=False* is passed so that the plot doesn't get
        # displayed and can be used as a basis for other plots.
        plot_chromaticity_diagram_CIE1931(standalone=False)
        
        # Plotting the *xy* chromaticity coordinates.
        x = cie_x
        y = cie_y
        xy = (x,y)

        # RGB coordinates Rec.2020
        redx = 0.708
        redy = 0.292
        grex = 0.170
        grey = 0.797
        blux = 0.131
        bluy = 0.046
        # RGB coordinates sRGB
        Rx = 0.64
        Ry = 0.33
        Gx = 0.30
        Gy = 0.60
        Bx = 0.15
        By = 0.06

        #CIE from spectrum
        plt.plot(x, y, 'o-', color='black')
        if (triangles):
            #points Rec2020
            plt.plot(redx, redy, 'ko')
            plt.plot(grex, grey, 'ko')
            plt.plot(blux, bluy, 'ko')
            #triangle Rec2020
            plt.plot([redx,grex],[redy,grey],'k--')
            plt.plot([grex,blux],[grey,bluy],'k--')
            plt.plot([blux,redx],[bluy,redy],'k--')
            
            #points sRGB
            plt.plot(Rx, Ry, 'o-', color='grey')
            plt.plot(Gx, Gy, 'o-', color='grey')
            plt.plot(Bx, By, 'o-', color='grey')
            #triangle sRGB
            plt.plot([Rx,Gx],[Ry,Gy],color='grey')
            plt.plot([Gx,Bx],[Gy,By],color='grey')
            plt.plot([Bx,Rx],[By,Ry],color='grey')

        # Annotating the plot.
        plt.annotate(extract_solvent_name(file_path),
                      xy=xy,
                      xytext=(-50, 30),
                      textcoords='offset points',
                      arrowprops=dict(arrowstyle='->', connectionstyle='arc3, rad=-0.2'))

        # Add CIE X and CIE Y labels
        plt.text(0.6, 0.9, f"CIE X: {cie_x:.3f}", ha='left', transform=plt.gca().transAxes)
        plt.text(0.6, 0.85, f"CIE Y: {cie_y:.3f}", ha='left', transform=plt.gca().transAxes)
        if (triangles):
            plt.text(0.6, 0.70, "Rec2020 ---", ha='left', transform=plt.gca().transAxes)
            plt.text(0.6, 0.65, "sRGB", ha='left', color='grey', transform=plt.gca().transAxes)

        # Set the x-axis range
        plt.xlim(-.1, 1)  # Replace x_min and x_max with your desired range
        plt.ylim(-.1, 1)  # Replace x_min and x_max with your desired range

        # Save the CIE diagram as a PNG
        plt.savefig("CIE_diagram.jpg", dpi=300, bbox_inches='tight', pad_inches=0)

        
        # Displaying the plot.
        render(
            standalone=True,
            limits=(-0.1, 0.9, -0.1, 0.9),
            x_tighten=True,
            y_tighten=True);


    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def extract_solvent_name(filename):
    # Remove the file extension if present
    filename_without_extension = filename.split('.')[0]
    # Extract the solvent name by removing any digits at the end
    solvent_name = ''.join([char for char in filename_without_extension if not char.isdigit()])
    return solvent_name


def process_folder(folder_path, cmf_file_path, showtri):
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(folder_path, filename)
            solvent_name = extract_solvent_name(filename)
            print(f"Processing file: {file_path}")
            read_emission_spectrum(file_path, cmf_file_path, showtri)



# Get the current working directory
current_directory = os.getcwd()

# Construct the path to the data folder
folder_name = "spectra"
folder_path = os.path.join(current_directory, folder_name)

cmf_file_path = "CMF.txt"
showtri = True

# Process all files in the folder
process_folder(folder_path, cmf_file_path, showtri)



# Specify the file paths
# file_path = "Toluene298.txt"
# cmf_file_path = "CMF.txt"
# showtri = True

# Read and interpolate the emission spectrum, and calculate the CIE coordinates
# read_emission_spectrum(file_path, cmf_file_path, showtri)



