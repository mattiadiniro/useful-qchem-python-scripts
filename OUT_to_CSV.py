"""
This script will automatically extract and calculate several properties of CMA complexes and create an organized list in a CSV file. 
It only works with CMA containing one metal center, Cu, Ag, or Au. 
The donor must be an amide (N must be bonded to M), and the acceptor must be a carbene (C must be bonded to M). 
When specifying the folder path, the script will automatically ignore input files, fchk files, and also the slurm-24652452.out files so don't worry about those.
"""


import re
import csv
import math
import os
import PySimpleGUI as sg


# Set the folder path containing the .out files
#folder_path = 'C:/Users/12089/Documents/USC/_COMPUTATIONS/QCHEM/12_diff_donors/PZI/PZI_outputs'


""" 
*** FUNCTIONS ***
"""

# Calculate the dihedral angle between four atoms
def calculate_dihedral_angle(atom1, atom2, atom3, atom4):
    x1, y1, z1 = map(float, atom1)
    x2, y2, z2 = map(float, atom2)
    x3, y3, z3 = map(float, atom3)
    x4, y4, z4 = map(float, atom4)

    v1 = [x2 - x1, y2 - y1, z2 - z1]
    v2 = [x3 - x2, y3 - y2, z3 - z2]
    v3 = [x4 - x3, y4 - y3, z4 - z3]

    n1 = [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]
    n2 = [v2[1]*v3[2] - v2[2]*v3[1], v2[2]*v3[0] - v2[0]*v3[2], v2[0]*v3[1] - v2[1]*v3[0]]

    dot = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
    norm_n1 = math.sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
    norm_n2 = math.sqrt(n2[0]**2 + n2[1]**2 + n2[2]**2)

    angle = math.degrees(math.acos(dot / (norm_n1 * norm_n2)))
    return angle

# Calculate the angle between three atoms
def calculate_angle(atom1, atom2, atom3):
    x1, y1, z1 = map(float, atom1)
    x2, y2, z2 = map(float, atom2)
    x3, y3, z3 = map(float, atom3)
    
    vector1 = [x1 - x2, y1 - y2, z1 - z2]
    vector2 = [x3 - x2, y3 - y2, z3 - z2]
    
    dot_product = sum(a * b for a, b in zip(vector1, vector2))
    magnitude1 = math.sqrt(sum(a * a for a in vector1))
    magnitude2 = math.sqrt(sum(a * a for a in vector2))
    
    angle_rad = math.acos(dot_product / (magnitude1 * magnitude2))
    angle_deg = math.degrees(angle_rad)
    
    return angle_deg

# Calculate the distance between two atoms
def calculate_distance(atom1, atom2):
    x1, y1, z1 = map(float, atom1)
    x2, y2, z2 = map(float, atom2)
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

# Extract the final energy from the file
def extract_final_energy(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        match = re.search(r'Final energy is\s+(-?\d+\.\d+)', content)
        if match:
            energy = float(match.group(1))
            return energy
        else:
            return None

# Calculate the final energy of the molecule
def calculate_final_energy_ev(final_energy_h):
    return final_energy_h * 27.211383971279

# Extract and Calculate the MOs
def extract_MOs(file_path):
    numbers = []
    optimization_converged = False
    occupied_found = False
    lumo = None

    with open(file_path, 'r') as file:
        for line in file:
            if not optimization_converged:
                if '**  OPTIMIZATION CONVERGED  **' in line:
                    optimization_converged = True
            elif not occupied_found:
                if '-- Occupied --' in line:
                    occupied_found = True
            else:
                if '-- Virtual --' in line:
                    lumo_line = next(file).strip()
                    lumo = float(lumo_line.split()[0])
                    break
                else:
                    line = line.strip()
                    numbers += line.split()
        if len(numbers) > 0:
            homo = float(numbers[-1])
    return homo, lumo

# Saving all the data to a CSV file
def save_to_csv(writer, name, dihedral, angle, cdist, ndist, cndist, homo, homoev, lumo, lumoev, gapev, gapnm, final_energy_h, singlet, osc, triplet, dest, singlet2, osc2, singlet3, osc3):
    final_energy_ev = calculate_final_energy_ev(final_energy_h)

    writer.writerow({
        'name': name,
        'dihedral': dihedral,
        'angle': angle,
        'C-M': cdist,
        'M-N': ndist,
        'C---N': cndist,
        'HOMO': homo,
        'HOMO(eV)': homoev,
        'LUMO': lumo,
        'LUMO(eV)': lumoev,
        'LU-HO(eV)': gapev,
        'LU-HO(nm)': gapnm,
        'Final E(H)': final_energy_h,
        'Final E(eV)': final_energy_ev,
        'S1': singlet, 
        'osc1': osc, 
        'T1': triplet, 
        'DEST(eV)': dest,
        'S2': singlet2,
        'osc2': osc2,
        'S3': singlet3,
        'osc3': osc3
    })


""" 
*** OPENING FILE ***
"""

layout = [
    [sg.Text("Analyze all your Q-Chem OUT files into a CSV file with all the important properties.")],
    [sg.Text("This will only work on mononuclear CMAs.")],
    [sg.Text("Enter folder path")],
    [sg.Input(key="-FOLDER-"), sg.Button("Convert")], 
    [sg.Text("", key="-OUTPUT-")]
]

window = sg.Window("Q-Chem OUT to CSV Converter", layout)

while True:
    event, values = window.read()
    
    if event == sg.WIN_CLOSED:
        break

    if event == "Convert":
        folder_path = values["-FOLDER-"]
        
        # Check if the specified folder path exists
        if not os.path.exists(folder_path):
            sg.popup_ok("Invalid folder path. Please enter a valid path.", title="Error")
            continue  # Skip the conversion if the path is invalid

        # Get the folder name from the folder path
        folder_name = os.path.basename(folder_path)

        # Create the CSV file path
        csv_file_path = os.path.join(folder_path, folder_name + '_properties.csv')

        # Open the CSV file for writing
        with open(csv_file_path, 'w', newline='') as csv_file:
            fieldnames = ['name', 'dihedral', 'angle', 'C-M', 'M-N', 'C---N', 'HOMO', 'HOMO(eV)', 'LUMO', 'LUMO(eV)', 'LU-HO(eV)', 'LU-HO(nm)', 'Final E(H)', 'Final E(eV)', 'S1', 'osc1', 'T1', 'DEST(eV)', 'S2', 'osc2', 'S3', 'osc3']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()

            # Iterate over the .out files in the folder
            for file_name in os.listdir(folder_path):
                if file_name.endswith('.out') and not file_name.startswith('slurm'):
                    file_path = os.path.join(folder_path, file_name)

                    #file name
                    name = file_name.split(".", 1)[0]

                    #constants
                    HtoeV = 27.211383971279
                    eVtonm = 1239.84193
                                
                    
                    """ 
                    *** COORDINATES ***
                    """
                    
                    # Read the file and extract the coordinates
                    coordinates = []
                    with open(file_path, 'r') as file:
                        lines = file.readlines()
                    
                        # Find the line with "OPTIMIZATION CONVERGED"
                        found_converged = False
                        for i, line in enumerate(lines):
                            if "OPTIMIZATION CONVERGED" in line:
                                start_index = i + 5  # Skip the next 4 lines
                                found_converged = True
                                break
            
                        # Skip the file if "OPTIMIZATION CONVERGED" is not found
                        if not found_converged:
                            save_to_csv(writer,name,'optimization did not converge',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
                            continue
                    
                        # Extract the coordinates
                        for j in range(start_index, len(lines)):
                            line = lines[j].strip().split()
                    
                            # Check if the line is empty
                            if not line:
                                break
                    
                            coordinates.append(line[1:])  # Append the coordinates
                    
                    # Separate C and N atoms
                    c_atoms = []
                    n_atoms = []
                    M_atom = []
                    
                    for atom in coordinates:
                        if atom[0] == "C":
                            c_atoms.append(atom[1:])
                        elif atom[0] == "N":
                            n_atoms.append(atom[1:])
                        elif (atom[0] == "Cu" or atom[0] == "Ag" or atom[0] == "Au") :
                            M_atom.append(atom[1:])
                    
                    
                    # Find the closest C atom to the M atom
                    closest_c_atom = None
                    closest_c_distance = float('inf')
                    
                    for c_atom in c_atoms:
                        distance = calculate_distance(M_atom[0], c_atom)
                        if distance < closest_c_distance:
                            closest_c_atom = c_atom
                            closest_c_distance = distance
                    
                    # Find the closest N atom to the M atom
                    closest_n_atom = None
                    closest_n_distance = float('inf')
                    
                    for n_atom in n_atoms:
                        distance = calculate_distance(M_atom[0], n_atom)
                        if distance < closest_n_distance:
                            closest_n_atom = n_atom
                            closest_n_distance = distance
                    
                    # Calculate the distance between the closest C and N atoms
                    closest_c_n_distance = calculate_distance(closest_c_atom, closest_n_atom)
                    
                    # Find the closest atom to the C atom (excluding itself and M atom)
                    closest_atom_to_c = None
                    closest_distance_to_c = float('inf')
                    
                    
                    for atom in coordinates:
                        atom_coordinates = atom[1:]  # Extract coordinates from the atom
                    
                        if atom_coordinates != closest_c_atom and atom_coordinates != M_atom[0]:
                            distance = calculate_distance(closest_c_atom, atom_coordinates)
                            if distance < closest_distance_to_c:
                                closest_atom_to_c = atom_coordinates
                                closest_distance_to_c = distance
                    
                    # Find the closest atom to the N atom (excluding itself and M atom)
                    closest_atom_to_n = None
                    closest_distance_to_n = float('inf')
                    
                    for atom in coordinates:
                        atom_coordinates = atom[1:]  # Extract coordinates from the atom
                    
                        if atom_coordinates != closest_n_atom and atom_coordinates != M_atom[0]:
                            distance = calculate_distance(closest_n_atom, atom_coordinates)
                            if distance < closest_distance_to_n:
                                closest_atom_to_n = atom_coordinates
                                closest_distance_to_n = distance
                    
                    # Calculate angles and bond lengths
                    dihedral_angle = calculate_dihedral_angle(closest_atom_to_c, closest_c_atom, closest_n_atom, closest_atom_to_n)
                    if (abs(dihedral_angle) > 90):
                        dihedral_angle = 180 - abs(dihedral_angle)
                    angle_c_Au_n = calculate_angle(closest_c_atom, M_atom[0], closest_n_atom)
                    closest_c_n_distance = calculate_distance(closest_c_atom, closest_n_atom)
                    
                    
                    """ 
                    *** MOLECULAR ORBITALS ***
                    """
                    
                    #MOs
                    homo, lumo = extract_MOs(file_path)
                    homoev = homo * HtoeV
                    lumoev = lumo * HtoeV
                    gapev = lumoev - homoev
                    gapnm = eVtonm / gapev
                    
                    """ 
                    *** EXCITED STATES ***
                    """
                    # Excitation energies in eV and osc075illator strengths
                    ecc = []
                    osc = []
                    spectrum = 0
                    
                    # Read Q-Chem output and store RPA excitations with stregths
                    o = open( file_path, "r")
                    rl = o.readlines()
                    count = 0
                    for line in rl:
                        count +=1
                        if "TDDFT Excitation " in line:
                            #print("here")
                            for i in range(count,len(rl)):
                                if "excitation energy" in rl[i]:
                                    ecc.append(float(rl[i].strip().split()[-1]))
                                    #print(i)
                                if "Strength   : " in rl[i]:
                                    osc.append(float(rl[i].strip().split()[-1]))
                    o.close()
                    
                    # Create singlets, triplets, and osc lists
                    singlets = []
                    triplets = []
                    osc_list = []
                    for i in range(len(osc)):
                        if osc[i] == 0:
                            triplets.append(ecc[i])
                            #osc_list.append("")
                        else:
                            singlets.append(ecc[i])
                            osc_list.append(osc[i])
                    
                    STgap = (singlets[0] - triplets[0])
                    
                    final_energy_h = extract_final_energy(file_path)

                    
                                # Write data to the CSV file
                    save_to_csv(writer, name, dihedral_angle, angle_c_Au_n, closest_c_distance, closest_n_distance, closest_c_n_distance, homo, homoev, lumo, lumoev, gapev, gapnm, final_energy_h, singlets[0], osc_list[0], triplets[0], STgap, singlets[1], osc_list[1], singlets[2], osc_list[2])

        print(f"Data saved to {csv_file_path}")
        output_text = "Conversion complete!"
        window['-OUTPUT-'].update(output_text)
        
window.close()
