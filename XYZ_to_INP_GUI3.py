import os
import PySimpleGUI as sg

def get_additional_text(solvent, dielectric, optical_dielectric, omega, convergence, cycles):
    if solvent == "No solvent":
        return f"""$end

$rem
   BASIS  =  LACVP
   ECP  =  fit-LACVP
   GUI  =  2
   JOB_TYPE  =  Optimization
   METHOD  =  B3LYP
   SCF_CONVERGENCE  =  {convergence}
   SCF_MAX_CYCLES  =  {cycles}
$end


@@@

$molecule
read
$end

$rem
   BASIS  =  LACVP
   CIS_MOMENTS  =  1
   CIS_MULLIKEN  =  1
   CIS_N_ROOTS  =  10
   ECP  =  fit-LACVP
   EXCHANGE  =  CAM-B3LYP
   GUI  =  2
   JOB_TYPE  =  SP
   NTO_PAIRS  =  2
   OMEGA  =  {omega}
   RPA  =  1
   SCF_CONVERGENCE  =  {convergence}
   SCF_MAX_CYCLES  =  {cycles}
   STS_MOM  =  1
$end"""

    else:
        return f"""$end

$rem
   BASIS  =  LACVP
   ECP  =  fit-LACVP
   GUI  =  2
   JOB_TYPE  =  Optimization
   METHOD  =  B3LYP
   SCF_CONVERGENCE  =  {convergence}
   SCF_MAX_CYCLES  =  {cycles}
   SOLVENT_METHOD  =  PCM
$end


$pcm
   THEORY  CPCM   
$end

$solvent
   DIELECTRIC  {dielectric}
   OPTICALDIELECTRIC  {optical_dielectric}   
$end



@@@

$molecule
read
$end

$rem
   BASIS  =  LACVP
   CIS_MOMENTS  =  1
   CIS_MULLIKEN  =  1
   CIS_N_ROOTS  =  10
   ECP  =  fit-LACVP
   EXCHANGE  =  CAM-B3LYP
   GUI  =  2
   JOB_TYPE  =  SP
   NTO_PAIRS  =  2
   OMEGA  =  {omega}
   RPA  =  1
   SCF_CONVERGENCE  =  {convergence} 
   SCF_MAX_CYCLES  =  {cycles}
   SOLVENT_METHOD  =  PCM
   STS_MOM  =  1
$end

$pcm
  THEORY  CPCM
$end

$solvent
  DIELECTRIC  {dielectric}
  OPTICALDIELECTRIC {optical_dielectric}
$end"""

solvent_options = ["No solvent", "Toluene", "THF", "Other"]

layout = [
    [sg.Text("Convert all your XYZ files into Q-Chem input files.\nThe input file will do a geometry optimization (B3LYP/LACVP) \nand then TD-DFT (CAM-B3LYP/LACVP)")],
    [sg.Text("Enter folder path")],
    [sg.Input(key="-FOLDER-"), sg.Button("Convert")], 
    [sg.Text("", key="-OUTPUT-")],
    [sg.Text("Convergence", size=15), sg.Input(default_text="8", key="-CONVERGENCE-", justification='left', size=5)],
    [sg.Text("Cycles",size=15), sg.Input(default_text="200", key="-CYCLES-", justification='left', size=5)],   
    [sg.Text("Omega",size=15), sg.Input(default_text="0.145",key="-OMEGA-",justification='left',size=5)],
    [sg.Text("Solvent",size=15), sg.Combo(solvent_options, default_value="No solvent", key="-SOLVENT-", enable_events=True)],
    [sg.Text("Dielectric",size=15, visible=False, key="-DIELECTRIC-TEXT-"), sg.Input(default_text="0",visible=False, key="-DIELECTRIC-",size=5)],
    [sg.Text("Optical Dielectric",size=15, visible=False, key="-OPTDIELECTRIC-TEXT-"), sg.Input(default_text="0",visible=False, key="-OPTDIELECTRIC-",size=5)]
]

window = sg.Window("XYZ to QChem INP Converter", layout)

while True:
    event, values = window.read()
    
    if event == sg.WIN_CLOSED:
        break

    if values["-SOLVENT-"] == "Other":
        window["-DIELECTRIC-TEXT-"].update(visible=True) 
        window["-OPTDIELECTRIC-TEXT-"].update(visible=True)
        window["-DIELECTRIC-"].update(visible=True)
        window["-OPTDIELECTRIC-"].update(visible=True)
    else:
        window["-DIELECTRIC-"].update(visible=False)
        window["-OPTDIELECTRIC-"].update(visible=False)
        window["-DIELECTRIC-TEXT-"].update(visible=False)
        window["-OPTDIELECTRIC-TEXT-"].update(visible=False)

    if event == "Convert":
        folder_path = values["-FOLDER-"]
        
        # Check if the specified folder path exists
        if not os.path.exists(folder_path):
            sg.popup_ok("Invalid folder path. Please enter a valid path.", title="Error")
            continue  # Skip the conversion if the path is invalid

        for file_name in os.listdir(folder_path):
            if file_name.endswith('.xyz'):
                file_path = os.path.join(folder_path, file_name)
                
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    
                lines[0] = "$molecule\n"
                lines[1] = "0 1\n"
                

                omega = values["-OMEGA-"].lstrip("0.")
                convergence = values["-CONVERGENCE-"]
                cycles = values["-CYCLES-"]
                solvent = values["-SOLVENT-"]
                if solvent == "Other":
                    dielectric = values["-DIELECTRIC-"]
                    optical_dielectric = values["-OPTDIELECTRIC-"]
                elif solvent == "Toluene":
                    dielectric = "2.3800"
                    optical_dielectric = "2.2410"
                elif solvent == "THF":
                    dielectric = "7.5800"
                    optical_dielectric = "2.0200"
                else:
                    dielectric = "0"
                    optical_dielectric = "0"    

                additional_text = get_additional_text(solvent, dielectric, optical_dielectric, omega, convergence, cycles)
                lines.append(additional_text)
                
                new_file_path = file_path.replace('.xyz', '.inp')
                with open(new_file_path, 'w') as f:
                    f.writelines(lines)
                    
                print(f"Modified file saved as: {new_file_path}")

        output_text = "Conversion complete!"
        window['-OUTPUT-'].update(output_text)
        
window.close()