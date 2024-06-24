import os

def find_pdb_file(input_folder):
    for file in os.listdir(input_folder):
        if file.endswith(".pdb"):
            return os.path.join(input_folder, file)
    raise FileNotFoundError("No PDB file found in the input folder")

def construct_output_paths(output_folder, base_name, timestamp):
    output_init_with_hydrogens_filename = f'{base_name}_with_hydrogens_{timestamp}.pdb'
    output_pdb_filename = f'{base_name}_out_{timestamp}.pdb'
    output_dcd_filename = f'{base_name}_traj_{timestamp}.dcd'
    output_csv_filename = f'{base_name}_{timestamp}.csv'

    return {
        'output_init_with_hydrogens_path': os.path.join(output_folder, output_init_with_hydrogens_filename),
        'output_pdb_path': os.path.join(output_folder, output_pdb_filename),
        'output_dcd_path': os.path.join(output_folder, output_dcd_filename),
        'output_csv_path': os.path.join(output_folder, output_csv_filename),
    }
