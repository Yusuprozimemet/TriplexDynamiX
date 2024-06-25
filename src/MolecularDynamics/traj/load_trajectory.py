import os
import mdtraj as md
import yaml

def load_trajectory(config):
    output_folder = config['output_folder']

    # Function to find the output PDB file in the output folder
    def find_output_pdb(output_folder):
        for file in os.listdir(output_folder):
            if file.endswith(".pdb") and "_out_" in file:
                return os.path.join(output_folder, file)
        raise FileNotFoundError("No suitable PDB file found in the output folder")

    # Find the output PDB file
    try:
        pdb_file = find_output_pdb(output_folder)
        print(f"Found PDB file: {pdb_file}")
    except FileNotFoundError as e:
        print(e)
        raise

    # Load the trajectory
    traj = md.load(pdb_file)

    return traj
