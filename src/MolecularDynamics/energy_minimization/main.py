import sys
import os
import yaml
from pathlib import Path
from datetime import datetime


# Add the parent directory to the Python path
current_dir = Path(__file__).resolve().parent
sys.path.append(str(current_dir.parent))

from energy_minimization.utils import construct_output_paths
from energy_minimization.simulation import MinimizeEnergy

# Function to find all PDB files in a folder
def find_pdb_files(input_folder):
    pdb_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".pdb"):
            pdb_files.append(os.path.join(input_folder, file))
    return pdb_files

def main():
    # Add the current directory to the Python path
    current_dir = Path(__file__).resolve().parent
    sys.path.append(str(current_dir.parent))

    # Load configuration from YAML file
    with open(current_dir / 'config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Find all PDB files in the input folder
    input_folder = config['input_folder']
    pdb_files = find_pdb_files(input_folder)

    # Process each PDB file
    for pdb_file in pdb_files:
        # Configure output paths for each PDB file
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_paths = construct_output_paths(config['output_folder'], base_name, timestamp)

        # Run MinimizeEnergy simulation for each PDB file
        md_simulation = MinimizeEnergy(config, pdb_file, output_paths)
        md_simulation.run()

if __name__ == "__main__":
    main()

