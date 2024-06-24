import sys
import os
import yaml
from pathlib import Path

# Add the current directory to the Python path
current_dir = Path(__file__).resolve().parent
sys.path.append(str(current_dir.parent))

# Now import your module
from molecular_dynamics.simulation import MolecularDynamicsSimulation

# Load configuration from YAML file
with open(current_dir / 'config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Run the simulation
md_simulation = MolecularDynamicsSimulation(config)
md_simulation.run()


