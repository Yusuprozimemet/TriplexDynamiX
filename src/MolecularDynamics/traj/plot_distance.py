import yaml
from save_plot_distance import save_plot
from load_trajectory import load_trajectory
from compute_distances import compute_distances

# Load YAML configuration
config_path = r'E:\MolecularDynamicsApp\src\MolecularDynamics\molecular_dynamics\config.yaml'
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

# Load trajectory
traj = load_trajectory(config)

# Compute distances
dist1, dist2, dist3, area = compute_distances(traj)

# Specify output folder for plots
output_folder = r'E:/MolecularDynamicsApp/src/MolecularDynamics/data/plot'

# Save plot results
save_plot(dist1, dist2, dist3, area, output_folder)
