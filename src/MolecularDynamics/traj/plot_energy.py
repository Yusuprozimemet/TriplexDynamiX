# plot_energy.py
import os
import pandas as pd
import matplotlib.pyplot as plt
from load_csv import load_config, find_output_csv

def plot_energy(config_path, output_folder):
    # Load YAML configuration
    config = load_config(config_path)
    output_folder = os.path.join(output_folder, 'plots')  # Ensure output folder for plots exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Find the output CSV file
    try:
        csv_file = find_output_csv(config['output_folder'])
        print(f"Found CSV file: {csv_file}")
    except FileNotFoundError as e:
        print(e)
        raise
    
    # Read CSV file
    df = pd.read_csv(csv_file)
    potential_energy = df['Potential Energy (kJ/mole)']
    
    # Plot potential energy
    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(potential_energy, color='orange', label="Potential Energy")
    ax.set_title('Potential Energy over Time')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Potential Energy (kJ/mole)')
    ax.legend(loc='best')
    
    # Save plot
    plot_filename = os.path.join(output_folder, 'potential_energy_plot.png')
    plt.savefig(plot_filename)
    print(f"Plot saved: {plot_filename}")
    
    # Show plot
    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    config_path = r'E:\MolecularDynamicsApp\src\MolecularDynamics\molecular_dynamics\config.yaml'
    output_folder = r'E:/MolecularDynamicsApp/src/MolecularDynamics/data/plot'
    plot_energy(config_path, output_folder)
