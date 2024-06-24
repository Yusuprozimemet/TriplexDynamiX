import os
import csv
import shutil
import yaml

def load_config(config_file):
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def combine_csv_files(output_folder, destination_folder):
    # Adjust to get the correct path for output_folder
    output_folder = os.path.abspath(output_folder)
    
    # Find all CSV files in the output folder
    csv_files = [filename for filename in os.listdir(output_folder) if filename.endswith('.csv')]
    
    combined_csv_path = os.path.join(destination_folder, 'combined_minimized_energies.csv')
    
    # Write header to the combined CSV file
    header_written = False
    
    with open(combined_csv_path, 'w', newline='') as combined_file:
        writer = csv.writer(combined_file)
        
        for csv_file in csv_files:
            csv_file_path = os.path.join(output_folder, csv_file)
            
            with open(csv_file_path, 'r', newline='') as individual_file:
                reader = csv.reader(individual_file)
                
                # Skip header in all but the first CSV file
                if not header_written:
                    header = next(reader)
                    writer.writerow(header)
                    header_written = True
                
                # Write remaining rows
                for row in reader:
                    writer.writerow(row)
    
    print(f"Combined CSV file saved to: {combined_csv_path}")
    
    # Move the combined CSV file to the destination folder
    destination_path = os.path.join(destination_folder, 'combined_minimized_energies.csv')
    shutil.move(combined_csv_path, destination_path)
    print(f"Moved combined CSV file to: {destination_path}")

def main():
    # Adjust these paths to match your folder structure
    config_file = 'E:\\MolecularDynamicsApp\\src\\MolecularDynamics\\data_preparation\\config.yaml'
    config = load_config(config_file)
    
    # Get output folder path from config
    output_folder = config.get('output_folder')
    if output_folder is None:
        print("Error: 'output_folder' not defined in config.yaml.")
        return
    
    # Adjust destination_folder path
    destination_folder = 'E:\\MolecularDynamicsApp\\src\\MolecularDynamics\\data\\preparedData'
    
    combine_csv_files(output_folder, destination_folder)

if __name__ == "__main__":
    main()
