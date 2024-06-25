
import os
import yaml

def load_config(config_path):
    with open(config_path, 'r') as config_file:
        config = yaml.safe_load(config_file)
    return config

def find_output_csv(output_folder):
    for file in os.listdir(output_folder):
        if file.endswith(".csv"):
            return os.path.join(output_folder, file)
    raise FileNotFoundError("No suitable CSV file found in the output folder")
