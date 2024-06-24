import csv
import yaml
import os

def load_config(config_file):
    # Load the configuration from the YAML file
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def clean_csv(input_file, output_file):
    try:
        # Open the input CSV file for reading
        with open(input_file, 'r') as infile:
            reader = csv.reader(infile)
            # Skip the initial headers (repeated in every line)
            next(reader)
            
            # Prepare the list to hold the cleaned data
            cleaned_data = []

            # Process each row in the input file
            for row in reader:
                # Skip rows that are just repeated headers
                if row[0] == 'pdb_file' and row[1] == 'minimized_energy':
                    continue
                # Add the cleaned row to the list
                cleaned_data.append(row)

        # Open the output CSV file for writing
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            # Write the single header row
            writer.writerow(['pdb_file', 'minimized_energy'])
            # Write the cleaned data rows
            writer.writerows(cleaned_data)

        print(f"Successfully processed and saved: {output_file}")

    except PermissionError:
        print(f"Permission denied: '{output_file}'. Unable to write the output file.")
    except Exception as e:
        print(f"An error occurred while processing '{input_file}': {e}")

def main():
    # Load the configuration
    config_file = 'config.yaml'
    config = load_config(config_file)

    # Get the input folder path from the config
    input_folder = config['rawData']
    # Define the output folder path (same as input folder or a subfolder)
    output_folder = os.path.join(input_folder, 'cleaned')

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Process each CSV file in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith('.csv'):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)  # Save cleaned files in the output folder
            clean_csv(input_file, output_file)

if __name__ == '__main__':
    main()
