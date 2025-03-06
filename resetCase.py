import os
import shutil

# Define the base directory where caseName is located
caseName = 'valveClose'

run_dir = f'{caseName}/run'  # Replace with the actual path to the caseName directory
for item in os.listdir(run_dir):
    # Get the full path of the item
    item_path = os.path.join(run_dir, item)

    # Check if the item is a directory and is not '0'
    if os.path.isdir(item_path) and item != '0':
        # Remove the directory and its contents
        shutil.rmtree(item_path)
        print(f"Deleted directory: {item_path}")