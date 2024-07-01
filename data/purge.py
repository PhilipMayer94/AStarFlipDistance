import os
import re

def rename_files(folder_path):
    # Get a list of all files in the folder
    file_list = os.listdir(folder_path)

    # Regular expression pattern for matching the old file name format
    pattern = r"t_(\d{4})_(\d{1,2})_(\d{1,2})\.csv"

    # Iterate through each file
    for file_name in file_list:
        # Check if the file name matches the pattern
        match = re.match(pattern, file_name)
        if match:
            # Extract the year, month, and day from the matched groups
            year = match.group(1)
            month = match.group(2)
            day = match.group(3)

            # Construct the new file name
            new_file_name = f"t_{year}_{month}_42.csv"

            # Get the current file's full path
            current_path = os.path.join(folder_path, file_name)

            # Get the new file's full path
            new_path = os.path.join(folder_path, new_file_name)

            # Rename the file
            os.rename(current_path, new_path)

def rename_folders(root_path, depth):
    if depth < 0 or not os.path.exists(root_path):
        return

    # Iterate through all items (files and folders) in the root path
    for item_name in os.listdir(root_path):
        item_path = os.path.join(root_path, item_name)

        # Check if the item is a folder
        if os.path.isdir(item_path):
            # Check if the folder name starts with "start"
            if item_name.startswith('start'):
                # Rename the folder
                new_name = 'triangulation' + item_name[5:]
                new_path = os.path.join(root_path, new_name)
                os.rename(item_path, new_path)

            # Recursively process subfolders with reduced depth
            rename_folders(item_path, depth - 1)

def delete_files_with_name(start_folder, file_name):
    for root, dirs, files in os.walk(start_folder):
        for file in files:
            if file == file_name:
                file_path = os.path.join(root, file)
                os.remove(file_path)
                print(f"Deleted file: {file_path}")



# Provide the folder path
folder_path = 'random_experiments_paper'
file_name = "results_eppstein.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_simple.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_heuristic.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_bfs.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_ilp.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_combined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_rootcombined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_decomposition.csv"
delete_files_with_name(folder_path, file_name)

folder_path = 'sealevel_experiments_paper'
file_name = "results_eppstein.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_simple.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_heuristic.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_bfs.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_ilp.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_combined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_rootcombined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_decomposition.csv"
delete_files_with_name(folder_path, file_name)

folder_path = 'random_hod_experiments_paper'
file_name = "results_eppstein.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_simple.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_heuristic.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_bfs.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_ilp.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_combined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_rootcombined.csv"
delete_files_with_name(folder_path, file_name)
file_name = "results_decomposition.csv"
delete_files_with_name(folder_path, file_name)


