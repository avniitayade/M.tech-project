#intersectional angle calculations between two adjacent residues in the protein structure in beta strands

import os
from Bio.PDB import PDBParser
import numpy as np

def read_sheet_data(sheet_file):
    # Read sheet data from the text file
    sheet_data = []
    with open(sheet_file, 'r') as file:
        next(file)  # Skip the first line (headers)
        for line in file:
            line = line.strip()
            if line:
                chain_id, residue_num, residue, sheet = line.split()
                sheet_data.append((chain_id, int(residue_num), residue, sheet))
    return sheet_data

def calculate_intersectional_angles(pdb_file, sheet_data):
    # Initialize a PDB parser
    parser = PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("protein", pdb_file)

    # Initialize a list to store intersectional angles and residue names
    intersectional_angles = []

    # Define a list of residue names to skip (glycine and proline)
    skip_residues = ["GLY", "PRO"]

    # Iterate through sheet data, excluding the last one
    for i in range(len(sheet_data) - 1):
        # Extract data from sheet data
        chain_id, residue_num, residue, sheet = sheet_data[i]
        next_chain_id, next_residue_num, next_residue, next_sheet = sheet_data[i + 1]

        # Check if the current residue is from the same sheet as the next one
        if sheet != next_sheet or chain_id != next_chain_id:
            continue

        # Check if the residue is in the skip list
        if residue in skip_residues or next_residue in skip_residues:
            continue

        # Get the chain object for the specified chain ID
        chain = structure[0][chain_id]

        # Select atoms for the current residue
        try:
            atom1 = chain[residue_num]['C']
            atom2 = chain[residue_num]['CA']
            atom3 = chain[residue_num]['CB']
        except KeyError:
            print(f"Skipping residue {residue_num} ({residue}) in chain {chain_id}: Atom 'CB' not found.")
            continue

        # Select atoms for the next residue
        try:
            atom4 = chain[next_residue_num]['N']
            atom5 = chain[next_residue_num]['CA']
            atom6 = chain[next_residue_num]['CB']
        except KeyError:
            print(f"Skipping residue {next_residue_num} ({next_residue}) in chain {next_chain_id}: Atom 'CB' not found.")
            continue

        # Check if all atoms are not None
        if atom1 is not None and atom2 is not None and atom3 is not None \
                and atom4 is not None and atom5 is not None and atom6 is not None:
            # Calculate vectors defining the planes
            vector1 = atom2.coord - atom1.coord
            vector2 = atom5.coord - atom4.coord

            # Calculate the intersectional angle between the two planes
            angle_rad = np.arccos(
                np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
            angle_deg = np.degrees(angle_rad)

            # Append the angle and residue names to the list
            intersectional_angles.append((residue, next_residue, angle_deg))

    return intersectional_angles

# Specify input and output folders
input_pdb_folder = '/home/pooja/project/different_cath_class/mainly_beta_architecture/lessorequalto2.7/2solenoidsheetinfo/2solenoid'  # Replace with your input PDB folder path
input_sheet_folder = '/home/pooja/project/different_cath_class/mainly_beta_architecture/lessorequalto2.7/2solenoidsheetinfo/extractedsheetinfo/'  # Replace with your input extracted sheet info folder path
output_folder = '/home/pooja/project/different_cath_class/mainly_beta_architecture/lessorequalto2.7/2solenoidsheetinfo/intersectional_angles/'  # Replace with your output folder path

# Ensure the output folder exists, create it if not
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Iterate through files in the input PDB folder
for sheet_filename in os.listdir(input_sheet_folder):
    if sheet_filename.endswith('_sheet_data_extracted_info.txt'):
        pdb_identifier = sheet_filename.split('_')[0]  # Extract PDB identifier from the filename

        # Find the corresponding PDB file
        pdb_filename = f"{pdb_identifier}.pdb"
        pdb_file = os.path.join(input_pdb_folder, pdb_filename)

        # Check if the corresponding PDB file exists
        if os.path.exists(pdb_file):
            sheet_file = os.path.join(input_sheet_folder, sheet_filename)
            output_file = os.path.join(output_folder, f"{pdb_identifier}_angle.txt")

            sheet_data = read_sheet_data(sheet_file)
            angles = calculate_intersectional_angles(pdb_file, sheet_data)

            # Write the calculated intersectional angles to a text file in table format
            with open(output_file, 'w') as file:
                file.write("Residue1\tResidue2\tAngle\n")
                for angle in angles:
                    file.write(f"{angle[0]}\t{angle[1]}\t{angle[2]}\n")

            print(f"Intersectional angles written to {output_file}")
        else:
            print(f"No corresponding PDB file found for {sheet_filename}. Skipping.")
