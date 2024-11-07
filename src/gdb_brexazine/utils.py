import pandas as pd
import ast
import re
import os
import openbabel as ob
from openbabel import pybel

#RDKit:
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem.rdmolops import *
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

FILE_PATH_READ = "file_path_read"
FILE_PATH_SAVE = "file_path_save"

class Utils:

    def __init__(self):
        # Initialize any instance variables here, if needed
        pass

    # TODO: append all the dataframes from a folder
    def append_dfs_in_folder(self, FOLDER_PATH, COLUMN_NAME_INPUT, SEPRATOR='\t') -> pd.DataFrame:

        iteration = 0
        row_count = 0

        for filename in os.listdir(FOLDER_PATH):
            
            if filename.endswith(".smi"):
                iteration += 1
                if iteration == 1:
                    df = pd.read_csv(FOLDER_PATH + filename, sep=SEPRATOR, names=COLUMN_NAME_INPUT)
                    #print("basis lenght = " + str(len(df)))
                    row_count += len(df)

                if iteration != 1:
                    df2 = pd.read_csv(FOLDER_PATH + filename, sep=SEPRATOR,names=COLUMN_NAME_INPUT)
                    #print("added lenght = " + str(len(df2)))
                    df = pd.concat([df, df2], axis=0, ignore_index=True)
                    row_count += len(df2)

        print("Total enteries that have been merged: \t" + str(row_count))

        return df
    

    # TODO: Load the data
    def load_data(self, FILE_PATH_READ, SEPRATOR='\t', file_type="csv", has_header=0, add_header=1, COLUMN_NAME_INPUT= ["SMILES"]):
        """Loads data from the specified file type"""
        if file_type == "csv":
            # Read the file twice: once with header and once without
            df_with_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, nrows=1)
            df_without_header = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR, header=None, nrows=1)

            if has_header==1:
                print("The file has a header row.")
                new_data = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR)

            else:
                print("The file does not have a header row.")
                if add_header==0:
                    new_data = pd.read_csv(FILE_PATH_READ, sep=SEPRATOR)
                    print("No header has been added.")
                else:
                    new_data = pd.read_csv(FILE_PATH_READ, names=COLUMN_NAME_INPUT, sep=SEPRATOR)
                    print("A header has been added.")
                
        elif file_type == "json":
            new_data = pd.read_json(FILE_PATH_READ)
        elif file_type == "excel":
            new_data = pd.read_excel(FILE_PATH_READ)
        else:
            raise ValueError("Unsupported file type. Supported types: csv, json, excel")
        return new_data
    
        
    # TODO: Save the output file
    def save_to_file(self, df, FILE_PATH_SAVE):
        """
        Saves the DataFrame to a csv file.
        """
        try:
            df.to_csv(FILE_PATH_SAVE, sep='\t', header=False, index=False)
            print(f"File saved successfully to {FILE_PATH_SAVE}")
        except Exception as e:
            print(f"An error occurred while saving the file: {e}")


    # TODO: Count the number of rings
    def number_of_ring(self, smiles):
        number_of_ring = 0
        for s in smiles:
            count_ring_indices = 0
            in_brackets = False
            two_digit_idx = False
            idx = ''
            for char in s:
                if char == '[':
                    in_brackets = True
                elif char == ']':
                    in_brackets = False
                elif char == '%':
                    two_digit_idx = True
                elif not in_brackets and char.isdigit():
                    if two_digit_idx:
                        if len(idx) == 0:
                            idx += char
                        elif len(idx) >= 1:
                            count_ring_indices += 1
                            idx = ''
                            two_digit_idx = False
                    else:
                        count_ring_indices += 1    
            ring_part = count_ring_indices / 2
            number_of_ring = number_of_ring + ring_part
        
        return number_of_ring


    def count_quaternary_centers_carbon(self, smiles):
        """
        Counts the number of quaternary carbons (carbon atoms bonded to four other carbons) in a molecule.
        
        Returns:
        - int: The number of quaternary carbons.
        """
        quaternary_count = 0
        mol = pybel.readstring("smi", smiles)
        obmol = mol.OBMol

        # Iterate through all the atoms in the OBMol object
        for atom in ob.OBMolAtomIter(obmol):
            # Check if the atom is a carbon (atomic number 6)
            if atom.GetAtomicNum() == 6:  # Carbon has atomic number 6
                # Count the number of bonded carbons
                bonded_carbon_count = 0

                # Iterate over all neighbors of the atom
                for neighbor in ob.OBAtomAtomIter(atom):
                    if neighbor.GetAtomicNum() == 6:  # Check if neighbor is also carbon
                        bonded_carbon_count += 1

                # If the carbon is bonded to exactly 4 other carbons, it's a quaternary carbon
                if bonded_carbon_count == 4:
                    quaternary_count += 1

        return quaternary_count
    
    # TODO: Calculate the FDV
    def divalent_nodes_fraction(self, smiles):
        mol=Chem.MolFromSmiles(smiles)
        atom_number = 0
        divalent_node = 0
        
        for atom in mol.GetAtoms():
            atom_number += 1
            degree = atom.GetDegree()

            if degree == 2:
                divalent_node += 1 
            else:
                continue

        divalent_ratio = round(divalent_node/atom_number,2)

        return divalent_ratio


    # TODO: Canonicalize molecules
    def canonicalize_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string provided")
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        return canonical_smiles
    

    # TODO: Check the aromaticity
    def is_aromatic(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        return any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    

    # TODO: To do a sub-structure search
    def substructure_search(self, df, smart):
        n = Chem.MolFromSmarts(smart)
        smiles_list = df["SMILES"].to_list()
        pattern = n
        true_list = []
        for idx,smiles in enumerate(smiles_list):
            m = Chem.MolFromSmiles(smiles)
            if m.HasSubstructMatch(pattern) == True:
                #print("Structure {}: pattern found {}".format(idx+1,m.HasSubstructMatch(pattern)))
                #print(idx)
                true_list.append(idx)

        df.loc[:, 'Substructures'] = "substructure"
        for i in true_list:
            df.loc[i, 'Substructures'] = True
        for j in range(len(df)):
            if j not in true_list:
                df.loc[j, 'Substructures'] = False 

        return df
    
    # TODO: To know how many repeating units of the sub-structures
    def substructure_units(self, df, smart):
        n = Chem.MolFromSmarts(smart)
        smiles_list = df["SMILES"].to_list()
        true_list = []
        substructure_appear_list = []
        for idx,smiles in enumerate(smiles_list):
            m = Chem.MolFromSmiles(smiles)
            # Get substructure matches
            matches = m.GetSubstructMatches(n)
            #print(matches)
            # Check if all atoms in the molecule belong to the substructure

            if matches:
                true_list.append(1)
                # Count the number of substructure matches   
                num_matches = len(matches)
                substructure_appear_list.append(num_matches)
            else:
                true_list.append(0)
                substructure_appear_list.append(0)


        df.loc[:, 'Substructures'] = "substructure"
        df.loc[:, 'Substructures appear times'] = "how many units"
        for i in range(len(df)):
            df.loc[i, 'Substructures'] = true_list[i]
            df.loc[i, 'Substructures appear times'] = substructure_appear_list[i]
        return df
    

    # TODO: Check the number of small rings
    def number_of_small_rings(self, smiles):
        number_of_small_ring = 0
        mol = pybel.readstring("smi", smiles)
        mol.OBMol.FindRingAtomsAndBonds()
        ring_data = mol.OBMol.GetSSSR()  # Get the smallest set of smallest rings (SSSR)
        
        for ring in pybel.ob.OBMolRingIter(mol.OBMol):
            if ring:
                # Check if the ring size is 3-4
                if ring.Size() == 3 or ring.Size() == 4:
                    number_of_small_ring += 1 
                else:
                    number_of_small_ring = number_of_small_ring

        return number_of_small_ring
    
