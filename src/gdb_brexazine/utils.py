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


    def count_quaternary_centers_carbon_openbabel(self, smiles):
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
    

    def count_quaternary_centers_carbon_rdkit(self, smiles):
        """
        Counts the number of quaternary carbons (carbon/general atoms bonded to four other carbons) in a molecule.
        
        Returns:
        - int: The number of quaternary carbons.
        """
        non_carbon_number = 0
        carbon_number = 0
        mol = Chem.MolFromSmiles(smiles)

        for x in mol.GetAtoms():
        
            if x.GetHybridization() == Chem.HybridizationType.SP3:
                i = x.GetIdx()
                neighbor_list_atom = []

                neighbor_list_RDKMol = mol.GetAtomWithIdx(i).GetNeighbors()

                #print(x.GetSymbol())
                #print(neighbor_list_RDKMol)

                if x.GetSymbol() == 'C' and x.GetTotalNumHs() == 0:

                    for n in neighbor_list_RDKMol:

                        neighbor_atom = n.GetSymbol()
                        neighbor_list_atom.append(neighbor_atom)
                    
                    #print(neighbor_list_atom)

                    if all([x == 'C' for x in neighbor_list_atom]) == True:
                        carbon_number += 1
                        #print('carbon:', carbon_number)

                    if all([x == 'C' for x in neighbor_list_atom]) == False:
                        non_carbon_number += 1
                        #print('general:',non_carbon_number)
            
            quaternary_count_carbon = carbon_number
            quaternary_count_general = carbon_number + non_carbon_number

        return quaternary_count_carbon, quaternary_count_general
    


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
    
    
    # TODO: Check if all atoms belong to the sub-structures
    def all_atoms_in_substructure(self, smart, smiles):
        n = Chem.MolFromSmarts(smart)
        m = Chem.MolFromSmiles(smiles)
        # Get substructure matches
        matches = m.GetSubstructMatches(n)
          # Check if all atoms in the molecule belong to the substructure
        if matches:
            # Find all atom indices involved in the matched substructure
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            # Check if all atoms in the molecule are part of the matched substructure
            if len(matched_atoms) == m.GetNumAtoms():
                return True
            else:
                #print("Not all atoms in the molecule belong to the substructure.")
                return False
        else:
            return False
    

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
    

    # TODO: Count the divalent and trivalent atoms
    def count_divalent_trivalent_atoms(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        divalent_count = 0
        trivalent_count = 0
        
        for atom in molecule.GetAtoms():
            # Count the number of bonds (valence) for each atom
            valence = len(atom.GetNeighbors())
            
            if valence == 2:
                divalent_count += 1
            elif valence == 3:
                trivalent_count += 1
        
        return [divalent_count, trivalent_count]
    

    # TODO: check a molecule contains how many 5/6/7/7+-membered rings
    def count_5to7_membered_rings(self, smiles):
        number_5 = 0
        number_6 = 0
        number_7 = 0
        number_large = 0


        mol = pybel.readstring("smi", smiles)
        mol.OBMol.FindRingAtomsAndBonds()
        ring_data = mol.OBMol.GetSSSR()  # Get the smallest set of smallest rings (SSSR)
        
        for ring in pybel.ob.OBMolRingIter(mol.OBMol):
            # Check the number of 5-membered-ring
            if ring.Size() == 5:
                number_5 += 1 
            
            elif ring.Size() == 6:
                number_6 += 1 
            elif ring.Size() == 7:
                number_7 += 1 

            elif ring.Size() > 7:
                number_large += 1


        return (number_5, number_6, number_7, number_large)
        

    # TODO: display the ring sizes of the rings in a molecules
    def show_ring_size(self, smiles):
        mol = pybel.readstring("smi", smiles)

        # Perform ring perception
        mol.OBMol.AddHydrogens()
        mol.OBMol.PerceiveBondOrders()

        # Get the number of rings and their sizes
        ring_list = []
        for ring in mol.OBMol.GetSSSR():
            print(ring)
            ring_list.append(ring)

        # Print the ring sizes
        print("Ring list:", ring_list)

        atom_ring_count = {}
        # Iterate over the rings
        for ring in ring_list:
            # Iterate over the atoms in the ring
            for atom_idx in ring._path:
                print(atom_idx, ring._path)
                if atom_idx in atom_ring_count:
                    
                    atom_ring_count[atom_idx] += 1
                else:
                    atom_ring_count[atom_idx] = 1
        atom_ring_count
        print(atom_ring_count)


    # TODO: HAC calculation
    def hac(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        size = mol.GetNumHeavyAtoms()

        return size


    # TODO: Make a counting table
    def count_table(self,length1, length2):
        # Create a dictionary with column names as keys and values for one row
        data = {
            'Size1': [length1],
            'Size2': [length2],

        }
        # Create a DataFrame from the dictionary
        df_count = pd.DataFrame(data)
        df_count

        return df_count


    # TODO: Count the number of acyclic bonds
    def count_acyclic_bonds(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        # Count the number of acyclic bonds
        num_acyclic_bonds = sum(1 for bond in mol.GetBonds() if not bond.IsInRing())
        return num_acyclic_bonds
