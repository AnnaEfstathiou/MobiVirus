import os
import pandas as pd
from typing import Dict

def process_csv(csv_file_path):

    """ Function that creates a dictionary (Dict[str, str]) with index as keys and sequences as values. """

    genomes_data = pd.read_csv(csv_file_path, header=None) # Read the CSV file into a DataFrame
    genomes_data = genomes_data.dropna()                   # Drop rows with NaN values
    genomes_data = genomes_data.astype(int)                # Convert all elements in the DataFrame to integers
    sequences = {str(idx): ''.join(map(str, row)) for idx, row in genomes_data.iterrows()} # Format > sequences: Dict[str, str]
    
    return sequences

def __filtering_sequences(sequences, coords_file):

    """ Extract rows from the coords_file corresponding to indexes of the sequences (keys of the 'sequences' dictionary). """
    
    info = pd.read_csv(coords_file)
    filtered_coords = info.loc[info.index.isin([int(key) for key in sequences.keys()])]
    return filtered_coords

"""
-------------------
POPULATION CREATION
-------------------
""" 

def __pop_coords(sequences, filtered_coords):

    """ Create 2 populations according to their coordinates (x,y). """

    # Initialize lists to store the populations
    population_1 = []
    population_2 = []

    extracted_values = filtered_coords.iloc[:, :2] # Extract rows from the coords_file corresponding to indexes of the sequences
    x_values = extracted_values['x'].astype(float) # Calculate the midpoint of x coordinates
    x_mid = (x_values.min() + x_values.max()) / 2  
    y_values = extracted_values['y'].astype(float) # Calculate the midpoint of y coordinates
    y_mid = (y_values.min() + y_values.max()) / 2

    # Iterate through extracted values and assign to populations based on coordinates
    for index, row in extracted_values.iterrows():
        x_coord = float(row['x'])
        y_coord = float(row['y'])
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequences:                               # Check if sequence is not None
            if x_coord < x_mid and y_coord < y_mid:
                population_1.append(sequence)
            else:
                population_2.append(sequence)

    return population_1, population_2

def __pop_infection_label(sequences, filtered_coords): 

    """ Create 2 populations according to the 'infection' label in the DataFrame. """
    """ Infection label=1 : no recombination | Infection label=2 : recombination  """

    # Initialize lists to store the populations
    population_1 = []
    population_2 = []

    label_values = filtered_coords['Infection label'] # Extract values including the 'label' column, which is column 2

    # Iterate through the DataFrame rows
    for index, label in label_values.items():
        sequence = sequences.get(str(index), None)    # Convert index to string to match sequence keys
        if sequence:                                  # Check if sequence is not None
            if label == 1.0:                          # Assume 'Type1' as a label type for Population 1 (aka no recombination)
                population_1.append(sequence)
            else:                                     # All other label types go to Population 2 (aka recombination)
                population_2.append(sequence)

    return population_1, population_2

def __pop_mutation_label(sequences, filtered_coords):

    """ Create 2 populations according to the 'mutation' label in the DataFrame. """
    """ Mutation label=1 : normal strain | Mutation label=2 : super strain       """
    
    # Initialize lists to store the populations
    population_1 = []
    population_2 = []

    mutation_values = filtered_coords['Mutation']   # Extract values including the 'mutation' column, which is column 5 

    # Iterate through the DataFrame rows
    for index, mutation in mutation_values.items():
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequence:                                # Check if sequence is not None
            if mutation == 1.0:                     # Assume 'Type1' as a mutation type for Population 1 (aka normal strain)
                population_1.append(sequence)
            else:                                   # All other mutation types go to Population 2 (aka super strain)
                population_2.append(sequence)

    return population_1, population_2

"""
------
CHECKS
------
""" 

def __check(sequences: Dict[str, str]) -> None:

    """ Check for valid sequence proportions. """
    """ Number of sequences > 1               """
    """ Sequences with the same length        """

    if len(sequences) < 2:
        raise ValueError("At least 2 sequences required!")

    sequences_list = list(sequences.values()) # list of sequences
    seqs_len = len(sequences_list[0])         # sequence's length 
    if not all(len(s) == seqs_len for s in sequences_list):
        raise ValueError("Sequences are required to have the same length!")
    
def __check_populations(population_1, population_2):
    
    """ Check for valid population types. """
    """ Non empty populations             """
    """ Equal genome lengths              """

    # Check if either population is empty
    if not population_1 or not population_2:
        raise ValueError(f"Both populations must contain at least one sequence. Population 1 size: {len(population_1)}, Population 2 size: {len(population_2)}")

    # Check if sequences in both populations have the same length
    pop1_length = len(population_1[0]) if population_1 else 0
    pop2_length = len(population_2[0]) if population_2 else 0

    for seq in population_1:
        if len(seq) != pop1_length:
            raise ValueError("All sequences in population 1 must have the same length. Found sequence of length: {}. Expected length: {}.".format(len(seq), pop1_length))

    for seq in population_2:
        if len(seq) != pop2_length:
            raise ValueError("All sequences in population 2 must have the same length. Found sequence of length: {}. Expected length: {}.".format(len(seq), pop2_length))

    if pop1_length != pop2_length:
        raise ValueError("Sequences in both populations must have the same length. Population 1 length: {}, Population 2 length: {}".format(pop1_length, pop2_length))
    
def validate_files(genomes_file, coords_file):

    """ Ensure that inpur files come from the same simulation and generation. """

    # Normalize paths
    genomes_path = os.path.normpath(genomes_file)
    coords_path = os.path.normpath(coords_file)
    
    # Split paths
    genomes_parts = genomes_path.split(os.sep)
    coords_parts = coords_path.split(os.sep)

    # Find indices for /genomes/ and /samples/
    genomes_index = genomes_parts.index('genomes')
    coords_index = coords_parts.index('samples')

    # Ensure that the directories before /genomes/ and /samples/ are the same
    if genomes_parts[:genomes_index] != coords_parts[:coords_index]:
        raise ValueError("The part of the path before /genomes/ and /samples/ must be the same.")

    genomes_filename = genomes_parts[-1]
    coords_filename = coords_parts[-1]

    # Extract the date and number part from the filenames
    genomes_date_num = genomes_filename.split('_')[1:3]
    coords_date_num = coords_filename.split('_')[1:3]

    # Compare the extracted parts
    if genomes_date_num != coords_date_num:
        raise ValueError("The input files must come from the same simulation_date file and have the same number.")