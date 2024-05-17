import pandas as pd
from typing import Dict
from Bio import SeqIO

def __csv_to_fasta(csv_file, output_fasta_file):

    """ Convert a CSV file to a FASTA file while excluding rows that are entirely NaN. """
    
    try:
        df = pd.read_csv(csv_file)  # Read the csv file
        with open(output_fasta_file, 'w') as fasta_file:
            for index, row in df.iterrows():
                if row.dropna().empty:
                    continue  # Skip rows that are all NaN
                # Convert each value to integer if it's a whole number, else keep original
                sequence = ''.join(str(int(val)) if isinstance(val, float) and val.is_integer() else str(val) for val in row)
                header = f'>{index + 1}'
                fasta_file.write(f'{header}\n{sequence}\n')
    except FileNotFoundError:
        print(f"Error: The file {csv_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __process_fasta(fasta_file, output_fasta):

    """ Processes a FASTA file to exclude rows that are entirely NaN. """

    try:
        # open input and output files
        with open(fasta_file, "r") as input_handle, open(output_fasta, "w") as output_handle:
            # parse the input FASTA file
            for line in SeqIO.parse(input_handle, "fasta"):
                # check if the sequence contains at least one position with 1
                if not 'nan' in line.seq:
                    # write the sequence to the output file
                    SeqIO.write(line, output_handle, "fasta")
    except FileNotFoundError:
        print(f"Error: The file {fasta_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __read_fasta(fasta_file):

    """ Function that creates a dictionary with headers as keys and sequences as values. """

    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences


def __filtering_sequences(sequences, coords_file):

    """ Extract rows from the coords_file corresponding to indexes of the sequences (keys of the 'sequences' dictionary). """
    info = pd.read_csv(coords_file)
    filtered_rows = info.loc[info.index.isin([int(key) for key in sequences.keys()])]
    return filtered_rows

"""
-------------------
POPULATION CREATION
-------------------
""" 

def __pop_coords(sequences, filtered_rows):

    """ Create 2 populations according to their coordinates (x,y). """

    # Create two empty populations
    population_1 = []
    population_2 = []

    extracted_values = filtered_rows.iloc[:, :2]
    # Calculate the midpoint of x and y coordinates
    x_values = extracted_values['x'].astype(float)
    y_values = extracted_values['y'].astype(float)
    x_mid = (x_values.min() + x_values.max()) / 2
    y_mid = (y_values.min() + y_values.max()) / 2

    # Iterate through extracted values and assign to populations based on coordinates
    for index, row in extracted_values.iterrows():
        x_coord = float(row['x'])
        y_coord = float(row['y'])
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequences:  # Check if sequence is not None
            if x_coord < x_mid and y_coord < y_mid:
                population_1.append(sequence)
            else:
                population_2.append(sequence)

    return population_1, population_2


def __pop_mutation_label(sequences, filtered_rows):

    """ Create 2 populations according to the 'mutation' label in the DataFrame. """
    population_1 = []
    population_2 = []

    # Extract values including the 'mutation' column, which is column 5 (index 4)
    mutation_values = filtered_rows['mutation']  # assuming 'mutation' is the name of the column

    # unique_mutations = mutation_values.unique()
    # if len(unique_mutations) != 2:
    #     raise ValueError("There must be exactly two unique mutation types. Found these mutations: {}".format(unique_mutations))

    # Iterate through the DataFrame rows
    for index, mutation in mutation_values.items():
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequence:  # Check if sequence is not None
            if mutation == 1.0:  # Assume 'Type1' as a mutation type for Population 1
                population_1.append(sequence)
            else:  # All other mutation types go to Population 2
                population_2.append(sequence)

    return population_1, population_2

"""
------
CHECKS
------
""" 

def __check(sequences: Dict[str, str]) -> None:

    """ Check for valid sequence proportions. """
    """ Number of sequences > 1 """
    """ Sequences with the same length """

    if len(sequences) < 2:
        raise ValueError("At least 2 sequences required!")

    sequences_list = list(sequences.values()) # list of sequences
    seqs_len = len(sequences_list[0]) # sequence's length 
    if not all(len(s) == seqs_len for s in sequences_list):
        raise ValueError("Sequences are required to have the same length!")
    
def __check_populations(population_1, population_2):
       
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