# import argparse
import os
import pandas as pd
from typing import Dict
from Bio import SeqIO

import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


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


def __read_fasta(file_path: str) -> Dict[str, str]:

    """ Function that creates a dictionary with headers as keys and sequences as values. """
    
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
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


# def __pop_coords(sequences, filtered_rows, n_clusters=2):
#     """
#     Create populations according to their coordinates (x,y) using KMeans clustering.
    
#     Parameters:
#     sequences (dict): A dictionary containing sequences keyed by index.
#     filtered_rows (pd.DataFrame): A DataFrame containing 'x' and 'y' coordinate columns.
#     n_clusters (int): Number of clusters for KMeans. Default is 2.
    
#     Returns:
#     tuple: Two lists representing the populations split by clusters.
#     """

#     # Create two empty populations
#     populations = [[] for _ in range(n_clusters)]

#     extracted_values = filtered_rows.iloc[:, :2].astype(float)
#     coordinates = extracted_values[['x', 'y']].values

#     # Perform KMeans clustering
#     kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(coordinates)
#     labels = kmeans.labels_

#     # Iterate through extracted values and assign to populations based on clusters
#     for index, label in zip(extracted_values.index, labels):
#         sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
#         if sequence is not None:
#             populations[label].append(sequence)
            
#     # Plotting the clusters
#     plt.figure(figsize=(10, 8))
#     plt.scatter(coordinates[:, 0], coordinates[:, 1], c=labels, cmap='viridis', marker='o')
#     centers = kmeans.cluster_centers_
#     plt.scatter(centers[:, 0], centers[:, 1], c='red', s=300, alpha=0.6, marker='x')
#     plt.title('KMeans Clustering of Coordinates')
#     plt.xlabel('X Coordinate')
#     plt.ylabel('Y Coordinate')
#     plt.show()

#     return populations



def __pop_infection_label(sequences, filtered_rows): # split the population according to recombination

    """ Create 2 populations according to the 'infection' label in the DataFrame. """
    population_1 = []
    population_2 = []

    # Extract values including the 'label' column, which is column 2
    label_values = filtered_rows['Infection label']  # assuming 'mutation' is the name of the column

    # Iterate through the DataFrame rows
    for index, label in label_values.items():
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequence:  # Check if sequence is not None
            if label == 1.0:  # Assume 'Type1' as a label type for Population 1
                population_1.append(sequence)
            else:  # All other label types go to Population 2
                population_2.append(sequence)

    return population_1, population_2

def __pop_mutation_label(sequences, filtered_rows):

    """ Create 2 populations according to the 'mutation' label in the DataFrame. """
    population_1 = []
    population_2 = []

    # Extract values including the 'mutation' column, which is column 5 
    mutation_values = filtered_rows['Mutation']  # assuming 'mutation' is the name of the column

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
    

def validate_files(genomes_file, coords_file):
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
        raise ValueError("The input files must come from the same results_date file and have the same number.")