""" 
INPUT: FASTA or CSV file with binary sequences
Calculating: 
- Fst
"""

import pandas as pd
from typing import Dict
import numpy as np
from Bio import SeqIO
from calc_stats import __csv_to_fasta, __process_fasta, __read_fasta
from configparser import ConfigParser
import os
import sys



"""
------------------
INTERNAL FUNCTIONS
------------------
""" 


def __filtering_sequences(sequences, coords_file):

    """ Extract rows from the coords_file corresponding to indexes of the sequences (keys of the 'sequences' dictionary). """
    info = pd.read_csv(coords_file)
    filtered_rows = info.loc[info.index.isin([int(key) for key in sequences.keys()])]
    return filtered_rows

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

    unique_mutations = mutation_values.unique()
    if len(unique_mutations) != 2:
        raise ValueError("There must be exactly two unique mutation types. Found these mutations: {}".format(unique_mutations))

    # Iterate through the DataFrame rows
    for index, mutation in mutation_values.items():
        sequence = sequences.get(str(index), None)  # Convert index to string to match sequence keys
        if sequence:  # Check if sequence is not None
            if mutation == 1.0:  # Assume 'Type1' as a mutation type for Population 1
                population_1.append(sequence)
            else:  # All other mutation types go to Population 2
                population_2.append(sequence)

    return population_1, population_2



file = '/home/anna/mobivirus/simulation/results_12_05_2024_12_14/genomes/genomes_210.csv'
new_fasta = 'lala.csv'
__csv_to_fasta(file, new_fasta)
sequences = __read_fasta(new_fasta)
# print(sequences)
# Extract rows corresponding to dictionary keys
filtered_rows = __filtering_sequences(sequences, '/home/anna/mobivirus/simulation/results_12_05_2024_12_14/samples/coords_210.csv')
# print(filtered_rows)
pop_coords_1, pop_coords_2 = __pop_coords(sequences, filtered_rows)
pop_label_1, pop_label_2 = __pop_mutation_label(sequences, filtered_rows)

# Print the populations
# print("COORDS")
# print("Population 1:",pop_coords_1)
# print("Population 2:", pop_coords_2)
# print("LABEL")
# print("Population 1:",pop_label_1)
# print("Population 2:", pop_label_2)




def allele_frequencies(population):
    """Calculate allele frequencies for a population."""
    allele_counts = np.sum(population, axis=0)
    total_alleles = 2 * len(population)
    return allele_counts / total_alleles

def expected_heterozygosity(allele_freq):
    """Calculate expected heterozygosity for a population."""
    return 1 - np.sum(allele_freq ** 2)

def observed_heterozygosity(population):
    """Calculate observed heterozygosity for a population."""
    heterozygous_individuals = np.any(population != population[0], axis=1)
    return np.mean(heterozygous_individuals)

# Convert populations to numpy arrays
population_1_array = np.array([list(map(int, seq)) for seq in pop_label_1])
population_2_array = np.array([list(map(int, seq)) for seq in pop_label_2])

# Calculate allele frequencies for each population
allele_freq_1 = allele_frequencies(population_1_array)
allele_freq_2 = allele_frequencies(population_2_array)

# Calculate expected heterozygosity for each population
He_1 = expected_heterozygosity(allele_freq_1)
He_2 = expected_heterozygosity(allele_freq_2)

# Calculate observed heterozygosity for each population
Ho_1 = observed_heterozygosity(population_1_array)
Ho_2 = observed_heterozygosity(population_2_array)

# Calculate Fst
Fst = (He_1 + He_2 - (Ho_1 + Ho_2)) / (He_1 + He_2)

print("Fst:", Fst)
