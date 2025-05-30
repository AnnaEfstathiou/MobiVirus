import pandas as pd
import libsequence
import argparse
from datetime import datetime
from collections import Counter
import re

#%% Data manipulation

def import_data(csv_file):
    """
    Load genome sequences from a CSV file and split into normal and super strain datasets.
    Args:
        csv_file (str): Path to the input CSV file containing genome sequences.
    Returns:
        tuple: (all_data, super_strains, normal_strains, genome_length, total_pop, super_pop, normal_pop)
    """

    all_seqs_data = pd.read_csv(csv_file, header=None)  # CSV file --> DataFrame
    all_seqs_data = all_seqs_data.dropna()              # Drop rows with NaN values
    all_seqs_data = all_seqs_data.astype(int)           # Turn values to integers
    l = len(all_seqs_data.columns)                      # Genome length

    ## Split the dataset ##
    ss_seqs_data = _retrieve_super_strains(all_seqs_data) # super strains
    ns_seqs_data = all_seqs_data.drop(ss_seqs_data.index)  # normal strains

    pop = len(all_seqs_data)   # Total number of all sequences
    pop_ss = len(ss_seqs_data) # Total number of super strains
    pop_ns = len(ns_seqs_data) # Total number of normal strains

    return all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns

def _retrieve_super_strains(df):  
    """
    Identify 'super strain' sequences with a SNP in the middle position(s).
    Args:
        df (DataFrame): Input DataFrame with binary genome sequences.
    Returns:
        DataFrame: Subset of input with only super strain sequences.
    """

    length = len(df.columns)
        
    if length % 2 == 1:
        # Odd length
        middle_column_index = length // 2
        middle_column_name = df.columns[middle_column_index]
        # Filter rows where the value in the middle column is 1
        ss_df = df[df[middle_column_name] == 1]
    else:
        # Even length
        middle_left = (length // 2) - 1
        middle_right = length // 2
        # Get the two middle column names
        middle_left_name = df.columns[middle_left]
        middle_right_name = df.columns[middle_right]
        # Filter rows where the value in EITHER middle column is 1
        ss_df = df[(df[middle_left_name] == 1) | (df[middle_right_name] == 1)]

    return ss_df

def event_number(csv_file):
    """
    Extract the numeric event identifier from the filename.
    Args:
        csv_file (str): Filename of the genome CSV.
    Returns:
        str: Extracted number from the filename (e.g., '100' from 'genomes_100.csv').
    Raises:
        ValueError: If no number is found in the filename.
    """

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

#%% Sampling

def _sampling(data, sample, pop):
    """
    Sample a subset of sequences from the population.
    Args:
        data (DataFrame): Input sequences.
        sample (int): Number of sequences to sample.
        pop (int): Total number of sequences available.
    Returns:
        DataFrame: Sampled subset.
    Raises:
        ValueError: If the input data is empty or sampling returns fewer than two entries.
    """

    if not data.empty:
        if sample <= pop:
            sampled_data = data.sample(n=sample, random_state=int(datetime.now().timestamp())) 
        else:
            sampled_data = data.sample(n=pop, random_state=int(datetime.now().timestamp())) 
            print(f"Sample size ({sample}) is larger than the population size ({pop}), so the entire population was used as the sample.")
        
        # Raise an error if the sampled data contains only one entry
        if len(sampled_data) == 1:
            raise ValueError("Sampled data contains only one entry, which is insufficient.")
    else:
        sampled_data = None
        raise ValueError("No data to sample from.") # Raise an error if sampled_data is None
    
    return sampled_data

#%% Simulate Data

def simulated_data(sampled_data, l):
    """
    Convert sampled sequence data into a format compatible with libsequence.
    Args:
        sampled_data (DataFrame): Sampled sequence data.
        l (int): Genome length.
    Returns:
        SimData: libsequence-compatible simulation object.
    """

    ## Polymorphic sites ##
    if sampled_data is not None:
        polym_sites = sampled_data.loc[:, (sampled_data != 0).any(axis=0)]
        pos = [(col / l) for col in list(polym_sites.columns)]
        seqs = [''.join(map(str, row)) for row in polym_sites.values]
    else:
        pos = None
        seqs = None
    ## Simulate Data for further analysis ##   
    if pos is not None and seqs is not None:
        sd = libsequence.SimData()
        sd.assign(pos, seqs)
    else:
        sd = None
    
    return sd

#%% Statistics Functions

def _summary_stats(simulated_data):
    """
    Compute key population genetics statistics using libsequence.
    Args:
        simulated_data (libsequence.SimData): A simulation object containing sequences and SNP positions.
    Returns:
        tuple: 
            - tajimasd (float): Tajima's D statistic, indicating neutrality deviations.
            - pi (float): Nucleotide diversity (Pi), measuring average pairwise differences.
            - thetaw (float): Watterson's Theta, estimating mutation rate based on segregating sites.
    """

    ps = libsequence.PolySIM(simulated_data) 
    tajimasd = ps.tajimasd()
    pi = ps.thetapi()
    thetaw = ps.thetaw()
    
    return tajimasd, pi, thetaw

def _haplotypes(sampled_data):
    """
    Calculate haplotype count and diversity from sampled sequence data.
    Args:
        sampled_data (DataFrame): Sampled binary sequence data.
    Returns:
        tuple: (number of unique haplotypes, haplotype diversity score)
    """

    # 1. DataFrame into dict (index:sequence)
    sample_dict = {index: ''.join(map(str, row)) for index, row in sampled_data.iterrows()}
    # 2. Get the total number of unique haplotypes
    unique_haplotypes = len(set(sample_dict.values()))
    # 3. Calculate haplotype diversity
    total_haplotypes = len(sample_dict)                                                                      # Total number of haplotypes (sequences)
    haplotype_frequencies = Counter(sample_dict.values())                                                    # Count the frequency of each unique haplotype
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # Compute haplotype diversity using the formula
    
    return unique_haplotypes, haplotype_diversity

def statistics(sequences, strain_type, sample):
    """
    Compute summary statistics from genome sequence data.
    Args:
        sequences (str): Path to genome CSV file.
        strain_type (str): Which strains to analyze ('mix_strains', 'ss_strains', 'ns_strains').
        sample (int): Sample size.
    Returns:
        tuple: (Tajima's D, Pi, Theta W, Number of unique haplotypes, Haplotype diversity)
    """

    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling
    sampled_data = None  # Initialize sampled_data to avoid unbound error
    if strain_type   == "mix_strains":
        sampled_data = _sampling(all_seqs_data, sample, pop)
    elif strain_type == "ss_strains":
        sampled_data = _sampling(ss_seqs_data, sample, pop_ss)
    elif strain_type == "ns_strains":
        sampled_data = _sampling(ns_seqs_data, sample, pop_ns)

    if sampled_data is None:
        raise ValueError(f"Invalid strain_type '{strain_type}' or sampling failed.") # Handle the case where strain_type is invalid or sampling fails

    # 3) Simulating sampled data
    sim_data = simulated_data(sampled_data, l)

    if sim_data is not None:

        # 4) Tajima's D, Pi, θw
        tajimasd, pi, thetaw = _summary_stats(sim_data)

        # 5) Number of unique seqs, haplotype diversity
        uniq_haplo, haplo_div = _haplotypes(sampled_data)

    return tajimasd, pi, thetaw, uniq_haplo, haplo_div 

#%% 
if __name__ == "__main__":

    '''Define parsers '''

    ## General Parsers
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics (Tajimas D, Pi-estimator, Theta Watterson, Number of unique haplotypes, Haplotype Diversity)")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-s', '--sample_size', required=True, type=int, help='Sample size.')
    parser.add_argument('-p', '--population_size', action="store_true", help='Population size.')

    ## Mutually exclusive group for strain statistics
    strain_group = parser.add_mutually_exclusive_group()
    strain_group.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    strain_group.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    
    args = parser.parse_args()
    
    ''' Results '''

    # Return the population size 
    if args.population_size:
        all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(args.genome_file)
        print(f"Number of sequences (normal & super strains): {pop}")
        print(f"Number of super sequences: {pop_ss}")
        print(f"Number of normal strains: {pop_ns}")  
        print(f"Genome length: {l}\n")

    # Determine strain_type based on command-line arguments
    # Default strain type: both normal and super strains
    strain_type = "mix_strains"
    if args.ss_strains:
        strain_type = "ss_strains"
    elif args.ns_strains:
        strain_type = "ns_strains"
    
    try:
        genome_file = args.genome_file
        sample_size = args.sample_size
        suffix = event_number(genome_file) # Determine suffix (e.g. File:genomes_100.csv -> suffix:100)
    except UnboundLocalError:
        suffix = None 
        genome_file = None
        sample_size = None

    ## Summary Statistics ##
    tajimasd, pi, thetaw, uniq_haplo, haplo_div = statistics(args.genome_file, strain_type, args.sample_size)

    # DataFrame with the results
    sumstats = {
        'Tajimas D': [tajimasd],
        'Pi-estimator': [pi],
        'Theta Watterson': [thetaw],
        'Number of unique haplotypes': [uniq_haplo],
        'Haplotype Diversity': [haplo_div]}
    sumstats_df = pd.DataFrame(sumstats)

    # Save the DataFrame to a CSV file
    sumstats_df.to_csv(f"sumstats_{suffix}.csv", index=False)
    print(f"The results are stored in the file: sumstats_{suffix}.csv")