import pandas as pd
import libsequence
import argparse
from datetime import datetime
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

def _filtering_sequences(sequences, coords_file):
    """
    Extract rows from the coordinates file corresponding to the indexes of the provided sequences DataFrame.
    Args:
        sequences (DataFrame): Genome sequences.
        coords_file (str): Path to the coordinates CSV file.
    Returns:
        DataFrame: Filtered coordinates matching the sequence indexes.
    """

    info = pd.read_csv(coords_file)
    filtered_coords = info.iloc[info.index.intersection(sequences.index)]
    return filtered_coords

def _validate_files(genomes_file, coords_file):
    """
    Ensure that the genome and coordinate files come from the same simulation and generation.
    Args:
        genomes_file (str): Path to the genome CSV file.
        coords_file (str): Path to the coordinates CSV file.
    Raises:
        ValueError: If the directories or filenames do not match as expected.
    """

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

def _simple_random_sampling(sampled_data):
    """
    Perform simple random sampling by separating sampled data into super and normal strains,
    then combining them into one dataset.
    Args:
        sampled_data (DataFrame): Randomly sampled genome sequences.
    Returns:
        tuple: (combined DataFrame, count of super strains, count of normal strains)
    Raises:
        ValueError: If any population is empty or has only one entry.
    """

    sampled_data_1 = _retrieve_super_strains(sampled_data)  # super strains
    sampled_data_2 = sampled_data.drop(sampled_data_1.index) # normal strains

    # Validation
    # Raise an error if one of the sampled data contains only one entry or no entries
    if len(sampled_data_1) == 1 or len(sampled_data_2) == 1:
        raise ValueError("One of the populations contain only one entry, which is insufficient.")
    if len(sampled_data_1) == 0 or len(sampled_data_2) == 0:
        raise ValueError("One of the populations contain no entries, which is insufficient.")

    simple_random_sample = pd.concat([sampled_data_1, sampled_data_2], axis=0).reset_index(drop=True) # Stack the 2 samples

    return simple_random_sample, len(sampled_data_1), len(sampled_data_2)

def _stratified_sampling(sampled_data_1, sampled_data_2):
    """
    Combine two sampled populations into a single DataFrame for stratified sampling.
    Args:
        sampled_data_1 (DataFrame): First sampled group (e.g., super strains).
        sampled_data_2 (DataFrame): Second sampled group (e.g., normal strains).
    Returns:
        tuple: (combined DataFrame, size of first group, size of second group)
    Raises:
        ValueError: If either group is empty or has only one entry.
    """

    # Validation
    # Raise an error if one of the sampled data contains only one entry or no entries
    if len(sampled_data_1) == 1 or len(sampled_data_2) == 1:
        raise ValueError("One of the populations contain only one entry, which is insufficient.")
    if len(sampled_data_1) == 0 or len(sampled_data_2) == 0:
        raise ValueError("One of the populations contain no entries, which is insufficient.")

    stratified_sample = pd.concat([sampled_data_1, sampled_data_2], axis=0).reset_index(drop=True) # Stack the 2 samples

    return stratified_sample, len(sampled_data_1), len(sampled_data_2)

def _space_sampling(sampled_data, filtered_coords):
    """
    Divide sampled data into two groups based on the x-axis coordinate value relative to the midpoint.
    Args:
        sampled_data (DataFrame): Genome sequences sampled from the population.
        filtered_coords (DataFrame): Coordinate information for the sampled data.
    Returns:
        tuple: (combined DataFrame, size of left group, size of right group)
    """

    common_indices = sampled_data.index.intersection(filtered_coords.index) # common indexes between sampled data and coordinates file
    # filtered_sampled_data = sampled_data.loc[common_indices]

    filtered_coords = filtered_coords.loc[common_indices] # filter the filtered_coords to only include entries with these common indices
    extracted_values = filtered_coords['x'] # Extract x coordinates from the coords file corresponding to indexes of the sequences
    
    x_values = extracted_values.astype(float) 
    x_mid = (x_values.min() + x_values.max()) / 2 # Calculate the midpoint of x coordinates

    # Initialize lists to store the populations
    population_1 = []  # For x < x_mid
    population_2 = []  # For x >= x_mid

    # Iterate through the filtered_coords DataFrame to assign sequences to populations
    for index, row in filtered_coords.iterrows():
        x_coord = row['x']
        sequence = sampled_data.loc[index]

        if x_coord < x_mid:
            population_1.append(sequence)  # Assign to population 1 if x < x_mid
        else:
            population_2.append(sequence)  # Assign to population 2 if x >= x_mid

    # Convert the populations into DataFrames
    population_1 = pd.DataFrame(population_1)
    population_2 = pd.DataFrame(population_2)

    coordinates_sample = pd.concat([population_1, population_2], axis=0).reset_index(drop=True) # Stack the 2 samples

    # Returning the two DataFrames and the population counts
    return coordinates_sample, len(population_1), len(population_2)

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

def Fst(sequences, sampling_type, strain_type, sample, coords_file = None):
    """
    Calculate the fixation index (Fst) between two populations using various sampling strategies.
    Args:
        sequences (str): Path to the input genome CSV file.
        sampling_type (str): Sampling method to use. Options: 'simple_random_sampling', 
                             'stratified_sampling', or 'coordinates_sampling'.
        strain_type (str): Type of strain to sample. Options: 'mix_strains', 'ss_strains', 'ns_strains'.
        sample (int): Number of sequences to sample.
        coords_file (str, optional): Path to the coordinates file. Required if using 'coordinates_sampling'.
    Returns:
        Fst: libsequence.Fst object containing calculated Fst statistics.
    Raises:
        ValueError: If strain_type is invalid or sampling fails.
    """

    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling  
    edited_samples = None 
    if sampling_type == "simple_random_sampling":
        sampled_data = _sampling(all_seqs_data, sample, pop)
        edited_samples, n_pop1, n_pop2 = _simple_random_sampling(sampled_data)

    elif sampling_type == "stratified_sampling":
        sampled_data_ss = _sampling(ss_seqs_data, sample, pop_ss)
        sampled_data_ns = _sampling(ns_seqs_data, sample, pop_ns)
        edited_samples, n_pop1, n_pop2 = _stratified_sampling(sampled_data_ss, sampled_data_ns)

    elif sampling_type == "coordinates_sampling": 
        if strain_type == "mix_strains":
            filtered_coords = _filtering_sequences(all_seqs_data, coords_file) # get x axis coordinates
            sampled_data = _sampling(all_seqs_data, sample, pop)
        elif strain_type == "ss_strains":
            filtered_coords = _filtering_sequences(ss_seqs_data, coords_file) # get x axis coordinates
            sampled_data = _sampling(ss_seqs_data, sample, pop_ss)
        elif strain_type == "ns_strains":
            filtered_coords = _filtering_sequences(ns_seqs_data, coords_file) # get x axis coordinates
            sampled_data = _sampling(ns_seqs_data, sample, pop_ns)
        edited_samples, n_pop1, n_pop2 = _space_sampling(sampled_data, filtered_coords)
    if edited_samples is None:
        raise ValueError(f"Invalid strain_type '{strain_type}' or sampling failed.") # Handle the case where strain_type is invalid or sampling fails
    
    # 3) Simulating sampled data
    sim_data = simulated_data(edited_samples, l)

    if sim_data is not None:

        # 4) Fst
        fst = libsequence.Fst(sim_data,[n_pop1,n_pop2])

    return fst  

#%% 
if __name__ == "__main__":

    '''Define parsers '''

    ## General Parsers
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-s', '--sample_size', required=True, type=int, help='Sample size.')
    parser.add_argument('-p', '--population_size', action="store_true", help='Population size.')

    # Parsers for Fst
    parser.add_argument('-sample_type', '--sampling_technique', type=str, required=True, help='Define the way the sampling of 2 populations will happen. 3 ways: "rdm","str","coords"')
    parser.add_argument('-c', '--coords_file', type=str, help='The path to the coordinates file. Required if -sampling_techning == "coords".')

    ## Mutually exclusive group for strain statistics
    strain_group = parser.add_mutually_exclusive_group()
    strain_group.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    strain_group.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    
    args = parser.parse_args()

    # Check if the file containing the coordinates is required
    if args.sampling_technique == "coords" and not args.coords_file:
        parser.error('--coords_file is required when --sampling_technique is set to "coords"')

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

    ## Summary Statistics & Fst ##
    if args.sampling_technique == "rdm":
        sampling_type = "simple_random_sampling"
        fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = None)
    elif args.sampling_technique == "str":
        sampling_type = "stratified_sampling"
        fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = None)
    elif args.sampling_technique == "coords":
        sampling_type = "coordinates_sampling"
        _validate_files(args.genome_file, args.coords_file)
        fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = args.coords_file)
        
    hsm_fst = fst.hsm()         # Hudson, Slatkin & Maddison Fst
    slatkin_fst = fst.slatkin() # Slatkin Fst
    hbk_fst = fst.hbk()         # Hudson, Boos & Kaplan Fst

    # DataFrame with the results
    sumstats = {
        'HSM Fst': [hsm_fst],
        'Slatkin Fst': [slatkin_fst],
        'HBK Fst': [hbk_fst]}
    sumstats_df = pd.DataFrame(sumstats)

    # Save the DataFrame to a CSV file
    sumstats_df.to_csv(f"fst_{suffix}.csv", index=False)
    print(f"The results are stored in the file: fst_{suffix}.csv")