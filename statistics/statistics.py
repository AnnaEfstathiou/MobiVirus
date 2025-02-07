import pandas as pd
import libsequence
import argparse
from datetime import datetime
from collections import Counter
import re
import os

## mix is refered to the whole pool of sequences aka super and normal strains ##
## ss is refered to super strains                                             ##
## ns is refered to normal strains                                            ##

#%% Data manipulation

def import_data(csv_file):
    ## Convert the CSV file (genomes) to a DataFrame (rows: sequence / columns: positions) ##
    ## Retrieve normal and super spreaders, store them into seperate DataFrames            ##

    all_seqs_data = pd.read_csv(csv_file, header=None)  # CSV file --> DataFrame
    all_seqs_data = all_seqs_data.dropna()              # Drop rows with NaN values
    all_seqs_data = all_seqs_data.astype(int)           # Turn values to integers
    l = len(all_seqs_data.columns)                      # Genome length

    ## Split the dataset ##
    ss_seqs_data = __retrieve_super_strains(all_seqs_data) # super strains
    ns_seqs_data = all_seqs_data.drop(ss_seqs_data.index)  # normal strains

    pop = len(all_seqs_data)   # Total number of all sequences
    pop_ss = len(ss_seqs_data) # Total number of super strains
    pop_ns = len(ns_seqs_data) # Total number of normal strains

    return all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns

def __retrieve_super_strains(df):  
    ## Create a dataframe with the sequences that have a SNP in the middle position (super strains) ##

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
    ## Isolate suffix number from the csv file ##
    ## Suffix represents the event number      ##

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

def __filtering_sequences(sequences, coords_file):
    ## Extract rows from the coords_file corresponding to indexes of the sequences DataFrame. ##
    
    info = pd.read_csv(coords_file)
    filtered_coords = info.iloc[info.index.intersection(sequences.index)]
    return filtered_coords

def __validate_files(genomes_file, coords_file):
    ## Ensure that inpur files come from the same simulation and generation. ##

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

def __sampling(data, sample, pop):
    ## Sampling n sequences from a pool                                             ##
    ## If sample size is bigger than the population size, take the whole population ##

    if not data.empty:
        if sample <= pop:
            sampled_data = data.sample(n=sample, random_state=int(datetime.now().timestamp())) 
        else:
            sampled_data = data.sample(n=pop, random_state=int(datetime.now().timestamp())) 
        
        # Raise an error if the sampled data contains only one entry
        if len(sampled_data) == 1:
            raise ValueError("Sampled data contains only one entry, which is insufficient.")
    else:
        sampled_data = None
        raise ValueError("No data to sample from.") # Raise an error if sampled_data is None
    
    return sampled_data

def __simple_random_sampling(sampled_data):
    ## 1. Split sampled sequences into normal and super strains       ##
    ## 2. Stack the 2 types of strains into a single pool (DataFrame) ##

    sampled_data_1 = __retrieve_super_strains(sampled_data)  # super strains
    sampled_data_2 = sampled_data.drop(sampled_data_1.index) # normal strains

    # Validation
    # Raise an error if one of the sampled data contains only one entry or no entries
    if len(sampled_data_1) == 1 or len(sampled_data_2) == 1:
        raise ValueError("One of the populations contain only one entry, which is insufficient.")
    if len(sampled_data_1) == 0 or len(sampled_data_2) == 0:
        raise ValueError("One of the populations contain no entries, which is insufficient.")

    simple_random_sample = pd.concat([sampled_data_1, sampled_data_2], axis=0).reset_index(drop=True) # Stack the 2 samples

    return simple_random_sample, len(sampled_data_1), len(sampled_data_2)

def __stratified_sampling(sampled_data_1, sampled_data_2):
    ## Stack the 2 types of strains into a single pool (DataFrame) ##    

    # Validation
    # Raise an error if one of the sampled data contains only one entry or no entries
    if len(sampled_data_1) == 1 or len(sampled_data_2) == 1:
        raise ValueError("One of the populations contain only one entry, which is insufficient.")
    if len(sampled_data_1) == 0 or len(sampled_data_2) == 0:
        raise ValueError("One of the populations contain no entries, which is insufficient.")

    stratified_sample = pd.concat([sampled_data_1, sampled_data_2], axis=0).reset_index(drop=True) # Stack the 2 samples

    return stratified_sample, len(sampled_data_1), len(sampled_data_2)

def __space_sampling(sampled_data, filtered_coords):
    ## Create 2 populations according to their coordinates in x axis. ##

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
    ## Create the necessary format of simulated data for the calculation of statistics ##

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

def __summary_stats(simulated_data):
    ## Calculate Tajima's D, Pi & θw ##

    ps = libsequence.PolySIM(simulated_data) 
    tajimasd = ps.tajimasd()
    pi = ps.thetapi()
    thetaw = ps.thetaw()
    
    return tajimasd, pi, thetaw

def __haplotypes(sampled_data):
    ## Calculates the number of unique sequences and the haplotype diversity ##

    # 1. DataFrame into dict (index:sequence)
    sample_dict = {index: ''.join(map(str, row)) for index, row in sampled_data.iterrows()}
    # 2. Get the total number of unique haplotypes
    unique_haplotypes = len(set(sample_dict.values()))
    # 3. Calculate haplotype diversity
    total_haplotypes = len(sample_dict)                                                                      # Total number of haplotypes (sequences)
    haplotype_frequencies = Counter(sample_dict.values())                                                    # Count the frequency of each unique haplotype
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # Compute haplotype diversity using the formula
    
    return unique_haplotypes, haplotype_diversity

def __selectice_sweep(simulated_data, ws, step):
    ## Calculate Tajima's D, Pi & θw for every sliding window ##

    tajimas_d = {}
    pi_est    = {}
    theta_w   = {}    
    w = libsequence.Windows(simulated_data, window_size=ws, step_len=step, starting_pos=0., ending_pos=1.0)
    for i in range(len(w)):
        wi = w[i]
        pswi = libsequence.PolySIM(wi)
        tajimas_d[i] = pswi.tajimasd()
        pi_est[i]    = pswi.thetapi()
        theta_w[i]   = pswi.thetaw()

    return tajimas_d, pi_est, theta_w

def statistics(sequences, strain_type, sample):
    ## Calculate Tajima's D, Pi, θw, the number of unique sequences and the haplotype diversity ##
    ## Calculations can be done to either super or normal strains, or both                      ##

    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling
    sampled_data = None  # Initialize sampled_data to avoid unbound error
    if strain_type   == "mix_strains":
        sampled_data = __sampling(all_seqs_data, sample, pop)
    elif strain_type == "ss_strains":
        sampled_data = __sampling(ss_seqs_data, sample, pop_ss)
    elif strain_type == "ns_strains":
        sampled_data = __sampling(ns_seqs_data, sample, pop_ns)

    if sampled_data is None:
        raise ValueError(f"Invalid strain_type '{strain_type}' or sampling failed.") # Handle the case where strain_type is invalid or sampling fails

    # 3) Simulating sampled data
    sim_data = simulated_data(sampled_data, l)

    if sim_data is not None:

        # 4) Tajima's D, Pi, θw
        tajimasd, pi, thetaw = __summary_stats(sim_data)

        # 5) Number of unique seqs, haplotype diversity
        uniq_haplo, haplo_div = __haplotypes(sampled_data)

    return tajimasd, pi, thetaw, uniq_haplo, haplo_div

def Fst(sequences, sampling_type, strain_type, sample, coords_file = None):
    ## Calculate Fst between 2 populations                                                                                                     ##
    ## Three Sampling techniques:                                                                                                              ##
    ##   1) Simple Random Sampling: collect n sequences from the pool without knowing how many are super or normal (len(super) != len(normal)) ##
    ##   2) Stratified Sampling: collect n sequences from normal strains and n from super strains (len(super) == len(normal))                  ##
    ##   3) Sampling according to space: divide the population into to 

    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling  
    edited_samples = None 
    if sampling_type == "simple_random_sampling":
        sampled_data = __sampling(all_seqs_data, sample, pop)
        edited_samples, n_pop1, n_pop2 = __simple_random_sampling(sampled_data)

    elif sampling_type == "stratified_sampling":
        sampled_data_ss = __sampling(ss_seqs_data, sample, pop_ss)
        sampled_data_ns = __sampling(ns_seqs_data, sample, pop_ns)
        edited_samples, n_pop1, n_pop2 = __stratified_sampling(sampled_data_ss, sampled_data_ns)

    elif sampling_type == "coordinates_sampling": 
        if strain_type == "mix_strains":
            filtered_coords = __filtering_sequences(all_seqs_data, coords_file) # get x axis coordinates
            sampled_data = __sampling(all_seqs_data, sample, pop)
        elif strain_type == "ss_strains":
            filtered_coords = __filtering_sequences(ss_seqs_data, coords_file) # get x axis coordinates
            sampled_data = __sampling(ss_seqs_data, sample, pop_ss)
        elif strain_type == "ns_strains":
            filtered_coords = __filtering_sequences(ns_seqs_data, coords_file) # get x axis coordinates
            sampled_data = __sampling(ns_seqs_data, sample, pop_ns)
        edited_samples, n_pop1, n_pop2 = __space_sampling(sampled_data, filtered_coords)
    if edited_samples is None:
        raise ValueError(f"Invalid strain_type '{strain_type}' or sampling failed.") # Handle the case where strain_type is invalid or sampling fails
    
    # 3) Simulating sampled data
    sim_data = simulated_data(edited_samples, l)

    if sim_data is not None:

        # 4) Fst
        fst = libsequence.Fst(sim_data,[n_pop1,n_pop2])

    return fst  

def LD(sequences, strain_type, sample):
    ## Calculate linkage_disequilibrium                                    ##
    ## Calculations can be done to either super or normal strains, or both ##

    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling
    if strain_type == "mix_strains":
        sampled_data = __sampling(all_seqs_data, sample, pop)
    elif strain_type == "ss_strains":
        sampled_data = __sampling(ss_seqs_data, sample, pop_ss)
    elif strain_type == "ns_strains":
        sampled_data = __sampling(ns_seqs_data, sample, pop_ns)
        
    # 3) Simulating sampled data
    sim_data = simulated_data(sampled_data, l)

    if sim_data is not None:

        # 4) LD values
        ld = pd.DataFrame(libsequence.ld(sim_data))

        return ld

def selective_sweep(sequences, window_size, step_size, strain_type, sample):
    ## Calculate Tajima's D, Pi & θw for every sliding window              ##
    ## Calculations can be done to either super or normal strains, or both ##

    # Validation for window_size and step_size
    if not isinstance(window_size, float):
        return print("Window size must be a float between 0 and 1.")
    if not isinstance(step_size, float):
        return print("Step size must be a float between 0 and 1.")
    if not (0 < window_size <= 1):
        return print("Window size must be a float between 0 and 1.")
    if not (0 < step_size <= 1):
        return print("Step size must be a float between 0 and 1.")
    
    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling
    if strain_type == "mix_strains":
        sampled_data = __sampling(all_seqs_data, sample, pop)
    elif strain_type == "ss_strains":
        sampled_data = __sampling(ss_seqs_data, sample, pop_ss)
    elif strain_type == "ns_strains":
        sampled_data = __sampling(ns_seqs_data, sample, pop_ns)
        
    # 3) Simulating sampled data
    sim_data = simulated_data(sampled_data, l)

    if sim_data is not None:

        # 4) Tajima's D, Pi & θw for every sliding window
        tajimas_d, pi_est, theta_w = __selectice_sweep(sim_data, window_size, step_size)
    
    # DataFrame with the results
    selsw_data = {
        "Window": tajimas_d.keys(),
        "Tajima's D": tajimas_d.values(),
        "Pi": pi_est.values(),
        "Theta W": theta_w.values()}
        
    selsw_df = pd.DataFrame(selsw_data)
    
    return selsw_df

def site_frequency_spectrum(sequences, strain_type, sample):
    
    # 1) Import data
    all_seqs_data, ss_seqs_data, ns_seqs_data, l, pop, pop_ss, pop_ns = import_data(sequences)

    # 2) Sampling
    sampled_data = None  # Initialize sampled_data to avoid unbound error
    if strain_type   == "mix_strains":
        sampled_data = __sampling(all_seqs_data, sample, pop)
    elif strain_type == "ss_strains":
        sampled_data = __sampling(ss_seqs_data, sample, pop_ss)
    elif strain_type == "ns_strains":
        sampled_data = __sampling(ns_seqs_data, sample, pop_ns)

    if sampled_data is None:
        raise ValueError(f"Invalid strain_type '{strain_type}' or sampling failed.") # Handle the case where strain_type is invalid or sampling fails

    # 3) Calculating Site Frequency Spectrum
    allele_counts = sampled_data.sum(axis=0)                       # Sum the number of SNP's for each genomic position 
    sfs = allele_counts.value_counts().sort_index()                # Count occurrences and sort by SNP count
    sfs_dataframe = sfs.reset_index()                              # Convert Series to DataFrame
    sfs_dataframe.columns = ['SNP Count', 'Frequency']             # Rename columns for clarity
    sfs_dataframe = sfs_dataframe[sfs_dataframe['SNP Count'] != 0] # Exclude rows where 'SNP Count' is 0

    return sfs_dataframe
#%% 
if __name__ == "__main__":

    '''1. Define parsers '''

    ## General Parsers
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-sample', '--sample_size', required=True, type=int, help='Sample size.')
    parser.add_argument('-p', '--population_size', action="store_true", help='Population size.')

    ## Parsers that define which statistics and metrics will be calculated
    # Parsers for summary statistics & Fst
    parser.add_argument('-sumstats', '--summary_statistics', action="store_true", help='Calculate summary statistics & Fst.')
    parser.add_argument('-sample_type', '--sampling_technique', type=str, help='Define the way the sampling of 2 populations will happen. 3 ways: "rdm","str","coords"')
    parser.add_argument('-c', '--coords_file', type=str, help='The path to the coordinates file. Required if -sampling_techning == "coords".')
    # Parsers for Linkage Disequilibrium
    parser.add_argument('-ld', '--linkage_disequilibrium', action="store_true", help='Calculate Linkage Disequilibrium.')
    # Parsers for Selective Sweep
    parser.add_argument('-selsw', '--selective_sweep', action="store_true", help='Search for Selective Sweep.')
    parser.add_argument('-winsize', '--window_size', type=float, help='Length of sliding window.')
    parser.add_argument('-stepsize', '--step_size', type=float, help='Step size.')
    # Parsers for site frequency spectrum
    parser.add_argument('-sfs', '--site_frequency_spectrum', action="store_true", help='Calculate Site Frequency Spectrum.')

    ## Mutually exclusive group for strain statistics
    strain_group = parser.add_mutually_exclusive_group(required=True)
    strain_group.add_argument('-mix', '--ss_and_ns_strains', action="store_true", help='Calculate statistics for both strains.')
    strain_group.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    strain_group.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    
    args = parser.parse_args()


    ''' 2. Parsers Validation'''

    # Check if selective sweep is specified, then window_size and step_size are mandatory
    if args.selective_sweep and (args.window_size is None or args.step_size is None):
        parser.error("--selective_sweep requires --window_size and --step_size to be specified.")
    # Check if a sampling technique is specified, when Fst is calculated
    if args.summary_statistics and args.sampling_technique is None:
        parser.error("--summary_statistics requires --sampling_technique to be specified.")
    # Check if the file containing the coordinates is required
    if (args.summary_statistics and args.sampling_technique == "coords") and not args.coords_file:
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
    if args.ss_and_ns_strains:
        strain_type = "mix_strains"
    elif args.ss_strains:
        strain_type = "ss_strains"
    elif args.ns_strains:
        strain_type = "ns_strains"
    else:
        strain_type = None

    # Determine suffix (e.g. File:genomes_100.csv -> suffix:100)
    try:
        genome_file = args.genome_file
        sample_size = args.sample_size
        suffix = event_number(genome_file)
    except UnboundLocalError:
        suffix = None 
        genome_file = None
        sample_size = None

    # Deternine names
    ld_filename = f"ld_{suffix}.csv"
    selsw_filename = f"selsw_{suffix}.csv"
    sfs_filename = f"sfs_{suffix}.csv"

    ## Summary Statistics & Fst ##
    if args.summary_statistics:
        tajimasd, pi, thetaw, uniq_haplo, haplo_div = statistics(args.genome_file, strain_type, args.sample_size)
        if args.sampling_technique == "rdm":
            sampling_type = "simple_random_sampling"
            fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = None)
        elif args.sampling_technique == "str":
            sampling_type = "stratified_sampling"
            fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = None)
        elif args.sampling_technique == "coords":
            sampling_type = "coordinates_sampling"
            __validate_files(args.genome_file, args.coords_file)
            fst = Fst(args.genome_file, sampling_type, strain_type, args.sample_size, coords_file = args.coords_file)
        
        hsm_fst = fst.hsm()         # Hudson, Slatkin & Maddison Fst
        slatkin_fst = fst.slatkin() # Slatkin Fst
        hbk_fst = fst.hbk()         # Hudson, Boos & Kaplan Fst

        # DataFrame with the results
        sumstats = {
            'Tajimas D': [tajimasd],
            'Pi-estimator': [pi],
            'Theta Watterson': [thetaw],
            'Number of unique haplotypes': [uniq_haplo],
            'Haplotype Diversity': [haplo_div],
            'HSM Fst': [hsm_fst],
            'Slatkin Fst': [slatkin_fst],
            'HBK Fst': [hbk_fst]}
        sumstats_df = pd.DataFrame(sumstats)

        # Save the DataFrame to a CSV file
        sumstats_df.to_csv(f"sumstats_{suffix}.csv", index=False)
        print(f"The results are stored in the file: sumstats_{suffix}.csv")

    ## Linkage Disequilibrium ##
    if args.linkage_disequilibrium:
        
        ld_df = LD(genome_file, strain_type, sample_size)
        ld_df.to_csv(ld_filename, index=False) 
        print(f"The results are stored in the file: ld_{suffix}.csv")

    ## Selective Sweep ##
    if args.selective_sweep:
        
        selsw_df = selective_sweep(args.genome_file, args.window_size, args.step_size, strain_type, args.sample_size)
        selsw_df.to_csv(selsw_filename, index=False)
        print(f"The results are stored in the file: selsw_{suffix}.csv")

    ## Site Frequency Spectrum ##
    if args.site_frequency_spectrum:

        sfs_df = site_frequency_spectrum(args.genome_file, strain_type, args.sample_size)
        sfs_df.to_csv(sfs_filename, index=False)
        print(f"The results are stored in the file: sfs_{suffix}.csv")