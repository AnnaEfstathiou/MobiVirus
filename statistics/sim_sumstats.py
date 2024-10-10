import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os

from statistics import statistics, Fst, __event_number

def sumstats_directory(directory, sample_size, strain_type, sampling_type):

    ## Define directory and file names
    genome_dir = os.path.join(directory, 'genomes/')
    samples_dir = os.path.join(directory, 'samples/')
    event_time_filename = os.path.join(samples_dir, 'event_type.csv')
    
    ## Get the genome, sample and time files
    genome_files = {__event_number(f): f for f in os.listdir(genome_dir) if f.startswith('genomes') and f.endswith('.csv')}
    simulation_times = pd.read_csv(event_time_filename, usecols=['Event', 'Simulation Time'], index_col='Event').squeeze().to_dict() # Convert simulation time and events into a dictionary
    if sampling_type == "coordinates_sampling":
         coords_files = {__event_number(f): f for f in os.listdir(samples_dir) if f.startswith('coords') and f.endswith('.csv')}
    
    ## Initialize directories to store statistics ##
    tajimasd_stats = {} # Initialize a dictionary to store Tajima's D values for every simulation time
    pi_stats = {}       # Initialize a dictionary to store Pi values for every simulation time
    thetaw_stats = {}   # Initialize a dictionary to store θw values for every simulation time
    uniq_stats = {}     # Initialize a dictionary to store number of unique seqs for every simulation time
    hp_div = {}         # Initialize a dictionary to store haplotype diversity for every simulation time
    fst_stats = {}      # Initialize a dictionary to store Fst values for every simulation time  
            
    ## Process matching files by suffix, to get corresponding generations
    for suffix, genome_file in genome_files.items():
        genome_data = os.path.join(genome_dir, genome_file)
        
        try:
            print(f"Processing file: {genome_file}")
            
            tajimasd, pi, thetaw, uniq_haplo, haplo_div = statistics(genome_data, strain_type, sample_size)
            time = simulation_times[int(suffix)] # Get the simulation time the current event happened

            if (sampling_type == "coordinates_sampling") and (suffix in coords_files):
                coords_file = coords_files[suffix]
                coords_data = os.path.join(samples_dir, coords_file)
                fst = Fst(genome_data, sampling_type, strain_type, sample_size, coords_file = coords_data)
            else:
                 fst = Fst(genome_data, sampling_type, strain_type, sample_size, coords_file = None)
            
            tajimasd_stats[time] = tajimasd
            pi_stats[time] = pi
            thetaw_stats[time] = thetaw
            uniq_stats[time] = uniq_haplo
            hp_div[time] = haplo_div
            fst_stats[time] = {'hsm': fst.hsm(), 'slatkin': fst.slatkin(), 'hbk': fst.hbk()}
                  
        except ValueError as e:
                    print(f"Error processing {genome_file}: {e}")
                    continue
    # Sort dictionaries according to simulation time (keys)
    tajimasd_stats = dict(sorted(tajimasd_stats.items()))
    pi_stats       = dict(sorted(pi_stats.items()))
    thetaw_stats   = dict(sorted(thetaw_stats.items()))
    uniq_stats     = dict(sorted(uniq_stats.items()))
    hp_div         = dict(sorted(hp_div.items()))
    fst_stats      = dict(sorted(fst_stats.items()))

    # Convert to DataFrame
    data = {
        'Time': list(tajimasd_stats.keys()),
        'TajimasD': list(tajimasd_stats.values()),
        'Pi': list(pi_stats.values()),
        'Theta_w': list(thetaw_stats.values()),
        'Unique_seqs': list(uniq_stats.values()),
        'Hp_Diversity': list(hp_div.values()),
        'Fst_hsm': [fst['hsm'] for fst in fst_stats.values()],
        'Fst_slatkin': [fst['slatkin'] for fst in fst_stats.values()],
        'Fst_hbk': [fst['hbk'] for fst in fst_stats.values()]
}
    stats = pd.DataFrame(data)

    return stats

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-d', '--directory', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-sample', '--sample_size', type=int, required=True, help='Sample size.')
    parser.add_argument('-sample_type', '--sampling_technique', type=str, required=True, help='Define the way the sampling of 2 populations will happen. 3 ways: "rdm","str","coords"')
    # Adding mutually exclusive group for strain type options
    strain_group = parser.add_mutually_exclusive_group(required=True)
    strain_group.add_argument('-mix', '--ss_and_ns_strains', action="store_true", help='Calculate statistics for both strains.')
    strain_group.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    strain_group.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    args = parser.parse_args()


    # Determine strain_type based on command-line arguments
    if args.ss_and_ns_strains:
        strain_type = "mix_strains"
    elif args.ss_strains:
        strain_type = "ss_strains"
    elif args.ns_strains:
        strain_type = "ns_strains"
    else:
        strain_type = None

    if args.sampling_technique == "rdm":
        sampling_type = "simple_random_sampling"
    elif args.sampling_technique == "str":
        sampling_type = "stratified_sampling"
    elif args.sampling_technique == "coords":
        sampling_type = "coordinates_sampling"
    else:
        sampling_type = None

    # Ensure strain_type is set properly
    if strain_type is None:
        print("Error: Please specify a strain type to calculate statistics (use -mix, -ss, or -ns).")
    else:
        stats = sumstats_directory(args.directory, args.sample_size, strain_type, sampling_type)
        # Saving statistics in a csv format
        directory = os.path.normpath(args.directory)
        directory_name = os.path.basename(directory)
        datetime_part = "_".join(directory_name.split('_')[1:])
        csv_filename = f"stats_{datetime_part}.csv"
        stats.to_csv(csv_filename)              