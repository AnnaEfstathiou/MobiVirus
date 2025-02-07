import pandas as pd
import libsequence
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import re
import os

#%% Data manipulation

def import_data(csv_file):
    ## Convert the CSV file (genomes) to a DataFrame (rows: sequence / columns: positions) ##
    ## Retrieve normal and super spreaders, store them into seperate DataFrames            ##

    seqs = pd.read_csv(csv_file, header=None)  # CSV file --> DataFrame
    seqs = seqs.dropna()                       # Drop rows with NaN values
    seqs = seqs.astype(int)                    # Turn values to integers
    l = len(seqs.columns)                      # Genome length
    num_seqs = len(seqs)                       # Total number of strains

    return seqs, l, num_seqs

def event_number(csv_file):
    ## Isolate suffix number from the csv file ##
    ## Suffix represents the event number      ##

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

#%% Sampling

def sampling(data, sample, num_seqs):
    ## Sampling n sequences from a pool                                             ##
    ## If sample size is bigger than the population size, take the whole population ##

    if not data.empty:
        if sample <= num_seqs:
            sampled_data = data.sample(n=sample, random_state=int(datetime.now().timestamp())) 
        else:
            sampled_data = data.sample(n=num_seqs, random_state=int(datetime.now().timestamp())) 
        
        # Raise an error if the sampled data contains only one entry
        if len(sampled_data) == 1:
            raise ValueError("Sampled data contains only one entry, which is insufficient.")
    else:
        sampled_data = None
        raise ValueError("No data to sample from.") # Raise an error if sampled_data is None
    
    return sampled_data

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

#%% Theta Watterson for a CSV file

def thetaW_calc(sequences, sample):
    ## Calculate θw ##

    # 1) Import data
    seqs, l, num_seqs = import_data(sequences)

    # 2) Sampling
    sampled_data = sampling(seqs, sample, num_seqs)

    # 3) Simulating sampled data
    sim_data = simulated_data(sampled_data, l)

    if sim_data is not None:

        # 4) Calculating Theta Watterson
        ps = libsequence.PolySIM(sim_data) 
        thetaw = ps.thetaw()

    return thetaw, num_seqs

#%% Theta Watterson for multiple CSV files 

def thetaW_dir(directory, sample):

    ## Define directory and file names
    genome_dir = os.path.join(directory, 'genomes/')
    samples_dir = os.path.join(directory, 'samples/')
    event_time_filename = os.path.join(samples_dir, 'event_type.csv')
    
    ## Get the genome, sample and time files
    genome_files = {event_number(f): f for f in os.listdir(genome_dir) if f.startswith('genomes') and f.endswith('.csv')}
    simulation_times = pd.read_csv(event_time_filename, usecols=['Event', 'Simulation Time'], index_col='Event').squeeze().to_dict() # Convert simulation time and events into a dictionary
    
    ## Initialize a dictionary to store θw values for every simulation time ##
    thetaw_dict = {} 
    
    ## Process matching files by suffix, to get corresponding generations
    for suffix, genome_file in genome_files.items():
        file_path = os.path.join(genome_dir, genome_file)
        
        try:            
            thetaw, num_seqs = thetaW_calc(file_path, sample)
            time = simulation_times[int(suffix)] # Get the simulation time the current event happened

            thetaw_dict[time] = (thetaw, num_seqs)
       
        except ValueError as e:
            print(f"Error processing {genome_file}: {e}")
            continue

    ## Sort dictionaries according to simulation time (keys)
    thetaw_dict = dict(sorted(thetaw_dict.items()))

    return thetaw_dict

#%% Plotting Ne (according to θw) 

def plot_ne(dictionary, sample, mutation_rate, l, png_name):

    ## Calculating Ne (effective population size based on θw)
    times = list(dictionary.keys())
    thetaW_values = [value[0] for value in thetaw_dict.values()]
    num_seqs = [value[1] for value in thetaw_dict.values()]
    ne_values = [theta_w / (2 * mutation_rate * l) for theta_w in thetaW_values]  # Calculate Ne
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [3, 2, 1]})
    
    # Plot 1: Ne
    axs[0].plot(times, ne_values, color="#003f5c", linewidth=1.0, label="Effective Population Size (Ne)")
    
    # Smooth line using LOESS
    smoothed = lowess(ne_values, times, frac=0.1)
    axs[0].plot(smoothed[:, 0], smoothed[:, 1], color="#cf405d", linewidth=2.0, label="Smoothed Ne")
    
    axs[0].set_ylabel("Effective Population Size (Ne)", fontsize=12)
    axs[0].set_title("Approximate Effective Population Size", fontsize=14, weight="bold")
    axs[0].legend(fontsize=10)
    
    # Plot 2: Numner of viral strains
    axs[1].plot(times, num_seqs, color="#7a5195", linewidth=1.0, label="Number of Sequences")

    # Plot 3: θw
    axs[2].plot(times, thetaW_values, color="#ffa600", linewidth=1.5, label="ThetaW (θw)")
    axs[2].set_ylabel("ThetaW (θw)", fontsize=12)
    axs[2].set_xlabel("Simulation Time", fontsize=12)
    axs[2].legend(fontsize=10)
    
    ## Add overall title and adjust layout
    fig.suptitle(
        f"Sample = {sample}, Genome Length = {l}, Mutation Rate = {mutation_rate:.2e}", 
        fontsize=10, color="gray", y=0.95
    )
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Adjust space for the overall title

    # Save the plot
    output_filename = f"plotNe_{png_name}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    # plt.show()



#%%

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Calculate approximate Effective Population size for a simulation.")
    parser.add_argument('-d', '--directory', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-s', '--sample_size', type=int, required=True, help='Sample size.')
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help='Sample size.')
    parser.add_argument('-l', '--genome_length', type=int, required=True, help='Sample size.')
    args = parser.parse_args()

    directory = os.path.normpath(args.directory)
    directory_name = os.path.basename(directory)
    datetime_part = "_".join(directory_name.split('_')[1:])
    thetaw_dict = thetaW_dir(args.directory, args.sample_size)
    plot_ne(thetaw_dict, args.sample_size, args.mutation_rate, args.genome_length, datetime_part) 