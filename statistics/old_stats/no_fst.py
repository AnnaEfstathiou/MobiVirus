import os
import pandas as pd
import numpy as np
import argparse
import re

import argparse
import random
from preprocessing import process_csv
from statistics import tajimas_d, pi_estimator, watterson_estimator, count_haplotypes, calculate_haplotype_diversity, Fst


def calculate_statistics(input_file, sample_size = None):

    """
    Calculates various statistics for a given genome sequence dataset.
    
    Args: input_file (str): Path to the input CSV file with binary sequences.
          coords_file (str): Path to the coordinates file for population information.
          sample_size (int, optional): Number of sequences to sample. If None, use all sequences.
    Returns: dict: A dictionary with calculated statistics.
    """
    
    infected_sequences = process_csv(input_file) # Process the input CSV file and extract binary sequences in a dictionary form

    '''Sampling'''
    # If sample_size is specified and valid, sample the sequences
    if sample_size and sample_size < len(infected_sequences):
        sampled_keys = random.sample(list(infected_sequences.keys()), sample_size)
        sampled_infected_sequences = {key: infected_sequences[key] for key in sampled_keys}
    else:
        sampled_infected_sequences = infected_sequences  # Use all sequences if sample size is invalid or not specified

    # Calculate the desired statistics
    statistics = {
        "tajimas_d_score": tajimas_d(sampled_infected_sequences),
        "pi_estimator_score": pi_estimator(sampled_infected_sequences),
        "watterson_estimator_score": watterson_estimator(sampled_infected_sequences),
        "unique_count": count_haplotypes(sampled_infected_sequences),
        "haplotype_diversity": calculate_haplotype_diversity(sampled_infected_sequences),
        "total_sequences": len(infected_sequences),
    }

    return statistics

def extract_suffix(filename):

    """
    Extract the suffix from the filename.
    e.g. 'genomes_150.csv' > '150'
    """

    match = re.search(r'_(\d+)\.csv$', filename)
    return match.group(1) if match else None

def calc_stats_for_dir(directory, sample_size=None):

    """
    Process multiple CSV files from the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.
    
    Args: directory (str): Directory containing the CSV or FASTA files to process.
    Returns: stats_df: DataFrame with the statistics.
    """

    ## Define directory and file names
    genome_dir = os.path.join(directory, 'genomes/')
    samples_dir = os.path.join(directory, 'samples/')
    event_time_filename = os.path.join(samples_dir, 'event_type.csv')
    
    results = {} # Initialize a dictionary to store results

    ## Get the genome, sample and time files
    genome_files = {extract_suffix(f): f for f in os.listdir(genome_dir) if f.startswith('genomes') and f.endswith('.csv')}
    time_file = pd.read_csv(event_time_filename)[['Event', 'Simulation Time']]
    selected_times = time_file.set_index('Event')['Simulation Time'].to_dict() # get the corresponding simulation time for each event

    ## Process matching files by suffix, to get corresponding generations
    for suffix, genome_file in genome_files.items():
        
        genome_path = os.path.join(genome_dir, genome_file)
   
        try:
            ## Calculate and handle statistics 
            stats = calculate_statistics(genome_path, sample_size) # dictionary containg all statistics
            total_sequences = stats['total_sequences']
            unique_count = stats['unique_count']
            # Handling data based on the sample size
            if total_sequences:
                if sample_size and total_sequences > sample_size:
                    num_unique_seqs = f"{unique_count / sample_size:.2f} ({unique_count}/{sample_size})"
                else:
                    num_unique_seqs = f"{unique_count / total_sequences:.2f} ({unique_count}/{total_sequences})"
            else:
                num_unique_seqs = "0.0 (0/0)"
                # Dictionary with results
            results[suffix] = [
                    stats['tajimas_d_score'],
                    stats['pi_estimator_score'],
                    stats['watterson_estimator_score'],
                    num_unique_seqs,
                    stats['haplotype_diversity']]
            ## Handle error cases ##
        except ValueError as e:
            if str(e) == "At least 2 sequences required!":
                results[suffix] = ["Not enough sequences"] * 8
            else:
                print(f"Error processing {genome_file}: {e}")
                continue

    ## DataFrame ##
    stats_df = pd.DataFrame.from_dict(results, orient='index', 
                                      columns=['Tajima\'s D', 'Pi-Estimator', 'Watterson-Estimator', 
                                               'Number of unique sequences', 'Haplotype Diversity'])   
    stats_df = stats_df.reset_index().rename(columns={'index': 'Events'})  # Include the suffix (index) as a column
    stats_df['Events'] = stats_df['Events'].astype(int)                    # Ensure the 'Events' column is treated as integers before sorting
    stats_df['Time'] = stats_df['Events'].map(lambda x: selected_times.get(int(x), np.nan) if str(x).isdigit() else np.nan) # Add 'Time' as a column
    stats_df = stats_df.sort_values(by='Events', ascending=True).reset_index(drop=True) # Sort the DataFrame based on the 'Events' column and modify the DataFrame in-place
    
    return stats_df

def main(directory, output_file, sample_size):
    
    """
    Main function to create the DataFrame with the statistics.
    
    Args: directory (str): Directory containing the CSV files to process.
          output_file (str): Name of the stored CSV file with the statistics (default name: statistics_summary.csv).
    Returns: standard output: DataFrame
             (OPTIONAL) saved CSV file with all the summary statistics for the files in the given directory.
    """

    statistical_df = calc_stats_for_dir(directory, sample_size)

    # Print the DataFrame to stdout
    print(f"Sample size: {sample_size}\n")
    print(statistical_df)

    # Save the DataFrame to a CSV file if specified
    if args.store_dataframe:
        statistical_df.to_csv(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute statistics.")
    parser.add_argument('directory', type=str, help='Directory containing CSV or FASTA files with genomes to process.')
    parser.add_argument('-o', '--output_filename', type=str, default='summary_statistics.csv', help='Output filename for the stored CSV file')
    parser.add_argument('-save', '--store_dataframe', action="store_true", help='Store the DataFrame in a CSV file')
    parser.add_argument('-s', '--sample_size', type=int, help='Sample size for the analysis')
    
    args = parser.parse_args()

    # Call main function with the specified directory
    main(args.directory, args.output_filename, args.sample_size)