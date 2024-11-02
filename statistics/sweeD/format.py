import pandas as pd
import argparse
from datetime import datetime
import numpy as np
import re

def allele_frequency_function(csv_file, sample_size):
    ## Data format consists of four columns, which describes four attributes for each SNP:    ##
    ## - location --> the location of a SNP                                                   ##
    ## - x --> the number of sequences carry the derived allele for a SNP                     ##
    ## - n --> the number of valid sequences at a SNP                                         ##
    ## - folded --> a binary character which denotes if the SNP is unfolded (0) or folded (1) ##

    genome_data = pd.read_csv(csv_file, header=None)  # CSV file --> DataFrame
    genome_data = genome_data.dropna()                # Drop rows with NaN values
    genome_data = genome_data.astype(int)             # Turn values to integers
    pop = len(genome_data)                            # Total number of all sequences
    
    # Sampling
    if sample_size <= pop:
        sampled_data = genome_data.sample(n=sample_size, random_state=int(datetime.now().timestamp())) 
    else:
        sampled_data = genome_data

    allele_frequency_data = [] # Initialize an empty list to store the results

    # Loop through each column (position) to count allele frequencies
    for position in sampled_data.columns:
        allele_counts = sampled_data[position].value_counts()               # Get the counts of alleles (0s and 1s) for the current position
        snp_count = allele_counts.get(1, 0)                                 # Count of allele 1
        allele_frequency_data.append([position, snp_count, sample_size, 0]) # Genome position , allele count (x), sample size (n), polarization (0)

    allele_frequency_df = pd.DataFrame(allele_frequency_data, columns=['position', 'x', 'n', 'folded']) # Convert the list to a DataFrame for easy export

    return allele_frequency_df

def event_number(csv_file):
    ## Isolate suffix number from the csv file ##
    ## Suffix represents the event number      ##

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

if __name__ == "__main__":

    ## General Parsers
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-sample', '--sample_size',type=int, required=True, help='Sample size.')    
    parser.add_argument('-o', '--output_filename',type=str, required=True, help='Output filenameeeeee.')   
    args = parser.parse_args()

    suffix = event_number(args.csv_file)
    filename = args.output_filename
    allele_frequency_data = allele_frequency_function(args.csv_file, args.sample_size)
    allele_frequency_data.to_csv(filename, index=False, sep='\t')  # Save to a CSV file