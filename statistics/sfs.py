import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from preprocessing import __csv_to_fasta, __read_fasta

def plot_sfs(num_snps, count_pos, pos, max_snps):

    """ Bar plot of the Site Frequency spectrum. """
    
    # Creating the bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(num_snps, count_pos, color='indigo')
    plt.legend([f'Pos of the {max_snps} SNPs: {pos}'])
    plt.xlabel('Count of Positions')
    plt.ylabel('Number of SNPs')
    plt.title('Site Frequency Spectrum')
    plt.grid(True) 
    plt.tight_layout()
    # Save the plot if save_png flag is provided, otherwise show the plot
    if args.save_png:
        plt.savefig(args.save_png, format="png")
    else:
        plt.show()
    plt.close()

def main(csv_file):

    output_fasta_file = 'temp_fasta.fa'
    __csv_to_fasta(csv_file, output_fasta_file)        # CSV to FASTA file
    sequences = __read_fasta(output_fasta_file)        # Read the FASTA file
    num_sequences = len(sequences)                     # number of sequences
    seq_length = len(sequences[next(iter(sequences))]) # length of sequences

    print(f"Loaded {num_sequences} sequences.")
        
    ## Checks ##
    if not sequences:
        raise ValueError("No sequences loaded. Please check the FASTA file.")
    for seq in sequences.values():
        if len(seq) != seq_length:
            raise ValueError("Sequences are not of the same length.")
    if num_sequences == 0 or seq_length == 0:
        raise ValueError("Number of sequences or sequence length is zero.")
    
    seq_array = np.array([list(seq) for seq in sequences.values()]) # Convert sequences to a numpy array
    seq_array = seq_array.astype(int)                               # Convert strings to integers
    count_snps = np.sum(seq_array, axis=0)                          # Count SNPs at each position (column)
    counter = Counter(count_snps)                                   # Use Counter to count occurrences
    snp_dict = dict(counter)                                        # Convert Counter object to dictionary
    if 0 in snp_dict.keys():                                        # Exclude the number of 0 occurrences 
        del snp_dict[0]                                             # No need to count how many pos have no SNP
    sorted_snp_dict = dict(sorted(snp_dict.items()))                # Sort by keys
    num_snps = list(sorted_snp_dict.keys())                         # Number of SNPs
    count_pos = list(sorted_snp_dict.values())                      # Count of Positions

    ## Find the position in the sequence where more SNPs exist ##
    max_index = -1                                                 
    max_value = float('-inf')                                   
    for i, j in enumerate(count_snps): # Iterate through count_snps using enumerate
        if j > max_value:              # Check if the current value (j) is greater than the current max_value
            max_value = j              # Update max_value
            max_index = i              # Update max_index

    # Create the bar plot
    plot_sfs(num_snps, count_pos, max_index, max_value)

    # Remove the temporarly FASTA file
    os.remove(output_fasta_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV file to compute the Site Frequency Spectrum (SFS) and plot it.")
    parser.add_argument('genome_file', type=str, help='The path to the input CSV file.')
    parser.add_argument('-s', '--save_png', type=str, help='Path to save the plot as a PNG file.')
    args = parser.parse_args()

    main(args.genome_file)