import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from preprocessing import process_csv

def plot_sfs(num_snps, count_pos, pos, max_snps, save_png=None):
    
    """ Bar plot of the Site Frequency Spectrum. """
    
    # Creating the bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(num_snps, count_pos, color='indigo')
    plt.legend([f'Pos of max SNPs ({max_snps}): {pos}'])
    plt.xlabel('Count of Positions')
    plt.ylabel('Number of SNPs')
    plt.title('Site Frequency Spectrum')
    plt.grid(True)
    plt.tight_layout()
    
    # Save the plot if save_png is provided, otherwise show the plot
    if save_png:
        plt.savefig(save_png, format="png")
    else:
        plt.show()
    plt.close()

def main(csv_file, save_png=None):
    
    sequences = process_csv(csv_file)                  # Process the input CSV file and extract binary sequences in a dictionary form
    num_sequences = len(sequences)                     # Number of sequences
    seq_length = len(sequences[next(iter(sequences))]) # Length of sequences

    print(f"Loaded {num_sequences} sequences.")

    ''' Checks '''
    if not sequences:
        raise ValueError("No sequences loaded. Please check the CSV file.")
    for seq in sequences.values():
        if len(seq) != seq_length:
            raise ValueError("Sequences are not of the same length.")
    if num_sequences == 0 or seq_length == 0:
        raise ValueError("Number of sequences or sequence length is zero.")
    
    # 1. Convert sequences to a numpy array and count SNPs at each position
    seq_array = np.array([list(seq) for seq in sequences.values()], dtype=int)
    count_snps = np.sum(seq_array, axis=0)
    
    # 2. Use Counter to count occurrences and sort by SNP count
    snp_dict = Counter(count_snps)
    snp_dict.pop(0, None)  # Remove zero occurrences if present
    sorted_snp_dict = dict(sorted(snp_dict.items()))
    
    # 3. Calculate the number of SNPs and the positions
    num_snps = list(sorted_snp_dict.keys())    # Number of SNPs
    count_pos = list(sorted_snp_dict.values()) # Count of Positions

    # 4. Find the position with the maximum number of SNPs
    max_index = np.argmax(count_snps)
    max_value = count_snps[max_index]

    # Create the bar plot
    plot_sfs(num_snps, count_pos, max_index, max_value, save_png)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV file to compute the Site Frequency Spectrum (SFS) and plot it.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-s', '--save_png', type=str, help='Path to save the plot as a PNG file.')
    args = parser.parse_args()

    main(args.genome_file, args.save_png)
