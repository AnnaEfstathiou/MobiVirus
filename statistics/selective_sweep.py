import argparse
import matplotlib.pyplot as plt
import random

from preprocessing import process_csv
from statistics import watterson_estimator, pi_estimator


def filter_super_strains(sequences):

    """ Filter all super strains. """

    filtered_dict = {}
    for index, sequence in sequences.items():
        length = len(sequence)
        
        ## Determine if the middle pos has a SNP
        if length % 2 == 1:
            # Odd length, single middle element
            middle_index = length // 2
            if sequence[middle_index] == '1':
                filtered_dict[index] = sequence
        else:
            # Even length, check two middle elements
            middle_left = (length // 2) - 1
            middle_right = length // 2
            if sequence[middle_left] == '1' or sequence[middle_right] == '1':
                filtered_dict[index] = sequence

    return filtered_dict


def window_slide(seqs, window_size, step_size):

    theta_w_dict = {}
    pi_est_dict = {}
    
    for i in range(0,len(next(iter(seqs.values()))),step_size):

        window_sequences = {key: seq[i:(i+window_size)] for key, seq in sequences.items()}
        theta_w = watterson_estimator(window_sequences)
        theta_w_dict[i] = theta_w
        pi_est = pi_estimator(window_sequences)
        pi_est_dict[i] = pi_est

    return theta_w_dict, pi_est_dict


def plot_selective_sweep(theta_w_dict, pi_est_dict, window_size, step_size, num_sequences, sequence_length, save_png = None):

    # Extract the keys (x values) and values (y values)
    window_pos = list(theta_w_dict.keys())
    theta_w_values = list(theta_w_dict.values())
    pi_est_values = list(pi_est_dict.values())

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(window_pos, theta_w_values, marker='.', linestyle='solid', color='navy', label='Watterson Estimator')
    plt.plot(window_pos, pi_est_values, marker='.', linestyle='solid', color='firebrick', label='Pi Estimator')

    # Add labels and title
    plt.xlabel("Window position")
    plt.ylabel("Estimator value")
    plt.title("Selective Sweep", loc='left')
    plt.title(f'Window size:{window_size}, Step size:{step_size}, Seqs length:{sequence_length}, Number of seqs: {num_sequences}', loc='right')
    plt.legend()

    # Show the plot
    plt.grid(True)
    ## Save or display the plot ##
    if save_png:
        plt.savefig(save_png, format="png")
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-w', '--window_size', type=int, required=True, help='Length of sliding window.')
    parser.add_argument('-s', '--step_size', type=int, required=True, help='Step size.')
    parser.add_argument('-sample', '--sample_size', type=int, help='Sample size.')
    parser.add_argument('-ss_sample', '--super_strains_sample', action='store_true', help='Sample size.')
    parser.add_argument('-save', '--save_png', type=str, help='Path to save the plot as a PNG file.')
    args = parser.parse_args()

    # Assign parsers to variables
    genomes = args.genome_file
    step_size = args.step_size
    window_size = args.window_size

    # Get all viral sequences
    sequences = process_csv(genomes)
    print("Total number sequences",len(sequences))

    if args.sample_size:
        # Plot pi estimator and theta watterson for a sample of viral sequences
        print(f"Taking randomly {args.sample_size} sequences.")
        sampled_seqs = dict(random.sample(list(sequences.items()), args.sample_size))
        print("Starting calculating statistics.")
        theta_w_dict, pi_est_dict =  window_slide(sampled_seqs, window_size, step_size)
        print("Done with calculations")
        print("Proceeding with plotting.")
        plot_selective_sweep(theta_w_dict, pi_est_dict, window_size, step_size, len(sampled_seqs), len(next(iter(sampled_seqs.values()))), args.save_png)
    elif args.super_strains_sample:
        # Plot pi estimator and theta watterson for all super strains
        ss_seqs = filter_super_strains(sequences)
        print("Number of super strains:",len(ss_seqs))
        print("Starting calculating statistics.")
        theta_w_dict, pi_est_dict =  window_slide(ss_seqs, window_size, step_size)
        print("Done with calculations")
        print("Proceeding with plotting.")
        plot_selective_sweep(theta_w_dict, pi_est_dict, window_size, step_size, len(ss_seqs), len(next(iter(ss_seqs.values()))), args.save_png)
    else:
        # Plot pi estimator and theta watterson for all viral sequences
        print("Starting calculating statistics.")
        theta_w_dict, pi_est_dict =  window_slide(sequences, window_size, step_size)
        print("Done with calculations")
        print("Proceeding with plotting.")
        plot_selective_sweep(theta_w_dict, pi_est_dict, window_size, step_size, len(sequences), len(next(iter(sequences.values()))), args.save_png)