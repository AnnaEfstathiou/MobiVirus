from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool

from preprocessing import process_csv
from statistics import watterson_estimator, pi_estimator

def process_window_chunk(chunk_info: Tuple[Dict[str, str], int, int, int]) -> Tuple[List[float], List[float]]:
    
    """
    Process a chunk of windows in parallel.
    
    Args: chunk_info: A tuple containing:
                         - sequences (Dict[str, str]): Dictionary of sequences.
                         - start (int): Start index of the window chunk.
                         - end (int): End index of the window chunk.
                         - window_size (int): Size of the sliding window.
    Returns: Tuple[List[float], List[float]]: Lists of Watterson and Pi estimator values.
    """
    
    sequences, start, end, window_size = chunk_info
    step_size = 20 # step_size: Step size of the sliding window.
    theta_w_list = []
    pi_est_list = []
    
    # Iterate over each window in the chunk
    for i in range(start, end, step_size):
        # Extract the subsequence for each sequence
        window_sequences = {key: seq[i:(i+window_size)] for key, seq in sequences.items()}
        # Compute estimators and append results
        theta_w = watterson_estimator(window_sequences)
        theta_w_list.append(theta_w)
        pi_est = pi_estimator(window_sequences)
        pi_est_list.append(pi_est)

    return theta_w_list, pi_est_list

def sliding_window(sequences: Dict[str, str], window_size: int, num_processes: int = 4) -> Tuple[List[float], List[float], int]:
    
    """
    Apply sliding window approach in parallel to compute Watterson and Pi estimators.
    
    Args: sequences (Dict[str, str]): Dictionary of sequences.
          window_size (int): Size of the sliding window.
          num_processes (int): Number of parallel processes to use.
    Returns: Tuple[List[float], List[float], int]: Lists of Watterson and Pi estimator values, and sequence length.
    """

    if not sequences:
        raise ValueError("No sequences found. Please check the input CSV file.")
    
    sequence_length = len(next(iter(sequences.values()))) # length of sequences
    num_windows = sequence_length - window_size + 1       # number of windows
    
    # Calculate chunk size for each process
    chunk_size = max(1, num_windows // num_processes)  # Ensure chunk_size is at least 1

    # Prepare arguments for each process
    chunk_info = [(sequences, i, min(i + chunk_size, num_windows), window_size) for i in range(0, num_windows, chunk_size)]

    # Use multiprocessing pool to process window chunks in parallel
    with Pool(num_processes) as pool:
        results = pool.map(process_window_chunk, chunk_info)

    # Flatten the results from each process
    theta_w_list = [item for sublist in (result[0] for result in results) for item in sublist]
    pi_est_list = [item for sublist in (result[1] for result in results) for item in sublist]

    return theta_w_list, pi_est_list, sequence_length

def plot_selective_sweep(theta_w: List[float], pi_est: List[float], window_size: int, num_sequences: int, sequence_length: int, save_png: str = None):
    
    """
    Plot the results of the sliding window estimations.
    
    Args: theta_w (List[float]): List of Watterson estimator values.
          pi_est (List[float]): List of Pi estimator values.
          window_size (int): Size of the sliding window.
          num_sequences (int): Number of sequences.
          sequence_length (int): Length of the sequences.
          save_png (str, optional): Path to save the plot as a PNG file.
    """
    
    step_size = 20 # step_size: Step size of the sliding window.
    
    plt.figure(figsize=(10, 6))
    ## Plot Watterson and Pi estimator values ##
    plt.plot(theta_w, marker='.', linestyle='solid', color='navy', label='Watterson Estimator')
    plt.plot(pi_est, marker='.', linestyle='solid', color='firebrick', label='Pi Estimator')
    ## Add a vertical line at a specific SNP position ##
    plt.axvline(x=(sequence_length/2), color='green', linestyle='dashed', label='SNP position')
    ## Add metadata to the plot ##
    # plt.plot([], [], ' ', label=f'Number of Sequences: {num_sequences}')
    # plt.plot([], [], ' ', label=f'Sequence Length: {sequence_length}')
    plt.xlabel('Window Position')
    plt.ylabel('Estimator Value')
    plt.title('Selective Sweep', loc='left')
    plt.title(f'Window size:{window_size}, Step size:{step_size}, Seqs length:{sequence_length}, Number of seqs: {num_sequences}', loc='right')
    plt.legend()
    ## Save or display the plot ##
    if save_png:
        plt.savefig(save_png, format="png")
    else:
        plt.show()
    plt.close()

def main(csv_file: str, window_size: int, save_png: str = None):

    """
    Main function to process CSV file, compute estimators, and plot results.
    
    Args: csv_file (str): Path to the CSV file containing sequences.
          window_size (int): Size of the sliding window.
          save_png (str, optional): Path to save the plot as a PNG file.
    """
    
    try:
        sequences = process_csv(csv_file) # Process CSV file to get sequences
        if not sequences:
            print("No sequences loaded. Please check the a CSV file.")
        else:
            print(f"Loaded {len(sequences)} sequences.")
            
            theta_w, pi_est, seq_len = sliding_window(sequences, window_size) # Compute estimators using sliding window approach
            plot_selective_sweep(theta_w, pi_est, window_size, len(sequences), seq_len, save_png) # Plot the results
        
    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-w', '--window_size', type=int, required=True, help='Length of sliding window.')
    parser.add_argument('-save', '--save_png', type=str, help='Path to save the plot as a PNG file.')
    args = parser.parse_args()

    # Execute main function with parsed arguments
    main(args.genome_file, args.window_size, args.save_png)
