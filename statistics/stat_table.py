import os
import pandas as pd
import numpy as np
import argparse
import re

# Check if 'calc_stats.py' exists in the current directory
if not os.path.exists('calc_stats.py'):
    raise ImportError("calc_stats.py is not in the current directory.")

from calc_stats import calculate_statistics 

def extract_suffix(filename):

    """
    Extract the suffix from the filename.
    For 'genomes_final.csv', return 'final'.
    For other files like 'genomes_150.csv', return '150'.
    """

    if 'final' in filename:
        return 'final'
    match = re.search(r'_(\d+)\.csv$', filename)
    return match.group(1) if match else None

def calc_stats_for_dir(directory, sample_size=None):

    """
    Process multiple CSV files from the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.
    
    Args: directory (str): Directory containing the CSV or FASTA files to process.
    Returns: stats_df: DataFrame with the statistics.
    """

    genome_dir = os.path.join(directory, 'genomes/')
    samples_dir = os.path.join(directory, 'samples/')
    results = {}

    # Get the genome and sample files
    genome_files = {extract_suffix(f): f for f in os.listdir(genome_dir) if f.startswith('genomes') and f.endswith('.csv')}
    coords_files = {extract_suffix(f): f for f in os.listdir(samples_dir) if f.startswith('coords') and f.endswith('.csv')}
    
    # Process matching files by suffix, to get corresponding generations
    for suffix, genome_file in genome_files.items():
        if suffix in coords_files:
            coords_file = coords_files[suffix]
            genome_path = os.path.join(genome_dir, genome_file)
            coords_path = os.path.join(samples_dir, coords_file)
            
            try:
                stats = calculate_statistics(genome_path, coords_path, sample_size) # dictionary containg all statistics
                total_sequences = stats['total_sequences']
                unique_count = stats['unique_count']
                
                if total_sequences:
                    if sample_size and total_sequences > sample_size:
                        num_unique_seqs = f"{unique_count / sample_size:.2f} ({unique_count}/{sample_size})"
                    else:
                        num_unique_seqs = f"{unique_count / total_sequences:.2f} ({unique_count}/{total_sequences})"
                else:
                    num_unique_seqs = "0.0 (0/0)"
                
                results[suffix] = [
                    stats['tajimas_d_score'],
                    stats['pi_estimator_score'],
                    stats['watterson_estimator_score'],
                    num_unique_seqs,
                    stats['haplotype_diversity'],
                    stats['Fst_coords'],
                    stats['Fst_inf_label'],
                    stats['Fst_mut_label']
                ]
            except ValueError as e:
                # Handle error cases
                if str(e) == "At least 2 sequences required!":
                    results[suffix] = ["Not enough sequences"] * 8
                else:
                    print(f"Error processing {genome_file}: {e}")
                    continue

    # Create DataFrame and sort by numeric suffix, with 'final' coming last
    stats_df = pd.DataFrame.from_dict(results, orient='index', 
                                      columns=['Tajima\'s D', 'Pi-Estimator', 'Watterson-Estimator', 'Number of unique sequences', 
                                               'Haplotype Diversity', 'Fst (coords)', 'Fst (label)', 'Fst (mutation)'])

    def sorting_key(suffix):

        """ Convert index to ensure correct sorting: 'final' comes last, numeric suffixes sorted naturally """
        
        if suffix == 'final':
            return float('inf') # Ensure 'final' comes last
        return int(suffix)      # Sort numerically for other suffixes

    stats_df.index = stats_df.index.map(sorting_key)
    stats_df.sort_index(inplace=True)

    # Convert the index back to string: remove decimals for numeric indices and retain 'final'
    stats_df.index = stats_df.index.map(lambda x: 'final' if x == float('inf') else str(int(x)))

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