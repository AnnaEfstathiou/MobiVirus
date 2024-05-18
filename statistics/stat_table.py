import os
import pandas as pd
import argparse
import re

# Check if 'calc_stats.py' exists in the current directory
if not os.path.exists('calc_stats.py'):
    raise ImportError("calc_stats.py is not in the current directory.")

from calc_stats import calculate_statistics 


def calc_stats_for_dir(directory, sample_size=None):

    """
    Process multiple CSV or FASTA files from the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.

    INPUT
    - directory (str): Directory containing the CSV or FASTA files to process.
    OUTPUT
    - stats_df: DataFrame with the statistics.
    """
    
    genome_dir = directory + 'genomes/'
    samples_dir = directory + 'samples/'
    results = {}


    def extract_suffix(filename):
    # Implement the logic to extract the suffix from the filename
    # For example, if the suffix is everything after the last underscore:
        parts = filename.split('_')
        if len(parts) > 1:
            return parts[-1].split('.')[0]
        return None

    # Create a dictionary to store genome files with their suffix as keys
    genome_files = {}
    for genome_filename in os.listdir(genome_dir):
        if genome_filename.startswith('genomes') and (genome_filename.endswith('.csv') or genome_filename.endswith('.fasta') or genome_filename.endswith('.fa')):
            suffix = extract_suffix(genome_filename)
            if suffix:
                genome_files[suffix] = genome_filename

    # Match genome files with coordinates files based on suffix and process them
    for coords_filename in os.listdir(samples_dir):
        if coords_filename.startswith('coords') and coords_filename.endswith('.csv'):
            suffix = extract_suffix(coords_filename)
            if suffix and suffix in genome_files:
                genome_file = genome_files[suffix]
            
                try:
                    stats = calculate_statistics(genome_dir+genome_file, samples_dir+coords_filename, sample_size)
                    print(stats)
                    if sample_size:
                        if stats['total_sequences'] != 0:
                            if stats['total_sequences'] > sample_size:
                                num_unique_seqs = stats['unique_count'] / sample_size
                                num_unique_seqs_formatted = f"{num_unique_seqs:.2f} ({stats['unique_count']}/{sample_size})"
                            else: 
                                num_unique_seqs = stats['unique_count'] / stats['total_sequences']
                                num_unique_seqs_formatted = f"{num_unique_seqs:.2f} ({stats['unique_count']}/{stats['total_sequences']})"
                        else:
                            num_unique_seqs_formatted = "0.0 (0/0)"
                    else:
                        if stats['total_sequences'] != 0:
                            num_unique_seqs = stats['unique_count'] / stats['total_sequences']
                            num_unique_seqs_formatted = f"{num_unique_seqs:.2f} ({stats['unique_count']}/{stats['total_sequences']})"
                        else:
                            num_unique_seqs_formatted = "0.0 (0/0)"
                    
                    results[genome_file] = [
                        stats['tajimas_d_score'],
                        stats['pi_estimator_score'],
                        stats['watterson_estimator_score'],
                        num_unique_seqs_formatted,
                        stats['haplotype_diversity'],
                        stats['Fst_coords'],
                        stats['Fst_label']
                    ]
                except ValueError as e:
                    print(f"Error processing {genome_filename}: {e} \n")
                    # Remove the created FASTA files that raise an error.
                    fasta_file_to_remove = os.path.splitext(genome_filename)[0] + "_processed.fa"
                    if os.path.exists(fasta_file_to_remove):
                        os.remove(fasta_file_to_remove)
                    continue  # Continue processing the next file
    
    stats_df = pd.DataFrame.from_dict(results, orient='index', columns=['Tajima\'s D', 'Pi-Estimator', 'Watterson-Estimator', 'Number of unique sequences', 'Haplotype Diversity', 'Fst (coords)', 'Fst (label)'])
    # extract numbers from the filenames and sort accordingly
    # the csv named as 'genomes_final.csv' will be the last one in the DataFrame
    stats_df['sort_key'] = stats_df.index.map(
        lambda x: float('inf') if x == 'genomes_final.csv' else (
            int(re.search(r'(\d+)', x).group(1)) if re.search(r'(\d+)', x) else 0))
    stats_df.sort_values('sort_key', inplace=True)
    stats_df.drop(columns=['sort_key'], inplace=True)  # remove the auxiliary column after sorting ('sort_key' column)
    return stats_df
    # return stats
  


def main(directory, output_file, sample_size):

    """
    Main function to create the DataFrame with the statistics.

    INPUT
    - directory (str): Directory containing the CSV files to process.
    - output_file (str): Name of the stored CSV file with the statistics (default name: statistics_summary.csv).
    OUTPUT
    - standard output: DataFrame
    - (OPTIONAL) saved CSV file with all the summary statistics for the files in the given directory.
    - (OPTIONAL) plot the scores stored in the saved CSV file with all the summary statistics.
    """

    ## ERROR HANDLING ##
    # in order to plot the summary statistics, the dataframe must be saved into a CSV file
    if not args.store_dataframe:
        if args.plot_statistics:
            raise ValueError("Store the DataFrame with the summary statistics as a CSV file in order to plot those values! In order to do that use the 'store' flag along with the 'plots' one")

    statistical_df = calc_stats_for_dir(directory, sample_size)
    
    # print the DataFrame to stdout
    print(f"Sample size: {sample_size}\n")
    print(statistical_df)
    
    # if args.store_dataframe:
    #     # save the DataFrame to a CSV file 
    #     statistical_df.to_csv(output_file)
        # if args.plot_statistics:
        #     plot_summary_statistics(output_file)
   
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('directory', type=str, help='Directory containing CSV or FASTA files with genomes to process.')
    parser.add_argument('-o', '--output_filename', type=str, default='summary_statistics.csv', help='Output filename for the stored CSV file')
    parser.add_argument('-save', '--store_dataframe', action="store_true", help='Store the DataFrame in a CSV file')
    parser.add_argument('-s', '--sample_size', type=int, help='Output filename for the stored CSV file')
    args = parser.parse_args()
    
    # call main function with the specified directory
    main(args.directory, args.output_filename, args.sample_size)