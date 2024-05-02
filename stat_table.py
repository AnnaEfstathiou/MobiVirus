import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

# Check if 'calc_stats.py' exists in the current directory
if not os.path.exists('calc_stats.py'):
    raise ImportError("calc_stats.py is not in the current directory.")

from calc_stats import statistics 


def calc_stats_for_dir(directory, sample_size = None):
    """
    Process multiple CSV or FASTA files from the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.

    INPUT
    - directory (str): Directory containing the CSV or FASTA files to process.
    OUTPUT
    - stats_df: DataFrame with the statistics.
    """

    results = {}
    for filename in os.listdir(directory):
        if filename.endswith('.csv') or filename.endswith('.fasta') or filename.endswith('.fa'):
            file_path = os.path.join(directory, filename)
            try:
                # Attempt to call statistics and process the file
                tajimas_d_score, pi_estimator_score, watterson_estimator_score, unique_count, haplotype_diversity, num_seqs = statistics(file_path, sample_size)
                if num_seqs != 0:
                    num_unique_seqs = unique_count / num_seqs
                    num_unique_seqs_formatted = f"{num_unique_seqs:.2f} ({unique_count}/{num_seqs})"
                else:
                    num_unique_seqs_formatted = "0.0 (0/0)"
                results[filename] = [tajimas_d_score, pi_estimator_score, watterson_estimator_score, num_unique_seqs_formatted, haplotype_diversity]
            except ValueError as e:
                # Handle errors by skipping the file and optionally logging the error
                print(f"Error processing {filename}: {e} \n")
                # Remove the created FASTA files that raise an error.
                fasta_file_to_remove = os.path.splitext(file_path)[0] + "_processed.fa"
                if os.path.exists(fasta_file_to_remove):
                    os.remove(fasta_file_to_remove)
                continue  # Continue processing the next file

    stats_df = pd.DataFrame.from_dict(results, orient='index', columns=['Tajima\'s D', 'Pi-Estimator', 'Watterson-Estimator', 'Number of unique sequences', 'Haplotype Diversity'])
    # extract numbers from the filenames and sort accordingly
    # the csv named as: 'genomes_final.csv' will be the last one in the DataFrame
    stats_df['sort_key'] = stats_df.index.map(
        lambda x: float('inf') if x == 'genomes_final.csv' else (
            int(re.search(r'(\d+)', x).group(1)) if re.search(r'(\d+)', x) else 0))
    stats_df.sort_values('sort_key', inplace=True)
    stats_df.drop(columns=['sort_key'], inplace=True)  # remove the auxiliary column after sorting ('sort_key' column)
    return stats_df

def plot_summary_statistics(csv_file):

    """
    Plot the summary statistics.
    1st plot: Tajima's D, Pi-Estimator score, Watterson-Estimator score
    2nd plot: number of unique sequences represented as a ratio (e.g 415/559 = 0.742)
    3rd plot: haplotype diversity score

    INPUT
    - csv_file: The path to the csv file for which to run the stats.py script.
    """
    
    # read the CSV file
    csv_df = pd.read_csv(csv_file, index_col=0)  # adjust the index_col parameter as needed

    # convert 'Number of unique sequences' column to just the numeric value before the parentheses
    csv_df['Unique Sequences Value'] = csv_df['Number of unique sequences'].apply(lambda x: float(x.split(' ')[0]))

    fig, ax = plt.subplots(3, 1, figsize=(12, 18))

    # Tajima's D score, Pi-Estimator score, and Watterson-Estimator score in the first plot
    csv_df[["Tajima\'s D", "Pi-Estimator", "Watterson-Estimator"]].plot(ax=ax[0], marker='o')
    ax[0].set_title("Tajima's D, Pi-Estimator, and Watterson-Estimator Scores")
    ax[0].set_ylabel('Summary Statistics')
    ax[0].grid(True)

    # Number of unique sequences in the second plot
    csv_df["Unique Sequences Value"].plot(ax=ax[1], marker='o', color='purple')
    ax[1].set_title("Number of unique sequences (Ratio)")
    ax[1].set_ylabel('Ratio')
    ax[1].grid(True)

    # Haplotype Diversity in the third plot
    csv_df["Haplotype Diversity"].plot(ax=ax[2], marker='o', color='brown')
    ax[2].set_title("Haplotype Diversity")
    ax[2].set_ylabel('Score')
    ax[2].grid(True)

    plt.tight_layout(pad=4.0)
    plt.show()
    plt.close(fig) 

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

    df = calc_stats_for_dir(directory, sample_size)
    
    # print the DataFrame to stdout
    print(f"Sample size: {sample_size}\n")
    print(df)
    
    if args.store_dataframe:
        # save the DataFrame to a CSV file 
        df.to_csv(output_file)
        if args.plot_statistics:
            plot_summary_statistics(output_file)
   
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('directory', type=str, help='Directory containing CSV or FASTA files with genomes to process.')
    parser.add_argument('-n', '--output_filename', type=str, default='summary_statistics.csv', help='Output filename for the stored CSV file')
    parser.add_argument('-s', '--store_dataframe', action="store_true", help='Store the DataFrame in a CSV file')
    parser.add_argument('-p', '--plot_statistics', action="store_true", help='Plot')
    parser.add_argument('-m', '--sample_size', type=int, help='Output filename for the stored CSV file')
    args = parser.parse_args()
    
    # call main function with the specified directory
    main(args.directory, args.output_filename, args.sample_size)
     