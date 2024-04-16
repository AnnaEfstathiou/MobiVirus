import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def run_stats_script(file_path, stats_script_path='calc_stats.py'):

    """
    Runs the stats.py script for a csv file and captures its output.
    
    INPUT
    - file_path (str): The path to the csv file for which to run the stats.py script.
    - stats_script_path (str): The path to the calc_stats.py script.
    OUTPUT
    - stats: A dictionary containing parsed statistical results from stats.py output.
    """

    # check if stats.py exists at the provided path
    if not os.path.exists(stats_script_path):
        raise FileNotFoundError(f"The calc_stats.py script was not found at the specified path: {stats_script_path}")
    
    # construct the command to execute stats.py with the given file path
    command = f'python3 "{stats_script_path}" "{file_path}"'
    # run the command and capture stdout (standard output)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    output = result.stdout

    stats = {}
    # iterate through each line of the captured output
    for line in output.split('\n'):
        if line.startswith("Tajima's D score:") or \
           line.startswith("Pi-Estimator score:") or \
           line.startswith("Watterson-Estimator score:") or \
           line.startswith("Number of unique sequences:") or \
           line.startswith("Haplotype diversity:"):
            # split the line into key and value parts and strip whitespace
            key, value = line.split(':', 1)
            stats[key.strip()] = value.strip()
    return stats

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

    # convert 'Number of unique sequences' column to ratio 
    # assume the 'Number of unique sequences' column is formatted as '415/537'
    # if 'Number of unique sequences' column is formatted as '0/0' then is set equal to 0 
    csv_df['Unique Sequences Ratio'] = csv_df['Number of unique sequences'].apply(lambda x: 0 if x == '0/0' else eval(x))

    fig, ax = plt.subplots(3, 1, figsize=(12, 18))

    # Tajima's D score, Pi-Estimator score, and Watterson-Estimator score in the first plot
    csv_df[["Tajima's D score", "Pi-Estimator score", "Watterson-Estimator score"]].plot(ax=ax[0], marker='o')
    ax[0].set_title("Tajima's D, Pi-Estimator, and Watterson-Estimator Scores")
    ax[0].set_ylabel('Scores')
    ax[0].grid(True)

    # Number of unique sequences ratio in the second plot
    csv_df["Unique Sequences Ratio"].plot(ax=ax[1], marker='o', color='purple')
    ax[1].set_title("Unique Sequences Ratio")
    ax[1].set_ylabel('Ratio')
    ax[1].grid(True)

    # Haplotype diversity in the third plot
    csv_df["Haplotype diversity"].plot(ax=ax[2], marker='o', color='brown')
    ax[2].set_title("Haplotype Diversity")
    ax[2].set_ylabel('Diversity')
    ax[2].grid(True)

    #plt.tight_layout()
    plt.tight_layout(pad=4.0)
    plt.show()
    plt.close(fig) 

def main(directory, output_file):

    """
    Main function to process multiple CSV files in the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.

    INPUT
    - directory (str): Directory containing the CSV files to process.
    OUTPUT
    - standard output: DataFrame
    - (OPTIONAL) statistics_summary.csv: saved CSV file with all the statistics for the files in the directory
    - (OPTIONAL) plot the scores stored in the statistics_summary.csv
    """
    ## ERROR HANDLING ##
    # in order to plot the summary statistics, the dataframe must be saved into a CSV file
    if not args.store_dataframe:
        if args.plot_statistics:
            raise ValueError("Store the DataFrame with the summary statistics as a CSV file in order to plot those values! In order to do that use the 'store' flag along with the 'plots' one")

    results = {}
    # iterate over all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            # construct the full path to the file
            file_path = os.path.join(directory, filename)
            # run statistical analysis on the file and capture the results
            stats = run_stats_script(file_path)
            # only store non-empty results 
            # theoretically there are no empty results...the number of unique sequences and the haplotype diversity will always be calculated 
            # the Tajima's D, Pi-Estimator score & Watterson-Estimator score might not be calculated if the number of sequences are less than 4
            if stats:  
                # store the results in the dictionary using the filename as the key
                results[filename] = stats

    # convert the collected results into a pandas DataFrame
    df = pd.DataFrame(results).T  # transpose so that each row represents a file
    df = df.sort_index() # sort the dataframe based on the index 
    
    # print the DataFrame to stdout
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
    args = parser.parse_args()
    
    # call main function with the specified directory
    main(args.directory, args.output_filename)
    