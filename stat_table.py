import os
import subprocess
import pandas as pd
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
        raise FileNotFoundError(f"The stats.py script was not found at the specified path: {stats_script_path}")
    
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


def main(directory):

    """
    Main function to process multiple CSV files in the specified directory,
    run statistical analysis on each, and compile results into a DataFrame.

    INPUT
    - directory (str): Directory containing the CSV files to process.
    OUTPUT
    - standard output: DataFrame
    - (OPTIONAL) statistics_summary.csv: saved CSV file with all the statistics for the files in the directory
    """

    results = {}
    # iterate over all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            # construct the full path to the file
            file_path = os.path.join(directory, filename)
            # run statistical analysis on the file and capture the results
            stats = run_stats_script(file_path)
            # store the results in the dictionary using the filename as the key
            results[filename] = stats

    # convert the collected results into a pandas DataFrame
    df = pd.DataFrame(results).T  # transpose so that each row represents a file
    
    # print the DataFrame to stdout
    print(df)
    ## OPTIONAL
    if args.store_dataframe:
        # save the DataFrame to a CSV file 
        df.to_csv("statistics_summary.csv")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    # arguments for directory path
    parser.add_argument('directory', type=str, help='Directory containing CSV or FASTA files with genomes to process.')
    parser.add_argument('-store', '--store_dataframe', action="store_true", help='Store the DataFrame in a CSV file')
    args = parser.parse_args()
    
    # call main function with the specified directory
    main(args.directory)
    