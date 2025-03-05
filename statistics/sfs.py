import pandas as pd
import argparse
from datetime import datetime
import re
import matplotlib.pyplot as plt

#%% Data manipulation
def import_data(csv_file):
    ## Convert the CSV file (genomes) to a DataFrame (rows: sequence / columns: positions) ##

    seqs_data = pd.read_csv(csv_file, header=None)  # CSV file --> DataFrame
    seqs_data = seqs_data.dropna()                  # Drop rows with NaN values
    seqs_data = seqs_data.astype(int)               # Turn values to integers
    l = len(seqs_data.columns)                      # Genome length
    pop = len(seqs_data)                            # Total number of all sequences

    return seqs_data, l, pop

def event_number(csv_file):
    ## Isolate suffix number from the csv file ##
    ## Suffix represents the event number      ##

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

#%% Sampling
def __sampling(data, sample, pop):
    ## Sampling n sequences from a pool                                             ##
    ## If sample size is bigger than the population size, take the whole population ##

    if not data.empty:
        if sample <= pop:
            sampled_data = data.sample(n=sample, random_state=int(datetime.now().timestamp())) 
        else:
            sampled_data = data.sample(n=pop, random_state=int(datetime.now().timestamp())) 
        
        # Raise an error if the sampled data contains only one entry
        if len(sampled_data) == 1:
            raise ValueError("Sampled data contains only one entry, which is insufficient.")
    else:
        sampled_data = None
        raise ValueError("No data to sample from.") # Raise an error if sampled_data is None
    
    return sampled_data

#%% Site Frequency Spectrum
def site_frequency_spectrum(sequences, sample):
    ## Calculates the Site Frequency Spectrum (SFS) from a given set of genetic sequences. ##
    
    # 1) Import data
    seqs_data, l, pop = import_data(sequences)

    # 2) Sampling
    sampled_data = __sampling(seqs_data, sample, pop)
    
    # 3) Calculating Site Frequency Spectrum 
    allele_counts = sampled_data.sum(axis=0)                       # Sum the number of SNP's for each genomic position 
    sfs = allele_counts.value_counts().sort_index()                # Count occurrences and sort by SNP count
    sfs_dataframe = sfs.reset_index()                              # Convert Series to DataFrame
    sfs_dataframe.columns = ['SNP Count', 'Frequency']             # Rename columns for clarity
    sfs_dataframe = sfs_dataframe[sfs_dataframe['SNP Count'] != 0] # Exclude rows where 'SNP Count' is 0

    return sfs_dataframe

def count_snps_per_position(sequences, sample):
    ## Counts the number of Single Nucleotide Polymorphisms (SNPs) at each genomic position. ##

    # 1) Import data
    seqs_data, l, pop = import_data(sequences)

    # 2) Sampling
    sampled_data = __sampling(seqs_data, sample, pop)

    # 3) Count SNPs (1) per position
    snp_df = pd.DataFrame(sampled_data.sum(axis=0)).reset_index()
    snp_df.columns = ['Position', 'SNP Count']

    snp_df = snp_df[snp_df["SNP Count"] > 0]

    return snp_df

#%% Plotting
def plot_bar_chart(df: pd.DataFrame, x_col: str, y_col: str, title: str = "Bar Chart", xlabel: str = None, ylabel: str = None, plot_type = False, save_path: str = None):
    
    ## Plots a bar chart from a given DataFrame. ##
    ## Parameters:
        # df (pd.DataFrame): The DataFrame containing the data.
        # x_col (str): Column name to be used for the x-axis.
        # y_col (str): Column name to be used for the y-axis.
        # title (str): Title of the bar chart. Default is "Bar Chart".
        # xlabel (str): Label for the x-axis. 
        # ylabel (str): Label for the y-axis. 
        # save_path (str): File path to save the plot as a PNG file. If None, the plot is not saved.

    plt.figure(figsize=(10, 6))
    if plot_type:
        plt.bar(df[x_col], df[y_col], color='#204035FF', alpha=0.7)
    else:
        plt.bar(df[x_col], df[y_col], width=80, color='#204035FF')
    plt.xlabel(xlabel if xlabel else x_col)
    plt.ylabel(ylabel if ylabel else y_col)
    plt.title(title)
    plt.xticks(rotation=0)
    if plot_type:
        plt.xticks(df[x_col])
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    # Save the plot if a save_path is provided
    if save_path:
        plt.savefig(save_path, format='png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


#%% 
if __name__ == "__main__":

    ## Parsers
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute Site Frequency Spectrum (SFS) and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-s', '--sample_size', type=int, required=True, help='Sample size.')
    parser.add_argument('-save', '--save_csv', action="store_true", help='Save Site Frequency Spectrum values in a CSV file.')
    parser.add_argument('-png', '--save_png', action="store_true", help='Save Site Frequency Spectrum into a png.')
    args = parser.parse_args()

    ## Deternine names
    suffix = event_number(args.genome_file)
    sfs_filename = f"sfs_{suffix}.csv"
    snp_filename = f"snp_{suffix}.csv"
    sfs_plot_name = f"plotSFS_{suffix}.png"
    snp_plot_name = f"plotSNPcount_{suffix}.png"

    ## Site Frequency Spectrum 
    sfs_df = site_frequency_spectrum(args.genome_file, args.sample_size) # Site Frequency Spectrum
    snp_df = count_snps_per_position(args.genome_file, args.sample_size) # Analytic Site Frequency Spectrum
    print(snp_df)
    # print(snp_df[snp_df["SNP Count"] >= 1]["Position"])
    ## Save results in a CSV file
    if args.save_csv:
        sfs_df.to_csv(sfs_filename, index=False)
        snp_df.to_csv(snp_filename, index=False)
        print(f"The results are stored in the files: sfs_{suffix}.csv & snp_{suffix}.csv")

    ## Save plots in png format
    if args.save_png:
        plot_bar_chart(sfs_df, "SNP Count", "Frequency", "Site Frequency Spectrum", "SNP Count", "Frequency", plot_type = True, save_path = sfs_plot_name)
        plot_bar_chart(snp_df, "Position", "SNP Count", "Number of SNPs in genome positions", "Genome Positions", "SNP Count", plot_type = False, save_path = snp_plot_name)
    else:
        plot_bar_chart(sfs_df, "SNP Count", "Frequency", "Site Frequency Spectrum", "SNP Count", "Frequency", plot_type = True)
        plot_bar_chart(snp_df, "Position", "SNP Count", "Number of SNPs in genome positions", "Genome Positions", "SNP Count", plot_type = False)