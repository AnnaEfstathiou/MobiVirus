import matplotlib.pyplot as plt
import pandas as pd
import argparse
import re

def extract_numeric_value(name):
    """Extract numeric value from a string using regex."""
    match = re.search(r'\d+', name)
    return int(match.group()) if match else None

def plot_summary_statistics(csv_file):
    """
    Plot the summary statistics.
    1st plot: Tajima's D, Pi-Estimator score, Watterson-Estimator score
    2nd plot: number of unique sequences represented as a ratio (e.g 415/559 = 0.742)
    3rd plot: haplotype diversity score

    INPUT
    - csv_file: The path to the csv file for which to run the stats.py script.
    """
    
    # Read the CSV file
    csv_df = pd.read_csv(csv_file, index_col=0)  # adjust the index_col parameter as needed

    # Extract numeric values from index (file names) for x-axis
    numeric_labels = [extract_numeric_value(label) for label in csv_df.index]

    # Convert 'Number of unique sequences' column to just the numeric value before the parentheses
    csv_df['Unique Sequences Value'] = csv_df['Number of unique sequences'].apply(lambda x: float(x.split(' ')[0]))

    fig, ax = plt.subplots(6, 1, figsize=(12, 22))  # Adjusted to 6 plots to accommodate new plot

    # Define a step size for x-axis labels (e.g., every 5th label)
    step = 5
    x_labels = [numeric_labels[i] if i % step == 0 else '' for i in range(len(numeric_labels))]

    # Plot Tajima's D score
    csv_df["Tajima's D"].plot(ax=ax[0], marker='o', color='orange')
    ax[0].set_title("Tajima's D")
    ax[0].set_ylabel('Score')
    ax[0].grid(True)
    ax[0].set_xticks(range(len(numeric_labels)))
    ax[0].set_xticklabels(x_labels)

    # Plot Pi-Estimator score
    csv_df["Pi-Estimator"].plot(ax=ax[1], marker='o', color='green')
    ax[1].set_title("Pi-Estimator")
    ax[1].set_ylabel('Score')
    ax[1].grid(True)
    ax[1].set_xticks(range(len(numeric_labels)))
    ax[1].set_xticklabels(x_labels)

    # Plot Watterson-Estimator score
    csv_df["Watterson-Estimator"].plot(ax=ax[2], marker='o', color='blue')
    ax[2].set_title("Watterson-Estimator")
    ax[2].set_ylabel('Score')
    ax[2].grid(True)
    ax[2].set_xticks(range(len(numeric_labels)))
    ax[2].set_xticklabels(x_labels)

    # Plot Haplotype Diversity
    csv_df["Haplotype Diversity"].plot(ax=ax[3], marker='o', color='brown')
    ax[3].set_title("Haplotype Diversity")
    ax[3].set_ylabel('Score')
    ax[3].grid(True)
    ax[3].set_xticks(range(len(numeric_labels)))
    ax[3].set_xticklabels(x_labels)

    # Plot Number of unique sequences
    csv_df["Unique Sequences Value"].plot(ax=ax[4], marker='o', color='purple')
    ax[4].set_title("Number of unique sequences (Ratio)")
    ax[4].set_ylabel('Ratio')
    ax[4].grid(True)
    ax[4].set_xticks(range(len(numeric_labels)))
    ax[4].set_xticklabels(x_labels)

    # Plot Fst types on the same subplot with different colors
    csv_df["Fst (coords)"].plot(ax=ax[5], marker='o', markersize=5, color='orangered', label='Fst (coords)')
    csv_df["Fst (label)"].plot(ax=ax[5], marker='o', markersize=5, color='limegreen', label='Fst (label)')
    csv_df["Fst (mutation)"].plot(ax=ax[5], marker='o', markersize=5, color='dimgrey', label='Fst (mutation)')
    ax[5].set_title("Fst (coords), Fst (label), Fst (mutation)")
    ax[5].set_ylabel('Scores')
    ax[5].grid(True)
    ax[5].set_xticks(range(len(numeric_labels)))
    ax[5].set_xticklabels(x_labels)
    ax[5].legend()  # Add legend to differentiate the lines

    plt.tight_layout(pad=4.0)
    if args.save_png:
        plt.savefig("plot_statistics.png", format="png")
    else:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('statistics_file', type=str, help='CSV file containing summary statistics.')
    parser.add_argument('-s','--save_png', action="store_true", help='Flag to save the plot as an PNG file.')    
    args = parser.parse_args()

    plot_summary_statistics(args.statistics_file)
