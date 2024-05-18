import matplotlib as plt   
import pandas as pd
import argparse


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process multiple genome CSV or FASTA files to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('statistics_csv', type=str, help='CSV file containing summary statistics.')
    args = parser.parse_args()

    plot_summary_statistics(args.statistics_csv)