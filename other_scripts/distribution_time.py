import pandas as pd
import matplotlib.pyplot as plt
import argparse

# # Define your rates
# rt_i = 0.5  # Example rate value for rt_i
# rt_m = 0.3  # Example rate value for rt_m

# # Calculate the scale parameter
# scale = 1 / (rt_i + rt_m)

# # Generate random samples from an exponential distribution
# samples = np.random.exponential(scale, size=1000)  # Adjust the size for more or fewer samples

# # Plot the histogram
# plt.hist(samples, bins=30, alpha=0.75, color='skyblue', edgecolor='black')

# # Add labels and title
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.title('Histogram of Exponential Distribution')

# # Show the plot
# plt.show()


def distribution_of_time(csv_file):
    
    ## Plotting the information from an instance of the simulation ##
    # all_inf csv file

    # Read the CSV file
    time_data = pd.read_csv(csv_file)

    # Extract columns
    time = time_data['Time']

    plt.hist(time, bins=30, alpha=0.75, color='skyblue', edgecolor='black')

    # # Add labels and title
    plt.xlabel('Time Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of Time')

    if args.save_png:
        plt.savefig("spreaders_plot.png", format="png")
    else:
        plt.show()
    plt.close()
  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze events over a specified time period from a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to a CSV file containing time-stamped events data.')
    parser.add_argument('-s','--save_png', action="store_true", help='Flag to save the plot as an PNG file.')
    args = parser.parse_args()

    distribution_of_time(args.csv_file)