"""
Run as:
python3 total_time_vs_rates.py -directories simulation_09_01_2025_14_01 simulation_09_01_2025_14_02 simulation_09_01_2025_14_03 -suffix 5000
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def process_directories(directories, suffix):
    results = []

    for directory in directories:
        samples_path = os.path.join(directory, "samples")

        if not os.path.isdir(samples_path):
            print(f"Directory not found: {samples_path}")
            continue

        all_inf_file = os.path.join(samples_path, f"all_inf_{suffix}.csv")
        coords_file = os.path.join(samples_path, f"coords_{suffix}.csv")

        try:
            if os.path.isfile(all_inf_file):
                all_inf_df = pd.read_csv(all_inf_file)
                time = all_inf_df['Time'].iloc[-1]
            else:
                print(f"File not found: {all_inf_file}")
                continue

            if os.path.isfile(coords_file):
                coords_df = pd.read_csv(coords_file)
                infection_rate = coords_df['Rate of infection'].iloc[0]
            else:
                print(f"File not found: {coords_file}")
                continue

            results.append((directory, time, infection_rate))

        except Exception as e:
            print(f"Error processing files in {samples_path}: {e}")

    return results

def generate_table_and_plot(results):
    
    # Create a plot
    times = [result[1] for result in results]
    infection_rates = [result[2] for result in results]

    plt.figure(figsize=(8, 6))
    plt.plot(times, infection_rates, marker='o', color='indigo', linestyle='-')
    plt.yticks(infection_rates)  # Ensure all rates are displayed on the y-axis
    # for i, result in enumerate(results):
    #     plt.annotate(result[0], (times[i], infection_rates[i]))

    plt.title("Total Simulation Time vs Infection Rate")
    plt.xlabel("Time")
    plt.ylabel("Rate of Infection")
    plt.grid(True)
    plt.savefig("total_time_VS_inf_rates.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Process directories and plot data.")
    parser.add_argument('-directories', nargs='+', required=True, help="List of directories to process")
    parser.add_argument('-suffix', required=True, type=str, help="Suffix for the files")

    args = parser.parse_args()

    directories = args.directories
    suffix = args.suffix

    results = process_directories(directories, suffix)

    if results:
        generate_table_and_plot(results)
    else:
        print("No valid data found to process.")

if __name__ == "__main__":
    main()
