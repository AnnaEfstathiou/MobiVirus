import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import imageio.v2 as imageio

def plot_spreaders(csv_file):
    
    ## Plotting the information from an instance of the simulation ##
    # all_inf csv file

    # Read the CSV file
    spreaders_data = pd.read_csv(csv_file)

    # Extract columns
    total_infected = spreaders_data['Total infected']
    super_spreaders = spreaders_data['Super spreaders']
    normal_spreaders = spreaders_data['Normal spreaders']
    time = spreaders_data['Time']
    
    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(time, total_infected, label='Total Infected')
    plt.plot(time, super_spreaders, label='Super Spreaders')
    plt.plot(time, normal_spreaders, label='Normal Spreaders')

    # Add labels and title
    plt.xlabel('Simulation Time')
    plt.ylabel('Number of Spreaders')
    plt.title('Simulation Results Over Time')
    plt.legend()

    if args.save_png:
        plt.savefig("spreaders_plot.png", format="png")
    else:
        plt.show()
    plt.close()


def plot_coordinates(csv_file):
    
    ## Create a scatter plot with the coordinates and mutation label for all individuals. ##
    # coords csv file

    coords_data = pd.read_csv(csv_file)

    plt.figure(figsize=(10, 6))

    # Assigning labels for every mutation value 
    mutation_labels = {
    0.0: 'Healthy individuals',
    1.0: 'Normal spreaders',
    2.0: 'Super spreaders'
    }

    # Assigning colors to different mutation values
    colors = {0.0: 'SpringGreen', 1.0: 'DarkOrange', 2.0: 'SlateBlue'}

    # Plotting each point with the corresponding color and legend
    for mutation_value in coords_data['Mutation'].unique():
        subset = coords_data[coords_data['Mutation'] == mutation_value]
        plt.scatter(subset['x'], subset['y'], color=colors[mutation_value], label=f'{mutation_labels[mutation_value]} = {len(subset)}')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Scatter Plot of x vs y with Different Mutations')
    plt.legend()
    if args.save_png:
        plt.savefig("scatter_plot.png", format="png")
    else:
        plt.show()
    plt.close()

def plot_time_events(csv_file):
    
    ## Plotting simulation time according to time of events. ##
    # all_inf csv file

    # Read the CSV file
    time_related_data = pd.read_csv(csv_file)

    # Extract columns
    time = time_related_data['Time']
    events = time_related_data['Events']
    
    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(events, time)

    # Add labels and title
    plt.xlabel('Events')
    plt.ylabel('Simulation Time')
    plt.title('Events over simulation time')
    plt.legend()

    if args.save_png:
        plt.savefig("time_over_events.png", format="png")
    else:
        plt.show()
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot simulation results over time from a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, help='Path to a CSV file. For the scatter plot use the coords.csv file. For the plots of the number of normal - super spreaders and time-events use an all_inf.csv file.')
    parser.add_argument('-spreaders', '--ns_ss_spreaders', action="store_true", help='Create a plot showing the number of normal and super spreaders.')
    parser.add_argument('-coords', '--coordinates_scatter_plot', action="store_true", help='Create a scatter plot with the coordinates and mutation label for all individuals.')
    parser.add_argument('-time', '--time_over_events_plot', action="store_true", help='Create a plot showing the relation between events and simulation time.')
    parser.add_argument('-s','--save_png', action="store_true", help='Flag to save the plot as an PNG file.')
    args = parser.parse_args()

    if args.ns_ss_spreaders:
        plot_spreaders(args.csv_file)
    elif args.coordinates_scatter_plot:
        plot_coordinates(args.csv_file)
    elif args.time_over_events_plot:
        plot_time_events(args.csv_file)
