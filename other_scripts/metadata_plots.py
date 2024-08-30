import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import imageio.v2 as imageio


"""
----------------------
MUTATION RELATED PLOTS
----------------------
"""

def plot_spreaders(csv_file):

    '''Plotting the number of individuals with each mutation over time.'''
    ## necessary argument: all_inf CSV file ##

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

    # Display or save the plot
    if args.save_png:
        plt.savefig("spreaders_plot.png", format="png")
    else:
        plt.show()
    plt.close()


"""
-------------------
SPACE RELATED PLOTS
-------------------
"""

def plot_coordinates(csv_file):
    
    '''Plotting (scatter plot) the coordinates and mutation label for all individuals.'''
    ## necessary argument: coords CSV file ##

    # Read the CSV file
    coords_data = pd.read_csv(csv_file)

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

    plt.figure(figsize=(10, 6))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Scatter Plot of xy coordinates')
    plt.legend()

    # Display or save the plot
    if args.save_png:
        plt.savefig("scatter_plot.png", format="png")
    else:
        plt.show()
    plt.close()


def plot_marginal_histograms(csv_file):

    '''Plotting the distribution of x and y coordinates (separately) for all individuals.'''
    ## necessary argument: coords CSV file ##
    
    # Read the data from the CSV file
    coords_data = pd.read_csv(csv_file)

    # Create figure and axes
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    # Histogram for x-coordinates
    ax[0].hist(coords_data['x'], bins=30, color='skyblue', edgecolor='black')
    ax[0].set_title('Histogram of X Coordinates')
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('Frequency')

    # Histogram for y-coordinates
    ax[1].hist(coords_data['y'], bins=30, color='salmon', edgecolor='black')
    ax[1].set_title('Histogram of Y Coordinates')
    ax[1].set_xlabel('y')
    ax[1].set_ylabel('Frequency')

    plt.tight_layout()

    # Display or save the plot
    if args.save_png:
        plt.savefig("marginal_histograms.png", format="png")
    else:
        plt.show()
    plt.close()


"""
------------------
TIME RELATED PLOTS
------------------
"""

def plot_time_events(csv_file):
    
    '''Plotting the number of events over the simulation time.'''
    ## necessary argument: all_inf CSV file ##

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

    # Display or save the plot
    if args.save_png:
        plt.savefig("time_over_events.png", format="png")
    else:
        plt.show()
    plt.close()


def distribution_of_time(csv_file):
    
    '''Plotting the distribution of simulation time (time vs frequency).'''    
    ## necessary argument: all_inf CSV file ##

    # Read the CSV file
    time_data = pd.read_csv(csv_file)

    # Extract columns
    time = time_data['Time']

    # Create a histogram
    plt.hist(time, bins=30, alpha=0.75, color='skyblue', edgecolor='black')

    # Add labels and title
    plt.xlabel('Time Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of Time')

    # Display or save the plot
    if args.save_png:
        plt.savefig("distribution_of_time.png", format="png")
    else:
        plt.show()
    plt.close()

"""
-------------------
EVENT RELATED PLOTS
-------------------
"""

def type_of_event(csv_file):

    '''Plotting the type of event (movement or infection).'''
    ## necessary argument: type_event CSV file ##
    
    # Read the CSV file
    data = pd.read_csv(csv_file)

    # Map 'Type of Event' to numeric values for plotting
    type_mapping = {'infection': 1, 'movement': 2}
    data['Type Value'] = data['Type of Event'].map(type_mapping)

    # Create a scatter plot of all data points
    plt.figure(figsize=(16, 6))  # Larger figure size for clarity
    scatter = plt.scatter(data['Event'], data['Type Value'], 
                          c=data['Type Value'], cmap='viridis', 
                          alpha=0.6, s=10)  # Use color mapping and transparency

    # Customize the plot
    plt.title('Type of Events (Infections VS Movements)', fontsize=14)
    plt.xlabel('Events', fontsize=12)
    plt.ylabel('Type of Event', fontsize=12)
    plt.yticks([1, 2], ['infection', 'movement'])
    plt.grid(True)

    # Display or save the plot
    if args.save_png:
        plt.savefig("type_of_event.png", format="png")
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot simulation results over time from a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, help='Path to a CSV file.')
    parser.add_argument('-spreaders', '--ns_ss_spreaders', action="store_true", help='Plotting the number of individuals with each mutation over time. (all_inf CSV)')
    parser.add_argument('-coords', '--coordinates_scatter_plot', action="store_true", help='Plotting (scatter plot) the coordinates and mutation label for all individuals. (coords CSV)')
    parser.add_argument('-hist', '--histogram', action="store_true", help='Plotting the distribution of x and y coordinates (separately) for all individuals. (coords CSV)')
    parser.add_argument('-time', '--time_over_events_plot', action="store_true", help='Plotting the number of events over the simulation time. (all_inf CSV)')
    parser.add_argument('-dist_time', '--distribution_of_time', action="store_true", help='Plotting the distribution of simulation time (time vs frequency). (all_inf CSV)')
    parser.add_argument('-type', '--type_of_event', action="store_true", help='Plotting the type of events. (type_event CSV)')
    parser.add_argument('-s','--save_png', action="store_true", help='Flag to save the plot as an PNG file.')
    args = parser.parse_args()

    if args.ns_ss_spreaders:
        plot_spreaders(args.csv_file)
    elif args.coordinates_scatter_plot:
        plot_coordinates(args.csv_file)
    elif args.histogram:
        plot_marginal_histograms(args.csv_file)
    elif args.time_over_events_plot:
        plot_time_events(args.csv_file)
    elif args.distribution_of_time:
        distribution_of_time(args.csv_file)
    elif args.type_of_event:
        type_of_event(args.csv_file)
    

