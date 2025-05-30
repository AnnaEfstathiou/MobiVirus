import os
import pandas as pd
import numpy as np
from configparser import ConfigParser
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess


from sweedReports_functions import process_SweeDReport, normalize_events, bar_chart_plot, distance_scatterplot, distance_histogram, distance_lineplot

#%% Parameters initialization (INI file)
"""
----------------------------------------------
PARSE THE .INI FILE & PARAMETER INITIALIZATION
----------------------------------------------
""" 

## Specify the directory and file name that contains the parameters ##
directory = './'
file_name = 'sweedReports_parameters.ini'
file_path = os.path.join(directory, file_name)

## Check if the file exists ##
if not os.path.exists(file_path):
    print(f"Error: {file_name} must be in the current directory.")
    sys.exit("Exiting the simulation.")

## File exists, proceed with parsing ##
config = ConfigParser()
config.read(file_path)

# Directory containing simulations (e.g. /home/directory)
directory = config.get('Parameters', 'directory').strip('"')
# Subdirectory containing the SweeD Reports (e.g. /home/directory).
reports_dir = config.get('Parameters', 'reports_dir').strip('"')
# Infection Rate of the super strains
infection_rate = config.getint('Parameters', 'infection_rate') 
# Total events of the simulations
total_events = config.getint('Parameters', 'total_events') 
# Filter the dataset to exclude all entries with lower Likelihood than the threshold
# If filter is <=0 then the dataset is not filtered
threshold = config.getfloat('Parameters', 'threshold')
# Beneficial positions
position = config.getint('Parameters', 'position') 
# Boolean parameter: True --> save pandas DataFrame into a CSV file
save_csv = config.getboolean('Parameters', 'save_csv')

# Boolean parameter: True --> create the plot
plot_bar_chart = config.getboolean('Plot Type', 'plot_bar_chart')
plot_distance_scatterplot = config.getboolean('Plot Type', 'plot_distance_scatterplot')
plot_histogram = config.getboolean('Plot Type', 'plot_histogram')
plot_lineplot = config.getboolean('Plot Type', 'plot_lineplot')

#%% Main
if __name__ == "__main__":
    
    """
    --------------------------
    Process SweeD Report files
    --------------------------
    """ 
    all_results = []
    for dir_name in os.listdir(directory):
        dir_path = os.path.join(directory, dir_name)
        if os.path.isdir(dir_path) and dir_name.startswith("simulation_"):
            result_dict = process_SweeDReport(infection_rate, position, dir_path, reports_dir)
            if result_dict:
                all_results.extend(result_dict.values()) 

    # Create DataFrame given a list
    data = pd.DataFrame(
    all_results, 
    columns=['Rate', 'Simulation', 'Event', 'Likelihood', 'Position', 'Distance'])
    data = data.sort_values(by=['Simulation', 'Event']).reset_index(drop=True) # Sort by simulation name and then by event number 

    # Filter the dataset to exclude all entries with lower Likelihood than the threshold
    if threshold > 0:
        data = data[data["Likelihood"] >= threshold]

    # Save or print the DataFrame
    if save_csv:
        csv_name = config.get('Parameters', 'csv_name').strip('"') # name of the CSV file
        data.to_csv(csv_name, index=False)
        print(f"Results saved to {csv_name}.")

    """
    --------
    Plotting
    --------
    """ 

    """Bar Chart"""

    if plot_bar_chart:

        # Parse Parameters
        figure_size = tuple(map(int, config.get('Bar Chart', 'figure_size').split(','))) # Parse figure size
        bins = list(map(float, config.get('Bar Chart', 'bins').split(',')))
        labels = config.get('Bar Chart', 'labels').split(',')                            # Parse labels
        colors = config.get('Bar Chart', 'colors').split(',')                            # Parse colors
        bar_chart_title = config.get('Bar Chart', 'bar_chart_title').strip('"')          # Parse plot title
        x_axis_label = config.get('Bar Chart', 'x_axis_label').strip('"')
        y_axis_label = config.get('Bar Chart', 'y_axis_label').strip('"')
        bar_chart_name = config.get('Bar Chart', 'bar_chart_name').strip('"')            # Parse name of the png
        
        # Extract, from every simulation, the fixation event (event where super strains > 90% population & no normal strains)
        fixation_events = data[data['Event'] % 5000 != 0] 

        bar_chart_plot(fixation_events, figure_size, bins, labels, colors, bar_chart_title, x_axis_label, y_axis_label, bar_chart_name)


    """Scatter Plot"""

    if plot_distance_scatterplot:

        # Parse Parameters
        figure_size = tuple(map(int, config.get('Scatter Plot', 'figure_size').split(','))) # Parse figure size
        alpha = config.getfloat('Scatter Plot', 'alpha') 
        frac = config.getfloat('Scatter Plot', 'frac')
        linewidth = config.getint('Scatter Plot', 'linewidth') 
        color_scatter = config.get('Scatter Plot', 'color_scatter')
        color_axvline = config.get('Scatter Plot', 'color_axvline')
        color_loess = config.get('Scatter Plot', 'color_loess')            
        scatter_plot_title = config.get('Scatter Plot', 'scatter_plot_title').strip('"')    # Parse plot title
        x_axis_label = config.get('Scatter Plot', 'x_axis_label').strip('"')
        y_axis_label = config.get('Scatter Plot', 'y_axis_label').strip('"')
        scatter_plot_name = config.get('Scatter Plot', 'scatter_plot_name').strip('"')      # Parse name of the png

        # Aligns event values in each simulation based on a unique fixation event.
        normalized_events_data = normalize_events(data)

        distance_scatterplot(normalized_events_data, figure_size, alpha, frac, linewidth, color_scatter, color_loess, color_axvline, scatter_plot_title, x_axis_label, y_axis_label, scatter_plot_name)


    """Histogram"""

    if plot_histogram:

        # Parse Parameters
        figure_size = tuple(map(int, config.get('Histogram', 'figure_size').split(','))) # Parse figure size  
        bins_range = tuple(map(int, config.get('Histogram', 'bins_range').split(',')))       
        color_bins = config.get('Histogram', 'color_bins')
        histogram_title = config.get('Histogram', 'histogram_title').strip('"')    # Parse plot title
        x_axis_label = config.get('Histogram', 'x_axis_label').strip('"')
        y_axis_label = config.get('Histogram', 'y_axis_label').strip('"')
        histogram_name = config.get('Histogram', 'histogram_name').strip('"')

        distance_histogram(data, figure_size, bins_range, color_bins, histogram_title, x_axis_label, y_axis_label, histogram_name)


    if plot_lineplot:

        figure_size = tuple(map(int, config.get('Line Plot', 'figure_size').split(','))) # Parse figure size  
        lineplot_title = config.get('Line Plot', 'lineplot_title').strip('"')            # Parse plot title
        x_axis_label = config.get('Line Plot', 'x_axis_label').strip('"')
        y_axis_label = config.get('Line Plot', 'y_axis_label').strip('"')
        lineplot_name = config.get('Line Plot', 'lineplot_name').strip('"')

        # Aligns event values in each simulation based on a unique fixation event.
        normalized_events_data = normalize_events(data)

        distance_lineplot(normalized_events_data, figure_size, lineplot_title, x_axis_label, y_axis_label, lineplot_name)