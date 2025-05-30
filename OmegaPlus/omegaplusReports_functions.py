import io
import pandas as pd
import numpy as np
import os
import subprocess
import re
import csv
import logging
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

#%% Process SweeD Report files 

def process_OmegaReport(infection_rate, position, directory, reports_dir):
    '''
    Processes OmegaPlus_Report files in a specified directory, identifying the genomic position with the highest 
    omega statistic. Computes the distance between this peak and a reference position, returning results
    as a dictionary.

    Parameters:
    - infection_rate (float): Infection rate associated with the dataset.
    - position (int): Reference genomic position.
    - directory (str): Root directory containing simulation results.
    - reports_dir (str): Subdirectory within the root directory containing OmegaPlus_Report files.

    Returns:
    - dict: A dictionary keyed by filenames, with values as lists containing:
        [infection_rate, directory, event_number, max_omega, max_position, distance_to_reference]
    '''

    full_reports_dir = os.path.join(directory, reports_dir)

    if not os.path.isdir(full_reports_dir):
        logging.error(f"Reports directory '{full_reports_dir}' does not exist.")
        return {}

    # Collect files starting with "OmegaPlus_Report"
    report_files = [f for f in os.listdir(full_reports_dir) if f.startswith("OmegaPlus_Report")]

    if not report_files:
        logging.warning(f"No OmegaPlus_Report files found in '{full_reports_dir}'.")
        return {}

    results_dict = {}

    for file in report_files:
        try:
            event_match = re.search(r'\.(\d+)$', file)
            
            if not event_match:
                logging.warning(f"Skipping '{file}': filename does not match expected pattern.")
                continue

            event_number = int(event_match.group(1))
            file_path = os.path.join(full_reports_dir, file)

            with open(file_path, 'r') as f:
                lines = f.readlines()

            # Remove metadata line if present
            data_lines = [line for line in lines if not line.strip().startswith("//")]

            # Read data into DataFrame
            data = pd.read_csv(
                io.StringIO("".join(data_lines)),  
                sep=r'\s+',
                header=None,
                names=["Position", "Omega"],
                engine='python'
            )

            if 'Omega' not in data.columns:
                logging.warning(f"Skipping '{file}': 'Omega' column missing.")
                continue

            # Out of positions with the same ω value, choose the one closest to the beneficial mutation
            max_omega = data['Omega'].max()
            max_rows = data[data['Omega'] == max_omega]
            # Compute distances to reference position
            max_rows['Distance'] = (max_rows['Position'] - position).abs()
            # Pick the row with the smallest distance
            closest_row = max_rows.loc[max_rows['Distance'].idxmin()]
            max_position = round(closest_row['Position'])
            distance = int(closest_row['Distance'])

            results_dict[file] = [
                infection_rate, directory, event_number,
                max_omega, max_position, distance
            ]

        except Exception as e:
            logging.error(f"Error processing '{file}': {e}")
            continue

    return results_dict

#%% Extra functions

def num_event(csv_file, num_normal, num_super):
    '''
    Identifies the first event in a CSV dataset where the number of super spreaders 
    reaches or exceeds a specified threshold, given a fixed number of normal spreaders.

    Parameters:
        csv_file (str): Path to the CSV file containing event data.
        num_normal (int): The required number of normal spreaders.
        num_super (int): The threshold number of super spreaders to look for.

    Returns:
        int or None: The event number where the condition is first met, or None if not found.
    '''

    df = pd.read_csv(csv_file) # Load the CSV file into a DataFrame
    df_filtered = df[df['Normal spreaders'] == num_normal]
    first_event = df_filtered[df_filtered['Super spreaders'] >= num_super].head(1)
    if not first_event.empty:
        event = int(first_event.iloc[0]['Event'])
    else:
        event = None
    
    return event


def normalize_events(data):
    '''
    Aligns event values in each simulation based on a unique fixation event.

    For each unique simulation in the input DataFrame, this function identifies a single
    fixation event (defined as an event where Event % 5000 != 0). It then normalizes all 
    event values in that simulation by subtracting the fixation event value, effectively 
    aligning the fixation event to zero. The normalized values are stored in a new column 
    called 'AlignedEvent'.

    Parameters:
        data (pd.DataFrame): A DataFrame containing 'Simulation' and 'Event' columns.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'AlignedEvent' column.
    '''
        
    # Create a new column in the DataFrame called "AlignedEvent"
    # Initialize it with NaN; this will store the normalized event values
    data["AlignedEvent"] = np.nan

    # Iterate over each unique Simulation ID in the DataFrame
    for sim in data["Simulation"].unique():

        # Filter the DataFrame to get only the rows for the current simulation
        sub_data = data[data["Simulation"] == sim]
            
        # Identify the "fixation event" within this simulation:
        # It's defined as an event where the event value is **not** a multiple of 5000
        fixation_event = sub_data[sub_data["Event"] % 5000 != 0]

        # If there is not exactly one fixation event, print a warning and skip this simulation
        if len(fixation_event) != 1:
            print(f"Warning: Simulation {sim} has {len(fixation_event)} fixation events.")
            continue
            
        # Get the numeric value of the fixation event
        fixation_event_value = fixation_event["Event"].values[0]
            
        # Normalize the Event values by subtracting the fixation event's value
        # This effectively aligns the fixation event to zero
        data.loc[data["Simulation"] == sim, "AlignedEvent"] = data.loc[data["Simulation"] == sim, "Event"] - fixation_event_value

    # Return the modified DataFrame with the new "AlignedEvent" column
    return data

#%% Plotting

def bar_chart_plot(data, figure_size, bins, labels, colors, title, x_label, y_label, png_name):
    '''
    Creates and saves a bar chart showing the distribution of distances from the 
    beneficial mutation, based on the position with the highest likelihood at the fixation event.

    The fixation event is defined as the moment when super strains constitute more than 90% 
    of the population and no normal strains remain. This function bins the 'Distance' values 
    into categories, counts how many values fall into each bin, and visualizes this distribution 
    in a bar chart.

    Parameters:
        data (pd.DataFrame): A DataFrame that must contain a 'Distance' column.
        figure_size (tuple): Size of the figure, e.g., (8, 6).
        bins (list): Bin edges used to group distances.
        labels (list): Labels for each bin (should match the number of bins - 1).
        colors (list): List of colors to use for each bar.
        title (str): Plot title.
        x_label (str): Label of the x axis
        y_label (str): Label of the y axis
        png_name (str): File name to save the plot (e.g., "bar_chart.png").

    Returns:
        None: The bar chart is saved to the specified file.
    '''

    # Bin the 'Distance' column into categorical ranges using the provided bins and labels
    data['Distance Bin'] = pd.cut(
        data['Distance'], 
        bins=bins, 
        labels=labels, 
        right=False  # Include left edge, exclude right edge
    )

    # Count occurrences in each bin and sort by bin order
    bin_counts = data['Distance Bin'].value_counts().sort_index()

    # Create the bar chart
    plt.figure(figsize=figure_size)
    plt.bar(
        bin_counts.index.astype(str),  # X-axis: bin labels as strings
        bin_counts.values,             # Y-axis: counts
        color=colors[:len(bin_counts)] # Limit colors to number of bars
    )

    # Add labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()

    # Save the figure
    # plt.savefig(png_name)
    plt.show()


def distance_scatterplot(data, figure_size, alpha, frac, linewidth, color_scatter, color_loess, color_axvline, title, x_label, y_label, png_name):
    '''
    Plots a scatter plot of Distance vs. Normalized Events with a LOESS-smoothed curve.

    This function visualizes the relationship between aligned events and distance by:
    - Plotting raw data as a scatter plot.
    - Applying LOESS (Locally Estimated Scatterplot Smoothing) to highlight trends.
    - Marking the fixation event (aligned at 0) with a vertical line.

    Parameters:
        data (pd.DataFrame): A DataFrame containing 'AlignedEvent' and 'Distance' columns.
        figure_size (tuple): Size of the figure, e.g., (8, 6).
        alpha (float): Alpha blending value, between 0 (transparent) and 1 (opaque)
        frac (float): The fraction of the data used when estimating each y-value in LOESS. 
                      Lower values make the curve more sensitive to local changes; 
                      higher values produce smoother curves.
        linewidth (int or float): Width of the LOESS smoothing line
        color_scatter (str): 
        color_axvline (str):
        color_loess (str):
        title (str): Plot title
        x_label (str): Label of the x axis
        y_label (str): Label of the y axis
        png_name (str): File name to save the plot (e.g., "scatter_plot.png").
        
    Returns:
        None: The plot is saved to the specified file.
    '''
    
    # Set up the figure size
    plt.figure(figsize=figure_size)

    # Scatter plot of raw data
    plt.scatter(data["AlignedEvent"], data["Distance"], alpha=alpha, color=color_scatter)

    # Apply LOESS smoothing to the data
    smoothed = lowess(data["Distance"], data["AlignedEvent"], frac=frac)

    # Plot the LOESS smoothed curve
    plt.plot(smoothed[:, 0], smoothed[:, 1], color=color_loess, linewidth=linewidth, label="LOESS Fit")

    # Mark the fixation event (aligned to zero)
    plt.axvline(0, linestyle='--', color=color_axvline, label="Fixation Event")

    # Label axes and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    # Add legend and improve layout
    plt.legend()
    plt.tight_layout()

    # Save the figure as a PNG file
    plt.savefig(png_name)


def distance_histogram(data, figure_size, bins_range, color_bins, title, x_label, y_label, png_name):
    '''
    Plots a histogram of Distance values.
    This function visualizes the distribution of the 'Distance' variable using a histogram.

    Parameters:
        data (pd.DataFrame): A DataFrame containing a 'Distance' column.
        figure_size (tuple): Size of the figure, e.g., (8, 6).
        bins_range (tuple): A tuple specifying the start, stop, and step for bin edges, e.g., (0, 100, 5).
        color_bins (str): Color used for the histogram bars.
        title (str): Title of the plot.
        x_label (str): Label for the x-axis.
        y_label (str): Label for the y-axis.
        png_name (str): File name to save the plot (e.g., "histogram.png").

    Returns:
        None: The plot is saved to the specified file.
    '''
    
    # Plot the histogram
    plt.figure(figsize=figure_size)
    plt.hist(data["Distance"], bins=range(*bins_range), color=color_bins, edgecolor="black")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    # Save the figure
    plt.savefig(png_name)


def distance_lineplot(data, figure_size, title, x_label, y_label, png_name):
    '''
    Plots line graphs of Distance vs AlignedEvent for each simulation group.
    This function generates a line plot where each line represents a different simulation, showing how 'Distance' changes over 'Event'.

    Parameters:
        data (pd.DataFrame): A DataFrame containing at least 'Simulation', 'AlignedEvent', and 'Distance' columns.
        figure_size (tuple): Size of the figure, e.g., (12, 6).
        title (str): Title of the plot.
        x_label (str): Label for the x-axis.
        y_label (str): Label for the y-axis.
        png_name (str): File name to save the plot (e.g., "lineplot.png").

    Returns:
        None: The plot is displayed but not saved unless savefig is uncommented.
    '''

    # Group the data by 'Simulation'
    grouped = data.groupby('Simulation')
    
    # Initialize the figure with specified size
    plt.figure(figsize=figure_size)

    # Plot a line for each simulation group
    for name, group in grouped:
        plt.plot(group['AlignedEvent'], group['Distance'], label=name)

    # Label axes and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    # plt.legend(title='Simulation', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    # plt.tight_layout()
    # plt.grid(True)

    # Save the figure
    plt.savefig(png_name)