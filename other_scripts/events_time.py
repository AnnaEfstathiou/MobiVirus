import argparse
import pandas as pd

'''This script checks how many events happend during the given period.'''

def events_per_time(csv_file, start_time, end_time):
    
    # Read the CSV file
    data = pd.read_csv(csv_file)
    
    # Ensure 'Time' column is string for pattern matching
    data['Time'] = data['Time'].astype(str)
    
    # Filter data for start and end times 
    filtered_start = data[data['Time'].str.startswith(start_time)]
    filtered_end = data[data['Time'].str.startswith(end_time)]

    # Handle cases where no matching events are found
    if filtered_start.empty:
        print(f"No events found starting at {start_time}.")
        return
    if filtered_end.empty:
        print(f"No events found starting at {end_time}.")
        return

    # Find the first event for both start and end time
    event_1 = filtered_start.iloc[0]
    event_2 = filtered_end.iloc[0]
    
    # Calculate how many events happened between the two times
    event_count = event_2['Event'] - event_1['Event']
    
    return print(f"\nBetween the period: {event_1['Time']}-{event_2['Time']}, {int(event_count)} events happened (events from {int(event_1['Event'])} to {int(event_2['Event'])}).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze events over a specified time period from a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to a CSV file containing time-stamped events data.')
    parser.add_argument('-s', '--start_time', type=str, required=True, help='Start time pattern to look for in the data.')
    parser.add_argument('-e', '--end_time', type=str, required=True, help='End time pattern to look for in the data.')
    args = parser.parse_args()

    events_per_time(args.csv_file, args.start_time, args.end_time)

