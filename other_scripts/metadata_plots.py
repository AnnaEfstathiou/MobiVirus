import argparse
import pandas as pd
import matplotlib.pyplot as plt
import imageio.v2 as imageio


def plot_simulation_results(csv_file):
    
    ## Plotting the information from an instance of the simulation ##
    # all_inf csv file

    # Read the CSV file
    df = pd.read_csv(csv_file)

    # Extract columns
    total_infected = df['Total infected']
    super_spreaders = df['Super spreaders']
    normal_spreaders = df['Normal spreaders']
    time = df['Time']
    

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

    # Show the plot
    plt.show()


def create_gif(directory):
    images = []
    names = []
    for filename in os.listdir(directory):
        if int(filename[:-4])==0 or int(filename[:-4])%100==0:
            names.append(os.path.join(directory, filename))
        
    names.sort()
    names.sort(key=len)
    
    for g in names:
        images.append(imageio.imread(g, format='jpg'))
    
    imageio.mimsave(directory+'figs.gif', images, duration=0.8)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot simulation results over time from a CSV file.')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file')
    args = parser.parse_args()

    plot_simulation_results(args.csv_file)
