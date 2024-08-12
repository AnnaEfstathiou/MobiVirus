import argparse
import matplotlib.pyplot as plt
import numpy as np

def middle_coords(lower_bound, higher_bound):
    # Calculate the midpoint
    x = (lower_bound + higher_bound) / 2
    y = (lower_bound + higher_bound) / 2
    return x, y

def distance(lower_bound, higher_bound, dist, extra_dist=None):
    x, y = middle_coords(lower_bound, higher_bound)

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Create the first circle
    circle1 = plt.Circle((x, y), dist, edgecolor='blue', facecolor='none')
    # Add the first circle to the plot
    ax.add_patch(circle1)

    # Plot the center of the circle
    ax.scatter(x, y, color='red', zorder=5)  # zorder is used to ensure the dot is on top

    # Check if the second radius is provided
    if extra_dist is not None:
        # Create the second circle
        circle2 = plt.Circle((x, y), extra_dist, edgecolor='green', facecolor='none')
        # Add the second circle to the plot
        ax.add_patch(circle2)

    # Set the limits of the plot to the bounds
    ax.set_xlim(lower_bound, higher_bound)
    ax.set_ylim(lower_bound, higher_bound)

    # Add grid, labels, and title
    ax.grid(True)
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Circle Visualization with Center Dot at Midpoint')

    # Show the plot
    plt.show()
    plt.close()

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description='Visualize a circle with a given center and radius.')
    parser.add_argument('-l', '--bound_l', type=float, required=True, help='Lower bound for the spatial axis')
    parser.add_argument('-u', '--bound_h', type=float, required=True, help='Upper bound for the spatial axis')
    parser.add_argument('-d', '--inf_dist', type=float, required=True, help='Radius of the circle')
    parser.add_argument('-ed', '--extra_inf_dist', type=float, help='Radius of the second circle (optional)')

    args = parser.parse_args()

    # Call the distance function with the parsed arguments
    distance(args.bound_l, args.bound_h, args.inf_dist, args.extra_inf_dist)

if __name__ == '__main__':
    main()
