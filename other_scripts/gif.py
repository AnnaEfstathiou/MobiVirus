import argparse
import os
import imageio.v2 as imageio

def create_gif(directory, frame_duration):
    images = []
    names = []
    
    # Check if the directory exists
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    # Get all .png files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.png'):
            names.append(os.path.join(directory, filename))
    
    # Check if any .png files were found
    if not names:
        print(f"No .png files found in directory '{directory}'.")
        return
    
    # Sort files first numerically and then by length
    try:
        names.sort(key=lambda x: (int(os.path.splitext(os.path.basename(x))[0]), len(x)))
    except ValueError:
        print("Error: Non-numeric filenames encountered. Ensure filenames are numeric for proper sorting.")
        return

    # Read images and append to the list
    for g in names:
        try:
            images.append(imageio.imread(g))
        except Exception as e:
            print(f"Error reading image {g}: {e}")

    # Check if images list is empty
    if not images:
        print("Error: No valid images to create a GIF.")
        return

    # Create output path for the GIF
    output_path = os.path.join(directory, 'plots.gif')
    
    # Save the GIF with the calculated duration per file
    try:
        imageio.mimsave(output_path, images, duration=frame_duration)
        print(f"GIF saved at {output_path}")
    except Exception as e:
        print(f"Error creating GIF: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a GIF from .png files in a directory.')
    parser.add_argument('-dir', '--info_dir', type=str, required=True, help='Path to a directory of .png files')
    parser.add_argument('-d', '--duration', type=float, required=True, help='Duration per frame in seconds')
    args = parser.parse_args()
    
    create_gif(args.info_dir, args.duration)
