# import argparse
# import os
# import imageio.v2 as imageio

# def create_gif(directory):
#     images = []
#     names = []
    
#     # Get all .png files in the directory
#     for filename in os.listdir(directory):
#         if filename.endswith('.png'):
#             names.append(os.path.join(directory, filename))
    
#     # Sort files first numerically and then by length
#     names.sort(key=lambda x: (int(os.path.splitext(os.path.basename(x))[0]), len(x)))
    
#     # Read images and append to the list
#     for g in names:
#         images.append(imageio.imread(g))
    
#     # Create output path for the GIF
#     output_path = os.path.join(directory, 'plots.gif')
    
#     # Save the GIF
#     imageio.mimsave(output_path, images, duration=0.4)
#     print(f"GIF saved at {output_path}")


import argparse
import os
import imageio.v2 as imageio

def create_gif(directory, frame_duration):
    images = []
    names = []
    
    # Get all .png files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.png'):
            names.append(os.path.join(directory, filename))
    
    # Sort files first numerically and then by length
    names.sort(key=lambda x: (int(os.path.splitext(os.path.basename(x))[0]), len(x)))
    
    # Read images and append to the list
    for g in names:
        images.append(imageio.imread(g))
    
    # Create output path for the GIF
    output_path = os.path.join(directory, 'plots.gif')
    
    # Save the GIF with the calculated duration per file
    imageio.mimsave(output_path, images, duration=frame_duration)
    print(f"GIF saved at {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a GIF from .png files in a directory.')
    parser.add_argument('-dir', '--info_dir', type=str, required=True, help='Path to a directory of .png files')
    parser.add_argument('-d', '--duration', type=float, required=True, help='Duration per frame in seconds')
    args = parser.parse_args()
    
    create_gif(args.info_dir, args.duration)
