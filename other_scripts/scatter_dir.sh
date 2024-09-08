#!/bin/bash

"""mko"""

# Check if a directory argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Set the directory path from the argument
DIR="$1"

# Check if the directory exists
if [ ! -d "$DIR/samples" ]; then
    echo "Directory '$DIR/samples' does not exist"
    exit 1
fi

# Loop through every file that matches coords_*.csv pattern
for file in "$DIR/samples/coords_"*.csv; do
    # Check if the file exists (for the case of no match)
    if [ -e "$file" ]; then
        # Execute the Python script with the required arguments
        python3 metadata_plots.py -csv "$file" -inf "$DIR/samples/event_type.csv" -coords -s
    else
        echo "No matching files found in '$DIR/samples/'"
        exit 1
    fi
done
