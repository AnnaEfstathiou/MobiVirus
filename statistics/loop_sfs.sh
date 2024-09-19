#!/bin/bash

# List of values
values=(0 1000 2000 3000 4000 5000 6000 7000 8000 9000 final)

# Loop over each value
for i in "${values[@]}"; do
    echo "genomes_${i}.csv"
    # Run the command with the current value
    python3 sfs.py -g "/home/anna/mobivirus/files/simulation_16_09_2024_19_15/genomes/genomes_${i}.csv" -s "sfs_${i}.png"
done
