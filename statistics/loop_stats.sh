#!/bin/bash

# Base path to the genomes directory
base_path="/home/anna/mobivirus/simulation/simulation_08_10_2024_13_52/genomes"

# List of values
values=(1000 25000 50000)

# Loop over each value
for i in "${values[@]}"; do
    genome_file="${base_path}/genomes_${i}.csv"
    echo "Processing the file: ${genome_file}"
    
    # Run the command with the current value
    python3 statistics.py -g "${genome_file}" -sample 100 -ld -save -mix
    python3 statistics.py -g "${genome_file}" -sample 100 -ld -save -ss
    python3 statistics.py -g "${genome_file}" -sample 100 -selsw -winsize 0.05 -stepsize 0.002 -save -mix
    python3 statistics.py -g "${genome_file}" -sample 100 -selsw -winsize 0.05 -stepsize 0.002 -save -ss
done
