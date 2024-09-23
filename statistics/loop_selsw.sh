#!/bin/bash

# List of values
values=(5000)

# Loop over each value
for i in "${values[@]}"; do
    echo "Processing file: genomes_${i}.csv"
    echo "Running 'selective_sweep.py' script for a sample of the viral strains."
    python3 selective_sweep.py -g "/home/anna/mobivirus/files/simulation_20_09_2024_13_30/genomes_${i}.csv" -w 500 -s 20 -sample 100 -save "selsw_${i}.png"
    echo "Running 'selective_sweep.py' script for all super strains."
    python3 selective_sweep.py -g "/home/anna/mobivirus/files/simulation_20_09_2024_13_30/genomes_${i}.csv" -w 500 -s 20 -ss_sample -save "selsw_ss_${i}.png"
done
