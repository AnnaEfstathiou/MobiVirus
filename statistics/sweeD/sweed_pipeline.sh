#!/bin/bash

# File parameters
path="/home/anna/mobivirus/simulation/selsw/simulation_01_11_2024_13_00/genomes"
suffix=15000 # genome file + SwweD report name
sample_size=100
# SweeD parameters
sweed_grid=1000
genome_length=10000
# Set monomorphic variable (TRUE or FALSE)
monomorphic=FALSE


# Conditional execution based on monomorphic value
if [ "$monomorphic" = TRUE ]; then
    output_filename="all_fr_${suffix}_mono.SF"

    # Create the SweepFinder format (from a CSV file containing binary genomes)
    python3 format.py -csv "${path}/genomes_${suffix}.csv" -sample ${sample_size} -o ${output_filename}
    
    # Run SweeD with monomorphic option
    SweeD -input ${output_filename} -name ${suffix}_mono -grid ${sweed_grid} -length ${genome_length} -monomorphic
    
    # Plot likelihood for Selective Sweep detection 
    Rscript sweed_plot.R -f "SweeD_Report.${suffix}_mono" -s ${suffix}_mono
else
    output_filename="all_fr_${suffix}.SF"

    # Create the SweepFinder format (from a CSV file containing binary genomes)
    python3 format.py -csv "${path}/genomes_${suffix}.csv" -sample ${sample_size} -o ${output_filename}
    
    # Run SweeD without monomorphic option
    SweeD -input ${output_filename} -name ${suffix} -grid ${sweed_grid} -length ${genome_length}
    
    # Plot likelihood for Selective Sweep detection
    Rscript sweed_plot.R -f "SweeD_Report.${suffix}" -s ${suffix}
fi
