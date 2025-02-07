#!/bin/bash

# File parameters
path="/home/anna/mobivirus/simulation/selsw/r2/simulation_07_02_2025_11_29/genomes"
sample_size=10
# SweeD parameters
sweed_grid=1000
genome_length=10000
# Set monomorphic variable (TRUE or FALSE)
monomorphic=TRUE

# List of suffixes to process
suffix_list=(5500 5600 5700 5800 5900 6000)

# Loop through each suffix
for suffix in "${suffix_list[@]}"; do
    echo "Processing suffix: ${suffix}"

    # Generate output filename and SweeD options based on monomorphic setting
    output_filename="all_fr_${suffix}$( [ "$monomorphic" = TRUE ] && echo '_mono' ).SF"
    sweed_name="${suffix}$( [ "$monomorphic" = TRUE ] && echo '_mono' )"
    sweed_mono_flag=$( [ "$monomorphic" = TRUE ] && echo '-monomorphic' )

    # Create the SweepFinder format (from a CSV file containing binary genomes)
    python3 format.py -csv "${path}/genomes_${suffix}.csv" -sample ${sample_size} -o ${output_filename}

    # Run SweeD with appropriate options
    SweeD -input ${output_filename} -name ${sweed_name} -grid ${sweed_grid} -length ${genome_length} ${sweed_mono_flag}

    # Plot likelihood for Selective Sweep detection
    Rscript sweed_plot.R -f "SweeD_Report.${sweed_name}" -s ${sweed_name} -sample ${sample_size} -event ${suffix}
    
done
