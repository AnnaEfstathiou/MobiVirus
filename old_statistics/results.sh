#!/bin/bash

# Load the config file
source results.cfg

# # Create the results subdirectory inside the directory
# mkdir -p "${directory}/results_${strain_type}"

# Assign the subdirectory path based on the directory
samples_dir="$directory/samples"
genomes_dir="$directory/genomes"

# Strains Plot
# Rscript /home/anna/mobivirus/statistics/plot_strains.R -csv "${samples_dir}/all_inf_${final_event}.csv" 

# Summary statistics, Fst, Haplotype Diversity, Number of unique haplotypes
python3 /home/anna/mobivirus/statistics/sim_sumstats.py -d ${directory} -sample ${sample_size} -sample_type ${sampling_tech} -${strain_type}
# Plotting
Rscript plot_sumstats.R -csv "stats_${datetime}.csv" -s ${sample_size} -l ${genome_length} -pop ${strain_type} -png TRUE -pdf TRUE

# for event in "${simulation_events[@]}"; do

#     genome_file="${genomes_dir}/genomes_${event}.csv"
#     echo "Processing the file: ${genome_file}"

#     # Deternine names
#     ld_filename="ld_${event}.csv"
#     selsw_filename="selsw_${event}.csv"
#     sfs_filename="sfs_${event}.csv"
    
#     # LD, Selective Sweep, Site Frequency Spectrum
#     python3 statistics.py -g "${genome_file}" -sample ${sample_size} -${strain_type} -ld 
#     python3 statistics.py -g "${genome_file}" -sample ${sample_size} -${strain_type} -selsw -winsize ${window_size} -stepsize ${step_size}
#     python3 statistics.py -g "${genome_file}" -sample ${sample_size} -${strain_type} -sfs

#     # Plotting
#     Rscript plot_ld.R -csv ld_${event}.csv -l ${genome_length} -png TRUE
#     Rscript plot_selsw.R -csv selsw_${event}.csv -s ${sample_size} -l ${genome_length} -pop ${population}
#     Rscript plot_sfs.R -csv sfs_${event}.csv 
# done

# Move all .png and .csv files from the current directory to the results subdirectory
mv *.png *.csv *.pdf "${directory}/results_${strain_type}"
# cp results.cfg "${directory}/results_${strain_type}"