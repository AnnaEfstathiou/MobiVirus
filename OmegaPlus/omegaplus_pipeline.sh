# !/bin/bash

# Base path
base_path="test"
results_dir="omega_plus"

# Loop through each directory that starts with "simulation_"
for sim_dir in "$base_path"/simulation_*/; do
    
    genomes_dir="${sim_dir}genomes/"
    # fasta_dir="${sim_dir}fasta_file/"
    
    # Check if genomes directory exists
    if [ -d "$genomes_dir" ]; then

        mkdir -p $results_dir
        
        # Loop through each file in the genomes directory
        for csv_file in "$genomes_dir"*; do

            filename=$(basename "$csv_file")        # Extract just the filename (e.g., genomes_100.csv)
            name_without_ext="${filename%%.*}"      # Remove everything after the first dot
            echo "${name_without_ext}.fa"                # Output: genomes_suffix (e.g. genomes_100)
            number="${name_without_ext##*_}"        # Extract suffix
            
            # Run the python script on the file
            python3 csv_to_ms.py -csv "$csv_file" -s 10
            rm "${name_without_ext}.fa"
            OmegaPlus -name $number -input "${name_without_ext}.ms" -minwin 50 -maxwin 10000 -grid 1000 -length 10000 -minsnps 5
            rm "${name_without_ext}.ms"
            # echo $genomes_dir
            # rm "${genomes_dir}genomes_*.fa"
            # rm "${name_without_ext}.fa"

            # if [ -f "${name_without_ext}.fa" ]; then
            #     rm "${name_without_ext}.fa"
            # fi
            
        done

        rm OmegaPlus_Warnings.*
        mv OmegaPlus_Info.* OmegaPlus_Report.* ${results_dir}/.
        # mv "$genomes_dir"*.fa "$fasta_dir" 
        mv ${results_dir} ${sim_dir}/.
    else
        echo "No genomes directory found in $sim_dir"
    fi
done
