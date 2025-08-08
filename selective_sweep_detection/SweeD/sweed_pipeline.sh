#!/bin/bash

## Default parameters (Used if not provided from the terminal)
DEFAULT_SAMPLE_SIZE=10
DEFAULT_SWEED_GRID=1000
DEFAULT_GENOME_LENGTH=10000
DEFAULT_MONOMORPHIC="TRUE"
DEFAULT_PLOT="FALSE"

## Check required argument: PATH
if [ -z "$1" ]; then
    echo "Error: Path argument is required."
    echo "Usage: $0 <path_to_genomes> [sample_size] [sweed_grid] [genome_length] [monomorphic] [plot]"
    exit 1
fi

## Assign parameters (Use provided values or defaults)
path="$1"
sample_size="${2:-$DEFAULT_SAMPLE_SIZE}"     # Use $2 if provided, otherwise use DEFAULT_SAMPLE_SIZE
sweed_grid="${3:-$DEFAULT_SWEED_GRID}"       # Use $3 if provided, otherwise use DEFAULT_SWEED_GRID
genome_length="${4:-$DEFAULT_GENOME_LENGTH}" # Use $4 if provided, otherwise use DEFAULT_GENOME_LENGTH
monomorphic="${5:-$DEFAULT_MONOMORPHIC}"     # Use $5 if provided, otherwise use DEFAULT_MONOMORPHIC
plot="${6:-$DEFAULT_PLOT}"                   # Use $6 if provided, otherwise use DEFAULT_PLOT

## Check if path exists
if [ ! -d "$path" ]; then
    echo "Error: Directory '$path' does not exist."
    exit 1
fi

## Validate parameters
if ! [[ "$sample_size" =~ ^[0-9]+$ ]]; then
    echo "Error: sample_size must be a positive integer."
    exit 1
fi

if ! [[ "$sweed_grid" =~ ^[0-9]+$ ]]; then
    echo "Error: sweed_grid must be a positive integer."
    exit 1
fi

if ! [[ "$genome_length" =~ ^[0-9]+$ ]]; then
    echo "Error: genome_length must be a positive integer."
    exit 1
fi

if [[ "$monomorphic" != "TRUE" && "$monomorphic" != "FALSE" ]]; then
    echo "Error: monomorphic must be either TRUE or FALSE."
    exit 1
fi

if [[ "$plot" != "TRUE" && "$plot" != "FALSE" ]]; then
    echo "Error: plot must be either TRUE or FALSE."
    exit 1
fi

## Suffix list
suffix_list=()
for file in "$path"/genomes/genomes_*.csv; do
    [[ -f "$file" ]] || continue  # Skip if no files match
    [[ $file =~ genomes_([0-9]+)\.csv ]] && suffix_list+=("${BASH_REMATCH[1]}")
done

## Create ouput directory
output_dir="${path}/sweed_output"
mkdir -p "$output_dir"

# Loop through each suffix
for suffix in "${suffix_list[@]}"; do
    echo "Processing suffix: ${suffix}"

    # Generate output filename and SweeD options based on monomorphic setting
    output_filename="all_fr_${suffix}$( [ "$monomorphic" = "TRUE" ] && echo '_mono' ).SF"
    sweed_name="${suffix}$( [ "$monomorphic" = "TRUE" ] && echo '_mono' )"
    sweed_mono_flag=$( [ "$monomorphic" = "TRUE" ] && echo '-monomorphic' )

    # Create the SweepFinder format (from a CSV file containing binary genomes)
    python3 format.py -csv "${path}/genomes/genomes_${suffix}.csv" -sample ${sample_size} -o ${output_filename}

    # Run SweeD with appropriate options
    SweeD -input ${output_filename} -name ${sweed_name} -grid ${sweed_grid} -length ${genome_length} ${sweed_mono_flag}

    # Conditionally Plot likelihood
    if [ "$plot" = "TRUE" ]; then
        Rscript sweed_plot.R -f "SweeD_Report.${sweed_name}" -s ${sweed_name} -sample ${sample_size} -event ${suffix}
    fi
done

mv all_fr_*.SF SweeD_Report.* SweeD_Info.* ${output_dir}/.
if [ "$plot" = "TRUE" ]; then
    mv plotSelSw_*.png ${output_dir}/.
fi