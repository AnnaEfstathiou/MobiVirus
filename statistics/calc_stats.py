import argparse
import random
from preprocessing import process_csv, __filtering_sequences, __pop_coords, __pop_mutation_label, __pop_infection_label, validate_files
from statistics import tajimas_d, pi_estimator, watterson_estimator, count_haplotypes, calculate_haplotype_diversity, Fst

"""
Main function for calculating genomic statistics from a given CSV file.

Calculates:
- Tajima's D 
- Pi-Estimator
- Watterson-Estimator
- Unique sequences count
- Haplotype diversity 
- Fst based on population, infection label, and mutation label
"""

"""
-------------
MAIN FUNCTION
-------------
""" 

def calculate_statistics(input_file, coords_file, sample_size = None):

    """
    Calculates various statistics for a given genome sequence dataset.
    
    Args: input_file (str): Path to the input CSV file with binary sequences.
          coords_file (str): Path to the coordinates file for population information.
          sample_size (int, optional): Number of sequences to sample. If None, use all sequences.
    Returns: dict: A dictionary with calculated statistics.
    """
    
    infected_sequences = process_csv(input_file) # Process the input CSV file and extract binary sequences in a dictionary form

    '''Sampling'''
    # If sample_size is specified and valid, sample the sequences
    if sample_size and sample_size < len(infected_sequences):
        sampled_keys = random.sample(list(infected_sequences.keys()), sample_size)
        sampled_infected_sequences = {key: infected_sequences[key] for key in sampled_keys}
    else:
        sampled_infected_sequences = infected_sequences  # Use all sequences if sample size is invalid or not specified

    '''Populations' division for Fst'''
    # Filter and split the population into two groups for Fst calculation
    filtered_rows = __filtering_sequences(sampled_infected_sequences, coords_file)
    pop_coords_1, pop_coords_2 = __pop_coords(sampled_infected_sequences, filtered_rows)
    pop_label_1, pop_label_2 = __pop_infection_label(sampled_infected_sequences, filtered_rows)
    pop_mut_1, pop_mut_2 = __pop_mutation_label(sampled_infected_sequences, filtered_rows)

    # Calculate the desired statistics
    statistics = {
        "tajimas_d_score": tajimas_d(sampled_infected_sequences),
        "pi_estimator_score": pi_estimator(sampled_infected_sequences),
        "watterson_estimator_score": watterson_estimator(sampled_infected_sequences),
        "unique_count": count_haplotypes(sampled_infected_sequences),
        "haplotype_diversity": calculate_haplotype_diversity(sampled_infected_sequences),
        "total_sequences": len(infected_sequences),
        "Fst_coords": Fst(pop_coords_1, pop_coords_2),
        "Fst_inf_label": Fst(pop_label_1, pop_label_2),
        "Fst_mut_label": Fst(pop_mut_1, pop_mut_2),
    }

    return statistics

def print_statistics(results, sample_size=None):
    """
    Prints the calculated statistics in a formatted manner.
    
    Args: results (dict): The dictionary containing calculated statistics.
          sample_size (int, optional): The sample size used. Default is None.
    """
    
    total_sequences = results['total_sequences']
    actual_sample_size = sample_size or total_sequences
    
    print(f"Sample size: {actual_sample_size}")
    print(f"Tajima's D score: {results['tajimas_d_score']}")
    print(f"Pi-Estimator score: {results['pi_estimator_score']}")
    print(f"Watterson-Estimator score: {results['watterson_estimator_score']}")
    print(f"Number of unique sequences: {results['unique_count']}/{actual_sample_size} ({results['unique_count']/actual_sample_size:.3f})")
    print(f"Haplotype diversity: {results['haplotype_diversity']}")
    print(f"Fst (2 populations divided by coordinates): {results['Fst_coords']}")
    print(f"Fst (2 populations divided by mutation label): {results['Fst_mut_label']}")
    print(f"Fst (2 populations divided by infection label): {results['Fst_inf_label']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate genomic statistics from a CSV file containing binary sequences.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input file (CSV).')
    parser.add_argument('-coords', '--coords_file', type=str, required=True, help='The path to the coordinates file.')
    parser.add_argument('-s', '--sample_size', type=int, default=None, help='Number of sequences to sample. Defaults to all.')
    parser.add_argument('-p', '--population_size', action="store_true", help='Show the total population size.')
    args = parser.parse_args()
    
    validate_files(args.genome_file, args.coords_file) # Validate the provided input files

    results = calculate_statistics(args.genome_file, args.coords_file, args.sample_size) # Calculate statistics based on input arguments

    # Optionally, display population size
    if args.population_size:
        print(f"Population size of infected individuals: {results['total_sequences']}")
    
    print_statistics(results, args.sample_size) # Print the statistics based on the calculated results