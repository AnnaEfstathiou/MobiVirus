""" 
INPUT: FASTA or CSV file with binary sequences
Calculating: 
- Tajima's D 
- Pi-Estimator score
- Watterson-Estimator score
- number of unique sequences
- haplotype diversity 
"""

import argparse
import os
import random
from preprocessing import __csv_to_fasta, __process_fasta, __read_fasta
from statistics import tajimas_d, pi_estimator, watterson_estimator, count_haplotypes, calculate_haplotype_diversity


"""
-------------
MAIN FUNCTION
-------------
""" 

def calculate_statistics(input_file, sample_size = None):
    
    file_extension = os.path.splitext(input_file)[1]  # get the file extension of the input file

    if file_extension.lower() == '.csv':  # check if the input file is a CSV file
        new_fasta = input_file.rsplit('.', 1)[0] + '_processed.fa'  # generate the processed FASTA file path by replacing the extension
        __csv_to_fasta(input_file, new_fasta)  # convert CSV to FASTA (exclude nan sequences)
        infected_sequences = __read_fasta(new_fasta)  # create a dictionary with sequences 
    elif file_extension.lower() == '.fasta' or file_extension.lower() == '.fa':  # check if the input file is a FASTA file
        new_fasta = input_file.rsplit('.', 1)[0] + '_processed.fa'  # generate the processed FASTA file path by replacing the extension
        __process_fasta(input_file, new_fasta)  # exclude nan sequences from FASTA
        infected_sequences = __read_fasta(new_fasta)  # create a dictionary with sequences 
    else:
        raise ValueError('The input file must be either a CSV or a FASTA file.')

    # Sampling sequences if a sample size is specified
    if sample_size and sample_size < len(infected_sequences):
        # print(f"Sample size: {sample_size}")
        # print(f"Population size: {len(sequences)}\n")
        sampled_keys = random.sample(list(infected_sequences.keys()), sample_size)
        sampled_infected_sequences = {key: infected_sequences[key] for key in sampled_keys}
    else:
        sampled_infected_sequences = infected_sequences  # Use all sequences if no valid sample size specified


    # Calculate statistics
    tajimas_d_score = tajimas_d(sampled_infected_sequences)
    pi_estimator_score = pi_estimator(sampled_infected_sequences)
    watterson_estimator_score = watterson_estimator(sampled_infected_sequences)
    unique_count = count_haplotypes(sampled_infected_sequences)
    haplotype_diversity = calculate_haplotype_diversity(sampled_infected_sequences)

    os.remove(new_fasta)  # delete the generated FASTA file

    # return tajimas_d_score, pi_estimator_score, watterson_estimator_score, unique_count, haplotype_diversity, len(sequences)
    return {
        "tajimas_d_score": tajimas_d_score,
        "pi_estimator_score": pi_estimator_score,
        "watterson_estimator_score": watterson_estimator_score,
        "unique_count": unique_count,
        "haplotype_diversity": haplotype_diversity,
        "total_sequences": len(infected_sequences) }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a genome CSV or FASTA file to compute the following statistics: Tajima's D score, Pi-Estimator score, Watterson-Estimator score, number of unique sequences and haplotype diversity.")
    parser.add_argument('input_file', type=str, help='The path to the input file (either CSV or FASTA).')
    parser.add_argument('-s','--sample_size', type=int, help='Number of sequences to sample for calculations. Use all if not specified or larger than available.', default=None)
    parser.add_argument('-p', '--population_size', action="store_true", help='See the population size, so it is easier to decide the sample size.')
    args = parser.parse_args()

    results = calculate_statistics(args.input_file, args.sample_size)

    if args.population_size:
        print(f"Population size of infected individuals: {results['total_sequences']}")

    if args.sample_size or not args.population_size:
        print(f"Sample size: {args.sample_size}")
        print(f"Tajima's D score: {results['tajimas_d_score']}")
        print(f"Pi-Estimator score: {results['pi_estimator_score']}")
        print(f"Watterson-Estimator score: {results['watterson_estimator_score']}")
        print(f"Number of unique sequences: {results['unique_count']}/{results['total_sequences']} ({results['unique_count']/results['total_sequences']:.3f})")
        print(f"Haplotype diversity: {results['haplotype_diversity']}")
    
