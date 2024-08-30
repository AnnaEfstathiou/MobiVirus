""" 
Calculating: 
- Tajima's D 
- Pi-Estimator score
- Watterson-Estimator score
- number of unique sequences
- haplotype diversity 
- Fst 
"""

from collections import Counter
from typing import Dict
import numpy as np
from itertools import combinations
from preprocessing import __check, __check_populations


"""
------------------
INTERNAL FUNCTIONS
------------------
""" 

# def __polymorphic_sites(sequences: Dict[str, str]) -> int:
    
#     """ Counts the number of positions which show differences (polymorphisms). """
#     """ Used in the calculation of tajimas'd and watterson estimator. """
    
#     sequences_list = list(sequences.values())
#     seg_sites = 0
#     for i in range(len(sequences_list[0])):
#         s = sequences_list[0][i]  # reference character
#         if any(seq[i] != s for seq in sequences_list):
#             seg_sites += 1
#     return seg_sites

def __polymorphic_sites(sequences: Dict[str, str]) -> int:
   
    """ Counts the number of positions which show differences (polymorphisms). """
    """ Used in the calculation of tajimas'd and watterson estimator. """

    sequences_list = list(sequences.values())
    seg_sites = sum(1 for i in range(len(sequences_list[0])) if len(set(seq[i] for seq in sequences_list)) > 1)
    return seg_sites

def __harmonic(n: int) -> float:

    """ Computes the n-1th harmonic number. """
    """ Used in the calculation of tajimas'd and watterson estimator. """

    return sum(1 / i for i in range(1, n))


"""
----------
STATISTICS
----------
""" 

def tajimas_d(sequences: Dict[str, str]) -> float:

    """ Computes Tajima's D. """

    
    __check(sequences)
    

    seg_sites = __polymorphic_sites(sequences)
    if seg_sites == 0:
        return 0

    # theta_pi = pi_estimator(sequences, safe=False)
    theta_pi = pi_estimator(sequences)
    num_seq = len(sequences)
    harmonic = __harmonic(num_seq)
    a2 = sum(1 / (i**2) for i in range(1, num_seq))

    b1 = (num_seq + 1) / (3 * (num_seq - 1))
    b2 = (2 * (num_seq**2 + num_seq + 3)) / (9 * num_seq * (num_seq - 1))

    c1 = b1 - 1 / harmonic
    c2 = b2 - ((num_seq + 2) / (harmonic * num_seq)) + (a2 / (harmonic**2))

    if c1 == 0 or c2 == 0:
        return float('nan')  # informing that calculation is not reliable

    e1 = c1 / harmonic
    e2 = c2 / (harmonic**2 + a2)

    delta_Theta = theta_pi - (seg_sites / harmonic)
    tD = delta_Theta / ((e1 * seg_sites + e2 * seg_sites * (seg_sites - 1)) ** 0.5)
    return float(tD)

# def pi_estimator(sequences: Dict[str, str]) -> float:

#     """ Computes Pi estimatorn (Θπ). """

#     __check(sequences)

#     sequences_list = list(sequences.values())
#     pairwise_combinations = combinations(sequences_list, 2)
#     cs = [sum(not char1 == char2 for char1, char2 in zip(seq1, seq2)) for seq1, seq2 in pairwise_combinations]
#     n = len(sequences_list)
#     binomial = ((n - 1) * n) / 2  

#     return sum(cs) / binomial

def pi_estimator(sequences: Dict[str, str]) -> float:
    
    """ Computes Pi estimator (Θπ). """
    
    __check(sequences)
    sequences_list = list(sequences.values())
    n = len(sequences_list)
    
    # Initialize the sum of pairwise differences
    pairwise_sum = 0
    
    # Iterate through each pair of sequences without using combinations
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences_list[i]
            seq2 = sequences_list[j]
            # Count differences between the two sequences
            pairwise_sum += sum(char1 != char2 for char1, char2 in zip(seq1, seq2))
    
    # Calculate the binomial coefficient
    binomial = n * (n - 1) / 2
    
    return pairwise_sum / binomial


# def pi_estimator(sequences: Dict[str, str]) -> float:
    
#     """ Computes Pi estimator (Θπ). """
    
#     __check(sequences)
    
#     sequences_list = np.array(list(sequences.values()))
#     n = len(sequences_list)
#     pairwise_sum = np.sum(np.array([np.sum(np.array(list(seq1)) != np.array(list(seq2))) for seq1, seq2 in combinations(sequences_list, 2)]))
#     binomial = n * (n - 1) / 2
    
#     return pairwise_sum / binomial


def watterson_estimator(sequences: Dict[str, str]) -> float:

    """ Computes Watterson estimator (Θw). """
    
    __check(sequences)


    seg_sites = __polymorphic_sites(sequences)
    harmonic = __harmonic(len(sequences))

    return seg_sites / harmonic


def count_haplotypes(sequences):

    """ Calculates the number of unique sequences. """

    return len(set(sequences.values()))


def calculate_haplotype_diversity(sequences):

    """ Calculates the diversity of haplotypes in a given set of sequences. """
    """ Formula: 1 - Σ(pi^2), pi is the frequency of the i-th haplotype divided by the total number of haplotypes. """

    total_haplotypes = len(sequences)
    haplotype_frequencies = Counter(sequences.values()) # count the frequencies of each unique haplotype 
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # calculates the haplotype diversity based on the frequencies of each haplotype.
    return haplotype_diversity

def Fst(population_1, population_2):

    try:
        __check_populations(population_1, population_2)
    except ValueError as e:
        # Return NaN if populations do not meet the criteria
        return float('nan')

    def allele_frequencies(population):

        """Calculate allele frequencies for a population."""

        allele_counts = np.sum(population, axis=0) # Count SNPs at each position (column)
        total_alleles = len(population)            # Each individual contributes one allele per locus (haplotypes)
        return allele_counts / total_alleles

    def expected_heterozygosity(allele_freq):

        """Calculate expected heterozygosity for a population."""

        return 1 - np.sum(allele_freq ** 2)

    def observed_heterozygosity(population):

        """Calculate observed heterozygosity for a population."""

        heterozygous_individuals = np.any(population != population[0], axis=1)
        return np.mean(heterozygous_individuals)
    
    # Convert populations to numpy arrays
    population_1_array = np.array([list(map(int, seq)) for seq in population_1])
    population_2_array = np.array([list(map(int, seq)) for seq in population_2])

    # Calculate allele frequencies for each population
    allele_freq_1 = allele_frequencies(population_1_array)
    allele_freq_2 = allele_frequencies(population_2_array)

    # Calculate expected heterozygosity for each population
    He_1 = expected_heterozygosity(allele_freq_1)
    He_2 = expected_heterozygosity(allele_freq_2)

    # Calculate observed heterozygosity for each population
    Ho_1 = observed_heterozygosity(population_1_array)
    Ho_2 = observed_heterozygosity(population_2_array)

    # Calculate Fst
    fst = (He_1 + He_2 - (Ho_1 + Ho_2)) / (He_1 + He_2)

    return fst
