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

def __polymorphic_sites(sequences: Dict[str, str]) -> int:
   
    """ 
    Counts the number of positions which show differences (polymorphisms). 
    Used in the calculation of tajimas'd and watterson estimator.          
    
    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: int: The number of segregating (polymorphic) sites.
    """

    sequences_list = list(sequences.values()) # Extract the sequences from the dictionary
    # Count positions where there is > 1 unique nucleotide (i.e., polymorphism exists)
    seg_sites = sum(1 for i in range(len(sequences_list[0])) if len(set(seq[i] for seq in sequences_list)) > 1)
    return seg_sites

def __harmonic(n: int) -> float:

    """ 
    Computes the n-1th harmonic number. 
    Used in the calculation of tajimas'd and watterson estimator. 

    Args: n (int): The number of sequences.
    Returns: float: The (n-1)th harmonic number.
    """

    return sum(1 / i for i in range(1, n))


"""
----------
STATISTICS
----------
""" 

def tajimas_d(sequences: Dict[str, str]) -> float:

    """ 
    Computes Tajima's D. 
    
    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: float: The value of Tajima's D statistic.
    """

    __check(sequences) # Check the validity of the input sequences
    
    seg_sites = __polymorphic_sites(sequences) # Get the number of polymorphic sites
    if seg_sites == 0:
        return 0                               # If no polymorphic sites, Tajima's D is zero

    theta_pi = pi_estimator(sequences)              # Compute Pi estimator (Θπ)
    num_seq = len(sequences)                        # Number of sequences
    harmonic = __harmonic(num_seq)                  # Compute harmonic number
    a2 = sum(1 / (i**2) for i in range(1, num_seq)) # Second harmonic component
    # Coefficients for Tajima's D formula
    b1 = (num_seq + 1) / (3 * (num_seq - 1))
    b2 = (2 * (num_seq**2 + num_seq + 3)) / (9 * num_seq * (num_seq - 1))
    # Constants for variance calculations
    c1 = b1 - 1 / harmonic
    c2 = b2 - ((num_seq + 2) / (harmonic * num_seq)) + (a2 / (harmonic**2))

    if c1 == 0 or c2 == 0:
        return float('nan')  # Return NaN if calculations are unreliable

    e1 = c1 / harmonic           # Scaling constant for e1
    e2 = c2 / (harmonic**2 + a2) # Scaling constant for e2

    delta_Theta = theta_pi - (seg_sites / harmonic)                                 # Difference between Pi estimator and Watterson estimator
    tD = delta_Theta / ((e1 * seg_sites + e2 * seg_sites * (seg_sites - 1)) ** 0.5) # Compute Tajima's D using the formula
    
    return float(tD)

def pi_estimator(sequences: Dict[str, str]) -> float:
    
    """ 
    Computes Pi estimator (Θπ). 
    
    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: float: The Pi estimator value.
    """
    
    __check(sequences) # Check the validity of the input sequences
    
    sequences_list = list(sequences.values()) # Convert dictionary to list of sequences
    n = len(sequences_list)                   # Number of sequences 
    pairwise_sum = 0 # Initialize sum of pairwise differences
    
    # Iterate through all pairs of sequences and count nucleotide differences
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences_list[i]
            seq2 = sequences_list[j]
            # Count the number of positions where the two sequences differ
            pairwise_sum += sum(char1 != char2 for char1, char2 in zip(seq1, seq2))
    
    binomial = n * (n - 1) / 2     # Calculate the binomial coefficient
    
    return pairwise_sum / binomial # Return average pairwise difference (π)

def watterson_estimator(sequences: Dict[str, str]) -> float:

    """ 
    Computes Watterson estimator (Θw). 
    
    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: float: Watterson's estimator value.
    """
    
    __check(sequences) # Check the validity of the input sequences

    seg_sites = __polymorphic_sites(sequences) # Get the number of segregating sites
    harmonic = __harmonic(len(sequences))      # Compute the harmonic number

    return seg_sites / harmonic # Return Watterson's estimator (Θw)

def count_haplotypes(sequences):

    """ 
    Calculates the number of unique sequences. 
    
    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: int: The number of unique haplotypes.
    """

    return len(set(sequences.values())) # Return the number of unique sequences (haplotypes)

def calculate_haplotype_diversity(sequences):

    """ 
    Calculates the diversity of haplotypes in a given set of sequences. 
    Formula: 1 - Σ(pi^2), where pi is the frequency of the i-th haplotype divided by the total number of haplotypes. 

    Args: sequences (Dict[str, str]): A dictionary where keys are sequence indices and values are the sequences themselves.  
    Returns: float: Haplotype diversity.
    """

    total_haplotypes = len(sequences)                   # Total number of haplotypes (sequences)
    haplotype_frequencies = Counter(sequences.values()) # Count the frequency of each unique haplotype
    # Compute haplotype diversity using the formula
    haplotype_diversity = 1 - sum((freq / total_haplotypes) ** 2 for freq in haplotype_frequencies.values()) # calculates the haplotype diversity based on the frequencies of each haplotype.
    
    return haplotype_diversity

def Fst(population_1, population_2):

    """
    Computes Fst, a measure of population differentiation due to genetic structure.
    
    Args: population_1 (List[str]): A list of sequences from population 1.
          population_2 (List[str]): A list of sequences from population 2.
    Returns: float: Fst value between the two populations.
    """

    try:
        __check_populations(population_1, population_2) # Check the validity of the input populations
    except ValueError as e:
        return float('nan') # Return NaN if populations do not meet the criteria

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
    
    # Convert populations to numpy arrays of integers
    population_1_array = np.array([list(map(int, seq)) for seq in population_1])
    population_2_array = np.array([list(map(int, seq)) for seq in population_2])

    # Calculate allele frequencies for each population
    allele_freq_1 = allele_frequencies(population_1_array)
    allele_freq_2 = allele_frequencies(population_2_array)

    # Calculate expected and observed heterozygosity for each population
    He_1 = expected_heterozygosity(allele_freq_1)
    He_2 = expected_heterozygosity(allele_freq_2)
    Ho_1 = observed_heterozygosity(population_1_array)
    Ho_2 = observed_heterozygosity(population_2_array)

    if (He_1 + He_2) == 0:
        return float('nan')  # Return NaN if heterozygosity is zero

    # Compute Fst using the formula
    fst = (He_1 + He_2 - (Ho_1 + Ho_2)) / (He_1 + He_2)

    return fst
