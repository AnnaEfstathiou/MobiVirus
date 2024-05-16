""" 
Calculating: 
- Tajima's D 
- Pi-Estimator score
- Watterson-Estimator score
- number of unique sequences
- haplotype diversity 
"""

from collections import Counter
from typing import Dict
from itertools import combinations
from preprocessing import __check


"""
------------------
INTERNAL FUNCTIONS
------------------
""" 

def __polymorphic_sites(sequences: Dict[str, str]) -> int:
    
    """ Counts the number of positions which show differences (polymorphisms). """
    """ Used in the calculation of tajimas'd and watterson estimator. """
    
    sequences_list = list(sequences.values())
    seg_sites = 0
    for i in range(len(sequences_list[0])):
        s = sequences_list[0][i]  # reference character
        if any(seq[i] != s for seq in sequences_list):
            seg_sites += 1
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

    theta_pi = pi_estimator(sequences, safe=False)
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

def pi_estimator(sequences: Dict[str, str], safe=True) -> float:

    """ Computes Pi estimatorn (Θπ). """

    if safe:
        __check(sequences)

    sequences_list = list(sequences.values())
    pairwise_combinations = combinations(sequences_list, 2)
    cs = [sum(not char1 == char2 for char1, char2 in zip(seq1, seq2)) for seq1, seq2 in pairwise_combinations]
    n = len(sequences_list)
    binomial = ((n - 1) * n) / 2  

    return sum(cs) / binomial


def watterson_estimator(sequences: Dict[str, str], safe=True) -> float:

    """ Computes Watterson estimator (Θw). """

    if safe:
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