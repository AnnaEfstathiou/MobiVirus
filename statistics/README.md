# Summary Statistics

## Abstract

The following scripts calculate different summary statistics for the population covered by the "MobiVirus Simulator". The statistics calculated are: Tajima's D, Pi-Estimator score, Watterson-Estimator score, number of unique sequences, haplotype diversity & Fst. 

## Requirements

Scripts for internal use:
- `preprocessing.py`:  Python script that does the necessary preparation of the files (CSV or FASTA) to calculate the statistics.
- `statistics.py`: Python scripts that contains all the functions that calculate the statistics.

Main scripts:
- `calc_stats.py`: Python script that calculates the summary statistics for the population at a point in time. (scripts required: `preprocessing.py`, `statistics.py`)
- `stat_table.py`: Python script that calculates the summary statistics for the population for all time points. (scripts required: `preprocessing.py`, `statistics.py`, `calc_stats.py`)

In order to calculate the summary statistics, all files must be in the same directory. 

