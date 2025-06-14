# MobiVirus Simulator

![MobiVirus Logo](MobiVirus_Logo.png)

## Table of contents

- [Abstract](#Abstract)
- [Requirements](#Requirements)
- [Main code](#Main-code)
  - [All about the strains](#All-about-the-strains)
  - [Output](#Output)
- [Examples](#Examples)

## Abstract

The MobiVirus Simulator simulates the spread of a virus within a population of individuals moving in a two-dimensional (2D) space. The purpose of the simulator is to study the evolution of viral genomes, diversity patterns, and adaptations of viruses as they evolve. The simulation resembles an epidemiological model, but the focus is solely on the evolution of the viral genome rather than its interaction with the human organism.

## Requirements

The MobiVirus Simulator works with the use of 3 files:

1. `MobiVirus_Simulator.py`: Main code where the simulation runs.
2. `mobifunctions.py`: Python script containing the functions used in the simulation.
3. `parameters.ini`: INI file containing the initial conditions of the simulation.

To the simulation, all files must be in the same directory. The chosen directory must, also, be added in the INI file.

All the required packages are mentioned in the `requirements.txt` file.

## Main code

The simulation starts with each individual located at a single position in a two-dimensional space and with a predefined number of infected individuals in the population carrying the viral genome. Initial sequences are either generated by msprime or all-0 genomes that differentiate slowly. At each simulated iteration, one of two events can occur <mark>movement</mark> or <mark>infection</mark>. If movement is chosen, it involves a single individual whose x and y coordinates in space change. If infection is selected then it is decided who will be the infector and who (can be more than one individual) will they infect by transferring their viral genome (infection is based on the distance between individuals). Also, sequences can optionally recombine. Alongside the events, infected individuals can recover (depending on recovery time). The simulation stops when everyone in the population is healthy (unless special conditions are laid down).
Note that the viral genome is binary. Positions with a value of 1.0 represent positions where a mutation has happened (SNPs).

### All about the strains... 

In default mode, the simulator creates only one strain. If '--manual_genomes' or '--msprime_genomes' flags is used then a second strain is created, with a different rate of infection, infection distance, probability to infect others and recovery time.

### Output

All the results are stored in one directory.
Path: same directory as the simulator’s code. 
Directory name: *simulation_timestamp* (e.g.*simulation_12_08_2024_21_33*)
1. *command_log.txt*: a text file containing the command, parsers, and parameters used for the simulation.
2. sub-directories:
    - *genome* directory:
        - genomes of all the individuals (healthy (NaN) and infected (binary format))
    - *samples* directory:
        - *all_inf*: information about viral strains
            This file is produced for every x number of events (x is determined in the INI file).
        - *coords*: coordinates and info for every individual
            This file is produced for every x number of events (x is determined in the INI file).
        - *event_type*:  info about the events
        - *infections*: info about who infected who
        - *recovery*: info about the recovery of the infected individuals
All the files are in a CSV format.

## Examples

A simulation that starts with genomes created by msprime and genetic recombinations.
```
python3 MobiVirus_Simulator.py -msprime -r 
```

A simulation that starts with all-0 genomes and genetic recombinations.
```
python3 MobiVirus_Simulator.py -manual -r 
```

A simulation that starts with genomes created by msprime and stops at a certain number of events.
```
python3 MobiVirus_Simulator.py -msprime -events 1000
```

A simulation that starts with genomes created by msprime and stops when super spreaders reach a certain percentage of the population.
```
python3 MobiVirus_Simulator.py -msprime -per_ss 0.95
```