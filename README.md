# MobiVirus Simulator

## Table of contents

- [Abstract](#Abstract)
- [Requirements](#Requirements)
- [Detailed Description](#Detailed-Description)
  - [Main code](#Main-code)
  - [Optional arguments](#Optional-arguments)
  - [Output](#Output)
- [Statistical calculations](#Statistical-calculations)

## Abstract

The MobiVirus Simulator simulates a virus infecting individuals in a population as they move around in a two-dimensional space.

## Requirements

The MobiVirus Simulator works with the use of 3 files:

1. `MobiVirus_Simulator.py`: Main code where the simulation runs.
2. `mobifunctions.py`: Python script containing the functions used in the simulation.
3. `parameters.ini`: INI file containing the initial conditions of the simulation.

In order to run the simulation, all files must be in the same directory. The chosen directory must, also, be added in the INI file.

All the required packages are mentioned in the `requirements.txt` file.

## Detailed Description

### Main code

The simulation starts with each individual located at a single position in a two-dimensional space and with a predefined number of infected individuals in the population carrying the viral genome. At each simulated iteration one of two events can occur: movement or infection. If movement is chosen, then it involves a single individual whose x and y coordinates in space change. If infection is chosen then it is decided who will be the infecting individual and who (can be more than one individual) will they infect by tranfering their viral genome (infection is based on the distance between individuals).Alongside the events occur uninfections in which, depending on the recovery time, some infected individuals are selected to become healthy. The simulation stops when everyone in the population is healthy.
Note that the viral genome consists of 0 and 1. Potitions with value 1 represent potitions that a mutation has happened.

### Optional arguments

The code contains several optional arguments that can be set from the command line.

Action arguments:

- **super_strain**: Create a 2nd strain in the simulation that has a different infectivity rate than the normal strain. The different strain has a specific number of important genome positions in the beginning. Both the infectivity rate and the important genome positions are defines in the INI file.
  If the flag is not used then the n_i parameter (number of important genome positions) in the INI file should be equal to 0.
- **initial_genomes**: Save the initial viral genomes. Healthy individuals have a viral genome that consists of 0. The number of non-zero genomes is equal to the ii parameter (initial number of infected individuals) in the INI file.
- **scatter_plots**: Create a scatter plot to visualize the positions (coordinates) along with the health status of the individuals.
- **visualize_data**: Visualize the data table as a dataframe in the console. the table contains information on the x,y cordinates, the health status of the individuals, the rate of movement and infection, the mutation label (which differs in case there is a super strain) and the susceptibility, of the individuals, to the virus.

String arguments - optional conditions under which the simulation stops:

- **percentage_infected**: The simulation stops if the percentage of infected individuals in the population reach a certain number. This number must be a float between 0 and 1.
- **max_infections**: The simulation stops if the infection are more than a certain number. This number must be a positive integer.
- **percentage_susceptibility**: The simulation stops if the susceptible, to the virus, individuals are less than a certain number. This number must be a float between 0 and 1. This argument is valid only when the super strain is present in the simulation. This is because it creates an imbalance between the susceptibility of the individuals regarding to the strain they have been infected with. In more detain, if an individual is infected with the normal strain, then when they get healthy they are again susceptible to the virus. In contrast, the healthy individuals, who had previously the super strain are no longer susceptible to any viral strain.
- **ratio_super_vs_normal**: The simulation stops if the number of individuals carrying the normal strain are less than a certain percentage of the individuals with the super strain. This number must be a float between 0 and 1. This argument is valid only when the super strain is present in the simulation.

### Output

The output directory contains 2 subdirectories: the genomes folder and the samples folder.

- The **"genomes" directory** contains multiple genome tables in CSV format. Ever table shows the genome of every infected individual in the sample. It has rows equal to the number of individuals and columns equal to the positions of the genome. Healthy individuals have nan values in the row that corresponds to their genome.
- The **"samples" directory** contains tables with information on infections. These tables keep the numbers of infected individuals in total, together with the number of individuals with each mutation (if there are more than one type) every time (in simulation time) that an infection happens. In addition, it contains tables with information on the position of the individuals in space (coordinates). These include the positions, together with labels on the health status of each individual, the mutation of the virus that they carry (if they are infected and if there are more than one type of mutation), and their susceptibility.
- Optional Output: The **"plots" directory** contains scatter plots to visualize the positions and the health status of the individuals.
  The samples of the simulation are taken every few generations (defined in the INI file).
  The collected data at the end of the simulation will be saved in the directory that is indicated in the INI file.