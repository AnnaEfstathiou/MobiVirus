# MobiVirus Simulator

## Table of contents

- [Abstract](## Abstract)
- Requirements
- Detailed Description
    - Main code 
    - Optional arguments
    - Output
- Statistical calculations

## Abstract

## Requirements

## Detailed Description
### Main code 
### Optional arguments
### Output

## Statistical calculations

## Description of the operating mode
The MobiVirus Simulator works using three files:
1. Main code where the simulation runs: MobiVirus_Simulator.py
2. Python file containing the functions used in the simulation: mobifunctions.py
3. Ini File containing the initial conditions of the simulation: parameters.ini

All files should be in the same directory in order to run the simulation. Your chosen directory must also be added in the "directory=" parameter both in Movement_inbox.py and in myfunctions.py

## Main Code partitions
The main code is run through three partitions. The partitioning was done for better control of the initialization before running the actual simulation. In order to avoid errors, it is best to run the partitions separately, but they can also be run as one. 

In the first partition, we state the main directory of the simulation and we also create three more directories where the output data will be saved. The output data and their corresponding directories are some samples of the simulation, plots to visualize the positions and the health status of the individuals, and the genome table. In addition, the output also contains some tables with extra data that might be useful for extra statistics. These tables have the time of infections and un-infections that happen, the number of total infections and total movements, the initial and final coordinate tables, and an index of who infected who every infection time. 

The samples of the simulation are taken every few generations. The sample consists of:
* The genome table, that shows the genome of every infected individual in the sample. The table has rows equal to the number of individuals and columns equal to the positions of the genome. Healthy individuals have nan values in the row that corresponds to their genome.
* A table containing information on infections. This table keeps the numbers of infected individuals in total, together with the number of individuals with each mutation every time (in simulation time) that an infection happens.
* A table containing information on the position of the individuals in space. This table contains the positions, together with labels on the health status of each individual, the mutation of the virus that they carry (if they are infected), and their susceptibility. 

The second partition produces the initial group of individuals and it initializes the important parameters for the simulation. Most parameters are imported in the code from the parameters.ini file. Most of the information on the individuals is in the coordinates table that is created in this step. The genome table is also created in this step. 

The last partition contains the loop where the actual simulation will run. 

## Simulation
The first loop dictates how long the simulation will run. As a default, we have initiated this loop to keep running while there are infected individuals. We also add the condition to stop the simulation if individuals carrying the Normal strain 

The collected data at the end of the simulation will be saved in the directory that is indicated in the directories stated in the beginning of the main code.
