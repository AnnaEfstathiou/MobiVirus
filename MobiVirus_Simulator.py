#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================
|      ------------------      |
|   | MobiVirus Simulator |    |
|      ------------------      |
================================
"""

#%% Import packages - functions

## Importing Packages ##
import argparse
import numpy as np
import pandas as pd
import time
import os
import sys
from datetime import datetime
from configparser import ConfigParser

initial=time.time()
func_time = time.time()-initial

"""
========================================
IMPORTING THE NECESSARY PYTHON FUNCTIONS
========================================
"""

from mobifunctions import log_command, coords_function, infection_label, msprime_genomes, ss_msprime_genomes, ss_mutation_position, mutation, rec_probi, recombine_genomes, movement, initial_distances, new_distances, ind_probi
from mobifunctions import sample_data, save_data

#%% Parameters initialization (INI file)
""" 
===================
PARSE THE .INI FILE
===================
""" 

## Specify the directory and file name that contains the parameters ##
directory = './'
file_name = 'parameters.ini'
file_path = os.path.join(directory, file_name)

## Check if the file exists ##
if not os.path.exists(file_path):
    print(f"Error: {file_name} must be in the current directory.")
    sys.exit("Exiting the simulation.")

## File exists, proceed with parsing ##
config = ConfigParser()
config.read(file_path)

""" 
========================
PARAMETER INITIALIZATION
========================
"""

"""
--------------------------------------------------
Read directory from the Initial_Parameters section
--------------------------------------------------
"""
directory = config.get('Directory', 'directory').strip('"')  # Directory where the scripts exist (.ini file & .py script with functions)

super_strain = config.getboolean('Type_of_strains', 'super_strain') # Boolean flag: True for two strains (super and normal), False for one (normal only)

n = config.getint('Initial_Parameters', 'n')                        # Number of individuals in the simulation
l = config.getint('Initial_Parameters', 'l')                        # Length of genome
bound_l = config.getfloat('Initial_Parameters', 'bound_l')          # Lower bound for the plot
bound_h = config.getfloat('Initial_Parameters', 'bound_h')          # Upper bound for the plot
ii = config.getint('Initial_Parameters', 'ii')                      # Number of infected individuals
r_m = config.getfloat('Initial_Parameters', 'r_m')                  # Mutation rate for each position in genome
n_i = config.getint('Initial_Parameters', 'n_i')                    # Important genome positions for the ss mutation; If there are no super strains, n_i = 0
region_pos = config.get('Initial_Parameters', 'region_pos')         # Area of the important positions in the viral genome that define the super strain
ri_n = config.getfloat('Initial_Parameters', 'ri_n')                # Rate of infection from normal strain 
ri_s = config.getfloat('Initial_Parameters', 'ri_s')                # Rate of infection from super strain 
rm_i = config.getfloat('Initial_Parameters', 'rm_i')                # Rate of movement for infected ind.
rm_h = config.getfloat('Initial_Parameters', 'rm_h')                # Rate of movement for healthy ind.
inf_dist_ns = config.getfloat('Initial_Parameters', 'inf_dist_ns')  # Infection distance (normal strain)
inf_dist_ss = config.getfloat('Initial_Parameters', 'inf_dist_ss')  # Infection distance (normal strain)
prob_inf_ns = config.getfloat('Initial_Parameters', 'prob_inf_ns')  # Probability that the infector infects an individual in their infection distance (normal strain)
prob_inf_ss = config.getfloat('Initial_Parameters', 'prob_inf_ss')  # Probability that the infector infects an individual in their infection distance (super strain)
r_rec = config.getfloat('Initial_Parameters', 'r_rec')              # Rate of recombination
rec_t_ns = config.getfloat('Initial_Parameters', 'rec_t_ns')        # Recovery time in simulation time (normal strain)
rec_t_ss = config.getfloat('Initial_Parameters', 'rec_t_ss')        # Recovery time in simulation time (super strain)
imm_t_ns = config.getfloat('Initial_Parameters', 'imm_t_ns')        # Immunity time in simulation time (normal strain)
imm_t_ss = config.getfloat('Initial_Parameters', 'imm_t_ss')        # Immunity time in simulation time (super strain)
rim_ns = config.getfloat('Initial_Parameters', 'rim_ns')            # Relative infection mobility (normal strain)
rim_ss = config.getfloat('Initial_Parameters', 'rim_ss')            # Relative infection mobility (super strain)
sample_times = config.getint('Initial_Parameters', 'sample_times')  # Generations where we take samples of our simulation's output
ss_form = config.getint('Super_strain_Parameters', 'ss_formation')  # Event that the super strain mutation is introduced for the 1st time
ss_events = config.getint('Super_strain_Parameters', 'ss_events')   # Period of events when the super strain is introduced (if previously there are no individuals with a super strain)

use_demography = config.getboolean('msprime_Parameters', 'use_demography') # Parameter defining whether to simulate with population demography (exponential growth) or not
initial_pop = config.getint('msprime_Parameters', 'initial_pop')           # Present-day population size
past_pop = config.getint('msprime_Parameters', 'past_pop')                 # Past population size (population t generation ago)
t_gen = config.getint('msprime_Parameters', 't_gen')                       # Time (in generations) at which population growth stops and becomes constant

#%% Parsed arguments
"""
===============
PARSE ARGUMENTS
===============
""" 

parser = argparse.ArgumentParser(description='MobiVirus Simulation Script')
## String parsers ##
parser.add_argument('-ratio', '--ratio_super_vs_normal', type=str, help='The simulation stops when the given ratio of super Strain individuals / normal Strain individuals becomes real.')
parser.add_argument('-per_inf', '--percentage_infected', type=str, help='The simulation stops when the infected individuals in the population reach the given percentage.')
parser.add_argument('-per_ss', '--percentage_super_strain', type=str, help='The simulation stops when the super spreaders in the population reach the given percentage.')
parser.add_argument('-per_ns', '--percentage_normal_strain', type=str, help='The simulation stops when the normal spreaders in the population reach the given percentage.')
parser.add_argument('-max_inf', '--max_infections', type=str, help='The simulation stops when the given number of infections happens.')
parser.add_argument('-max_mv', '--max_movements', type=str, help='The simulation stops when the given number of movements happens.')
parser.add_argument('-time', '--end_time', type=str, help='The simulation stops at the given simulation time.')
parser.add_argument('-events', '--end_events', type=str, help='The simulation stops when the given number of events (movements+infections) happpens.')
## Action parsers ##
# group = parser.add_mutually_exclusive_group(required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-manual', '--manual_genomes', action="store_true", help='Begin the simulation with all-0 genomes. Introduce the super strain mutation at the "ss_form" event.')
group.add_argument('-msprime', '--msprime_genomes', action="store_true", help='Begin the simulation with genomes created by msprime. Introduce the super strain mutation at the "ss_form" event.')
parser.add_argument('-all_inf', '--all_infected_once', action="store_true", help='The simulation stops when all the individuals are infected at least once.')
parser.add_argument('-r', '--recombination', action="store_true", help='Provide the ability to recombinate genomes, if 2 infections happen.')
parser.add_argument('-sample_inf', '--sample_infection', action="store_true", help='Ignore sample times from the INI file and sample every time an infection happens.')
parser.add_argument('-vis_data', '--visualize_data', action="store_true", help='Visualize the data table as a dataframe in the console.')
parser.add_argument('-g0', '--initial_genomes', action="store_true", help='Save the initial genomes of the population in a CSV.')
args = parser.parse_args()

"""
-------------------------
Validate string arguments
-------------------------
"""

if args.ratio_super_vs_normal:
    ratio_super_vs_normal = float(args.ratio_super_vs_normal)
    if not 0 <= ratio_super_vs_normal <= 1:
        raise ValueError("The ratio of number of Super Strain individuals/number of Normal Strain individuals must be between 0 and 1!")
else:
    ratio_super_vs_normal = None

if args.percentage_infected:
    percentage_infected = float(args.percentage_infected)
    if not 0 <= percentage_infected <= 1:
        raise ValueError("The percentage of infected individuals in the population must be between 0 and 1!")
else:
    percentage_infected = None

if args.percentage_super_strain:
    percentage_super_strain = float(args.percentage_super_strain)
    if not 0 <= percentage_super_strain <= 1:
        raise ValueError("The percentage of super spreaders in the population must be between 0 and 1!")
else:
    percentage_super_strain = None

if args.percentage_normal_strain:
    percentage_normal_strain = float(args.percentage_normal_strain)
    if not 0 <= percentage_normal_strain <= 1:
        raise ValueError("The percentage of normal spreaders in the population must be between 0 and 1!")
else:
    percentage_normal_strain = None

if args.max_infections:
    max_infections = int(args.max_infections)
    if not max_infections >= 1:
        raise ValueError("The number of maximum infections must be a positive integer!")
else:
    max_infections = None

if args.max_movements:
    max_movements = int(args.max_movements)
    if not max_movements >= 1:
        raise ValueError("The number of maximum movements must be a positive integer!")
else:
    max_movements = None

if args.end_time:
    end_time = float(args.end_time)
    if end_time <= 0:
        raise ValueError("The simulation time cannot be negative.")
else:
    end_time = None

if args.end_events:
    end_events = int(args.end_events)
    if not end_events >= 1:
        raise ValueError("The number of maximum events (movements+infections) must be a positive integer!")
else:
    end_events = None

"""
-------------------------------------
Validate arguments from the .INI file
-------------------------------------
"""

if (super_strain == False and n_i > 0) or (super_strain == True and n_i == 0):
    raise ValueError("For super strain to exist parameter 'super_strain = true' and 'n_i > 0', otherwise 'super_strain = false' and 'n_i = 0'")

if region_pos not in {'start', 'middle', 'end'}:
    raise ValueError(f"Invalid region_pos: '{region_pos}'. Must be one of {'start', 'middle', 'end'}.")

if ss_form < 0 or ss_events < 0:
    raise ValueError("The events regarding the formation of the super strain mutation must be positive integers!")

if n_i > l:
    raise ValueError("The number of important genome positions (n_i) in the INI file must be postitive intiger, smaller than the genome length!")

if ii > n:
    raise ValueError("The initial number of infected individuals (ii) can't be bigger than the number of individuals (n) in the simulation")

if not 0 <= prob_inf_ns <= 1 and not 0 <= prob_inf_ss <= 1:
    raise ValueError("The probability that the infector infects an individual in their infection distance must be between 0 and 1!")

#%% Results directory
"""
=================================
RESULTS' DIRECTORY INITIALIZATION 
=================================
"""

## Create the directories where the output will be saved ##

## Append current timestamp to directorie's name ##
timestamp = datetime.now().strftime("%d_%m_%Y_%H_%M")
results_directory = directory + f'simulation_{timestamp}/'
## Create the names for the subdirectories ##
samples = results_directory + 'samples'  # Directory to store genomes per sampled generation (CSV format) 
genomes =  results_directory + 'genomes' # Directory to store general information about infections and movements per sampled generation (CSV format)  

if not os.path.exists(results_directory):
    os.mkdir(results_directory)
    os.mkdir(samples)
    os.mkdir(genomes)
else:
    print(f"The directory {results_directory} already exists")

#%% Command output
"""
=============================================================
CREATING A TXT FILE WITH THE COMMAND USED FOR THIS SIMULATION
=============================================================
"""

command = ' '.join(sys.argv) # Construct the command string from sys.argv
## Define what each flag means ##
flag_explanations = {
    '-r': 'Provide the ability to recombinate genomes, if 2 infections happen during a certain time period.',
    '-all_inf': 'The simulation stops when all the individuals are infected at least once.',
    '-ratio': 'The simulation stops when the given ratio of super Strain individuals / normal Strain individuals becomes real.',
    '-per_inf': 'The simulation stops when the infected individuals in the population reach the given percentage.',
    '-per_ss': 'The simulation stops when the super spreaders in the population reach the given percentage.',
    '-per_ns': 'The simulation stops when the normal spreaders in the population reach the given percentage.',
    '-max_inf': 'The simulation stops when the given number of infections happens.',
    '-max_mv': 'The simulation stops when the given number of movements happens.',
    '-time': 'The simulation stops at the given simulation time.',
    '-events': 'The simulation stops when the given number of events (movements+infections) happpens.',
    '-vis_data': 'Visualize the data table as a dataframe in the console.',
    '-g0': 'Save the initial genomes of the population in a CSV.',
    '-ratio': 'The simulation stops when the given ratio of super Strain individuals / normal Strain individuals becomes real.',
    '-sample_inf': 'Ignore sample times from the INI file and sample every time an infection happens.',
    '-manual': 'Begin the simulation with all-0 genomes. Introduce the super strain mutation at the "ss_form" event.',
    '-msprime': 'Begin the simulation with genomes created by msprime. Introduce the super strain mutation at the "ss_form" event.'
    }
used_flags = {flag: explanation for flag, explanation in flag_explanations.items() if flag in command} # Filter the used flags and their explanations
log_command(results_directory, command, used_flags)                                                    # Log the command and flags

#%% Data table initialization
"""
=======================================================
CREATING A DATA TABLE WITH THE INFO FOR EACH INDIVIDUAL
=======================================================
"""

"""
The columns in Coords_(2 or t)
------------------------------
# 0 = x coords
# 1 = y coords
# 2 = label (not infected = 0, one infection = 1, two infections (recombination) = 2)
# 3 = movement rate
# 4 = infection rate
# 5 = mutation label (normal strain = 1, super strain = 2)
# 6 = susceptibility (susceptible to the virus = 1, not susceptible to the virus = 0)
------------------------------
"""

'''Genomes'''
label_i = infection_label(n,ii)                                                     # Divide the people in infected and healthy 
g = pd.DataFrame(index=range(n), columns=range(l))                                  # Initialization of an empty data frame where the rows correspond to the genomes of different individuals and the columns to the genome positions
total_rate = r_m * l                                                                # Total rate of mutation of genome
if super_strain == True:
    ss_pos = ss_mutation_position(n_i, l, region_pos)                               # For Super Strain(s): generate the random position in their important genome area where mutation will happen
initially_infected_ind = np.where(label_i == 1)[0]                                  # Indices of infected individuals
healthy_ind = np.where(label_i == 0)[0]                                             # Indices of healthy individuals
if args.msprime_genomes:
    infected_g = msprime_genomes(n, ii, l, r_rec, r_m, use_demography, initial_pop, past_pop, t_gen) # Generate genomes using msprime simulator with or without population expansion
    if super_strain == True:
        infected_g = ss_msprime_genomes(infected_g, ss_pos)                         # Set super strain position all to 0.0 to reassure initial genomes are all normal viruses
    for i, idx in enumerate(initially_infected_ind):                                # For each infected individual, copy the corresponding genome from msprime_genomes to g
        g.iloc[idx] = infected_g.iloc[i]
else:
    for idx in initially_infected_ind:
        g.iloc[idx] = 0.0                                                           # Infected individuals have all-0 genomes

'''Data Table'''
coords_2 = np.concatenate([coords_function(n, bound_l, bound_h), label_i], axis =1) # Initiate the array data table with the x,y coordinates and the infection label
# Rate of movement #
probm = np.zeros(n)                                                                 # Initialization of movement rate array
probm[initially_infected_ind] = rm_i                                                # All infected individuals have the corresponding rate of movement 
probm[healthy_ind] = rm_h                                                           # All healthy individuals have the corresponding rate of movement 
# Rate of infection #
probi = np.zeros(n)                                                                 # Initialization of infection rate array
probi[initially_infected_ind] = ri_n                                                # All the infected individuals have the corresponding rate of infection (initially there are only normal strains)
# Susceptibility #
sus = np.ones(n)                                                                    # Initialization of susceptibility rate array
sus[initially_infected_ind] = 0                                                     # The infected individuals are not susceptible to the virus
# Mutation label #
mut = np.zeros(n)                                                                   # Initialization of mutations array
mut[initially_infected_ind] = 1                                                     # For those with Normal Strain, they are labeled as mutation 1

coords_2 = np.concatenate([coords_2, np.column_stack(probm).T, np.column_stack(probi).T, np.column_stack(mut).T, np.column_stack(sus).T], axis=1) # Gather all the information in the final array table of info
coords_t = coords_2.copy() # Create a copy of the data table, to use during the simulation

# OPTIONAL: Save the initial genomes of the sample in a csv
if args.initial_genomes:
    g.to_csv(samples+'/initial_genomes.csv', header=False, index=False) 

# OPTIONAL: Print the data table as a dataframe (for visualization in the console)
if args.visualize_data:
    print(pd.DataFrame(data=coords_2.T, index=("x","y", "label", "rate of movement", "rate of infection", "mutation", "susceptibility")).T)

#%% Variables initialization
"""
======================================================================================
INITIALIZATION OF THE SIMULATION VARIABLES AND OF THE INFORMATION THAT WE WANT TO SAVE
======================================================================================
"""

df_i = initial_distances(coords_t)                  # Initial distance matrix 
rt_m = sum(coords_t[:, 3])                          # Total rate of movement 
rt_i = sum(coords_t[:, 4])                          # Total rate of infection 
total_mv = 0                                        # Total number of movements
total_inf = 0                                       # Total number of infections
totat_rec = 0                                       # Total number of recoveries
hah = np.zeros((1, 5))                              # Empty array to store who infected who
recovered_ind = []                                  # List that keeps who got recovered 
all_inf = np.zeros((1,5))                           # List with the number of Total infected, Super spreaders, Normal spreaders, Infection times for each event and total events (up to that point)
if args.all_infected_once:
    all_infected_once = np.empty((n, 1))            # Initiate the array to store the individuals that got infected at least once
    all_infected_once[:] = np.nan                   # Using NaN to signify empty positions
    for j in initially_infected_ind:
        all_infected_once[j] = j                    # Add the initially infected
event_type = np.zeros((1,4))                        # List with the type of each event
'''Time'''
coord_time = time.time()-func_time                  # Time to run the functions
distm_time = time.time()-coord_time                 # Time to calculate the distance matrix
t_s = 0                                             # Event time (initialization)
tt = -1                                             # Variable that corresponds to the event number (1st event = 0)
t_ss = []                                           # Array that saves the time for every time step in the simulation
t_rec = []                                          # Times of recovery
t_i = np.full((n,1), 999999999, dtype=float)        # Array that keeps the time when individuals get infected (initialization with a big value)
t_recomb = np.full((n,1), 999999999, dtype=float)   # Array that keeps the time when individuals get infected for the 2nd time - time of recombination (initialization with a big value)
t_immunity = np.full((n,1), 999999999, dtype=float) # Array that keeps the time an individual is supposed to loose their immunity after recovery (initialization with a big value)
t_recovery = np.full((n,1), 999999999, dtype=float) # Array that keeps the time an individual is supposed to recover based on their infection and recovery time (initialization with a big value)

## Initially infected ##
for j in initially_infected_ind:
    t_i[j] = 0                 # Add 0 as the infection time for those initially infected   
    t_recovery[j] = rec_t_ns   # Add the time of recovery for those initially infected (0 + recovery time = recovery time)
t_minrec = np.min(t_recovery)  # The minimum of the list of the times that individuals are supposed to recover

#%% Beggining of simulation - Break Scenarios
"""
======================
RUNNING THE SIMULATION
======================
"""

## Run the simulation until everyone is healthy ##
while sum(coords_t[:,2])!= 0: 
    
    #%% Recovery
    """
    --------
    RECOVERY
    --------
    """

    ## If there are any infected individuals proceed with the recovery process ##
    if (coords_t[:,2] == 1).any() or (coords_t[:,2] == 2).any(): 
        
        ## If there are infected people and the simulation time is bigger than the minimum recovery time ...                 ##
        ## Recovery takes place for specific individuals according to their recovery time and the time of the 1st infection  ##
        ## The time of recombination is not taken into consideration                                                         ##

        t_minrec = np.min(t_recovery)                      # Minimum recovery time
        recovery_idx = np.where(t_recovery <= t_minrec)[0] # Individuals who are about to recover (indices where t_recovery equals t_minrec)
        
        if (len(recovery_idx) > 0) and (t_s >= t_minrec):

            print("Recovered individuals (normal spreaders): " + ", ".join(map(str, [int(x) for x in recovery_idx])))

            ## Create a table with the individuals that recover ##                                                    
            recovered_ind = np.concatenate([recovered_ind, recovery_idx]) # List that keeps who got recovered
            totat_rec += len(list(recovery_idx))                          # Update the counter totat_rec to count the number or recoveries

            normal_idx = recovery_idx[coords_t[recovery_idx, 5] == 1]  # Normal spreaders to recover
            super_idx = recovery_idx[coords_t[recovery_idx, 5] == 2]   # Super spreaders to recover

            ## Updating time related variables ##
            for i in range(len(list(recovery_idx))):
                t_rec.append(t_s)                     # Save the simulation time of recovery 
            t_immunity[normal_idx] = t_s + imm_t_ns   # Time normal spreaders will stop having immunity
            t_immunity[super_idx] = t_s + imm_t_ss    # Time super spreaders will stop having immunity
            t_i[recovery_idx] = 999999999             # Re-initialize the infection times for those that recover 
            t_recomb[recovery_idx] = 999999999        # Re-initialize the recombination times for those that recover 
            t_recovery[recovery_idx] = 999999999      # Re-initialize the times individuals are supposed to recover         

            ## Updating the data table ##                                                                       
            coords_t[recovery_idx,2] = 0              # Change the infection label to healthy (0)
            coords_t[recovery_idx,3] = rm_h           # Change the rate of movement back to the one for healthy individuals
            coords_t[recovery_idx,4] = 0              # Change the rate of infection back to 0
            coords_t[recovery_idx,5] = 0              # Change the label of mutation back to 0
            if args.recombination:
                coords_t[recovery_idx,6] = 0          # Change the susceptibility label to 0, explicitly, for those that their genome didn't recombine,and still have sus = 1
            
            ## Updating genomes ##
            g.iloc[recovery_idx] = np.nan             # Remove the genome of the recovered individual, meanining make all of the positions nan again
            
            ## Updating rates (movement & infection) ##
            rt_m = sum(coords_t[:,3])                 # New rate of movement
            rt_i = sum(coords_t[:,4])                 # New rate of infection    
        
            ## Explicitly go back to while ... ## 
            continue 

    #%% Susceptibility
    """
    --------------
    SUSCEPTIBILITY
    --------------
    """ 

    ## If there are any healthy individuals that are not susceptible to the virus proceed with the susceptibility process ##
    if ((coords_t[:, 2] == 0) & (coords_t[:, 6] == 0)).any():

        ## Find the indices of the individuals that are healthy and not susceptible ##
        h_not_sus = np.where((coords_t[:, 2] == 0) & (coords_t[:, 6] == 0))[0]
        
        ## Go through the individuals that are healthy and not susceptible (indexing them with j) ... ##
        for j in h_not_sus:
            
            ## If the current simulation time is bigger than the recovery + immunity time ##
            if t_s >= t_immunity[j]:
                coords_t[j,6] = 1          # Change the susceptibility label to 1 (susceptible to the virus) 
                t_immunity[j] = 999999999  # Re-initialize the time an individual is supposed to loose their immunity after recovery

    #%% Break scenarios
    """
    ---------------
    BREAK SCENARIOS
    ---------------
    """    
    
    ## OPTIONAL ##

    ## If the number of individuals with Normal Strain are less than a certain % (ratio_super_vs_normal) of the individuals with Super Strain, stop the simulation! ##
    if ratio_super_vs_normal and (sum(coords_t[:, 5]==1) < (ratio_super_vs_normal * sum(coords_t[:, 5]==2))): # Ratio of # Super Strain ind / # Normal Strain ind
        print(f"The simulation ended because individuals with Normal Strain are less than {ratio_super_vs_normal*100}% of the individuals with Super Strain.") 
        break
    
    ## If the number of the infected individuals is more than a certain % (percentage_infected) of the population, stop the simulation! ##
    if percentage_infected and (sum(coords_t[:, 2]==1) >= (percentage_infected * n)): 
        print(f"The simulation ended because the {percentage_infected*100}% of the population is infected.")
        break
    
    ## If the number of individuals with Super Strain are more than a certain % (percentage_super_strain) of the individuals, stop the simulation! ##
    if percentage_super_strain and (sum(coords_t[:, 5]==2) >= (percentage_super_strain * n)): 
        print(f"The simulation ended because individuals with Super Strain are more than {percentage_super_strain*100}% of the individuals in the population.") 
        break

    ## If the number of individuals with Normal Strain are more than a certain % (percentage_normal_strain) of the individuals, stop the simulation! ##
    if percentage_normal_strain and (sum(coords_t[:, 5]==1) >= (percentage_normal_strain * n)): 
        print(f"The simulation ended because individuals with Normal Strain are more than {percentage_normal_strain*100}% of the individuals in the population.") 
        break
    
    ## If the total infections are more than a certain number (max_infections), stop the simulation! ##
    if max_infections and (total_inf >= max_infections): # Maximum infections
        print(f"The simulation ended because {max_infections} infections happended, totally, in the population.")
        break
       
    ## If the total movemets are more than a certain number (max_movements), stop the simulation! ##
    if max_movements and (total_mv >= max_movements): 
        print(f"The simulation ended because {max_movements} movements happended, totally, in the population.")
        break
    
    ## If simulation time is equal or higher than the one given, stop the simulation! ##
    if args.end_time and (t_s >= end_time):
        print(f"The simulation ended because it was running for {end_time} time.")
        break

    ## If events (movements+infections) is equal or higher than the one given, stop the simulation! ##
    if args.end_events and (tt >= end_events):
        print(f"\nThe simulation ended because {end_events} (movements+infections) happened during the simulation.")
        break

    ## If all the individuals got infected at least once, stop the simulation! ##
    if args.all_infected_once and np.all(~np.isnan(all_infected_once)):
        print("The simulation ended because all the individuals got infected at least once!")
        break

    #%% Simulation Time
    """
    ----------
    EVENT TIME
    ----------
    """ 
    print("----------------------------------")
    ## Calculate the time that an event will happen ##
    t_s += np.random.exponential(scale=(1/(rt_i+rt_m))) # time when an event will happen (the scale is 1/rate, because of how the exponential function is defined!)
    t_ss.append(t_s)                                    # Keep the event times in a list
    tt += 1                                             # Event number
    print(f"Event:{tt}, Simulation Time:{t_s}")

    #%% Choose event
    """
    --------------------------------
    CHOOSING WHICH EVENT WILL HAPPEN
    --------------------------------
    """

    time_bfloop = time.time()- distm_time # Time before event-loop

    ## Calculating the probabilities for the infection and the movement event ##
    p_i = rt_i/(rt_i+rt_m) # Probability of infection (rt_i = sum(coords_t[:, 4]))
    p_m = rt_m/(rt_i+rt_m) # Probability of movement  (rt_m = sum(coords_t[:, 3]))
  
    s1 = np.random.random() # Random number in [0,1]
    
    #%% Movement Event
    ## If s1 is smaller than p_m (probability of movement), the event that will happen is movement ##
    if s1 <= p_m: 
        
        ## Calculate the cumulative sum of the movement rates from all the individuals in order to create the probability axis of the movement event ##
        cum_sum = np.cumsum(coords_t[:, 3]) 

        ## Select a random number s2 in [0, maximum value of the cumulative sum] to choose which individual might move ##
        s2 = np.random.uniform(0, max(cum_sum))
        
        ## Look for any part of the cumsum array where s2 is smaller (or equal) to the closest (to the right) value of the cumsum ##
        change = np.where(s2<=cum_sum, True, False)
        
        """
        --------
        MOVEMENT
        --------
        """
        
        ## If the condition is met for one value of cumcum then proceed to move the selected individual ##
        if np.any(change):
            
            ## Find the index of the first True value in the "change" array, which corresponds to the selected individual's index "c" ##
            c = np.amin(np.where(change)) 
            
            ## Calculate the new position (in coordinates) for the selected individual ##
            if coords_t[c,5] == 2: 
                coords_f = movement(coords_t, bound_l, bound_h, c, rim_ss) # Individual with a Super Strain
            else: 
                coords_f = movement(coords_t, bound_l, bound_h, c, rim_ns) # Individual with a Normal Strain or healthy individuals

            ## Calculate the new distance matrix with the updated position of the selected individual ##
            df_f = new_distances(coords_t, coords_f, df_i, c)
            
            ## Store the newly computed distance matrix (df_f) into a variable df_i for further calculations ##
            df_i = df_f.copy()
            
            ## Update the position of the selected individual in the coords_t array with the newly computed position coords_f ##
            coords_t[c, :] = coords_f.copy()
            
            ## Update the total_mv counter of movements ##
            total_mv += 1

            ## Store the type of event (movement) and the individual that moved ##
            event_type = np.concatenate([event_type, np.column_stack((np.array([tt]), np.array([t_s]), np.array(['movement']), np.array([c])))])
            
            ## Optional: Print the event that just happened ##
            print(f"Movement of:{c}")

        ## Updating rates (movement & infection) ##
        rt_m = sum(coords_t[:,3]) # New rate of movement
        rt_i = sum(coords_t[:,4]) # New rate of infection

    #%% Infection Event
    ## In the other case, where s1 is bigger than p_m (probability of movement), infection will be the event happening ##
    else:                                       
        
        ## Calclulate the cumulative sum of the infected individuals' infection rates in order to make the probability axis of the infection event ##
        cum_sum = np.cumsum(coords_t[:, 4]) 

        ## Select a random number s3 in [0, maximum value of the cumulative sum] to choose which individual will infect ##
        s3 = np.random.uniform(0, max(cum_sum)) 

        ## Create an array "change" that has False in every cell that s3<=cum_sum, and True in every cell that s3>=cum_sum ## 
        change = np.where((s3<=cum_sum), True, False)
             
        ## Find the index of the first True value in the change array, which corresponds to the infected individual who will be the source of the infection (c) ##
        c = np.amin(np.where(change)) # First one to have cum_sum > s3, therefore the one that infects

        ## Store the type of event (infection) and the infector  ##
        event_type = np.concatenate([event_type, np.column_stack((np.array([tt]), np.array([t_s]), np.array(['infection']), np.array([c])))])

        """
        ---------
        INFECTION
        ---------
        """
        
        ## Calculate the probability of the selected individual to infect other individuals depending on the distance between them ##
        if coords_t[c,5] == 2: 
            ipi = ind_probi(df_i, c, inf_dist_ss, prob_inf_ss) # If the infector is a Super Spreader
        else: 
            ipi = ind_probi(df_i, c, inf_dist_ns, prob_inf_ns) # If the infector is a Normal Spreader
        
        ## Get the number of infected individuals (label:1) before the new infection process ##
        ib = len(coords_t[coords_t[:,2] == 1])
        
        ## Generate a new  ##

        ## List to store the infected individuals (label:1) that their genome will get recombined in the currect event ##
        recombined_g = []

        #%% Super Strain Formation
        """
        ----------------------
        Super Strain Formation
        ----------------------
        """ 
        
        ## If 1. super strains are enabled in the simulation                              ##
        ##    2. time variable that runs the simulation (tt) is between the given numbers ## 
        ##    3. there are no infected individuals with the super strain                  ##
        if super_strain == True and (ss_form <= tt <= ss_form + ss_events) and len(coords_t[coords_t[:,5] == 2]) == 0:
           
            ## Go through all the individuals (indexing them with j) ... ##
            for j in range(n): 

                ## Find those who have a non-zero probability of getting infected due to their distance ##           
                if ipi[j] != 0:
                  
                    ## Pick a random number s4 ##
                    s4 = np.random.random() 

                    ## Find those who:  1. Are not already infected                                 ## 
                    ##                  2. Are susceptible to the virus                             ##
                    ##                  3. The random number s4 is smaller or equal to their ipi[j] ##
                    if coords_t[j,2] == 0 and coords_t[j,6] == 1 and s4 <= ipi[j]: 

                        ## Updating genomes ##
                        g.iloc[j] = g.iloc[c].copy()     # The individual's genome is passed from the infector (c) to the newly infected individual (j)
                        g.loc[j, ss_pos] = 1.0           # Set the value at 1.0 in the selected positions in the beneficial area
                        
                        ## Updating the data table ##                                                                       
                        coords_t[j,2] = 1                # Change the infection label to infected (1)
                        coords_t[j,3] = rm_i             # Change the rate of movement for the infected individuals
                        coords_t[j,4] = ri_s             # Change the rate of infection for the infected individuals with a super strain
                        coords_t[j,5] = 2                # Change the label of mutation (2:infected individuals with a super strain)
                        if not args.recombination:
                            ## If the parser for recombination is used, the sus remains equal to 1 because they are prone to 2nd infection ##
                            coords_t[j,6] = 0            # Change the susceptibility label to 0 (not susceptible)             

                        ## Updating time-related variables ##
                        t_i[j] = t_s                     # Update the event time of infection for the infected individual in the list t_i 
                        t_recovery[j] = t_s + rec_t_ss   # Update the recovery time for the infected individual in the list t_recovery 

                        ## Update the total_inf counter of infections ##
                        total_inf += 1

                        ## Collect the data for the infected and infector, the times of infection (1st infection), the event time and the infection rate ##
                        hah = np.concatenate([hah, np.column_stack(np.array((c, j, 1.0, float(t_s), float(coords_t[j,4])), dtype=float))], axis=0)

                        if args.all_infected_once:
                            ## Add the individual who got infected to the all_infected_once array to keep track of the infividuals that got infected at least once ##
                            all_infected_once[j] = j         
                        
                        ## Explicitly continue to the next iteration (next individual) ##           
                        continue  

        #%% Infection  
        else:

            ## Go through all the individuals (indexing them with j) ... ##
            for j in range(n): 

                ## Find those who have a non-zero probability to get infected due to their distance ##           
                if ipi[j] != 0:

                    ## Pick a random number s4 ##
                    s4 = np.random.random() 

                    ## Find those who:  1. Are not already infected                                 ## 
                    ##                  2. Are susceptible to the virus                             ##
                    ##                  3. The random number s4 is smaller or equal to their ipi[j] ##
                    if coords_t[j,2] == 0 and coords_t[j,6] == 1 and s4 <= ipi[j]: 

                        ## Updating genomes ##
                        g.iloc[j] = g.iloc[c].copy()                # The individual's genome is passed from the infector (c) to the newly infected individual (j)
                        g.iloc[j] = mutation(g.iloc[j], total_rate) # The individual's (j) genome goes through the mutation procedure
                        
                        ## Updating the data table and recovery time ##                                                                       
                        coords_t[j,2] = 1                           # Change the infection label to infected (1)
                        coords_t[j,3] = rm_i                        # Change the rate of movement for the infected individuals
                        if super_strain == True and n_i > 0 and (g.iloc[j, ss_pos] == 1.0).any(): # Check if the genome is a super strain
                            coords_t[j,4] = ri_s                    # Change the rate of infection for the infected individuals with a super strain
                            coords_t[j,5] = 2                       # Change the label of mutation (2:infected individuals with a super strain)
                            t_recovery[j] = t_s + rec_t_ss          # Change the recovery time for the infected individual in the list t_recovery 
                        else:
                            coords_t[j,4] = ri_n                    # Change the rate of infection for the infected individuals with a normal strain
                            coords_t[j,5] = 1                       # Change the label of mutation (1:infected individuals with a normal strain)
                            t_recovery[j] = t_s + rec_t_ns          # Change the recovery time for the infected individual in the list t_recovery 
                        if not args.recombination:
                            ## If the parser for recombination is used, the sus remains equal to 1 because they are prone to 2nd infection ##
                            coords_t[j,6] = 0                       # Change the susceptibility label to 0 (not susceptible)               

                        ## Updating time-related variables ##
                        t_i[j] = t_s                                # Change the event time of infection for the infected individual in the list t_i 
                            
                        ## Update the total_inf counter of infections ##
                        total_inf += 1
                        
                        ## Collect the data for the infected and infector, the times of infection (1st infection) the event time and the infection rate ##
                        hah = np.concatenate([hah, np.column_stack(np.array((c, j, 1.0, float(t_s), float(coords_t[j,4])), dtype=float))], axis=0)

                        if args.all_infected_once:
                            ## Add the individual who got infected to the all_infected_once array to keep track of the infividuals that got infected at least once ##
                            all_infected_once[j] = j      
                            
                        ## Explicitly continue to the next iteration (next individual) ##
                        continue  
                    
                    #%% Genetic recombination
                    ## Find those who:  1. Are infected once                                        ## 
                    ##                  2. Are susceptible to the virus                             ##
                    ##                  3. The random number s4 is smaller or equal to their ipi[j] ##
                    elif args.recombination and (coords_t[j,2] == 1 and coords_t[j,6] == 1 and s4 <= ipi[j]) :

                        p = rec_probi(r_rec, l) # Number to determine whether recombination will take place 

                        ## If p = 0 : No genetic recombination ##
                        if p == 0:
                            continue # Continue to the next iteration (next j)

                        ## If p > 0 : Genetic recombination ##
                        else:
                            
                            """
                            ---------------------
                            Genetic recombination
                            ---------------------
                            """ 

                            ## If the individual carries a strain and their time for recovery hasn't come yet ... ##
                            if coords_t[j,2] == 1 and (t_s <= t_recovery[j]):

                                recombined_g.append(j) # List of the infected individuals (label:1) that their genome will get recombined

                                ## Updating genomes ##
                                g.iloc[j] = recombine_genomes(g.iloc[c], g.iloc[j], p) # The individual's genome is recombined with the infector's (c) 
                                g.iloc[j] = mutation(g.iloc[j], total_rate)            # The individual's (j) genome goes through the mutation procedure 
                                
                                ## Updating the data table ## 
                                # the rate of movement for the infected individuals remains the same (rm_i)
                                coords_t[j,2] = 2                       # Change the infection label to infected with a recombined strain (2)
                                if super_strain == True and n_i > 0 and (g.iloc[j, ss_pos] == 1.0).any(): # Check if the genome, after recombination, is a super strain
                                    coords_t[j,4] = ri_s                # Change the rate of infection for the infected individuals with a super strain
                                    coords_t[j,5] = 2                   # Change the label of mutation (2:infected individuals with a super strain)
                                else:                                   # Check if the genome, after recombination, is a normal strain
                                    coords_t[j,4] = ri_n                # Change the rate of infection for the infected individuals with a normal strain
                                    coords_t[j,5] = 1                   # Change the label of mutation (1:infected individuals with a normal strain)
                                coords_t[j,6] = 0                       # Change the susceptibility label to 0 (not susceptible) 

                                ## Updating time related variables ##
                                t_recomb[j] = t_s  # Update the time of recombination for the infected individual in the list t_recomb

                                ## Collect the data for the infected and infector, the event time and the infection rate ##
                                hah = np.concatenate([hah, np.column_stack(np.array((c, j, 2.0, float(t_s), float(coords_t[j,4])), dtype=float))], axis=0)                      

                                ## Explicitly continue to the next iteration (next individual) ##
                                continue 

                            else:
                                ## Explicitly continue to the next iteration (next individual) ##
                                continue 
                                
                    ## If neither condition is met, continue to the next iteration (next individual) ##
                    else:
                        continue

        #%% ... Infection

        ## Updating rates (movement & infection) ##
        rt_m = sum(coords_t[:,3])  # New rate of movement
        rt_i = sum(coords_t[:,4])  # New rate of infection

        ## Check if the number of infected individuals before the infection event (ib) is greater than or equal to the number of newly infected individuals (ia) after the infection ##
        ## If so, it means no new infections occurred, and "Nobody got infected" is printed.                                                                                         ##
        ## Otherwise, the indices of the newly infected individuals are printed as "These got infected:" using the hah data collected during the loop                                ##
        ia = len(coords_t[coords_t[:,2] == 1])

        if ib >= ia:
            if not recombined_g or not args.recombination:
                print("Nobody got infected.")
            else:
                print("Nobody got newly infected but this individual's genome got recombined:", ", ".join(map(str, recombined_g)))
                if args.sample_infection: # Sample when an infection happens
                    ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
                    sample_data(samples, genomes, g, tt, coords_t, all_inf, tt)
        else:
            if not recombined_g or not args.recombination:
                print("This/These individuals got newly infected: " + ", ".join(map(str, [int(x) for x in hah[-int(ia-ib):,1]]))) # array hah[-int(ia-ib):,1]: contains the newly infected individuals 
                if args.sample_infection: # Sample when an infection happens
                    ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
                    sample_data(samples, genomes, g, tt, coords_t, all_inf, tt)
            else:
                print("This/These individuals got newly infected: " + ", ".join(map(str, [int(x) for x in hah[-int(ia-ib):,1]]))) # array hah[-int(ia-ib):,1]: contains the newly infected individuals 
                print("This/These individuals' genome got recombined: ", ", ".join(map(str, recombined_g)))
                if args.sample_infection: # Sample when an infection happens
                    ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
                    sample_data(samples, genomes, g, tt, coords_t, all_inf, tt)
        
        ## Collect the data for the number of Total infected, Super spreaders, Normal spreaders, Infection times for each event and total events (up to that point) ##
        all_inf = np.concatenate([all_inf, np.column_stack(np.array((sum(coords_t[:,2] != 0), sum(coords_t[:,5]==2), sum(coords_t[:,5]==1), float(t_s), int(tt)), dtype=float))], axis=0)    

    #%% Save data - Info per event

    ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
    if not args.sample_infection:
        sample_data(samples, genomes, g, tt, coords_t, all_inf, sample_times)
    else:
        pass
    ## Time to run the simulation loop ##
    loop_t = time.time()-time_bfloop 

    print(f"Totally infected:{len(coords_t[coords_t[:,2]==1])+len(coords_t[coords_t[:,2]==2])} (ns:{len(coords_t[coords_t[:,5]==1])},ss:{len(coords_t[coords_t[:,5]==2])})")
    
## Save the data from the simulation in the correct directory ##
save_data(samples, genomes, coords_2, coords_t, g, all_inf, recovered_ind, hah, t_rec, event_type, tt)

## Time for the whole simulation to run ##
end=time.time()-initial

print(f"\nEvents: {tt}")
print("Total simulation time (minutes):", end/60)   