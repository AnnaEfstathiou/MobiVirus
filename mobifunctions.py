#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
---------------------------------
Functions for MobiVirus Simulator
---------------------------------
"""

# Packages
import numpy as np
import random
from datetime import datetime
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import statistics
import math
import os
import sys
import msprime
from configparser import ConfigParser

"""
===============================
CREATE A "command_log.txt" FILE
===============================
"""

def log_command(directory, command, flags):

    ## Creates a txt containing the command used for the simulation. ##

    # Create the directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Create a log file in the results directory
    log_file_path = os.path.join(directory, 'command_log.txt')

    # Open the log file and append the command and flags information
    with open(log_file_path, 'a') as log_file:
        log_file.write(f"Command: {command}\n")
        for flag, explanation in flags.items():
            log_file.write(f"{flag}: {explanation}\n")

        # Checking why the simulation stopped
        if not any(flag in flags for flag in ['-ratio', '-per_inf', '-per_ss', '-per_ns', '-max_inf', '-max_mv', '-all_inf', '-time', '-events']):
            log_file.write("\nThe simulation stopped because everyone is healthy.\n")
        
        # Append the initial parameters from the config
        log_file.write("\nInitial Parameters:\n")
        log_file.write(f"Parameter that defines if there will be super strains in the simulation (super_strain): {config.getboolean('Type_of_strains', 'super_strain')}\n")
        log_file.write(f"Number of individuals (n): {config.getint('Initial_Parameters', 'n')}\n")
        log_file.write(f"Length of genome (l): {config.getint('Initial_Parameters', 'l')}\n")
        log_file.write(f"Initial number of infected individuals (ii): {config.getint('Initial_Parameters', 'ii')}\n")
        log_file.write(f"Lower bound for the spacial axis (bound_l): {config.getfloat('Initial_Parameters', 'bound_l')}\n")
        log_file.write(f"Upper bound for the spacial axis (bound_h): {config.getfloat('Initial_Parameters', 'bound_h')}\n")
        log_file.write(f"Mutation rate for each genome position (r_m): {config.getfloat('Initial_Parameters', 'r_m')}\n")
        log_file.write(f"Number of important genome positions (n_i): {config.getint('Initial_Parameters', 'n_i')}\n")
        log_file.write(f"Area of the important positions in the viral genome that define the super strain (region_pos): {config.get('Initial_Parameters', 'region_pos')}\n")
        log_file.write(f"Infection rate (ri): {config.getfloat('Initial_Parameters', 'ri_n')} (ns), {config.getfloat('Initial_Parameters', 'ri_s')} (ss)\n")
        log_file.write(f"Movement rate (rm): {config.getfloat('Initial_Parameters', 'rm_i')} (infected), {config.getfloat('Initial_Parameters', 'rm_h')} (healthy)\n")
        log_file.write(f"Infection distance (inf_dist): {config.getfloat('Initial_Parameters', 'inf_dist_ns')} (ns), {config.getfloat('Initial_Parameters', 'inf_dist_ss')} (ss)\n")
        log_file.write(f"Probability of infection (prob_inf): {config.getfloat('Initial_Parameters', 'prob_inf_ns')} (ns), {config.getfloat('Initial_Parameters', 'prob_inf_ss')} (ss)\n")
        log_file.write(f"Recombination rate (r_rec): {config.getfloat('Initial_Parameters', 'r_rec')}\n")
        log_file.write(f"Recovery time (rec_t): {config.getfloat('Initial_Parameters', 'rec_t_ns')} (ns), {config.getfloat('Initial_Parameters', 'rec_t_ss')} (ss)\n")
        log_file.write(f"Immunity time (imm_t): {config.getfloat('Initial_Parameters', 'imm_t_ns')} (ns), {config.getfloat('Initial_Parameters', 'imm_t_ss')} (ss)\n")
        log_file.write(f"Relative infected mobility (rim): {config.getfloat('Initial_Parameters', 'rim_ns')} (ns), {config.getfloat('Initial_Parameters', 'rim_ss')} (ss)\n")
        log_file.write(f"Generations to get a sample (sample_times): {config.getint('Initial_Parameters', 'sample_times')}\n")        
        log_file.write(f"Event that the super strain mutation is introduced for the 1st time: {config.getint('Super_strain_Parameters', 'ss_formation')}\n")   
        log_file.write(f"Period of events when the super strain is introduced (if previously there are no individuals with a super strain): {config.getint('Super_strain_Parameters', 'ss_events')}\n")       
        log_file.write(f"msprime parameter: Parameter defining whether to simulate with population demography or not (use_demography): {config.getboolean('msprime_Parameters', 'use_demography')}\n") 
        log_file.write(f"msprime parameter: Present-day population size (initial_pop): {config.getint('msprime_Parameters', 'initial_pop')}\n")    
        log_file.write(f"msprime parameter: Past population size (past_pop): {config.getint('msprime_Parameters', 'past_pop')}\n") 
        log_file.write(f"msprime parameter: Time (in generations) at which population growth stops and becomes constant (t_gen): {config.getint('msprime_Parameters', 't_gen')}\n")

"""
==================================
READ THE PARAMETERS FROM .INI FILE
==================================
"""

# Specify the directory and file name that contains the parameters
directory = './'
file_name = 'parameters.ini'
file_path = os.path.join(directory, file_name)

# Check if the file exists
if not os.path.exists(file_path):
    print(f"Error: {file_name} must be in the current directory.")
    sys.exit("Exiting the simulation.")

# File exists, proceed with parsing
config = ConfigParser()
config.read(file_path)


"""
-----------------------------------------------
Read directory from the Initial_Parameters section
-----------------------------------------------
"""

# Directory where the scripts exist
directory = config.get('Directory', 'directory').strip('"')                 

"""
===================================
BUILDING THE DATASET OF INDIVIDUALS 
===================================
"""

''' 
--------------------
Coordinates Function 
--------------------
'''

def coords_function(n, bound_l, bound_h):
    
    ## Creates a stacked column array of the (x,y) 2D coordinates of each individual in the dataset. ##
    ## The coordinates are assigned randomly using a uniform distribution.                           ##
    
    # n = total number of individuals in the sample
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # x = random coordinate in the x axis 
    # y = random coordinate in the y axis
        
    x = np.random.uniform(low=bound_l, high=bound_h, size=n)
    y = np.random.uniform(low=bound_l, high=bound_h, size=n)

    return np.column_stack((x, y))

'''
------------------------
Infection Label Function
------------------------
'''

def infection_label(n, ii):
    
    ## Creates a stacked column array where every cell corresponds to an individual, with a label depending on if they are infected with the virus or not. ##
    ## In the initial stage of the simulation, there is a sample of infected individuals among the healthy ones.                                           ##
    
    # n = total number of individuals in the sample
    # ii = number of infectd individuals
    # label 0 = healthy individual
    # label 1 = infected individual (only 1 infection)
    
    infected_array = np.array([1]*ii)
    healthy_array = np.array([0]*(n-ii))
    all = np.concatenate((infected_array, healthy_array))
    random.shuffle(all)
    labels = np.column_stack(all)

    return labels.T

'''
---------------
Genome Function
---------------
'''

def msprime_genomes(n,ii,l,r_rec,r_m,use_demography,initial_pop,past_pop,t_gen):

    ## Generate genomes using msprime simulator with or without population expansion ##

    # n = population size
    # ii = number of genomes (infected individuals)
    # l = genomes length
    # r_rec = recombination rate
    # r_m = mutation rate for each genome position
    # use_demography = Parameter defining whether to simulate with population demography (exponential growth) or not
    #                  true = use demography, false = constant population size
    # initial_pop = Present-day population size (used as the final size after growth if demography is enabled)
    # past_pop = Past population size (population t generation ago)
    # t_gen = Time (in generations) at which population growth stops and becomes constant

    ## Note in case population demography is used
    # If pop_initial_size < pop_past_size → negative growth rate
    # If pop_initial_size > pop_past_size → positive growth rate
    # If pop_initial_size = pop_past_size → zero growth rate (no growth)

    # Function to apply the transformation: turn even numbers to 0.0 and odd (>1) to 1.0 to make the genome binary
    def transform_value(x):
        if x > 1.0:
            if x % 2 == 0:  # even
                return 0.0
            else:           # odd
                return 1.0
        return x
    
    if use_demography == True:

        """ Calculate Growth Rate"""
        a = (1 / t_gen) * math.log(initial_pop / past_pop) # Note: math.log(x) returns the natural logarithm, i.e., ln(x) (base e)

        """ Define Demography with Expansion """
        demography = msprime.Demography()
        
        # Assume population expanded from 100 to n over 500 generations
        demography.add_population(
            name = "source", 
            initial_size = initial_pop, # present-day size after growth
            growth_rate = a)            # exponential growth from the past

        # The past population size gets smaller and smaller, approaching zero as 𝑡 → ∞
        # Else...
        demography.add_population_parameters_change(
            time = t_gen,                                   # 'time' generations ago
            growth_rate = 0.0,                              # stop growth at 'time' generations behind
            population = "source")

        """ Simulate ancestry with demography """
        tree_sequence = msprime.sim_ancestry(
            samples = {"source": ii},                       # number of genomes
            demography = demography,                        # specified demography
            sequence_length = l,                            # genomes length
            recombination_rate = r_rec,                     # recombination rate
            ploidy = 1,                                     # haplotypes
            random_seed = int(datetime.now().timestamp()))  # seed for reproducibility

    else:

        """ Simulate ancestry without demography """
        tree_sequence = msprime.sim_ancestry( 
            population_size = n,                            # population size
            samples = ii,                                   # number of genomes
            sequence_length = l,                            # genomes length
            recombination_rate = r_rec,                     # recombination rate
            ploidy = 1,                                     # haplotypes
            random_seed = int(datetime.now().timestamp()))  # seed for reproducibility
    
    """ Introduce mutations on the simulated ancestry """
    mutated_tree_sequence = msprime.sim_mutations(
        tree_sequence, 
        rate = r_m, 
        random_seed = int(datetime.now().timestamp()))

    """ Process the simulated sequences """
    initial_genomes = np.zeros((ii, l)) # initialize an array for the initial viral genomes 

    # Loop through each variant in the mutated tree sequence
    for var in mutated_tree_sequence.variants():
        pos_index = int(var.site.position)                # convert the site position to an integer index (column in initial_genomes)
        if pos_index < l:                                 # safeguard for floating-point position errors
            initial_genomes[:, pos_index] = var.genotypes # insert the variants into the corresponding column of initial_genomes

    viral_genomes = pd.DataFrame(initial_genomes)         # convert the genomes into a DataFrame
    viral_genomes = viral_genomes.map(transform_value)    # munipulate the mutations
    
    return viral_genomes

def ss_msprime_genomes(viral_genomes, positions):

    ## Set super strain position all to 0.0 to reassure initial genomes are all normal viruses ##

    # viral_genomes = genomes generated from msprime
    # positions = position in the region of determinant genome positions for the infected with mutation 2 (Super Strain)
    
    viral_genomes[positions] = 0.0 # set super strain position all to 0.0
    
    return viral_genomes

def ss_mutation_position(n_i, l, position):
    
    ## Assign one random mutation (1.0) in the region of important genome positions for the infected with mutation 2 (Super Strain) ##
    
    # n_i = important positions in the genome
    # l = total length of the genome
    # position ('start', 'middle', or 'end') = area of the important positions in the viral genome
    
    if position == 'middle':
        middle_start = (l - n_i) // 2
        middle_end = middle_start + n_i
        positions = list(range(middle_start, middle_end))
    elif position == 'start':
        positions = list(range(n_i))
    elif position == 'end':
        positions = list(range(l - n_i, l))
    else:
        raise ValueError("Invalid position type. Choose from 'middle', 'start', or 'end'.")
    
    return positions
    # return random.choice(positions)

'''
--------------------------------
Initial Distance Matrix Function
--------------------------------
'''

def initial_distances(coords):

    ## Creates a dataframe containing the distance between all the individuals in the sample ##
    ## The rows and the columns of the dataframe both correspond to the individuals          ##

    # coords: numpy array of the individual's spatial coordinates

    distances = pdist(coords[:, :2], metric='euclidean') # Compute the pairwise distances using scipy's pdist, which is optimized
    distance_matrix = squareform(distances) # Convert the condensed distance matrix into a square form

    return distance_matrix

"""
================================================
SIMULATING THE SPACIAL MOVEMENT OF AN INDIVIDUAL 
================================================
"""

'''
-----------------
Movement Function
-----------------
'''

def movement(coords, bound_l, bound_h, c, rim): 

    ## Calculates the new spatial coordinates of the individual that will move                                                                                                                                    ##
    ## Returns the full data table but with new coordinates                                                                                                                                                       ##
    ## The new coordinates are randomly selected from two different Gaussian distributions (one for each coordinate)                                                                                              ##
    ## The mean value of the distribution is the initial coordinate (x or y accordingly) and the stdev is taken from the stdev of the distribution of all of the coordinates from the sample (x or y accordingly) ##
    ## If the individual is infected, the stdev of the Gaussian distribution is rims times smaller than the one for a healthy individual                                                                          ##
    ## This is because that an infected individual is more likely to stay in the same place or move closer to their initial position                                                                              ##
    ## Individuals are not allowed to move outside the box (or space)                                                                                                                                             ##

    # coords = the full coordinate table in order to get the standard deviations
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # c = the individual that will move
    # rim = relative infected mobility. The parameter that differentiates between the mobility of a healthy individual and an infected one
    
    x_c = 0                 # new x coordinate
    y_c = 0                 # new y coordinate
    label_c = coords[c, 2]  # label of individual (healthy or infected)
    probm_c = coords[c, 3]  # rate of movement
    probi_c = coords[c, 4]  # rate of infection 
    mut_c = coords[c, 5]    # mutation label
    sus_c = coords[c, 6]    # susceptibility

    stdev_x = statistics.stdev(coords[:, 0]) # std from the distribution of x coordinate from our sample
    stdev_y = statistics.stdev(coords[:, 1]) # std from the distribution of y coordinate from our sample

    # if the individual is healthy
    if coords[c, 2] == 0: 
        # x coordinates
        a = np.random.normal(coords[c,0], stdev_x)
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], stdev_x)
        # y coordinates        
        b = np.random.normal(coords[c,1], stdev_y)
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], stdev_y)
        x_c = a
        y_c = b
    
    # if the individual is infected
    else:
        # x coordinates
        a = np.random.normal(coords[c,0], stdev_x*rim)
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], stdev_x*rim)
        # y coordinates
        b = np.random.normal(coords[c,1], stdev_y*rim)
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], stdev_y*rim)          
        x_c = a
        y_c = b

    return np.column_stack((x_c, y_c, label_c, probm_c, probi_c, mut_c, sus_c))


'''
---------------------------------------
Distance Matrix After Movement Function
---------------------------------------
'''

def new_distances(coords, coords_f, distances, c):

    ## Changes only the column of the initial distance matrix that coresponds to the individual that moved       ##
    ## For that column only, it calculates the new distances between the individual c and the rest in the sample ##

    # coords = the full coordinate table in order to get the standard deviations
    # coords_f = full data table with coordinates after movement
    # distances = distance matrix
    # c = the individual that will move
    
    x1, y1 = coords_f[0, :2] # Extract moved individual's new coordinates
    x2, y2 = coords[:, :2].T # Calculate distances from the moved individual to all others
    
    new_dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2) # Euclidean distance

    # Update distance matrix
    df = distances.copy()
    df[c, :] = new_dist  # Update row for individual c
    df[:, c] = new_dist  # Update column for individual c to maintain symmetry
    df[c, c] = 0.0       # Set distance of the infector to themselves to 0

    return df

'''
---------
Infection
---------
'''

def ind_probi(df, c, inf_dist, prob_inf):
    
    ## Gives an array of length as many as the individuals and shape (n,1)                                  ##
    ## In the cells, non-zero values are for those who are in the "infection distance" and 0 for the others ##
        
    n = len(df) # number of individuals
    ipi = np.zeros((n)) # initialize an array to store the probabilities of someone to get infected
    for i in range(n):
        ipi[i] = np.where((df[i,c] < inf_dist) & (df[i,c] != 0), prob_inf, 0) # 1st condition: infection distance
                                                                              # 2nd condition: exclude the individual that will infect
    
    return ipi

    
def mutation(g_i, r_tot): 
    
    ## Mutation of an individuals' genome ##
    
    # g_i = genome of an individual 
    # r_tot = total rate of mutation of genome
    
    p = np.random.poisson(r_tot)                 # Random Poisson number with lambda = r_tot to give me how many mutations will happen
    N = np.random.randint(0, len(g_i), size = p) # Positions in the genome that mutations will happen (as many as the poisson)
    for i in range(len(N)):
        g_i[N[i]] = np.where(g_i[N[i]] == 0.0, 1.0, 0.0) # Mutate the selected positions: if the original value is 0, set it to 1; otherwise, set it to 0
    
    return g_i

'''
--------------------
Genome Recombination
--------------------
'''

def rec_probi(r_rec, l):

    ## Number to determine whether recombination will take place ##
    ## If p = 0 >> No genetic recombination                      ##
    ## If p > 0 >> Genetic recombination                         ##
    ## Genetic recombination is less possible for small r_rec    ##

    rec_total = r_rec * (l-1)        # recombination rate per genome
    p = np.random.poisson(rec_total) # Random Poisson number with lambda = rec_total to determine if recombination will happen

    return p


def recombine_genomes(genome1: pd.Series, genome2: pd.Series, p: int) -> pd.Series:

    ## Creating a recombined genome out of two genomes ##
    ## e.g.

        # INPUT
        # genome1: 0.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0 (pd.Series format)
        # genome2: 1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,1.0 (pd.Series format)
        # p = 3 (positive integer)

        # STEPS
        # 1. Choose 3 (p) random positions e.g. 2, 10, 5 > sort them: 2, 5, 10
        # 2. Choose randomly the genome to start e.g. genome2
        # 3. a. Take the 1st part of genome2 until position 2,
        #    b. the 2nd part of genome1 from position 2 until 5, 
        #    c. the 3rd part of genome2 from position 5 until 10, 
        #    d. the 4th part of genome1 from position 10 until the end.

        # OUTPUT
        # Recombined genome will be: 1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0 (pd.Series format)
   
    if len(genome1) != len(genome2):
        raise ValueError("Genomes must be of the same length")
    
    if p <= 0:
        raise ValueError("p must be a positive number")
    
    positions = sorted(np.random.choice(range(1, len(genome1)), p, replace=False)) # Generate p random unique positions
    start_genome = np.random.choice([1, 2])                                        # Choose the initial genome to start
    
    new_genome = []  # Initialize the new genome
    previous_position = 0
    
    # Perform recombination
    for i, pos in enumerate(positions + [len(genome1)]):
        if start_genome == 1:
            new_genome.extend(genome1[previous_position:pos])
            start_genome = 2
        else:
            new_genome.extend(genome2[previous_position:pos])
            start_genome = 1
        previous_position = pos
    
    return pd.Series(new_genome)

'''
--------------------------
Output - Results Functions 
--------------------------
'''

''' DATA '''

def sample_data(samples_directory, genomes_directory, g, tt, coords_t, all_inf, sample_times):

    if tt==0 or tt%sample_times==0:
                
        g.to_csv(genomes_directory+'/genomes_'+str(tt)+'.csv',header=False, index=False)
        
        coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
        coords_t.to_csv(samples_directory+'/coords_'+str(tt)+'.csv', header=True, index=False)

        all_inf = all_inf[1:, :] # Remove the first row of zeros using numpy slicing
        all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time', 'Event'])
        all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Event']] = all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Event']].astype(int) # Convert some columns to integer while keeping 'Time' as float
        all_inf.to_csv(samples_directory+'/all_inf_'+str(tt)+'.csv', header=True, index=False)


def save_data(samples_directory, genomes_directory, coords_2, coords_t, g, all_inf, unin, hah, t_un, event_type, tt):
    
    g.to_csv(genomes_directory+'/genomes_'+str(tt)+'.csv',header=False, index=False)
    # g.to_csv(genomes_directory+'/genomes_'+'final'+'.csv',header=False, index=False)
    
    unin = np.concatenate([np.column_stack(np.array((unin, t_un), dtype=float))], axis=0)
    unin = pd.DataFrame(data=unin, columns=['Recovered individual', 'Recovery Time'])
    unin.to_csv(samples_directory+'/recovery.csv', header=True, index=False)
    
    hah = hah[1:, :] # Remove the first row of zeros using numpy slicing
    hah = pd.DataFrame(data=hah, columns=['Infecting', 'Infected', 'Infection Time', 'Simulation Time', 'Infection Rate'])
    hah.to_csv(samples_directory+'/infections.csv', header=True, index=False)
    
    coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
    coords_t.to_csv(samples_directory+'/coords_'+str(tt)+'.csv', header=True, index=False)
    # coords_t.to_csv(samples_directory+'/coords_final.csv', header=True, index=False)
    
    coords_2 = pd.DataFrame(data=coords_2, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
    coords_2.to_csv(samples_directory+'/initial_coords.csv', header=True, index=False)

    all_inf = all_inf[1:, :] # Remove the first row of zeros using numpy slicing
    all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time', 'Event'])
    all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Event']] = all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Event']].astype(int) # Convert some columns to integer while keeping 'Time' as float
    all_inf.to_csv(samples_directory+'/all_inf_'+str(tt)+'.csv', header=True, index=False)
    # all_inf.to_csv(samples_directory+'/all_inf_'+'final'+'.csv', header=True, index=False)

    event_type = event_type[1:, :] # Remove the first row of zeros using numpy slicing
    event_type = pd.DataFrame(data=event_type, columns=['Event', 'Simulation Time', 'Event Type','Individual'])
    event_type.to_csv(samples_directory+'/event_type.csv', header=True, index=False)