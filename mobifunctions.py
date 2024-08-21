#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions for MobiVirus Simulator
"""

#Packages
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import statistics
import os
import sys
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
        if not any(flag in flags for flag in ['-ratio', '-per_inf', '-per_ss', '-max_inf', '-max_mv', '-sus', '-all_inf', '-time', '-events']):
            log_file.write("\nThe simulation stopped because everyone is healthy.\n")
        # elif 'all_inf' in flags.keys():
        #     log_file.write("\nThe simulation might stop when all the individuals are infected at least once.\n")
        # elif 'ratio' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the given ratio of super Strain individuals/normal Strain individuals becomes real.\n")      
        # elif 'per_inf' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the infected individuals in the population reach the given percentage.\n")
        # elif 'per_ss' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the super spreaders in the population reach the given percentage.\n")
        # elif 'max_inf' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the given number of infections happen.\n")
        # elif 'max_mv' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the given number of movements happen.\n")
        # elif 'sus' in flags.keys():
        #     log_file.write("\nThe simulation might stop when only the given number of individuals are susceptible to the virus.\n")
        # elif 'time' in flags.keys():
        #     print("FOUND")
        #     log_file.write("\nThe simulation might stop at the given simulation time.\n")
        # elif 'events' in flags.keys():
        #     log_file.write("\nThe simulation might stop when the given number of events (movements+infections) happpen.\n")
        
        # Append the initial parameters from the config
        log_file.write("\nInitial Parameters:\n")
        log_file.write(f"Number of individuals (n): {config.getint('Initial_Parameters', 'n')}\n")
        log_file.write(f"Length of genome (l): {config.getint('Initial_Parameters', 'l')}\n")
        log_file.write(f"Initial number of infected individuals (ii): {config.getint('Initial_Parameters', 'ii')}\n")
        log_file.write(f"Lower bound for the spacial axis (bound_l): {config.getfloat('Initial_Parameters', 'bound_l')}\n")
        log_file.write(f"Upper bound for the spacial axis (bound_h): {config.getfloat('Initial_Parameters', 'bound_h')}\n")
        log_file.write(f"Mutation rate for each genome position (r_m): {config.getfloat('Initial_Parameters', 'r_m')}\n")
        log_file.write(f"Number of important genome positions (n_i): {config.getint('Initial_Parameters', 'n_i')}\n")
        log_file.write(f"Rate of infection of 1st strain (ri_n): {config.getfloat('Initial_Parameters', 'ri_n')}\n")
        log_file.write(f"Rate of infection of 2nd strain (ri_s): {config.getfloat('Initial_Parameters', 'ri_s')}\n")
        log_file.write(f"Rate of movement of infected individuals (rm_i): {config.getfloat('Initial_Parameters', 'rm_i')}\n")
        log_file.write(f"Rate of movement of healthy individuals (rm_h): {config.getfloat('Initial_Parameters', 'rm_h')}\n")
        log_file.write(f"Infection distance (inf_dist): {config.getfloat('Initial_Parameters', 'inf_dist_ns')} (ns), {config.getfloat('Initial_Parameters', 'inf_dist_ss')} (ss)\n")
        log_file.write(f"Probability that the infector infects an individual in their infection distance (prob_inf): {config.getfloat('Initial_Parameters', 'prob_inf_ns')} (ns), {config.getfloat('Initial_Parameters', 'prob_inf_ss')} (ss)\n")
        log_file.write(f"Rate of recombination (r_rec): {config.getfloat('Initial_Parameters', 'r_rec')}\n")
        log_file.write(f"Recovery time (rec_t): {config.getfloat('Initial_Parameters', 'rec_t_ns')} (ns), {config.getfloat('Initial_Parameters', 'rec_t_ss')} (ss)\n")
        log_file.write(f"Relative infected mobility (rim): {config.getfloat('Initial_Parameters', 'rim')}\n")
        log_file.write(f"Generations to get a sample (sample_times): {config.getint('Initial_Parameters', 'sample_times')}\n")

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

def coords_func(n, bound_l, bound_h):
    
    ## Creates a stacked column array of the (x,y) 2D coordinates of each individual in the dataset ##
    ## The coordinates are assingned randomly using a uniform distribution ##
    
    # n = total number of individuals in the sample
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # x = random coordinate in the x axis 
    # y = random coordinate in the y axis
    
    # Default values:  The defined space is a rectangle with the same bounds in the x and y axis 
    
    x = np.random.uniform(low=bound_l, high=bound_h, size=n)
    y = np.random.uniform(low=bound_l, high=bound_h, size=n)

    return np.column_stack((x, y))

'''
------------------------
Infection Label Function
------------------------
'''

def infection_label(n, i):
    
    ## Creates a stacked column array where every cell corresponds to an individual, with a label depending on if they are infected with the virus or not ##
    ## In the initial stage of the simulation, there is a sample of infected individuals among the healthy ones ##
    
    # n = total number of individuals in the sample
    # i = number of infectd individuals
    # label 0 = healthy individual
    # label 1 = infected individual (only 1 infection)
    
    z = np.array([1]*i)
    o = np.array([0]*(n-i))
    s = np.concatenate((z, o))
    random.shuffle(s)
    l = np.column_stack(s)
    
    return l.T

'''
---------------
Genome Function
---------------
'''

def genome(n, l, r_m):
    
    ## Creates an empty data frame where the rows correspond to the genomes of different individuals and the columns to the genome positions ##
    ## Returns the dataframe together with the total mutation rate of a genome ##
    
    # n = total number of individuals in the sample
    # l = total length of the genome
    # r_m = mutation rate of a position in the genome
    
    # Default result: If the genome of an individual has a specific value in at least one of the important positions, the individual is infected by the Super Strain
    
    g =  pd.DataFrame(index=range(n), columns=range(l))
    r_tot = r_m * l # Total rate of mutation of genome
    
    return g, r_tot

def ss_mutation_position(n_i, l, position):
    
    ## Assign one random mutation (1.0) in the region of important genome positions for the infected with mutation 2 (Super Strain) ##
    
    # n_i = important positions in the genome
    # l = total length of the genome
    # position ('start', 'middle', or 'end') = area of the important positions in the individual's genome
    
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
    
    return random.choice(positions)

'''
----------------------
Movement Rate Function
----------------------
'''

def probs(coords, rm_i, rm_h):
    
    ## Creates a stacked column array where every cell corresponds to an individual with a value for the rate of movement for that individual ##
    ## The rate of movement depends on the infection label (see parameters.ini) ##
    
    # coords = the full coordinate table
    # rm_i = rate of movement for infected individual
    # rm_h = rate of movement for healthy individual
    
    probm = np.zeros(len(coords)) # random rate of movement for each individual

    mask = coords[:, 2] != 0.0    # for those who are infected (label 1 (or 2))
    probm[mask] = rm_i            # they get rm_i in the corresponding positions in probm

    mask = coords[:, 2] == 0.0    # for those that are healthy
    probm[mask] = rm_h            # they get rm_h in the corresponding positions in probm
    
    return np.column_stack(probm).T 


'''
-----------------------
Infection Rate Function
-----------------------
'''

def infectivity(in_probi, g, n_i, ri_s, position):

    ## Gives the specific value of the infection rate of an individual, depending on the strain of the virus ##
    ## The strain is decided by the n_i important positions in the individual's genome ##
    
    # g = genome of individual
    # n_i = important positions in the genome
    # in_probi = infection rate of the infector
    # ri_s = infectivity rate of super strain
    # position ('start', 'middle', or 'end') = area of the important positions in the individual's genome
    
    if position == "start":

        if np.any(np.array(g[:n_i])==1):
            probi = ri_s
        else:
            probi = in_probi

    elif position == "middle":

        genome_length = len(g)
        middle_start = (genome_length - n_i) // 2
        middle_end = middle_start + n_i
        if np.any(np.array(g[middle_start:middle_end]) == 1):
            probi = ri_s
        else:
            probi = in_probi
    
    elif position == "end":

        genome_length = len(g)
        end_start = genome_length - n_i
        if np.any(np.array(g[end_start:]) == 1):
            probi = ri_s
        else:
            probi = in_probi
    
    else:
        raise ValueError("Invalid position value. Expected 'start', 'middle', or 'end'.")
    
    return probi

'''
--------------------------------
Initial Distance Matrix Function
--------------------------------
'''

def ini_distm(coords):
    
    ## Creates a dataframe containing the distance between all the individuals in the sample ##
    ## The rows and the columns of the dataframe both correspond to the individuals ##
    
    # coords = dataframne of the individual's spacial coordinates
    
    n = len(coords)
    df = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            x1, y1 = coords[i, :2]
            x2, y2 = coords[j, :2]
            df[j,i] = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    
    return df


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

    ## Calculates the new spatial coordinates of the individual that will move ##
    ## Returns the full data table but with new coordinates ##
    ## The new coordinates are randomly selected from two different Gaussian distributions (one for each coordinate) ##
    ## The mean value of the distribution is the initial coordinate (x or y accordingly) and the standard deviation is taken from the standard deviation of the distribution of all of the coordinates from our sample (x or y accordingly) ##
    ## If the individual is infected, the standard deviation of the Gaussian distribution is rim*10 times smaller than the one for a healthy individual ##
    ## In that way we add the condition that an infected individual is more likely to stay in the same place or move closer to its initial position ##
    
    
    # coords = the full coordinate table in order to get the standard deviations
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # c = the specific individual that will move
    # If these standard deviations are defined, we use them and not the standard deviation from all of the coordinates in the sample
    # rim = relative infected mobility. The parameter that differentiates between the mobility of a healthy individual and an infected one
    ## Individuals are not allowed to move outside the box (or space) ##
    
    x_f = 0                 # new x coordinate
    y_f = 0                 # new y coordinate
    label_f = coords[c, 2]  # label of individual (healthy or infected)
    probm_f = coords[c, 3]  # random rate of movement
    probi_f = coords[c, 4]  # rate of infection of individual
    mut_f = coords[c, 5]    # mutation label
    sus_f = coords[c, 6]    # susceptibility

    r_x = statistics.stdev(coords[:, 0]) # std from the distribution of x coordinate from our sample
    r_y = statistics.stdev(coords[:, 1]) # std from the distribution of y coordinate from our sample

    # if the individual is healthy
    if coords[c, 2] == 0: 

        a = np.random.normal(coords[c,0], r_x )
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], r_x )
                
        b = np.random.normal(coords[c,1], r_y )
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], r_y )
                 
        x_f = a
        y_f = b
    
    # if the individual is infected
    else:
            
        a = np.random.normal(coords[c,0], r_x*rim )
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], r_x*rim )
                
        b = np.random.normal(coords[c,1], r_y*rim )
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], r_y*rim )
                 
        x_f = a
        y_f = b

    return np.column_stack((x_f, y_f, label_f, probm_f, probi_f, mut_f, sus_f))


'''
---------------------------------------
Distance Matrix After Movement Function
---------------------------------------
'''

def new_dist(coords, coords_f, initial_dist, c): # where c in the index of the coords that corresponds the individual that moved
    
    ## Changes only the column of the initial distance matrix that coresponds to the individual that moved ##
    ## For that column only, it calculates the new distances between the individual c and the rest in the sample ##

    n = len(coords)
    df = initial_dist.copy()
    for i in range(n):
        for j in range(i+1, n):
            x1, y1 = coords_f[0,:2]
            x2, y2 = coords[j, :2]
            df[j,i] = np.where((i==c or j==c), np.sqrt((x2-x1)**2 + (y2-y1)**2), df[j,i])  # Change only the distance matrix cells of the individul that moved
    
    return df

'''
---------
Infection
---------
'''

def ind_probi(df, c, inf_dist, prob_inf):
    
    ## Gives an array of length as many as the individuals and shape (n,1), only for the one individual that infects in that moment ##
    ## In the cells, 1 is for those who are in the "infection distance" and 0 for the others ##
    
    # 1 is also the probability of infection
    
    n = len(df) # the number of individuals
    ipi = np.zeros((n))
    for i in range(n):
        ipi[i] = np.where((df[i,c] < inf_dist) & (df[i,c] != 0), prob_inf, 0) # The first condition defines the infection distance
    
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

# def mutation(g_i, r_tot):
    
#     # Mutation of an individual's genome #
    
#     # g_i = genome of an individual (Pandas Series or DataFrame)
#     # r_tot = total rate of mutation of genome
    
#     p = np.random.poisson(r_tot) # Random Poisson number with lambda = r_tot to give how many mutations will happen
#     N = np.random.randint(0, len(g_i), size=p) # Positions in the genome where mutations will happen (as many as the Poisson number)
    
#     # Create a copy to avoid modifying the original DataFrame/Series in place
#     g_i_copy = g_i.copy()
    
#     for i in range(len(N)):
#         g_i_copy[N[i]] = np.where(g_i_copy[N[i]]==0, 1, 0) # Mutate the selected positions: if the original value is 0, set it to 1; otherwise, set it to 0
    
#     return g_i_copy

'''
--------------------
Genome Recombination
--------------------
'''

def rec_probi(r_rec, l):

    ## Number to determine whether recombination will take place ##
    ## If p = 0 >> No genetic recombination ##
    ## If p > 0 >> Genetic recombination ##
    ## Genetic recombination is less possible for small r_rec ##

    rec_total = r_rec * (l-1) # recombination rate per genome
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
    
    # Generate p random unique positions
    positions = sorted(np.random.choice(range(1, len(genome1)), p, replace=False))
    
    # Choose the initial genome to start
    start_genome = np.random.choice([1, 2])
    
    # Initialize the new genome
    new_genome = []
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

# ''' PLOTS '''

# def plot_map(plots_directory, coords, tt, ts, mv, un, super_strain = False):
    
#     ## Create plot with all of the regions ##

#     fig, ax = plt.subplots()
#     ax.scatter(200, 200, c='k', s=0)  # Add a point at coordinates (200, 200) with no size (s=0) for scaling purposes
#     # Set the x and y axis labels
#     ax.set_xlabel('x')
#     ax.set_ylabel('y')
#     # Set the x and y axis limits
#     ax.set_xlim(-2, 8)
#     ax.set_ylim(-0.5, 10)

#     ax.set_title("Simulation time =""{0:.4f}".format(ts)) # Set the title of the plot with the current simulation time

#     # Scatter plot infected individuals with the normal strain in red
#     ax.scatter(list(coords[(coords[:, 2] != 0) & (coords[:, 5] == 1)][:, 0]),
#                list(coords[(coords[:, 2] != 0) & (coords[:, 5] == 1)][:, 1]), s=3, c='indianred',
#                label='Normal infections = %i' % (len(coords[coords[:, 5] == 1][:, 0])))  
#     if super_strain == True:
#         # Scatter plot infected individuals with the super strain in blue
#         ax.scatter(list(coords[(coords[:, 2] != 0) & (coords[:, 5] == 2)][:, 0]),
#                 list(coords[(coords[:, 2] != 0) & (coords[:, 5] == 2)][:, 1]), s=3, c='blue',
#                 label='Super spreaders = %i' % (len(coords[coords[:, 5] == 2][:, 0])))  
#     # Scatter plot uninfected individuals in gold
#     # ax.scatter(list(coords[coords[:, 2] == 0][:, 0]), list(coords[coords[:, 2] == 0][:, 1]), s=3, c='gold')
#     ax.scatter(list(coords[coords[:, 2] == 0][:, 0]), list(coords[coords[:, 2] == 0][:, 1]), s=3, c='gold',
#                label='Not infected individuals = %i' % len(coords[coords[:, 2] == 0][:, 0]))
    
#     ax.scatter([], [], label="Total Un-infections = %i" % un, s=0) # empty scatter plot for displaying legend information

#     plt.legend(title="Cumulative Movements = %i" % mv, loc='upper center') # Add legend with title for cumulative movements
#     plt.savefig(plots_directory + '/' + str(tt) + '.jpg', dpi=1200)
#     plt.close(fig)

# def do_plots(plots_directory, coords_t, tt, t_s, mv, un, super_strain):
    
#     ## Make plots according to the following conditions:   ##
#     ## if its the 0th generation of the simulation,        ##
#     ## if the number of generation can be devided by 1000, ##
#     ## if the "Super Strain" disappears.                   ##
    
#     if tt==0 or tt%1000==0 or sum(coords_t[:,5]==1)==0: 
#         plot_map(plots_directory, coords_t, tt, t_s, mv, un, super_strain)


''' DATA '''

def sample_data(samples_directory, genomes_directory, g, tt, coords_t, all_inf, sample_times):

    if tt==0 or tt%sample_times==0:
                
        g.to_csv(genomes_directory+'/genomes_'+str(tt)+'.csv',header=False, index=False)
        
        coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
        coords_t.to_csv(samples_directory+'/coords_'+str(tt)+'.csv', header=True, index=False)

        all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time', 'Events'])
        all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Events']] = all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Events']].astype(int) # Convert some columns to integer while keeping 'Time' as float
        all_inf.to_csv(samples_directory+'/all_inf_'+str(tt)+'.csv', header=True, index=False)


def save_data(samples_directory, genomes_directory, coords_2, coords_t, g, all_inf, unin, hah, t_un):
    
    g.to_csv(genomes_directory+'/genomes_'+'final'+'.csv',header=False, index=False)
    
    unin = np.concatenate([np.column_stack(np.array((unin, t_un), dtype=float))], axis=0)
    unin = pd.DataFrame(data=unin, columns=['Uninfected individual', 'Time of uninfection'])
    unin.to_csv(samples_directory+'/uninfections.csv', header=True, index=False)
    #np.savetxt(directory+'uninfections.txt', pd.DataFrame(data = unin), fmt=['%3d','%.3f'])
    
    hah = pd.DataFrame(data=hah, columns=['Infecting', 'Infected', 'Times of Infections', 'Time', 'Infection Rate'])
    hah.to_csv(samples_directory+'/infections.csv', header=True, index=False)
    #np.savetxt(directory+'infections.txt', pd.DataFrame(data = hah, columns=['Infects', 'Infected', 'Infection time', 'Infection Probability']), fmt=['%3d','%3d','%.3f', '%.1f'])
    
    coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
    coords_t.to_csv(samples_directory+'/coords_final.csv', header=True, index=False)
    
    coords_2 = pd.DataFrame(data=coords_2, columns=["x","y", "Infection label", "Rate of movement", "Rate of infection", "Mutation", "Susceptibility"])
    coords_2.to_csv(samples_directory+'/initial_coords.csv', header=True, index=False)

    all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time', 'Events'])
    all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Events']] = all_inf[['Total infected', 'Super spreaders', 'Normal spreaders', 'Events']].astype(int) # Convert some columns to integer while keeping 'Time' as float
    all_inf.to_csv(samples_directory+'/all_inf_'+'final'+'.csv', header=True, index=False)

    # extra_data = np.column_stack(np.array((ss, ns, mv), dtype=int))
    # extra_data = pd.DataFrame(data=extra_data, columns=['Total Super spreaders', 'Total Normal spreaders', 'Total movements'])
    # extra_data.to_csv(samples_directory+'/extra_data.csv', header=True, index=False)
    #np.savetxt(directory+'extra_data.txt', pd.DataFrame(data=np.column_stack(np.array((ss, ns, mv),dtype=int)), columns=['Total Super spreaders', 'Total Normal spreaders', 'Total movements']), fmt=['%1d', '%1d', '%1d'], header='Total Super infections, Total Normal infections , Total movements', comments='')
