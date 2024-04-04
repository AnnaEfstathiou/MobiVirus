#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
"""
================================
|      ------------------     |
| + | MobiVirus Simulator | + |
|      ------------------     |
================================

"""
#Packages
import argparse
import numpy as np
import random
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

from mobifunctions import coords_func, label, genome, mutation, infectivity, probs, movement, ini_distm, new_dist, ind_probi, plot_map_normal, plot_map_normal_super, sample_data_normal, sample_data_normal_super, do_plots_normal, do_plots_normal_super, save_data_normal, save_data_normal_super

"""
===================
PARSE THE .INI FILE
===================
""" 

# Specify the directory and file name that contains the parameters
directory = './'
file_name = 'parameters.ini'
file_path = os.path.join(directory, file_name)

# Check if the file exists
if not os.path.exists(file_path):
    print(f"Error: {file_name} must be in the current directory.")
    sys.exit("Exiting the program.")

# File exists, proceed with parsing
config = ConfigParser()
config.read(file_path)

""" 
========================
PARAMETER INITIALIZATION
========================
"""

"""
-----------------------------------------------
Read directory from the Initial_Parameters section
-----------------------------------------------
"""
directory = config.get('Directory', 'directory').strip('"')  # Directory where the scripts exist (.ini file & .py script with functions)

"""
-----------------------------------------------
Read values from the Initial_Parameters section
-----------------------------------------------
coords_2 = The characteristics of the individuals, meaning coordinates, labels, probabilities of movement/infection
g = genome of the individuals, with the shape of a matrix where each line is each individuals' genome, in the form of an array
r_tot = total mutation rate 
n_i = number of important positions in the genome
n = total number of individuals
l = length of genome


The columns in Coords_(2 or t)
------------------------------
# 0 = x coords
# 1 = y coords
# 2 = label (infected = 1, not infected = 0 )
# 3 = rate of movement
# 4 = rate of infection
# 5 = mutation label
# 6 = susceptibility
------------------------------
"""
r_m = config.getfloat('Initial_Parameters', 'r_m')                  # Mutation rate for each position in genome
n_i = config.getint('Initial_Parameters', 'n_i')                    # Important genome positions for the ss mutation
                                                                    # If there are no super strains, n_i = 0
n = config.getint('Initial_Parameters', 'n')                        # Number of individuals in the simulation
l = config.getint('Initial_Parameters', 'l')                        # Length of genome
bound_l = config.getfloat('Initial_Parameters', 'bound_l')          # Lower bound for the plot
bound_h = config.getfloat('Initial_Parameters', 'bound_h')          # Upper bound for the plot
ii = config.getint('Initial_Parameters', 'ii')                      # Number of infected individuals
ri_n = config.getfloat('Initial_Parameters', 'ri_n')                # Rate of infection from normal/strain 1
rm_i = config.getfloat('Initial_Parameters', 'rm_i')                # Rate of movement for infected ind.
rm_h = config.getfloat('Initial_Parameters', 'rm_h')                # Rate of movement for healthy ind.
rec_t = config.getfloat('Initial_Parameters', 'rec_t')              # Recovery time in simulation time
label_i = label(n,ii)                                               # Divide the people in infected and not infected
std_x = config.getfloat('Initial_Parameters', 'std_x')              # Standard deviation of x coordinate
std_y = config.getfloat('Initial_Parameters', 'std_y')              # Standard deviation of y coordinate
rim = config.getfloat('Initial_Parameters', 'rim')                  # Relative infection mobility
sample_times = config.getint('Initial_Parameters', 'sample_times')  # Generations where we take samples of our simulation's output
# s_n = config.getfloat('Break parameters', 's_n')                    # Ratio of #Super Strain ind/#Normal Strain ind (break argument)
# in_per = config.getfloat('Break parameters', 'in_per')              # Percentage of infected people (break argument)
# l_in = config.getfloat('Break parameters', 'l_in')                  # Limit of infection (break argument) 

"""
===============
PARSE ARGUMENTS
===============
""" 

parser = argparse.ArgumentParser(description='Creating other possible senarios (e.g. super strain), break arguments, visualization options. All the arguments are optional!')
## string parsers (for break scenarions)
parser.add_argument('-ratio', '--ratio_super_vs_normal', type=str, help='Ratio of number of Super Strain individuals/number of Normal Strain individuals.')
parser.add_argument('-per_inf', '--percentage_infected', type=str, help='Percentage of infected individuals in the population.')
parser.add_argument('-max_inf', '--max_infections', type=str, help='Maximum infections (if the infection are more, stop).')
parser.add_argument('-sus', '--percentage_susceptibility', type=str, help='Minimum susceptible individuals (if the susceptible individuals are less, stop).')
## action parsers
parser.add_argument('-s', '--super_strain', action="store_true", help='Create a super strain with higher infectivity that the normal one.')
parser.add_argument('-vis_data', '--visualize_data', action="store_true", help='Visualize the data table as a dataframe in the console.')
parser.add_argument('-g0', '--initial_genomes', action="store_true", help='Save the initial genomes of the sample in a csv.')
parser.add_argument('-plots', '--scatter_plots', action="store_true", help='Create scatter plots of the coordinates of the individuals.')
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

if args.max_infections:
    max_infections = int(args.max_infections)
    if not max_infections >= 1 or max_infections > n:
        raise ValueError("The number of maximum infections must be a positive integer and less than the total number of individuals!")
else:
    max_infections = None

if args.percentage_susceptibility:
    percentage_susceptibility = float(args.percentage_susceptibility)
    if not 0 <= percentage_susceptibility <= 1:
        raise ValueError("The percentage of susceptable individuals in the population must be between 0 and 1!")   
else:
    percentage_susceptibility = None

# OPTIONAL: if the strains are devided into normal and super
if args.super_strain:
    ri_s = config.getfloat('Initial_Parameters', 'ri_s')            # Rate of infection from super/strain 2
    L = [ri_n, ri_s]                                                # Infection rates in order to divide the infected population in different strains
else:
    # Only normal strains
    ri_s = ri_n                                                     # Rate of infection for normal strains equal to the rate of infection for super strains
    L = [ri_n, ri_s]                                                # Infection rates in order to divide the infected population in different strains
                                                                    # BUT since there is only one strain, 

#%%

"""
=================================
RESULTS' DIRECTORY INITIALIZATION 
=================================
"""
# Make the directories where the output will be saved

# Append current timestamp to subdirectory names
timestamp = datetime.now().strftime("%d_%m_%Y_%H_%M")
results_directory = directory + f'results_{timestamp}/'

samples = results_directory + 'samples'
plots = results_directory +'plots'
genomes =  results_directory + 'genomes'

if not os.path.exists(results_directory):
    os.mkdir(results_directory)
    os.mkdir(samples)
    os.mkdir(genomes)
    if args.scatter_plots: # create the plots directory only if scatter plots will be created
        os.mkdir(plots)
else:
    print(f"The directory {results_directory} already exists")

#%%

"""
=======================================================
CREATING A DATA TABLE WITH THE INFO FOR EACH INDIVIDUAL
=======================================================
"""

coords_2 = np.concatenate([coords_func(n, bound_l, bound_h), label_i], axis =1) # Initiate the array data table with 
g, r_tot = genome(n, l, r_m, n_i)                                               # Initialization of genome table (as nan with shape n*l)
probm = probs(coords_2,rm_i, rm_h)                                              # Initialization of movement rate array
probi = np.zeros(n)                                                             # Initialization of infection rate array
mut = np.zeros(n)                                                               # Initialization of mutations array
sus = np.ones(n)                                                                # Initialization of susceptibility rate array
for i in range(n):
    probi[i] = np.where(coords_2[:,2][i]==1, random.choice(L), probi[i])        # Randomly assign infected individuals with each strain
                                                                                # If there is only one strain (normal) all the individuals will have that one
    g.iloc[i] = np.where(coords_2[:,2][i]==1, np.zeros((1,l)), g.iloc[i])       # For the infected individuals, there is a row in g with the genome that they carry
    sus[i] = np.where(coords_2[:,2][i]==1, 0, sus[i])                           # For the healthy individuals, they have a Susceptibility to the virus
    if args.super_strain:
        mut[i] = np.where(probi[i]==ri_s, 2, mut[i])                            # For those with Super Strain, they are labeled as mutation 2
        mut[i] = np.where(probi[i]==ri_n, 1, mut[i])                            # For those with Normal Strain, they are labeled as mutation 1
        g.loc[i,0] = np.where(mut[i]==2, 1, g.loc[i,0])                         # For the infected with mutation 2 (Super Strain), they have the value 1 in the first position of their genome    
    else:
        mut[i] = np.where(probi[i]==ri_n, 1, mut[i])                            # For those with Normal Strain, they are labeled as mutation 1

coords_2 = np.concatenate([coords_2, probm, np.column_stack(probi).T, np.column_stack(mut).T, np.column_stack(sus).T], axis=1) # Gather all the information in the final array table of info


# OPTIONAL: Save the initial genomes of the sample in a csv
if args.initial_genomes:
    g.to_csv(samples+'/initial_genomes.csv', header=False, index=False) 

# OPTIONAL: Print the data table as a dataframe (for visualization in the console)
if args.visualize_data:
    print(pd.DataFrame(data=coords_2.T, index=("x","y", "label", "rate of movement", "rate of infection", "mutation", "susceptibility")).T)

"""
======================================================================================
INITIALIZATION OF THE SIMULATION VARIABLES AND OF THE INFORMATION THAT WE WANT TO SAVE
======================================================================================
"""

mut_all = np.column_stack(mut).T                        # Create a perpendicular array with the mutations in each infected individual
coord_time = time.time()-func_time                      # Time to run the functions
coords_t = coords_2.copy()                              # Create a copy of the data table, to use during the simulation
df_i = ini_distm(coords_t)                              # Initial Distance Matrix 
#df_ff = np.zeros((len(coords_t), len(coords_t))) 
distm_time = time.time()-coord_time                     # Time to calculate the distance matrix
#print("time after distance matrix:", distm_time)
rt_m = sum(coords_t[:, 3])                              # Total rate of movement 
rt_i = sum(coords_t[:, 4])                              # Total rate of infection 
t_s = 0                                                 # Event time (initialization)
tt = 0                                                  # Time variable that runs the simulation (or number of generation)
t_ss = []                                               # Array that saves the event time for every time step in the simulation
t_i = np.full((n,1), 999999999, dtype=float)            # Array that keeps the event time when an individual (or more) get infected 
                                                        # We initialise it with a really big value
for i in range(n):
    t_i[i] = np.where(coords_t[:,2][i]==1, 0, t_i[i])   # Add 0 as the infection time for those initially infected  
                                  
t_im = np.min(t_i)                                      # The minimum of the list of the specific times that individuals got infected
t_un = []                                               # Uninfection time
mv = 0                                                  # Count for the movements
if args.super_strain:
    ss = sum(coords_t[:, 5]==2)                             # Number of Super-spreaders
ns = sum(coords_t[:, 5]==1)                             # Number of Normal-spreades
un = 0                                                  # Total number of uninfections
hah = np.zeros((1, 4))                                  # Empty array that will be the array of who infected who
unin = []                                               # List that keeps who got uninfected 
if args.super_strain:
    all_inf = np.zeros((1,4))                           # List with the number of Total infected, Super spreaders, Normal spreaders and Infextion times for each generation
else:
    all_inf = np.zeros((1,3))                           # List with the number of Total infected, Normal spreaders and Infextion times for each generation
#equ=[]                                                 # Times when infected individuals are ~equal to uninfected 

#%%


"""
======================
RUNNING THE SIMULATION
======================
"""

if args.super_strain:
    print(f"The simulation contains 2 types of strains, a normal strain with {ri_n} rate of infection and a super strain with {ri_s} rate of infection.")
else:
    print(f"The simulation contains one type of strain, a normal strain with {ri_n} rate of infection.")

# Run the simulation until everyone becomes uninfected 
while sum(coords_t[:,2])!=0: # activate lines 236, 240

# Run the simulation until everyone becomes uninfected or until everyone is infected
# while sum(coords_t[:,2]) != 0 or sum(coords_t[:,2]) != n: 

# Run the simulation for i events
# while tt <= 10000:
    
    # OPTIONAL: Print some plots during the run
    if args.scatter_plots:
        if args.super_strain:
            do_plots_normal_super(plots, coords_t, tt, t_s, mv, ss, un) 
        else: 
            do_plots_normal(plots, coords_t, tt, t_s, mv, un)

    """
    ---------------
    BREAK SCENARIOS
    ---------------
    """    
    
    ## OPTIONAL ##

    # If there is a super strain in the population 
    if args.super_strain:
        ## If the number of individuals with Normal Strain are less than a certain % (ratio_super_vs_normal) of the individuals with Super Strain, stop the simulation! ##
        if ratio_super_vs_normal: # Ratio of # Super Strain ind / # Normal Strain ind
            if sum(coords_t[:, 5]==1) < ratio_super_vs_normal * sum(coords_t[:, 5]==2): 
                print(f"The simulation ended because individuals with Normal Strain are less than {ratio_super_vs_normal*100}% of the individuals with Super Strain.") 
                break
    else:
        if ratio_super_vs_normal:
           raise ValueError("In order to calculate the ratio of super vs normal strains in the population, the super strain must exist. Use the appropriate argument to create the super strain.")

    ## If the number of the infected individuals is more than a certain % (percentage_infected) of the population, stop the simulation! ##
    if percentage_infected: # Percentage of infected individuals in the population
        if sum(coords_t[:, 2]==1) > percentage_infected * n:
            print(f"The simulation ended because the {percentage_infected*100}% of the population is infected.")
            break
    
    ## If the total infections are more than a certain number (max_infections), stop the simulation! ##
    if max_infections: # Maximum infections
        if args.super_strain:
            if (ss + ns) > max_infections:
                print(f"The simulation ended because {max_infections} infections happended, totally, in the population.")
                break
        else:
            if ns > max_infections:
                print(f"The simulation ended because {max_infections} infections happended, totally, in the population.")
                break
    
    ## If the number of the susceptable individuals is less than a certain % (percentage_susceptibility) of the population, stop the simulation! ##
    if percentage_susceptibility:
        if sum(coords_t[:, 6]==1) < percentage_susceptibility * n:
            print(f"The simulation ended because less than {percentage_susceptibility*100}% of the population is susceptable to the virus.")
            break
    

    # Calculate the time that an event will happen
    t_s += np.random.exponential(scale=(1/(rt_i+rt_m))) #The scale is 1/rate, because of how the exponential function is defined! t_s is the time when an event will happen
    
    # Print the event time
    
    # print(t_s)
    
    t_ss.append(t_s) # Keep the event times in a list
    
    """
    ------------
    UNINFECTIONS
    ------------
    """

    if t_s >= rec_t+t_im and coords_t[:,2].any() == 1: # If the event time is bigger than the minimum recovery time and if there are infected people,
                                                       # uninfect specific individuals (according to their recovery time)
        

        uninfection_idx = np.where(t_i == t_im)[0]  # Find indices where t_i equals t_im    

        ## Optional: Prints the index of people that get uninfected (from their line in the coords data table), AKA those that their time of infection is the "oldest" ##
        print("It's un-infection time! These people get un-infected: " + ", ".join(map(str, [int(x) for x in uninfection_idx])))
        
        # Create a table with the individuals that get uninfected 
                                                                  
        unin = np.concatenate([unin, uninfection_idx]) 
        for i in range(len(list(uninfection_idx))):
            t_un.append(t_s) # Save the simulation time of uninfection   
            
        un+=len(list(uninfection_idx))     # Keep count of the uninfections                                                                        
        coords_t[uninfection_idx,2] = 0    # Change the infection label to uninfected (0)
        coords_t[uninfection_idx,3] = rm_h # Change the rate of movement back to the one for uninfected individuals
        coords_t[uninfection_idx,4] = 0    # Change the rate of infection back to 0

        if args.super_strain:
            if percentage_susceptibility:
                # For individuals with a mutation label = 1 (Normal Spreaders)
                normal_idx = uninfection_idx[coords_t[uninfection_idx, 5] == 1] # find the indexes of the normal spreders
                coords_t[normal_idx, 5] = 0  # Change the label of mutation back to 0
                coords_t[normal_idx, 6] = 1  # Change the susceptibility label back to 1

                # For individuals with a mutation label of 2 (Super Spreaders)
                super_idx = uninfection_idx[coords_t[uninfection_idx, 5] == 2] # find the indexes of the super spreders
                coords_t[super_idx, 5] = 0  # hange the label of mutation back to 0
                                            # For those with mutation label 2, no change at the susceptibility label
            else:
                coords_t[uninfection_idx, 5] = 0  # Change the label of mutation back to 0
                coords_t[uninfection_idx, 6] = 1  # Change the susceptibility label back to 1 to all unifected individuals 
        else:
            coords_t[uninfection_idx, 5] = 0  # Change the label of mutation back to 0
            coords_t[uninfection_idx, 6] = 1  # Change the susceptibility label back to 1 to all unifected individuals

        g.iloc[uninfection_idx] = np.nan   # Remove the genome of the recovered individual, meanining make all of the positions nan again

        ## Optional: Print the index of the uninfected individuals ##
        
        # print(list(t_i[np.where(t_i == t_im)[0]])) 
        
        ## Optional: Print the minimum infection time (the "oldest" one) ##
        
        #print(t_im) 
        
        ## Re-initialize the infection times for those individuals that got uninfected so that there is a new minimum ##
        t_i[uninfection_idx] = 999999999 
        
        ## Optional: Check the minimum infection time (the "oldest" one) after ##
        
        #print(np.min(t_i))
        
        ## Optional: Check that the individuals got uninfected ##
        
        # print(list(t_i[np.where(t_i == t_im)[0]])) 
        
        ## New minimum infection time ##
        
        t_im = np.min(t_i) 
            
        ## New rate of movement ##
        rt_m = sum(coords_t[:, 3]) 
        
        ## New rate of infection ##
        rt_i = sum(coords_t[:,4]) 
            
        ## Fix the event time ##
        t_s = t_ss[tt-1] # t_s is the time when an event will happen
        t_ss.pop(-1)
        #t_ss[tt] = t_s # Keep the times of events in a list
        continue
    
    
    """
    --------------------------------
    CHOOSING WHICH EVENT WILL HAPPEN
    --------------------------------
    """

    ## Calculating the probabilities for the infection event and the movement event ##
    
    p_i = rt_i/(rt_i+rt_m) # Probability of infection  
    p_m = rt_m/(rt_i+rt_m) # Probability of movement 
    
    ## Optional: Print the two probabilities ##
    
    #print("Movement over infection:", p_m/p_i)
    #print("Infection over movement:", p_i/p_m)
    
    ## Select a random number s1 to see which of the two events will happen ##
    
    s1 = np.random.random() # Random number in  [0,1]
    
    time_bfloop = time.time()-distm_time # Time before event-loop

    ## If s1 is smaller than p_m, the event that will happen is movement ##
    
    if s1<=p_m: 
        
        ## Calculate the cumulative sum of the rate of movement from all the individuals in order to create the probability axis of the movement event ##
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
            coords_f = movement(coords_t, bound_l, bound_h, c, std_x, std_y, rim)

            ## Calculate the new distance matrix now with the updated position of the selected individual ##
            df_f = new_dist(coords_t, coords_f, df_i, c)
            
            ## Store the newly computed distance matrix (df_f) into a variable df_i for further calculations ##
            df_i = df_f.copy()
            
            ## Update the position of the selected individual in the coords_t array with the newly computed position coords_f ##
            coords_t[c, :] = coords_f.copy()
            
            ## Add one movement in the mv counter of movements that have occurred ##
            mv += 1
            
            ## Optional: Print the event that just happened ##
            # print("movement")
        
        ## Update the total rate of movement ##
        rt_m = sum(coords_t[:,3]) 
        
        ## Update total rate of infection ##
        rt_i = sum(coords_t[:,4])
        
    ## In the other case, where s1 is bigger than p_m, infection will be the event happening ##
    
    else:
        
        ## Calclulate the cumulative sum of the infected individuals' infection rates (but taking into account only those with non-zero infection rates) in order to make the probability axis of the infection event ##
        cum_sum = np.where(coords_t[:, 4]!= 0 , np.cumsum(coords_t[:, 4]), 0) 
        
        ## Select a random number s3 in [0, maximum value of the cumulative sum] to choose which individual will infect ##
        s3 = np.random.uniform(0, max(cum_sum)) 

        ## Create an array "change" with non-zero cells of cum_sum that have a larger value than s3 ## 
        change = np.where((s3<=cum_sum)& (cum_sum!=0), True, False) # This array contains the indices of infected individuals whose cumulative infection rates are >=s3

        ## If the condition is met for one value in the change array, then proceed with the infection event ##
        if np.any(change):
            
            ## Find the index of the first True value in the change array, which corresponds to the infected individual who will be the source of the infection (c) ##
            c = np.amin(np.where(change)) #First one to have >=s3, therefore the one that infects
            
            ## Optional: Print the event that just occured ##
            # print("infection")
            
        """
        ---------
        INFECTION
        ---------
        """
        
        ## Calculate the probability of the selected individual to infect each other individual depending on the distance between them ##
        ipi = ind_probi(df_i, c) 
        
        ## Get the number of infected individuals before the new infection process ##
        ib = np.sum(coords_t[:,2])
        
        ## Go through all the individuals (indexing them with j) and find those who: 1. Have a non-zero probability to get infected due to their distance ##
        ##                                                                           2. Are not already infected                                          ## 
        ##                                                                           3. Are susceptable to the virus
        for j in range(n): 
            if ipi[j] != 0 and coords_t[j,2] == 0 and coords_t[j,6] == 1: 
                
                ## Pick a random number s4 ##
                ## The conditions are: If the random number s4 is less than or equal to ipi[j], it means that individual j can get infected by individual c AND the individual j must be susceptible to the virus (coords_t[6,j]==1) ##
                s4 = np.random.random() 
                
                ## The individual's genome is passed from the infector (c) to the newly infected individual (j) ##
                g.iloc[j] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), g.iloc[c], g.iloc[j]) 
                
                ## The individual's (j) genome goes through the mutation procedure ##
                g.iloc[j] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), mutation(g.iloc[j], r_tot), g.iloc[j])
                
                ## The label of the individual (j) is updated to infected after the infection ##
                coords_t[j,2] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), 1, 0) 
                
                ## The  individual's (j) rate of movement is updated ##
                coords_t[j,3] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), coords_t[c,3], coords_t[j,3]) 
                
                ## The  individual's (j) rate of infection is updated and is the same as the infectivity of the infector (c) ##
                coords_t[j,4] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), infectivity(coords_t[c, 4], g.iloc[j], n_i, ri_s), coords_t[j, 4]) 
                
                ## Add the mutation label depending on the strain of the infection + the mutation produced ##
                if args.super_strain:
                    coords_t[j,5] = np.where((s4<=ipi[j])&(coords_t[j,6]==1)&(coords_t[j,4]==ri_s), 2, np.where((s4<=ipi[j])&(coords_t[j,6]==1)&(coords_t[j,4]==ri_n), 1, coords_t[j,5]))
                else:
                    coords_t[j, 5] = np.where((s4 <= ipi[j]) & (coords_t[j, 6] == 1) & (coords_t[j, 4] == ri_n), 1, coords_t[j, 5])

                ## Update the counters to track the number of Super Strain infections (if are existed) and Normal Strain infections ##
                if args.super_strain:
                    ss = np.where((s4<=ipi[j])&(coords_t[j,6]==1)&(coords_t[j,4]==ri_s), ss+1, ss) #Count the super infections
                ns = np.where((s4<=ipi[j])&(coords_t[j,6]==1)&(coords_t[j,4]==ri_n), ns+1, ns) #Count the normal infections
                
                ## Update the event time of infection for the infected individual in the list t_i ##
                t_i[j] = np.where(s4<=ipi[j], t_s, t_i[j]) 
                
                ## Update the susceptibility label for the newly infected individual, as he is not susceptible anymore ##
                coords_t[j,6] = np.where((s4<=ipi[j])&(coords_t[j,6]==1), 0, coords_t[j,6]) 
                
                ## Update minimum infection time parameter ##
                t_im = np.min(t_i)
                
                ## Collect the data for the infected and infector, the event time and the infection rate ##
                hah = np.concatenate([hah, np.column_stack(np.array((c, j, float(t_s), float(coords_t[j,4])), dtype=float))], axis=0)
                #t_un = t_im+100
                
        ## Check if the number of infected individuals before the infection event (ib) is greater than or equal to the number of infected individuals (ia) after the infection ##
        ## If so, it means no new infections occurred, and "Nobody got infected" is printed. Otherwise, the indices of the newly infected individuals are printed as "These got infected:" using the hah data collected during the loop ##
        ia= np.sum(coords_t[:,2])
        if ib>=ia:
            print("Nobody got infected.")
        else:
            # The array hah[-int(ia-ib):,1] contains all the indices of the newly infected individuals 
            print("These individuals got infected: " + ", ".join(map(str, [int(x) for x in hah[-int(ia-ib):,1]])))
            
        ## Optional: print the number of infected individuals in each strain ##       
        # print("Super spreaders:",sum(coords_t[:, 5]==2) )
        # print("Normal infected:", sum(coords_t[:, 5]==1))
                
        ## Update the total rate of movement ##
        rt_m = sum(coords_t[:,3])   
        
        ## Update the total rate of infection ##
        rt_i = sum(coords_t[:,4])
        
        ## Remove the lines in the genome table that correspond to the genomes of recovered individuals ##
        g.iloc[coords_t[:, 2] == 0] = np.zeros((1, l)) 
        
        if args.super_strain:
            ## Collect the data for the number of Total infected, Super spreaders, Normal spreaders and infection time for each generation ##
            all_inf = np.concatenate([all_inf, np.column_stack(np.array((sum(coords_t[:,2]), sum(coords_t[:,5]==2), sum(coords_t[:,5]==1), float(t_s)), dtype=float))], axis=0)    
        else:
            ## Collect the data for the number of Total infected, Normal spreaders and infection time for each generation ##
            all_inf = np.concatenate([all_inf, np.column_stack(np.array((sum(coords_t[:,2]), sum(coords_t[:,5]==1), float(t_s)), dtype=float))], axis=0)    

    if args.super_strain:
        ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
        sample_data_normal_super(samples, genomes, g, tt, coords_t, all_inf, sample_times)
    else:
        ## Save a sample of data from the simulation, according to the conditions of the sample_data() function##
        sample_data_normal(samples, genomes, g, tt, coords_t, all_inf, sample_times)
    
    ## Time to run the simulation loop ##
    loop_t = time.time()-time_bfloop 
    
    tt += 1 #Moving on in the simulation loop
        
    #print("time for this loop:", loop_t)
    
## OPTIONAL: Plot a visual map with the final positions of the sample ##
if args.scatter_plots:
    if args.super_strain:
        plot_map_normal_super(plots, coords_t, tt, t_s, mv, ss, un)
    else:
        plot_map_normal(plots, coords_t, tt, t_s, mv, un)

if args.super_strain:
    ## Save the data from the simulation in the correct directory ##
    save_data_normal_super(samples, genomes, coords_2, coords_t, tt, g, all_inf, unin, hah, ss, ns, mv, t_un)
else:    ## Save the data from the simulation in the correct directory ##
    save_data_normal(samples, genomes, coords_2, coords_t, tt, g, all_inf, unin, hah, ns, mv, t_un)

#coords_t = pd.DataFrame(data=coords_t.T, index=("x","y", "label", "probability of movement", "probability of infection", "mutation", "susceptibility")).T
#coords_t.to_csv('/Users/katiagim/Desktop/MobiVirus/attempt_2/'+'final_data'+'.csv',header=True, index=False)
print("Loops:", tt)
## Time for the whole code to run ##
end=time.time()-initial
print("Total simulation time (minutes):", end/60)   

## Create a gif from all the plots of positions in order to visualise the movement (takes a lot of time) ##
#create_gif("/Users/katiagim/Desktop/MobiVirus/attempt_2/figs/")

# %%