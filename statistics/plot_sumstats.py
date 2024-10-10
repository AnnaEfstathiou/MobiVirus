import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os

from statistics import statistics, Fst, __event_number

def sumstats_directory(directory, sample_size, strain_type, sampling_type):

    ## Define directory and file names
    genome_dir = os.path.join(directory, 'genomes/')
    samples_dir = os.path.join(directory, 'samples/')
    event_time_filename = os.path.join(samples_dir, 'event_type.csv')
    
    ## Get the genome, sample and time files
    genome_files = {__event_number(f): f for f in os.listdir(genome_dir) if f.startswith('genomes') and f.endswith('.csv')}
    simulation_times = pd.read_csv(event_time_filename, usecols=['Event', 'Simulation Time'], index_col='Event').squeeze().to_dict() # Convert simulation time and events into a dictionary
    if sampling_type == "coordinates_sampling":
         coords_files = {__event_number(f): f for f in os.listdir(samples_dir) if f.startswith('coords') and f.endswith('.csv')}
    
    ## Initialize directories to store statistics ##
    tajimasd_stats = {} # Initialize a dictionary to store Tajima's D values for every simulation time
    pi_stats = {}       # Initialize a dictionary to store Pi values for every simulation time
    thetaw_stats = {}   # Initialize a dictionary to store θw values for every simulation time
    uniq_stats = {}     # Initialize a dictionary to store number of unique seqs for every simulation time
    hp_div = {}         # Initialize a dictionary to store haplotype diversity for every simulation time
    fst_stats = {}      # Initialize a dictionary to store Fst values for every simulation time  
            
    ## Process matching files by suffix, to get corresponding generations
    for suffix, genome_file in genome_files.items():
        genome_data = os.path.join(genome_dir, genome_file)
        
        try:
            print(f"Processing file: {genome_file}")
            
            tajimasd, pi, thetaw, uniq_haplo, haplo_div = statistics(genome_data, strain_type, sample_size)
            time = simulation_times[int(suffix)] # Get the simulation time the current event happened

            if (sampling_type == "coordinates_sampling") and (suffix in coords_files):
                coords_file = coords_files[suffix]
                coords_data = os.path.join(samples_dir, coords_file)
                fst = Fst(genome_data, sampling_type, strain_type, sample_size, coords_file = coords_data)
            else:
                 fst = Fst(genome_data, sampling_type, strain_type, sample_size, coords_file = None)
            
            tajimasd_stats[time] = tajimasd
            pi_stats[time] = pi
            thetaw_stats[time] = thetaw
            uniq_stats[time] = uniq_haplo
            hp_div[time] = haplo_div
            fst_stats[time] = {'hsm': fst.hsm(), 'slatkin': fst.slatkin(), 'hbk': fst.hbk()}
                  
        except ValueError as e:
                    print(f"Error processing {genome_file}: {e}")
                    continue
    # Sort dictionaries according to simulation time (keys)
    tajimasd_stats = dict(sorted(tajimasd_stats.items()))
    pi_stats       = dict(sorted(pi_stats.items()))
    thetaw_stats   = dict(sorted(thetaw_stats.items()))
    uniq_stats     = dict(sorted(uniq_stats.items()))
    hp_div         = dict(sorted(hp_div.items()))
    fst_stats      = dict(sorted(fst_stats.items()))

    # Convert to DataFrame
    data = {
        'Time': list(tajimasd_stats.keys()),
        'TajimasD': list(tajimasd_stats.values()),
        'Pi': list(pi_stats.values()),
        'Theta_w': list(thetaw_stats.values()),
        'Unique_seqs': list(uniq_stats.values()),
        'Hp_Diversity': list(hp_div.values()),
        'Fst_hsm': [fst['hsm'] for fst in fst_stats.values()],
        'Fst_slatkin': [fst['slatkin'] for fst in fst_stats.values()],
        'Fst_hbk': [fst['hbk'] for fst in fst_stats.values()]
}

    stats = pd.DataFrame(data)

    # return tajimasd_stats, pi_stats, thetaw_stats, uniq_stats, hp_div, fst_stats
    return stats

# def plotting(dict_tajimasd_stats, dict_pi_stats, dict_thetaw_stats, dict_fst_stats, uniq_stats, hp_div, save_png = None):

#     # Create 6 subplots
#     fig, ax = plt.subplots(5, 1, figsize=(12, 28)) 

#     # Plot Tajima's D
#     ax[0].plot(list(dict_tajimasd_stats.keys()), list(dict_tajimasd_stats.values()), marker='.', color='green')
#     ax[0].set_title("Tajima's D", fontsize=12)
#     ax[0].set_ylabel('Scores')
#     ax[0].set_xlabel('Simulation Time')
#     ax[0].grid(True) 

#     # Plot Pi & θw scores
#     ax[1].plot(list(dict_pi_stats.keys()), list(dict_pi_stats.values()), marker='.', color='orange', label = "Pi-Estimator")
#     ax[1].plot(list(dict_thetaw_stats.keys()), list(dict_thetaw_stats.values()), marker='.', color='blue', label = "Watterson-Estimator")
#     ax[1].set_title("Pi and Theta Waterson Estimators", fontsize=12)
#     ax[1].set_ylabel('Scores')
#     ax[1].set_xlabel('Simulation Time')
#     ax[1].grid(True)
#     ax[1].legend()

#     # Plot Fst scores
#     simulation_times = list(dict_fst_stats.keys())
#     hsm_values = [dict_fst_stats[time]['hsm'] for time in simulation_times]
#     slatkin_values = [dict_fst_stats[time]['slatkin'] for time in simulation_times]
#     hbk_values = [dict_fst_stats[time]['hbk'] for time in simulation_times]
#     ax[2].plot(simulation_times, hsm_values, marker='.', color='orangered', label='HSM')
#     ax[2].plot(simulation_times, slatkin_values, marker='.', color='limegreen', label='Slatkin')
#     ax[2].plot(simulation_times, hbk_values, marker='.', color='dimgrey', label='HBK')
#     ax[2].set_title("Fst", fontsize=12)
#     ax[2].set_ylabel('Scores')
#     ax[2].set_xlabel('Simulation Time')
#     ax[2].grid(True)
#     ax[2].legend() 

#     # Plot Haplotype Diversity
#     ax[3].plot(list(hp_div.keys()), list(hp_div.values()), marker='.', color='steelblue')
#     ax[3].set_title("Haplotype Diversity", fontsize=12)
#     ax[3].set_ylabel('Scores')
#     ax[3].set_xlabel('Simulation Time')
#     ax[3].grid(True)

#     # Plot number of unique sequences
#     ax[4].plot(list(uniq_stats.keys()), list(uniq_stats.values()), marker='.', color='rebeccapurple')
#     ax[4].set_title("Unique Sequences", fontsize=12)
#     ax[4].set_ylabel('Number')
#     ax[4].set_xlabel('Simulation Time')
#     ax[4].grid(True)

#     plt.tight_layout(pad=5.0)
#     plt.subplots_adjust(hspace=0.5)  # Manually set height space between plots
#     if save_png:
#         plt.savefig(save_png, format="png")
#     else:
#         plt.show()
#     plt.close(fig)

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-d', '--directory', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-sample', '--sample_size', type=int, required=True, help='Sample size.')
    parser.add_argument('-sample_type', '--sampling_technique', type=str, required=True, help='Define the way the sampling of 2 populations will happen. 3 ways: "rdm","str","coords"')
    parser.add_argument('-save', '--save_png', action="store_true", help='Save the plot in a PNG format.')
    # Adding mutually exclusive group for strain type options
    strain_group = parser.add_mutually_exclusive_group(required=True)
    strain_group.add_argument('-mix', '--ss_and_ns_strains', action="store_true", help='Calculate statistics for both strains.')
    strain_group.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    strain_group.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    args = parser.parse_args()

    

    # Determine strain_type based on command-line arguments
    if args.ss_and_ns_strains:
        strain_type = "mix_strains"
    elif args.ss_strains:
        strain_type = "ss_strains"
    elif args.ns_strains:
        strain_type = "ns_strains"
    else:
        strain_type = None

    if args.sampling_technique == "rdm":
        sampling_type = "simple_random_sampling"
    elif args.sampling_technique == "str":
        sampling_type = "stratified_sampling"
    elif args.sampling_technique == "coords":
        sampling_type = "coordinates_sampling"
    else:
        sampling_type = None

    # Ensure strain_type is set properly
    if strain_type is None:
        print("Error: Please specify a strain type to calculate statistics (use -mix, -ss, or -ns).")
    else:
        # tajimasd_stats, pi_stats, thetaw_stats, uniq_stats, hp_div, fst_stats = sumstats_directory(args.directory, args.sample_size, strain_type, sampling_type)
        stats = sumstats_directory(args.directory, args.sample_size, strain_type, sampling_type)

        directory = os.path.normpath(args.directory)
        directory_name = os.path.basename(directory)
        datetime_part = "_".join(directory_name.split('_')[1:])
        csv_filename = f"stats_{datetime_part}.csv" # Form the CSV filename using the extracted datetime part of the simulation directory
        stats.to_csv(csv_filename) # Save the stats to a CSV file
        
        # if args.save_png:
        #     # Determine png title based on command-line arguments
        #     if args.ss_and_ns_strains:
        #         png_title = f'sim_stats.png'
        #     elif args.ss_strains:
        #         png_title = f'sim_stats_ss.png'
        #     elif args.ns_strains:
        #         png_title = f'sim_stats_ns.png'
        #     plotting(tajimasd_stats, pi_stats, thetaw_stats, fst_stats, uniq_stats, hp_div, save_png = png_title)
        # else:
        #     plotting(tajimasd_stats, pi_stats, thetaw_stats, fst_stats, uniq_stats, hp_div, save_png = None)