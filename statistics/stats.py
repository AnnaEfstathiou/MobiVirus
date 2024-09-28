import pandas as pd
import libsequence
import matplotlib.pyplot as plt
import argparse
import pandas as pd
from datetime import datetime
import re
from ld_plot.ld_plot import ld_plot

## all is refered to mix sequences aka super and normal strains  ##
## ss is refered to super strains                                ##
## ns is refered to normal strains                               ##

#%% Internal Functions

def event_number(csv_file):

    ## Isolate suffix number from the csv file ##
    ## Suffix represents the event number      ##

    match = re.search(r'_(\d+)\.csv$', csv_file)
    if match:
        number = match.group(1)
    else:
        print("No valid number found in the file name.")

    return number

def filter_middle_column(df):
    
    ## Create a dataframe with the sequences that have a SNP in the middle position ##

    length = len(df.columns)
    
    if length % 2 == 1:
        # Odd length
        middle_column_index = length // 2
        middle_column_name = df.columns[middle_column_index]
        # Filter rows where the value in the middle column is 1
        filtered_df = df[df[middle_column_name] == 1]
    else:
        # Even length
        middle_left = (length // 2) - 1
        middle_right = length // 2
        # Get the two middle column names
        middle_left_name = df.columns[middle_left]
        middle_right_name = df.columns[middle_right]
        # Filter rows where the value in EITHER middle column is 1
        filtered_df = df[(df[middle_left_name] == 1) | (df[middle_right_name] == 1)]

    return filtered_df

def simulated_data(csv_file, sample):

    ## Read the dataset ##
    all_seqs_data = pd.read_csv(csv_file, header=None) # CSV file --> DataFrame
    all_seqs_data = all_seqs_data.dropna()             # Drop rows with NaN values
    all_seqs_data = all_seqs_data.astype(int)          # Turn values to integers
    l = len(all_seqs_data.columns)                     # Genome length

    ## Split the dataset ##
    ss_seqs_data = filter_middle_column(all_seqs_data)    # super strains
    ns_seqs_data = all_seqs_data.drop(ss_seqs_data.index) # normal strains

    n = len(all_seqs_data)   # Total number of all sequences
    n_ss = len(ss_seqs_data) # Total number of super strains
    n_ns = len(ns_seqs_data) # Total number of normal strains

    ## Sampling ##
    all_sampled = all_seqs_data.sample(n=sample, random_state=int(datetime.now().timestamp())) # sampling from all sequences
    ss_sampled  = ss_seqs_data.sample(n=sample, random_state=int(datetime.now().timestamp()))  # sampling from all super strains
    ns_sampled  = ns_seqs_data.sample(n=sample, random_state=int(datetime.now().timestamp()))  # sampling from all normal strains

    ## Polymorphic sites ##
    # All sequences
    all_polym_sites = all_sampled.loc[:, (all_sampled != 0).any(axis=0)]
    all_pos  = [(col / l) for col in list(all_polym_sites.columns)] 
    all_seqs = [''.join(map(str, row)) for row in all_polym_sites.values]
    # All super strains
    ss_polym_sites = ss_sampled.loc[:, (ss_sampled != 0).any(axis=0)]
    ss_pos  = [(col / l) for col in list(ss_polym_sites.columns)] 
    ss_seqs = [''.join(map(str, row)) for row in ss_polym_sites.values]
    # All normal strains
    ns_polym_sites = ns_sampled.loc[:, (ns_sampled != 0).any(axis=0)]
    ns_pos  = [(col / l) for col in list(ns_polym_sites.columns)] 
    ns_seqs = [''.join(map(str, row)) for row in ns_polym_sites.values]

    ## Simulate Data for further analysis ##
    # All sequences
    all_sd = libsequence.SimData()
    all_sd.assign(all_pos, all_seqs)
    # All super strains
    ss_sd = libsequence.SimData()
    ss_sd.assign(ss_pos, ss_seqs)
    # All normal strains
    ns_sd = libsequence.SimData()
    ns_sd.assign(ns_pos, ns_seqs)
    
    return all_sd, ss_sd, ns_sd, l, n, n_ss, n_ns

#%% Statistics Functions

def summary_stats(sd):

    ## Calculate Tajima's D, Pi & θw ##

    ps = libsequence.PolySIM(sd) 
    tajimasd = ps.tajimasd()
    pi = ps.thetapi()
    thetaw = ps.thetaw()
    
    return tajimasd, pi, thetaw

def selectice_sweep(sd, ws, step):

    ## Calculate Tajima's D, Pi & θw for every sliding window ##

    tajimas_d = {}
    pi_est    = {}
    theta_w   = {}    
    w = libsequence.Windows(sd, window_size=ws, step_len=step, starting_pos=0., ending_pos=1.0)
    for i in range(len(w)):
        wi = w[i]
        pswi = libsequence.PolySIM(wi)
        tajimas_d[i] = pswi.tajimasd()
        pi_est[i]    = pswi.thetapi()
        theta_w[i]   = pswi.thetaw()

    return tajimas_d, pi_est, theta_w

#%% Plotting Functions

def plot_selectice_sweep(tajimas_d_dict, pi_est_dict, theta_w_dict, window_size, step_size, seq_length, num_seqs, custom_title=None, save_png = None):

    # Extract the keys (x values) and values (y values)
    window_pos = list(theta_w_dict.keys())
    tajimas_d_values = list(tajimas_d_dict.values())
    pi_est_values = list(pi_est_dict.values())
    theta_w_values = list(theta_w_dict.values())
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(window_pos, tajimas_d_values, marker='.', linestyle='solid', color='darkmagenta', label='Tajima s D')
    plt.plot(window_pos, pi_est_values, marker='.', linestyle='solid', color='darkslategrey', label='Pi Estimator')
    plt.plot(window_pos, theta_w_values, marker='.', linestyle='solid', color='orangered', label='Watterson Estimator')

    plt.suptitle('Selective Sweep', fontsize=15, fontweight='bold')
    plt.xlabel("Window position")
    plt.ylabel("Estimator value")
    if custom_title:
        plt.title(custom_title, loc='left')
    plt.title(f'Window size:{int(window_size*seq_length)}, Step size:{int(step_size*seq_length)}, l:{seq_length}, seqs:{num_seqs}', loc='right')
    plt.legend()
    plt.grid(True)
    ## Save or display the plot ##
    if save_png:
        plt.savefig(save_png, format="png")
    else:
        plt.show()
    plt.close()
  
def plot_ld(LD, custom_title=None, save_png = None):

    # Creating the pivot tables
    rsq_ld_table = LD.pivot(index='i', columns='j', values='rsq')
    d_table = LD.pivot(index='i', columns='j', values='D')
    dprime_table = LD.pivot(index='i', columns='j', values='Dprime')

    # Calculating means
    meanrsq = LD['rsq'].mean()
    meanD = LD['D'].mean()
    meanDprime = LD['Dprime'].mean()

    # Creating subplots
    fig, ax = plt.subplots(3, 1, figsize=(16, 24)) 
    
    # Adding a general title for the figure
    fig.suptitle('Linkage Disequilibrium Heatmaps', fontsize=15, fontweight='bold')

    # Plotting r^2 heatmap
    cax0 = ax[0].imshow(rsq_ld_table, cmap="viridis", aspect='auto', origin='lower')
    fig.colorbar(cax0, ax=ax[0], label='$r^2$') 
    ax[0].set_title('Heatmap of $r^2$ Values')
    if custom_title:
        ax[0].set_title(custom_title, loc='left')
    ax[0].set_xlabel('position j')
    ax[0].set_ylabel('position i')
    ax[0].text(0.02, 0.98, f'Mean $r^2$: {meanrsq:.2f}', 
               verticalalignment='top', horizontalalignment='left',
               transform=ax[0].transAxes, color='black', fontsize=12)

    # Plotting D heatmap
    cax1 = ax[1].imshow(d_table, cmap="viridis", aspect='auto', origin='lower')
    fig.colorbar(cax1, ax=ax[1], label='D') 
    ax[1].set_title('Heatmap of D Values')
    if custom_title:
        ax[1].set_title(custom_title, loc='left')
    ax[1].set_xlabel('position j')
    ax[1].set_ylabel('position i')
    ax[1].text(0.02, 0.98, f'Mean D: {meanD:.2f}', 
               verticalalignment='top', horizontalalignment='left',
               transform=ax[1].transAxes, color='black', fontsize=12)

    # Plotting D' heatmap
    cax2 = ax[2].imshow(dprime_table, cmap="viridis", aspect='auto', origin='lower')
    fig.colorbar(cax2, ax=ax[2], label="D'")  
    ax[2].set_title("Heatmap of D' Values")
    if custom_title:
        ax[2].set_title(custom_title, loc='left')
    ax[2].set_xlabel('position j')
    ax[2].set_ylabel('position i')
    ax[2].text(0.02, 0.98, f"Mean D': {meanDprime:.2f}", 
               verticalalignment='top', horizontalalignment='left',
               transform=ax[2].transAxes, color='black', fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.94])  
    plt.subplots_adjust(top=0.92, hspace=0.4)
    ## Save or display the plot ##
    if save_png:
        plt.savefig(save_png, format="png")
    else:
        plt.show()
    plt.close()

#%% Main
if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a genome file (csv format) to compute various statistics and plot results.")
    parser.add_argument('-g', '--genome_file', type=str, required=True, help='The path to the input CSV file.')
    parser.add_argument('-sample', '--sample_size', type=int, help='Sample size.')
    parser.add_argument('-sumstats', '--summary_statistics', action="store_true", help='Calculate summary statistics.')
    parser.add_argument('-ld', '--linkage_disequilibrium', action="store_true", help='Calculate Linkage Disequilibrium.')
    parser.add_argument('-fst', '--fst_statistics', action="store_true", help='Calculate Fst.')
    parser.add_argument('-selsw', '--selective_sweep', action="store_true", help='Search for Selective Sweep.')
    parser.add_argument('-winsize', '--window_size', type=float, help='Length of sliding window.')
    parser.add_argument('-stepsize', '--step_size', type=float, help='Step size.')
    parser.add_argument('-mix', '--ss_and_ns_strains', action="store_true", help='Calculate statistics for both strains.')
    parser.add_argument('-ss', '--ss_strains', action="store_true", help='Calculate statistics for super strains.')
    parser.add_argument('-ns', '--ns_strains', action="store_true", help='Calculate statistics for normal strains.')
    parser.add_argument('-p', '--population_size', action="store_true", help='Population size.')
    parser.add_argument('-save', '--save_png', action="store_true", help='Save the plot in a PNG format.')
    args = parser.parse_args()


    """ Data manipulation """
    suffix = event_number(args.genome_file)
    all_sd, ss_sd, ns_sd, l, n, n_ss, n_ns = simulated_data(args.genome_file, args.sample_size)
    
    if args.population_size:
        print(f"Number of sequences (normal & super strains): {n}")
        print(f"Number of super sequences: {n_ss}")
        print(f"Number of normal strains: {n_ns}")    
    
    """ Statistics """
    ## Summary Statistics ##
    if args.summary_statistics:

        if args.ss_and_ns_strains:
            all_tajimasd_stats, all_pi_stats, all_thetaw_stats = summary_stats(all_sd) # Total sequences
            print(f'Normal & Super Strains --> Tajimas D: {all_tajimasd_stats}, Pi-estimator: {all_pi_stats}, Theta Watterson: {all_thetaw_stats}')
        if args.ss_strains:
            ss_tajimasd_stats,  ss_pi_stats,  ss_thetaw_stats  = summary_stats(ss_sd)  # Super strains
            print(f'Super Strains --> Tajimas D: {ss_tajimasd_stats}, Pi-estimator: {ss_pi_stats}, Theta Watterson: {ss_thetaw_stats}')
        if args.ns_strains:
            ns_tajimasd_stats,  ns_pi_stats,  ns_thetaw_stats  = summary_stats(ns_sd)  # Normal strains
            print(f'Normal Strains --> Tajimas D: {ns_tajimasd_stats}, Pi-estimator: {ns_pi_stats}, Theta Watterson: {ns_thetaw_stats}')

    ## LD ##
    if args.linkage_disequilibrium:

        if args.ss_and_ns_strains:
            all_LDstats = pd.DataFrame(libsequence.ld(all_sd)) # Total sequences
            if args.save_png:
                plot_ld(all_LDstats, custom_title='Normal & Super strains', save_png = f'ld_{suffix}.png')
            else:
                plot_ld(all_LDstats, custom_title='Normal & Super strains', save_png = None)
        if args.ss_strains:
            ss_LDstats  = pd.DataFrame(libsequence.ld(ss_sd))  # Super strains
            if args.save_png:
                plot_ld(ss_LDstats,  custom_title='Super strains',  save_png = f'ld_ss_{suffix}.png')
            else:
                plot_ld(ss_LDstats,  custom_title='Super strains',  save_png = None)
        if args.ns_strains:
            ns_LDstats  = pd.DataFrame(libsequence.ld(ns_sd))  # Normal strains
            if args.save_png:
                plot_ld(ns_LDstats,  custom_title='Normal strains', save_png = f'ld_ns_{suffix}.png')
            else:
                plot_ld(ns_LDstats,  custom_title='Normal strains', save_png = None)

    ## Selective Sweep ##
    if args.selective_sweep:

        """ Checks """
        if args.window_size is None or args.step_size is None:
            parser.error("The arguments --window_size and --step_size are required when --selective_sweep is used.")
        if not (0 < args.window_size <= 1):
            raise ValueError("Window size must be a float between 0 and 1.")
        if not (0 < args.step_size <= 1):
            raise ValueError("Step size must be a float between 0 and 1.")
        
        if args.ss_and_ns_strains:
            all_tajimasd_selsw, all_pi_selsw, all_thetaw_selsw = selectice_sweep(all_sd, args.window_size, args.step_size) # Total sequences
            if args.save_png:
                plot_selectice_sweep(all_tajimasd_selsw, all_pi_selsw, all_thetaw_selsw, args.window_size, args.step_size, l, len(all_sd), custom_title='Normal & Super strains', save_png = f'selsw_{suffix}.png')
            else:
                plot_selectice_sweep(all_tajimasd_selsw, all_pi_selsw, all_thetaw_selsw, args.window_size, args.step_size, l, len(all_sd), custom_title='Normal & Super strains', save_png = None) 
        if args.ss_strains:
            ss_tajimasd_selsw,  ss_pi_selsw,  ss_thetaw_selsw  = selectice_sweep(ss_sd, args.window_size, args.step_size)  # Super strains
            if args.save_png:
                plot_selectice_sweep(ss_tajimasd_selsw,  ss_pi_selsw,  ss_thetaw_selsw,  args.window_size, args.step_size, l, len(ss_sd),  custom_title='Super strains',  save_png = f'selsw_ss_{suffix}.png')
            else:
                plot_selectice_sweep(ss_tajimasd_selsw,  ss_pi_selsw,  ss_thetaw_selsw,  args.window_size, args.step_size, l, len(ss_sd),  custom_title='Super strains',  save_png = None)
        if args.ns_strains:
            ns_tajimasd_selsw,  ns_pi_selsw,  ns_thetaw_selsw  = selectice_sweep(ns_sd, args.window_size, args.step_size)  # Normal strains
            if args.save_png:
                plot_selectice_sweep(ns_tajimasd_selsw,  ns_pi_selsw,  ns_thetaw_selsw,  args.window_size, args.step_size, l, len(ns_sd),  custom_title='Normal strains', save_png = f'selsw_ns_{suffix}.png')
            else:
                plot_selectice_sweep(ns_tajimasd_selsw,  ns_pi_selsw,  ns_thetaw_selsw,  args.window_size, args.step_size, l, len(ns_sd),  custom_title='Normal strains', save_png = None)

    ## Fst ##
    if args.fst_statistics:
        if args.ss_and_ns_strains:
            all_fst = libsequence.Fst(all_sd,[500,500])
            print(f'Normal & Super Strains --> {all_fst.hsm()}, {all_fst.slatkin()}, {all_fst.hbk()}')
        if args.ss_strains:
            ss_fst = libsequence.Fst(ss_sd,[500,500])
            print(f'Super Strains --> {ss_fst.hsm()}, {ss_fst.slatkin()}, {ss_fst.hbk()}')
        if args.ns_strains:
            ns_fst = libsequence.Fst(ns_sd,[500,500])
            print(f'Normal Strains --> {ns_fst.hsm()}, {ns_fst.slatkin()}, {ns_fst.hbk()}')