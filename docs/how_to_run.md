# Parameter Overview 

## INI file Parameters

| Parameter      | Type    | Description                                                           | Constraints        |
| -------------- | ------- | --------------------------------------------------------------------- | ------------------ |
| `directory`    | string  | Directory where the simulation results are stored                     | No constraints     |
| `super_strain` | boolean | True: super strains exist (otherwise only normal strains exist)       | true / false       |
| `n`            | int     | Number of individuals in the simulation                               | `n > 0`            |
| `ii`           | int     | Initial number of infected individuals                                | `0 <= ii <= n`     |
| `l`            | int     | Length of genome                                                      | `l > 0`            |
| `bound_l`      | float   | Lower bound for spatial axis                                          | No constraints     |
| `bound_h`      | float   | Upper bound for spatial axis                                          | No constraints     |
| `rm_i`         | float   | Movement rate of infected individuals                                 | `rm_i >= 0`        |
| `rm_h`         | float   | Movement rate of healthy individuals                                  | `rm_h >= 0`        |
| `r_m`          | float   | Mutation rate per site                                                | `0 <= r_m <= 1`    |
| `r_rec`        | float   | Recombination rate per site                                           | `0 <= r_rec <= 1`  |
| `sample_times` | int     | Period of events for simulation output sampling                       | `sample_times > 0` |
| `ri_n`         | float   | Infection rate (normal strain)                                        | `ri_n > 0`         |
| `inf_dist_ns`  | float   | Infection distance (normal strain)                                    | `inf_dist_ns > 0`  |
| `prob_inf_ns`  | float   | Probability of infection (normal strain)                              | `0 <= prob_inf_ns <= 1` |
| `rec_t_ns`     | float   | Recovery time (normal strain)                                         | `rec_t_ns >= 0`    |
| `imm_t_ns`     | float   | Immunity time (normal strain)                                         | `imm_t_ns >= 0`    |
| `rim_ns`       | float   | Relative infected mobility (normal strain)                            | `rim_ns >= 0`      |
| `n_i`          | int     | Number of sites that form the beneficial region of the super strain   | `0 <= n_i <= l` & if `super_strain = false` then  `n_i = 0` |
| `region_pos`   | string  | Location of the beneficial region in the viral genome                 | `start`, `middle`, `end` |
| `ri_s`         | float   | Infection rate (super strain)                                         | `ri_s > 0`         |
| `inf_dist_ss`  | float   | Infection distance (super strain)                                     | `inf_dist_ss > 0`  |
| `prob_inf_ss`  | float   | Probability of infection (super strain)                               | `0 <= prob_inf_ss <= 1` |
| `rec_t_ss`     | float   | Recovery time (super strain)                                          | `rec_t_ss >= 0`    |
| `imm_t_ss`     | float   | Immunity time (super strain)                                          | `imm_t_ss >= 0`    |
| `rim_ss`       | float   | Relative infected mobility (super strain)                             | `rim_ss >= 0`      |
| `ss_formation` | int     | Event when super strain mutation is introduced                        | `ss_formation >= 0`|
| `ss_events`    | int     | Period of events for super strain introduction                        | `ss_events >= 0`   |
| `use_demography` | boolean | True: simulate with population demography using msprime             | true / false       |
| `initial_pop`  | int     | Present-day population size                                           | `initial_pop >= 0` |
| `past_pop`     | int     | Past population size                                                  | `past_pop >= 0`    |
| `t_gen`        | int     | Time (in generations)                                                 | `t_gen >= 0`       |

## Parsers

#### String Parsers

| Argument                        | Type  | Description                                                                | Constraints       |
| ------------------------------- | ----- | -------------------------------------------------------------------------- | ----------------- |
| `--ratio_super_vs_normal`       | float | Stops simulation when ratio of super strain to normal strain reaches value | `0 <= value <= 1` |
| `--percentage_infected`         | float | Stops simulation when infected population reaches given percentage         | `0 <= value <= 1` |
| `--percentage_super_strain`     | float | Stops simulation when super spreaders reach given percentage               | `0 <= value <= 1` |
| `--percentage_normal_strain`    | float | Stops simulation when normal spreaders reach given percentage              | `0 <= value <= 1` |
| `--max_infections`              | int   | Stops simulation after given number of infections                          | `value >= 1`      |
| `--max_movements`               | int   | Stops simulation after given number of movements                           | `value >= 1`      |
| `--end_time`                    | float | Stops simulation at given simulation time                                  | `value > 0`       |
| `--end_events`                  | int   | Stops simulation after given number of events                              | `value >= 1`      |
| `--fixation_event`              | float | Sample when the normal and super strains will have reached the assigned percentages. 1st value: normal & 2nd value: super strains.| `0 <= value <= 1`      |

#### Action Parsers (Boolean Flags)

| Argument              | Type | Description                                                                 |
| --------------------- | ---- | --------------------------------------------------------------------------- |
| `--all_infected_once` | bool | Stops simulation when all individuals are infected at least once            |
| `--manual_genomes`    | bool | Starts simulation with all-0 genomes, introducing super strain at `ss_form` event if `super_strain = true` |
| `--msprime_genomes`   | bool | Starts simulation with genomes created by `msprime`, introducing super strain at `ss_form` event if `super_strain = true` using the msprime parameters |
| `--recombination`     | bool | Enables genome recombination if two infections occur                         |
| `--sample_infection`  | bool | Ignores `.INI` `sample_times` and samples every infection event               |
| `--visualize_data`    | bool | Displays data table as a dataframe in console                                |
| `--initial_genomes`   | bool | Saves the initial genomes of the population to a CSV format in the sample directory  |


> If the user does not specify any of the parsed arguments that stop the simulation under certain conditions, it will stop, by default, when all individuals become healthy. However, this can be problematic, as many parameter combinations may prevent this outcome, causing the simulation to run indefinitely!

# Simulation’s Output Data

The results are stored in one directory, which is given in the INI file.  
Directory name: *simulation_timestamp* (e.g., *simulation_13_03_2026_21_33*)

In addition to generating files containing simulation data, a **log file** (command_log.txt) is created to record all parameters and the exact command used.  
This ensures that users can track the conditions under which the results were produced and easily reproduce the simulation with the same parameters later.

## Sampled Data
Generated every `sample_times` events & at the end of simulation.

| File | Format | Description | Info (columns) |
| ---- | ------ | ----------- | --------------- |
| `genomes_sample_times.csv` <br /> (e.g., `genomes_100.csv`) | CSV | Viral genome sequences | - Dimensions: `(n x l)`<br>- `n`: number of individuals<br>- `l`: genome length<br>- Binary genomes: 0 = ancestral state, 1 = mutation<br>- Healthy individuals appear as "NaN" |
| `coords_sample_times.csv` <br /> (e.g., `coords_100.csv`) | CSV | Spatial coordinates and infection-related properties of all individuals | - `x`, `y`: spatial coordinates<br>- `Infection label`: 0 = healthy, 1 = infected, 2 = recombination-infected (if allowed)<br>- `Rate of movement`: `rm_i` if infected, `rm_h` if healthy<br>- `Rate of infection`: `ri_n` if normal, `ri_s` if super, 0 if healthy<br>- `Mutation`: 0 = healthy, 1 = normal, 2 = super<br>- `Susceptibility`: 0 = non-susceptible, 1 = susceptible<br>- Rows correspond to individuals |

## Final Data
Generated at the end of simulation.

| File | Format | Description | Info (columns) |
| ---- | ------ | ----------- | --------------- |
| `all_inf.csv` | CSV | Number of infected individuals during the simulation | - `Total infected`: Sum of super and normal spreaders<br>- `Super spreaders`: Sum of individuals carriying a super strain<br>- `Normal spreaders`: Sum of individuals carriying a normal strain<br>- `Time`: Corresponding simulation time<br>- `Event`: Corresponding event number |
| `recovery.csv` | CSV | List of recovered individuals and their recovery times | - `Recovered individual`: Index of the corresponding individual<br>- `Recovery Time`: Simulation time when recovery happened |
| `infections.csv` | CSV | Log of all infection events | - `Infecting`: Index of the infector<br>- `Infected`: Index of the infected<br>- `Infection Time`: 1st infection = 1.0, 2nd via recombination = 2.0<br>- `Simulation Time`: Simulation time of infection<br>- `Infection Rate`: `ri_n` if normal, `ri_s` if super spreader |
| `initial_coords.csv` | CSV | Initial spatial coordinates and infection-related properties of all individuals before simulation starts | - `x`, `y`: Spatial coordinates<br>- `Infection label`: 0 = healthy, 1 = infected, 2 = recombination-infected (if allowed)<br>- `Rate of movement`: `rm_i` if infected, `rm_h` if healthy<br>- `Rate of infection`: `ri_n` if normal, `ri_s` if super, 0 if healthy<br>- `Mutation`: 0 = healthy, 1 = normal, 2 = super<br>- `Susceptibility`: 0 = non-susceptible, 1 = susceptible<br>- Rows correspond to individuals |
| `event_type.csv` | CSV | Log of all simulation events | - `Event`: Event number/index (0 to final event)<br>- `Simulation Time`: Time of event<br>- `Event Type`: “Movement” or “Infection”<br>- `Individual`: Index of the moving individual (if movement) or index of infector (if infection) |

# Code Execution – Examples

The MobiVirus Simulator requires three files in the same directory:  
1. `MobiVirus_Simulator.py` – main simulation code  
2. `mobifunctions.py` – functions used by the simulator  
3. `parameters.ini` – simulation initial conditions  

Programming Language: Python  
All required packages are listed in `requirements.txt`.

**Examples**

A simulation that starts with viral genomes created by msprime, allows viral recombination, and stops at 100 events.

```python
python3 MobiVirus_Simulator.py -msprime -r -events 1000
```

A simulation that starts with viral genomes created by msprime and stops when super spreaders reach 95% of the population.

```python
python3 MobiVirus_Simulator.py -msprime -per_ss 0.95
```

A simulation that contains that starts with all-0 viral genomes, allows viral recombination, stops when all individuals are infected at least once and saves the initial genomes of the population.

```python
python3 MobiVirus_Simulator.py -manual -r -all_inf -g0
```

A simulation that starts with viral genomes created by msprime, allows viral recombination, stops when 500 infections take place and samples, extra, when normal strains disappear and super strains reach 90% of the population.

```python
python3 MobiVirus_Simulator.py -msprime -r -max_inf 500 -fixation 0 0.9
```