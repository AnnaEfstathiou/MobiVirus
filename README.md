# MobiVirus Simulator
MobiVirus is a <ins>Python-based simulator</ins> that allows users to configure epidemiological parameters and simulate disease dynamics over time.

It is a forward-in-time, stochastic, discrete-event simulator (DES) designed to study the evolution of airborne viral genomes and transmission dynamics in populations that span along spatial coordinates.
The model consists of the host and the viral populations.
Hosts live and move within a 2D rectangular space.
The interaction between hosts and viral strains occurs through infection mechanisms and follows a Susceptible-Infected-Recovered-Susceptible (SIRS) epidemiological model, where recovered hosts lose immunity over time and become susceptible again.
The viral population may consist of one or two distinct strains, a normal strain and a super strain with increasing fitness. 

## Repository Structure

```
MobiVirus/
│
├── src/
│   ├── MobiVirus_Simulator.py
│   └── mobifunctions.py
│   └── parameters.ini
|
├── docs/
|   ├── how_to_run.md
│
├── example/
|   ├── example_run.sh
|   ├── example_parameters.ini
|   ├── example_terminal_output.txt 
│   └── example_simulation/
|       └── command_log.txt
|       └── genomes/
|       └── samples/
|
├── README.md
├── .gitignore
├── requirements.txt
└── LICENSE
```

## Installation

Clone the repository:

```
git clone https://github.com/AnnaEfstathiou/MobiVirus.git
cd MobiVirus
```

Install required Python packages using conda:

```
conda create --name mobivirus --file requirements.txt   # The requirements.txt is a conda environment export.
conda activate mobivirus
```

## How to Run

The MobiVirus Simulator works with the use of 3 files, all located in the `src/` directory:
1. `MobiVirus_Simulator.py`: Main code where the simulation runs.
2. `mobifunctions.py`: Python script containing the functions used in the simulation.
3. `parameters.ini`: INI file containing the initial conditions and parameters of the simulation.

All files must be placed in the same directory for the simulation to run correctly. The chosen directory must also be added to the INI file.

Run the simulator with the default parameters in the INI file:

```
# General format to run: python3 MobiVirus_Simulator.py -flags
# e.g.
# Run the simulation for 100 events and generate the initial genomes using the msprime simulator
python3 MobiVirus_Simulator.py -events 100 -msprime 
```

## Output

An example simulation output is shown below:

```
example/example_simulation/
```

The simulator periodically generates output every E events, where E is a user-defined interval in the INI file under the variable name `sample_times`. 

Each sampling produces two major categories of data:
1. **Viral Genomes**: Complete viral genomes for all currently infected individuals. 
2. **Host Population Data**: For every host in the simulation, the following details are recorded:
    - Epidemiological state (e.g., susceptible, infected, recovered)
    - Movement and infection rate
    - Spatial coordinates
    - Viral strain type (if infected)

When the simulation ends, a final dataset is produced containing all of the above plus additional summary information:
- **Recovery Details**: which individuals recovered and when
- **Infection Details**: transmission history
- **Viral Strain Summary**: population size of each viral strain throughout the simulation
- **Event Breakdown**: event types throughout the simulation

## Author

Lead developer: Efstathiou Anna <br />
Original code base: Katia Gkimisi, Stefanos Papadantonakis <br />
Scientific supervision: Prof. Pavlos Pavlidis <br />
Contact: [efstathiouanna@gmail.com](mailto:efstathiouanna@gmail.com)

## License

Copyright (c) 2026 <br />
Foundation for Research and Technology – Hellas <br />

The repository will be updated after the associated publication.