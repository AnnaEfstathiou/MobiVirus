import pandas as pd
from typing import Dict
from Bio import SeqIO

def __csv_to_fasta(csv_file, output_fasta_file):

    """ Convert a CSV file to a FASTA file while excluding rows that are entirely NaN. """
    
    try:
        df = pd.read_csv(csv_file)  # Read the csv file
        with open(output_fasta_file, 'w') as fasta_file:
            for index, row in df.iterrows():
                if row.dropna().empty:
                    continue  # Skip rows that are all NaN
                # Convert each value to integer if it's a whole number, else keep original
                sequence = ''.join(str(int(val)) if isinstance(val, float) and val.is_integer() else str(val) for val in row)
                header = f'>{index + 1}'
                fasta_file.write(f'{header}\n{sequence}\n')
    except FileNotFoundError:
        print(f"Error: The file {csv_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __process_fasta(fasta_file, output_fasta):

    """ Processes a FASTA file to exclude rows that are entirely NaN. """

    try:
        # open input and output files
        with open(fasta_file, "r") as input_handle, open(output_fasta, "w") as output_handle:
            # parse the input FASTA file
            for line in SeqIO.parse(input_handle, "fasta"):
                # check if the sequence contains at least one position with 1
                if not 'nan' in line.seq:
                    # write the sequence to the output file
                    SeqIO.write(line, output_handle, "fasta")
    except FileNotFoundError:
        print(f"Error: The file {fasta_file} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def __read_fasta(fasta_file):

    """ Function that creates a dictionary with headers as keys and sequences as values. """

    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()
    return sequences


def __check(sequences: Dict[str, str]) -> None:

    """ Check for valid sequence proportions. """
    """ Number of sequences > 1 """
    """ Sequences with the same length """

    if len(sequences) < 2:
        raise ValueError("At least 2 sequences required!")

    sequences_list = list(sequences.values()) # list of sequences
    seqs_len = len(sequences_list[0]) # sequence's length 
    if not all(len(s) == seqs_len for s in sequences_list):
        raise ValueError("Sequences are required to have the same length!")