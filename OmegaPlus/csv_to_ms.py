import argparse
import os
import pandas as pd
from Bio import SeqIO

def csv_to_fasta(csv_file, sample_size, output_fasta_file):

    """ Function to convert a CSV file to a FASTA file """
    # The FASTA file is stored at the same directoey as the csv file
    # The FASTA file has the same name (different extension) as the csv file

    df = pd.read_csv(csv_file)                           # Read the CSV file into a DataFrame
    df = df.dropna()                                     # Drop rows that contain any NaN values (non existent genomes)
    df = df.apply(pd.to_numeric, errors='coerce')        # Ensures all values are numeric
    df = df.dropna()                                     # Drop any rows that had non-numeric values
    sampled_df = df.sample(n=sample_size, replace=False) # Sample sequences from the population
                                                         # Cannot take a larger sample than population when 'replace = False'
    
    if len(df) < sample_size:
        raise ValueError(f"Not enough valid rows to sample {sample_size} from {csv_file}. Only {len(df)} rows available.")

    with open(output_fasta_file, 'w') as fasta_file:
        for index, row in sampled_df.iterrows(): # Iterate through each row in the DataFrame          
            header = f'>{index + 1}' # Create the FASTA header using the row index
            # sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values]) # Concatenate the values in the row to form the sequence
            sequence = ''.join([
                str(int(val)) if isinstance(val, float) and val.is_integer() else str(val)
                for val in row.values])
            fasta_file.write(f'{header}\n{sequence}\n') # Write the header and sequence to the FASTA file


def fasta_to_ms_format(fasta_file, output_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError("No sequences found in FASTA file.")
    
    num_sites = len(records[0].seq)
    num_individuals = len(records)

    # Check all sequences are the same length
    for rec in records:
        if len(rec.seq) != num_sites:
            raise ValueError(f"Sequence length mismatch in {rec.id}")
    
    # Generate evenly spaced positions from 0.1 to 1.0
    positions = [round((i + 1) / num_sites, 5) for i in range(num_sites)]

    with open(output_file, "w") as out:
        out.write("//\n")
        out.write(f"segsites: {num_sites}\n")
        out.write("positions: " + " ".join(map(str, positions)) + "\n")
        for rec in records:
            out.write(str(rec.seq) + "\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Path to a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, help='The path to the input file (CSV).')
    parser.add_argument('-s', '--sample_size', type=int, help='The path to the input file (CSV).')
    args = parser.parse_args()

    csv_file =  args.csv_file
    sample_size = args.sample_size

    script_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the script
    base_name = os.path.splitext(os.path.basename(csv_file))[0]  # e.g., genomes_5000

    output_fasta_file_path = os.path.join(script_dir, base_name + '.fa')
    output_ms_file_path = os.path.join(script_dir, base_name + '.ms')

    ## CSV --> FASTA file
    # output_fasta_file_path = csv_file.rsplit('.', 1)[0] + '.fa' # Generate the output FASTA file path by replacing the extension
    csv_to_fasta(csv_file, sample_size, output_fasta_file_path) # Convert the CSV file to a FASTA file
    
    ## FASTA --> MS-like file
    # output_ms_file_path = csv_file.rsplit('.', 1)[0] + '.ms'        # Generate the output MS-like file path by replacing the extension
    fasta_to_ms_format(output_fasta_file_path, output_ms_file_path) # Convert the FASTA file to a MS-like file