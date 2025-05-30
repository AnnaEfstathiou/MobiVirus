# 🧬 Genome Sequence Statistics Pipeline

A Python script for processing genomic sequence data (in CSV format) to compute population genetics statistics, including **Tajima's D**, **nucleotide diversity (π)**, **Watterson’s theta (θw)**, **haplotype count**, and **haplotype diversity**. The script supports analysis on either **super strains**, **normal strains**, or a **mix** of both.

---

## 📌 Features

- Efficiently reads and processes binary-encoded genome sequences.
- Identifies "super spreader" strains based on mutations in central genome positions.
- Supports flexible sampling for different strain groups.
- Computes core population genetics statistics using `libsequence`.
- Outputs results to a structured CSV file.

---

## 📁 Input Format

The input must be a `.csv` file with:

- **No headers**
- **Each row**: a binary sequence (e.g., `0,1,0,0,...`)
- Each value represents a **SNP (Single Nucleotide Polymorphism)** position.

Example (`genomes_100.csv`):

```python
0,0,0,1,0
0,1,0,1,1
1,0,1,0,1

# OR 

0.0,0.0,0.0,1.0,0.0
0.0,1.0,0.0,1.0,1.0
1.0,0.0,1.0,0.0,1.0
```

---

## 🚀 How It Works

1. **Data Import**  
   Reads genome data into a DataFrame and separates:
   - **Super Strains**: contain a `1` in the middle SNP position(s)
   - **Normal Strains**: all others

2. **Sampling**  
   Randomly samples a number of strains from a group (or entire group if sample size exceeds group size).

3. **Simulation Preparation**  
   Filters out monomorphic positions and prepares input for `libsequence`.

4. **Statistical Analysis**  
   Calculates:
   - Tajima’s D
   - Pi (nucleotide diversity)
   - Theta Watterson (θw)
   - Number of unique haplotypes
   - Haplotype diversity

5. **Output**  
   Saves results to `sumstats_<event_number>.csv` (e.g., `sumstats_100.csv`)

---

## 📊 Output Example

Output CSV (`sumstats_100.csv`):

| Tajimas D | Pi-estimator | Theta Watterson | Number of unique haplotypes | Haplotype Diversity |
|-----------|--------------|------------------|------------------------------|----------------------|
| -1.23     | 0.034        | 0.045            | 17                           | 0.89                 |
| -1.84     | 0.094        | 0.055            | 20                           | 0.92                 |

---

## 🧪 Command-Line Interface

```bash
python genome_stats.py -g genomes_100.csv -s 30
```

### Arguments
Flag | Description|
|-----------|--------------|
|-g, --genome_file | Path to the genome CSV file (required) |
|-s, --sample_size | Number of sequences to sample (required) |
|-p, --population_size | Print population sizes (no sampling or stats) |
|-ss, --ss_strains | Use only super strains |
|-ns, --ns_strains | Use only normal strains |

If neither --ss_strains nor --ns_strains is specified, the script uses both strains (mix).

---

## 🧬 Super Strains Detection Logic
The script defines super strains as sequences where at least one of the middle SNP positions equals 1. It handles both even- and odd-length sequences.

---

## 🛠 Dependencies
Make sure you have the following Python packages installed:
```python
pip install pandas libsequence
```

---

## 📂 File Structure
```bash
genome_stats.py       # Main script
genomes_100.csv       # Example input file
sumstats_100.csv      # Output file (generated)
```

---

## 🔐 Error Handling
Raises ValueError if:
- Sample size is invalid
- No data is available
- Sampling returns only one sequence
- Warns if sample size > population (uses full dataset instead)

---

## 👨‍🔬 Example Use Case
You have simulated or real SNP data for a viral population. You want to:
- Identify potential super-spreader variants
- Compare diversity between normal vs. super strains
- Track how diversity metrics change across simulated events (genomes_100.csv, genomes_200.csv, etc.)
- Generate CSV-ready summary stats for plotting or meta-analysis

---

## 📞 Contact
For questions or improvements, feel free to open an issue or contact the developer. Contributions are welcome!