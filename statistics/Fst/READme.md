# 🧬 Genome Sequence Fst Analysis Pipeline

A Python script for processing binary-encoded genomic sequence data (CSV format) to compute **fixation index (Fst)** between populations using various **sampling strategies**. The script supports analysis on either **super strains**, **normal strains**, or a **mix** of both, and allows for stratified or coordinate-based sampling.

---

## 📌 Features

- Efficiently loads and processes binary genome sequence data.
- Identifies "super spreader" (super) strains via SNPs in central genome positions.
- Supports three sampling modes: simple random, stratified, and spatial (by x-coordinate).
- Computes pairwise Fst using `libsequence` with three estimators:
  - Hudson, Slatkin & Maddison (HSM)
  - Slatkin
  - Hudson, Boos & Kaplan (HBK)
- Saves results to a CSV file with a suffix corresponding to the input event number.

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
    Samples sequences from the population based on the chosen strategy:
    - simple_random: sample across all strains
    - stratified: equal sampling from super and normal strains
    - coordinates: splits sequences based on x-axis midpoint in a separate coordinates file

3. **Simulation Preparation**  
    Transforms the sampled data into a format compatible with libsequence.

4. **Fst Calculation**
    Uses libsequence.Fst() to compute:
    - Hudson, Slatkin & Maddison (HSM)
    - Slatkin
    - Hudson, Boos & Kaplan (HBK)

5. **Output**  
    Results saved to a file: `fst_<event_number>.csv` (e.g., fst_100.csv)

---

## 📊 Output Example

Output CSV (`fst_100.csv`):

|HSM Fst | Slatkin Fst | HBK Fst |
|0.103 | 0.129 | 0.091 |

---

## 🧪 Command-Line Interface
```bash
python genome_fst.py -g genomes_100.csv -s 50 -sample_type coords -c coords_100.csv
``` 

### Arguments
Flag | Description|
|-----------|--------------|
|-g, --genome_file | Path to the genome CSV file (required) |
|-s, --sample_size | Number of sequences to sample (required) |
|-p, --population_size | Print population sizes (no sampling or stats) |
|-sample_type, --sampling_technique | One of: rdm, str, coords (required) |
|-c, --coords_file | Required for coords sampling mode |
|-ss, --ss_strains | Use only super strains |
|-ns, --ns_strains | Use only normal strains |

If neither --ss_strains nor --ns_strains is specified, the script uses both strains (mix).

---

## 🧬 Super Strains Detection Logic
The script defines super strains as sequences where at least one of the middle SNP positions equals 1. It handles both even- and odd-length sequences.

---

## 🛠 Dependencies
Install required libraries:

```bash
pip install pandas libsequence
```

---

## 📂 File Structure
```bash
genome_fst.py          # Main script
genomes_100.csv        # Example genome input
coords_100.csv         # Coordinates file (for spatial sampling)
fst_100.csv            # Output file
```

---

## 🔐 Error Handling
Raises exceptions if:
- Sampling fails or returns < 2 entries
- Sample size > population (warns and uses full data)
- Coordinate/genome file mismatch (in coordinate-based sampling)

---

## 👨‍🔬 Example Use Case
You're simulating the spread of mutations in a population of genomes. This script lets you:
- Identify highly mutated super strains
- Compare genetic divergence (Fst) between strain types or spatial regions
- Generate CSV-ready summary stats for plotting or meta-analysis

---

## 📞 Contact
For questions or improvements, feel free to open an issue or contact the developer. Contributions are welcome!