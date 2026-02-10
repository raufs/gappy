# gappy

A genome assembly gap simulator that splits contigs/scaffolds at random positions to simulate gaps with gap lengths sampled from a beta distribution. It can also simulate sequencing errors/point mutations.

> [!NOTE]
> This is a pretty naive implementation for the task largely written using AI. A more realistic simulation of worsening assembly quality can be achieved by downsampling sequencing reads and performing re-assembly. This however requires sequencing reads to be available and is more involved computationally.

## Overview

`gappy` is a Python tool designed to simulate gaps and sequencing errors in genome assemblies for testing and benchmarking purposes of other bioinformatics software. It takes a genome file in FASTA format and splits sequences at random positions to create fragmented assemblies, with gap lengths controlled by a beta distribution.

Instead of inserting gap characters (N's), gappy splits your sequences into multiple smaller contigs, simulating real assembly fragmentation.

## Features

- Read and write FASTA files using BioPython
- **Sequencing Error simulation:**: Optional mutation rate to introduce "sequencing errors" or "mutations" before fragmenting/introducing gaps.
  - **Ts/Tv ratio control**: Adjustable transition/transversion ratio (default: 2.0)
  - Biologically realistic mutation patterns
- Sample gap lengths from a beta distribution with customizable parameters
- Split sequences at random positions to simulate gaps
- Automatic contig renaming (e.g., chr1 → chr1_0, chr1_1, chr1_2, etc.)
- Reproducible results with random seed option
- Detailed statistics on mutations (including Ts/Tv) and introduced gaps

## Installation

The project uses modern `pyproject.toml` packaging (PEP 517/518 compliant).

### Quick install (requirements only)

```bash
pip install -r requirements.txt
```

### Development install (editable mode)

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On macOS/Linux

# Install in editable mode
pip install -e .

# Now 'gappy' command is available
gappy --help
```

### Build from source

```bash
# Install build tools
pip install build

# Build distribution
python -m build

# Install the wheel
pip install dist/gappy-1.2.0-py3-none-any.whl
```

## Usage

### Basic usage

```bash
python gappy.py input_genome.fasta -o fragmented_output.fasta -p 80
```

This will remove 20% of the genome by splitting sequences, leaving 80% intact. Uses smart defaults (peak at ~10kb, gaps up to 1Mbp).

### Command-line options

```
positional arguments:
  input_fasta           Input genome file in FASTA format

required arguments:
  -o, --output          Output FASTA file with split sequences
  -p, --percent-remaining  Percentage of genome to remain (0-100)

optional arguments:
  -h, --help            Show help message and exit
  -a, --alpha           Alpha parameter for beta distribution (default: 1.09)
  -b, --beta            Beta parameter for beta distribution (default: 10.0)
  --min-gap-size        Minimum gap size in base pairs (default: 1)
  --max-gap-size        Maximum gap size in base pairs (default: 1000000)
  --mutation-rate       Point mutation rate, 0.0-1.0 (default: 0.0)
  --ts-tv-ratio         Transition/transversion ratio (default: 2.0)
  --seed                Random seed for reproducibility
```

### Beta Distribution Parameters

The beta distribution controls the distribution of gap lengths:

- **Alpha (α) and Beta (β)**: Shape parameters of the beta distribution
  - α > β: Distribution skewed towards larger gaps
  - α < β: Distribution skewed towards smaller gaps
  - α = β: Symmetric distribution
  - Both = 1: Uniform distribution
  - Both > 1: Bell-shaped distribution

The sampled values from Beta(α, β) are scaled to the range [min-gap-size, max-gap-size].

#### Default Gap Distribution Visualization

The default parameters (α=1.09, β=10.0, range=[1, 1M bp]) create a heavily right-skewed distribution with a mode around 10kb and a long tail up to 1 Mbp, mimicking realistic gap size distributions in draft assemblies:

![Gap Length Distribution](gap_distribution.png)

*Figure: Gap length distribution with default parameters. Left panel shows linear scale with peak near 10kb; right panel shows log scale revealing the full range up to 1M bp.*

### Examples

**Example 1: High quality assembly (90% remains)**

```bash
python gappy.py genome.fasta -o high_quality.fasta -p 90
```

Simulates a high-quality assembly with 10% genome removed as gaps.

**Example 2: Moderate quality (70% remains)**

```bash
python gappy.py genome.fasta -o moderate.fasta -p 70
```

Simulates moderate fragmentation with 30% genome removed.

**Example 3: Adding point mutations**

```bash
python gappy.py genome.fasta -o mutated.fasta -p 80 --mutation-rate 0.001
```

This removes 20% of genome as gaps AND adds ~0.1% point mutations (with realistic Ts/Tv ratio of 2.0).

**Example 4: Highly fragmented with custom gap sizes**

```bash
python gappy.py genome.fasta -o fragmented.fasta -p 50 --max-gap-size 50000
```

Removes 50% of genome, with gaps limited to 50kb max.

**Example 5: Reproducible results**

```bash
python gappy.py genome.fasta -o output.fasta -p 80 --seed 42
```

Using the same seed ensures reproducible mutation and split patterns.

## Requirements

- Python >= 3.7
- BioPython >= 1.79
- NumPy >= 1.21.0

**Note**: Project uses modern `pyproject.toml` packaging (PEP 517/518 compliant).

## How it works

1. **Read Input**: Parses the input FASTA file using BioPython
2. **Introduce Mutations**: (Optional) Applies point mutations at specified rate
3. **Sample Gap Lengths**: Draws gap lengths from Beta(α, β) distribution and scales them to [min-gap-size, max-gap-size]
4. **Distribute Gaps**: Allocates split points to sequences proportionally based on sequence length
5. **Split Sequences**: Randomly selects positions within each sequence and splits them into fragments
6. **Rename Contigs**: Assigns new names to fragments (original_name_0, original_name_1, etc.)
7. **Write Output**: Saves all fragment sequences to a new FASTA file

### Point Mutations

Point mutations are applied **before** gap simulation if `--mutation-rate` is set > 0. For each base:
- Probability of mutation = mutation_rate
- Only standard nucleotides (A, T, C, G) are mutated
- Mutations are substitutions to a different nucleotide
- Case is preserved (lowercase → lowercase, uppercase → uppercase)

### Transition/Transversion Ratio

The Ts/Tv ratio controls the relative frequency of transitions vs transversions:

**Transitions** (purine ↔ purine, pyrimidine ↔ pyrimidine):
- A ↔ G (purines)
- C ↔ T (pyrimidines)

**Transversions** (purine ↔ pyrimidine):
- A ↔ C, A ↔ T
- G ↔ C, G ↔ T

**Default ratio: 2.0** (transitions 2x more common than transversions)
- Biologically realistic for many organisms
- Can be adjusted with `--ts-tv-ratio`

## Output

The program outputs:
- FASTA file with split sequences (more contigs than input)
- Statistics about the gaps (min, max, mean, median, total length)
- Information about each split sequence showing:
  - Original sequence name
  - Number of fragments created
  - Total sequence length preserved

## License

This project is licensed under the G License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue.
