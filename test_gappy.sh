#!/bin/bash
# Test script to demonstrate gappy functionality

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "======================================"
echo "Testing gappy - Gap Simulator"
echo "======================================"
echo

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna example_genome.fasta

# Test 1: High quality (90% remaining)
echo "Test 1: High quality assembly (90% genome remains)"
gappy example_genome.fasta -o test1_high_quality.fasta -p 90 --seed 42
echo

# Test 2: Moderate quality (70% remaining)
echo "Test 2: Moderate quality (70% genome remains)"
gappy example_genome.fasta -o test2_moderate.fasta -p 70 --seed 42
echo

# Test 3: Highly fragmented (50% remaining)
echo "Test 3: Highly fragmented (50% genome remains)"
gappy example_genome.fasta -o test3_fragmented.fasta -p 50 --seed 42
echo

# Test 4: Poor quality (30% remaining)
echo "Test 4: Poor quality (30% genome remains)"
gappy example_genome.fasta -o test4_poor.fasta -p 30 --seed 42
echo

# Test 5: With default parameters (using smart defaults)
echo "Test 5: Using all default parameters (80% remaining)"
gappy example_genome.fasta -o test5_defaults.fasta -p 80 --seed 42
echo

# Test 6: With point mutations (default Ts/Tv = 2.0)
echo "Test 6: With point mutations (80% remaining, 1% mutation rate)"
gappy example_genome.fasta -o test6_mutations.fasta -p 80 --mutation-rate 0.01 --seed 42
echo

# Test 7: With custom gap distribution
echo "Test 7: Custom gap sizes (70% remaining, small gaps only)"
gappy example_genome.fasta -o test7_custom.fasta -p 70 --max-gap-size 10000 --seed 42
echo

echo "======================================"
echo "All tests completed!"
echo "======================================"
echo
echo "Output files created:"
ls -lh test*.fasta
echo
echo "To view gap statistics, check the output above."
echo "To examine the sequences, use: head -20 test1_small_gaps.fasta"
