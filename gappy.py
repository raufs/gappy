#!/usr/bin/env python3
"""
gappy - A tool for simulating gaps in genome assemblies

This tool takes a genome FASTA file and splits contigs/scaffolds at random positions
to simulate gaps, with gap lengths sampled from a beta distribution. Each input sequence
is split into multiple smaller contigs with modified names (e.g., contig_0, contig_1, etc.).
"""

import argparse
import bisect
import sys
import random
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Simulate a worse quality genomic assembly.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "input_fasta",
        type=str,
        help="Input genome file in FASTA format"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output FASTA file with simulated gaps. If not specified, writes to stdout"
    )
    
    parser.add_argument(
        "-p", "--percent-remaining",
        type=float,
        required=True,
        help="Percentage of genome to remain after introducing gaps (0-100). E.g., 90 = 90%% remains, 10%% becomes gaps"
    )
    
    parser.add_argument(
        "-a", "--alpha",
        type=float,
        default=1.09,
        help="Alpha parameter for beta distribution. Default: 1.09 (with beta=10 creates peak around 10kb)"
    )
    
    parser.add_argument(
        "-b", "--beta",
        type=float,
        default=10.0,
        help="Beta parameter for beta distribution. Default: 10.0 (with alpha=1.09 creates peak around 10kb)"
    )
    
    parser.add_argument(
        "--min-gap-size",
        type=int,
        default=1,
        help="Minimum gap size in base pairs. Default: 1"
    )
    
    parser.add_argument(
        "--max-gap-size",
        type=int,
        default=1000000,
        help="Maximum gap size in base pairs. Default: 1000000 (1 Mbp)"
    )
    
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=0.0,
        help="Point mutation rate (probability per base, 0.0-1.0). Applied before gap simulation. Default: 0.0 (no mutations)"
    )
    
    parser.add_argument(
        "--ts-tv-ratio",
        type=float,
        default=2.0,
        help="Transition/transversion ratio. Transitions (A↔G, C↔T) vs transversions (A↔C, A↔T, G↔C, G↔T). Default: 2.0"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    
    return parser.parse_args()


def sample_gap_lengths(num_gaps, alpha, beta, min_size, max_size):
    """
    Sample gap lengths from a beta distribution.
    
    Args:
        num_gaps: Number of gaps to generate
        alpha: Alpha parameter for beta distribution
        beta: Beta parameter for beta distribution
        min_size: Minimum gap size
        max_size: Maximum gap size
    
    Returns:
        List of gap lengths (integers)
    """
    # Sample from beta distribution (values between 0 and 1)
    beta_samples = np.random.beta(alpha, beta, num_gaps)
    
    # Scale to desired range [min_size, max_size]
    gap_lengths = min_size + beta_samples * (max_size - min_size)
    
    # Convert to integers and ensure they're at least min_size
    gap_lengths = np.maximum(np.round(gap_lengths).astype(int), min_size)
    
    return gap_lengths.tolist()


def introduce_mutations(records, mutation_rate, ts_tv_ratio=2.0):
    """
    Introduce point mutations into sequences with transition/transversion bias.
    
    Args:
        records: List of SeqRecord objects
        mutation_rate: Probability of mutation per base (0.0 to 1.0)
        ts_tv_ratio: Transition/transversion ratio (default: 2.0)
    
    Returns:
        Tuple of (mutated_records, total_mutations, transitions, transversions)
    """
    if mutation_rate <= 0:
        return records, 0, 0, 0
    
    # Define transitions and transversions for each base
    # Transitions: A↔G (purines), C↔T (pyrimidines)
    # Transversions: all other changes
    mutations_map = {
        'A': {'transitions': ['G'], 'transversions': ['C', 'T']},
        'G': {'transitions': ['A'], 'transversions': ['C', 'T']},
        'C': {'transitions': ['T'], 'transversions': ['A', 'G']},
        'T': {'transitions': ['C'], 'transversions': ['A', 'G']},
        'a': {'transitions': ['g'], 'transversions': ['c', 't']},
        'g': {'transitions': ['a'], 'transversions': ['c', 't']},
        'c': {'transitions': ['t'], 'transversions': ['a', 'g']},
        't': {'transitions': ['c'], 'transversions': ['a', 'g']},
    }
    
    # Calculate weights for transitions vs transversions
    # If ts_tv_ratio = 2.0, transitions should be 2x more common than transversions overall
    # Since there are 2 transversion options and 1 transition option:
    # P(transition) / P(all transversions) = ts_tv_ratio
    # If each transversion has weight w, transition has weight t:
    # t / (2w) = ts_tv_ratio, so t = 2 * ts_tv_ratio * w
    transversion_weight = 1.0
    transition_weight = 2.0 * ts_tv_ratio * transversion_weight
    
    mutated_records = []
    total_mutations = 0
    total_transitions = 0
    total_transversions = 0
    
    for record in records:
        seq_str = str(record.seq)
        seq_list = list(seq_str)
        seq_length = len(seq_list)
        
        # Pre-calculate number of mutations to introduce based on mutation rate
        # This ensures proper random distribution across the genome
        expected_mutations = int(seq_length * mutation_rate)
        
        # Add probabilistic component for fractional part
        if random.random() < (seq_length * mutation_rate - expected_mutations):
            expected_mutations += 1
        
        if expected_mutations == 0:
            mutated_records.append(record)
            continue
        
        # Randomly select positions to mutate (uniformly distributed).
        # Filter to only mutable positions (standard nucleotides) using a
        # vectorized numpy mask instead of a per-base Python loop, since the
        # latter is a memory/CPU bottleneck on chromosome/genome-scale inputs
        seq_bytes = np.frombuffer(seq_str.encode('ascii'), dtype='S1')
        mutable_bases = np.array(list(mutations_map.keys()), dtype='S1')
        mutable_positions = np.flatnonzero(np.isin(seq_bytes, mutable_bases))

        if mutable_positions.size == 0:
            mutated_records.append(record)
            continue

        # Randomly sample positions without replacement
        num_mutations = min(expected_mutations, mutable_positions.size)
        mutation_positions = np.random.choice(mutable_positions, size=num_mutations, replace=False)
        
        mutations_in_seq = 0
        transitions_in_seq = 0
        transversions_in_seq = 0
        
        # Apply mutations at randomly selected positions
        for i in mutation_positions:
            base = seq_list[i]
            
            # Get possible transitions and transversions
            transitions = mutations_map[base]['transitions']
            transversions = mutations_map[base]['transversions']
            
            # Build choices and weights
            choices = transitions + transversions
            weights = [transition_weight] * len(transitions) + [transversion_weight] * len(transversions)
            
            # Select mutation based on weights
            new_base = random.choices(choices, weights=weights)[0]
            seq_list[i] = new_base
            mutations_in_seq += 1
            
            # Track transition vs transversion
            if new_base in transitions:
                transitions_in_seq += 1
            else:
                transversions_in_seq += 1
        
        total_mutations += mutations_in_seq
        total_transitions += transitions_in_seq
        total_transversions += transversions_in_seq
        
        # Create new record with mutated sequence
        mutated_record = SeqRecord(
            Seq(''.join(seq_list)),
            id=record.id,
            description=record.description
        )
        mutated_records.append(mutated_record)
    
    return mutated_records, total_mutations, total_transitions, total_transversions


def split_sequence_at_gaps(sequence, gap_positions, gap_lengths):
    """
    Split a sequence at specified gap positions, removing the gap regions.

    Args:
        sequence: Original sequence string
        gap_positions: List of positions where gaps should start (sorted)
        gap_lengths: List of gap lengths corresponding to each position

    Returns:
        List of (fragment, original_start) tuples, where original_start is the
        1-based position of the fragment's first base in the original sequence.
    """
    seq_str = str(sequence)

    if not gap_positions:
        return [(seq_str, 1)]

    # Sort positions and corresponding lengths together
    gaps = sorted(zip(gap_positions, gap_lengths))

    # Create fragments by removing gap regions, tracking each fragment's
    # original 1-based start position alongside it so the two can never
    # drift out of sync (e.g. if a gap produces a zero-length fragment)
    fragments = []
    start = 0
    orig_start = 1

    for pos, gap_len in gaps:
        if pos > start:
            # Add fragment from start to gap position
            fragments.append((seq_str[start:pos], orig_start))

        # Skip over the gap region - next fragment starts after the gap.
        # Overlapping/nested gaps must merge into one continuous removed
        # region, so only extend start/orig_start if this gap's end goes
        # further than what's already been removed - otherwise a gap fully
        # nested inside an earlier one would wrongly un-delete bases.
        gap_end = pos + gap_len
        if gap_end > start:
            start = gap_end
            orig_start = gap_end + 1

    # Add final fragment if there's sequence remaining
    if start < len(seq_str):
        fragments.append((seq_str[start:], orig_start))

    return fragments


def insert_merged_interval(intervals, start, end):
    """
    Insert half-open interval [start, end) into a sorted list of disjoint,
    merged (start, end) intervals, merging any overlaps in place.

    This lets a caller track the true, overlap-aware union of many gaps
    incrementally -- one gap at a time -- without ever materializing the
    actual sequence fragments those gaps would produce.

    Args:
        intervals: Sorted list of disjoint (start, end) tuples, mutated in place
        start, end: The new interval to insert (end > start)

    Returns:
        Number of previously-uncovered bases newly covered by this insertion
    """
    if end <= start:
        return 0

    # Find the first existing interval that could overlap or touch the new
    # one, then absorb every subsequent interval that overlaps it too
    i = bisect.bisect_left(intervals, (start,))
    if i > 0 and intervals[i - 1][1] >= start:
        i -= 1

    j = i
    merged_start, merged_end = start, end
    while j < len(intervals) and intervals[j][0] <= merged_end:
        merged_start = min(merged_start, intervals[j][0])
        merged_end = max(merged_end, intervals[j][1])
        j += 1

    old_covered = sum(e - s for s, e in intervals[i:j])
    intervals[i:j] = [(merged_start, merged_end)]
    return (merged_end - merged_start) - old_covered


def build_fragments_for_record(record, gaps_here):
    """
    Split a single record at its assigned gaps (if any) into fragment
    SeqRecords, plus the corresponding mapping-data entries.

    Args:
        record: SeqRecord to split
        gaps_here: List of (local_position, gap_length) for this record, or
            None/empty if the record has no gaps

    Returns:
        Tuple of (fragment_records, mapping_entries)
    """
    if not gaps_here:
        return [record], [(record.id, 1, len(record.seq), record.id, 1, len(record.seq))]

    positions, gap_lengths = zip(*gaps_here)
    fragments = split_sequence_at_gaps(record.seq, positions, gap_lengths)

    # Preserve any extra metadata from the original description
    original_desc = record.description
    extra_desc = original_desc[len(record.id):].strip() if original_desc.startswith(record.id) else original_desc

    fragment_records = []
    mapping_entries = []
    for i, (fragment, frag_orig_start) in enumerate(fragments):
        fragment_id = f"{record.id}_{i}"
        split_note = f"Split from {record.id} (fragment {i+1}/{len(fragments)})"
        description = f"{split_note} {extra_desc}" if extra_desc else split_note
        fragment_records.append(SeqRecord(Seq(fragment), id=fragment_id, description=description))

        fragment_len = len(fragment)
        orig_end = frag_orig_start + fragment_len - 1
        # New coordinates are always relative to the fragment itself, since
        # each fragment is its own output record
        mapping_entries.append((record.id, frag_orig_start, orig_end, fragment_id, 1, fragment_len))

    return fragment_records, mapping_entries


def introduce_gaps(records, percent_remaining, alpha, beta, min_size, max_size):
    """
    Introduce gaps into genome sequences by splitting them.
    
    Args:
        records: List of SeqRecord objects
        percent_remaining: Percentage of genome that should remain (0-100)
        alpha: Alpha parameter for beta distribution
        beta: Beta parameter for beta distribution
        min_size: Minimum gap size
        max_size: Maximum gap size
    
    Returns:
        Tuple of (split_records, mapping_data)
        - split_records: List of split SeqRecord objects (more records than input)
        - mapping_data: List of tuples (orig_seq_id, orig_start, orig_end, new_seq_id, new_start, new_end)
    """
    # Calculate total genome length
    total_length = sum(len(record.seq) for record in records)
    
    if total_length == 0:
        print("Error: Input sequences have zero total length", file=sys.stderr)
        sys.exit(1)
    
    # Calculate target gap size based on percent remaining
    target_gap_size = int(total_length * (1 - percent_remaining / 100.0))
    target_remaining_length = total_length - target_gap_size
    
    if target_gap_size <= 0:
        print("Error: Percent remaining is too high (>= 100%). No gaps to introduce.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Total genome length: {total_length:,} bp", file=sys.stderr)
    print(f"Target: {percent_remaining:.1f}% remaining ({target_remaining_length:,} bp)", file=sys.stderr)
    print(f"To remove: {target_gap_size:,} bp in gaps", file=sys.stderr)

    # Sample and place gaps iteratively, batch by batch, tracking each
    # record's merged gap intervals as we go so we always know the true,
    # overlap-aware amount of sequence removed so far -- without ever
    # building a single fragment/string until the loop is done and gap
    # placement is final.
    seq_lengths = np.array([len(r.seq) for r in records], dtype=float)
    probs = seq_lengths / seq_lengths.sum()
    mean_gap_size = max(min_size + (alpha / (alpha + beta)) * (max_size - min_size), 1)

    assignments = defaultdict(list)
    intervals_by_idx = defaultdict(list)
    gap_lengths = []
    actual_removed = 0

    max_rounds = 100
    round_num = 0
    actual_percent = 0.0
    while (total_length - actual_removed > target_remaining_length
           and abs(actual_percent - percent_remaining) >= 1.0
           and round_num < max_rounds):
        round_num += 1

        needed_removal = (total_length - actual_removed) - target_remaining_length
        batch_size = max(int(needed_removal / mean_gap_size) + 1, 1)
        batch_lengths = sample_gap_lengths(batch_size, alpha, beta, min_size, max_size)
        target_indices = np.random.choice(len(records), size=len(batch_lengths), p=probs)

        for idx, gap_len in zip(target_indices, batch_lengths):
            idx = int(idx)
            seq_len = len(records[idx].seq)
            if seq_len < 2:
                continue
            position = random.randint(1, seq_len - 1)
            assignments[idx].append((position, gap_len))
            gap_lengths.append(gap_len)

            gap_end = min(position + gap_len, seq_len)
            actual_removed += insert_merged_interval(intervals_by_idx[idx], position, gap_end)

        actual_percent = (total_length - actual_removed) / total_length * 100

    if abs(actual_percent - percent_remaining) >= 1.0:
        print(f"Warning: stopped after {round_num} round(s) without reaching within 1% of "
              f"target (actual: {actual_percent:.2f}%, target: {percent_remaining:.1f}%)", file=sys.stderr)

    # Gap placement is final -- materialize the actual split fragments and
    # mapping data exactly once
    fragments_by_idx = {}
    mapping_by_idx = {}
    for idx, record in enumerate(records):
        fragments_by_idx[idx], mapping_by_idx[idx] = build_fragments_for_record(record, assignments.get(idx))

    split_records = [f for idx in range(len(records)) for f in fragments_by_idx[idx]]
    mapping_data = [m for idx in range(len(records)) for m in mapping_by_idx[idx]]

    # Final statistics
    num_gaps = len(gap_lengths)
    actual_percent_remaining = actual_percent

    print(f"Introduced {num_gaps} gaps across {round_num} round(s)...", file=sys.stderr)
    if gap_lengths:
        print(f"Gap length statistics:", file=sys.stderr)
        print(f"  Min: {min(gap_lengths)} bp", file=sys.stderr)
        print(f"  Max: {max(gap_lengths)} bp", file=sys.stderr)
        print(f"  Mean: {np.mean(gap_lengths):.1f} bp", file=sys.stderr)
        print(f"  Median: {np.median(gap_lengths):.1f} bp", file=sys.stderr)
        print(f"  Sampled total: {sum(gap_lengths):,} bp", file=sys.stderr)
    print(f"  Actually removed: {actual_removed:,} bp", file=sys.stderr)
    print(f"  Actual genome remaining: {actual_percent_remaining:.2f}%", file=sys.stderr)

    # Print per-sequence statistics, grouped by the ground-truth orig_id
    # recorded in mapping_data (rather than guessing parentage from fragment
    # ID suffixes, which breaks on original IDs that already end in digits)
    seq_fragment_lengths = defaultdict(list)
    for orig_id, _, _, _, new_start, new_end in mapping_data:
        seq_fragment_lengths[orig_id].append(new_end - new_start + 1)

    for orig_id, fragment_lengths in seq_fragment_lengths.items():
        if len(fragment_lengths) == 1:
            print(f"  {orig_id}: {fragment_lengths[0]:,} bp -> 1 contig (no gaps)", file=sys.stderr)
        else:
            orig_len = next((len(r.seq) for r in records if r.id == orig_id), 0)
            total_len = sum(fragment_lengths)
            num_gaps_in_seq = len(fragment_lengths) - 1
            print(f"  {orig_id}: {orig_len:,} bp -> {len(fragment_lengths)} contigs, {total_len:,} bp total ({num_gaps_in_seq} gaps)", file=sys.stderr)
    
    return split_records, mapping_data


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Random seed set to: {args.seed}", file=sys.stderr)
    
    # Validate input file
    input_path = Path(args.input_fasta)
    if not input_path.exists():
        print(f"Error: Input file '{args.input_fasta}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Validate parameters
    if args.alpha <= 0 or args.beta <= 0:
        print("Error: Alpha and beta parameters must be positive", file=sys.stderr)
        sys.exit(1)
    
    if args.min_gap_size < 1:
        print("Error: Minimum gap size must be at least 1", file=sys.stderr)
        sys.exit(1)
    
    if args.max_gap_size < args.min_gap_size:
        print("Error: Maximum gap size must be >= minimum gap size", file=sys.stderr)
        sys.exit(1)
    
    if args.percent_remaining < 0 or args.percent_remaining >= 100:
        print("Error: Percent remaining must be between 0 and 100 (exclusive)", file=sys.stderr)
        sys.exit(1)
    
    if args.mutation_rate < 0 or args.mutation_rate > 1:
        print("Error: Mutation rate must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)
    
    if args.ts_tv_ratio < 0:
        print("Error: Ts/Tv ratio must be non-negative", file=sys.stderr)
        sys.exit(1)
    
    output_dest = args.output if args.output else "stdout"
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"GAPPY - Genome Assembly Gap Simulator", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)
    
    print(f"Input file: {args.input_fasta}", file=sys.stderr)
    print(f"Output: {output_dest}", file=sys.stderr)
    print(f"\nMutation parameters:", file=sys.stderr)
    print(f"  Mutation rate: {args.mutation_rate}", file=sys.stderr)
    if args.mutation_rate > 0:
        print(f"  (Expected ~{args.mutation_rate * 100:.4f}% of bases mutated)", file=sys.stderr)
        print(f"  Ts/Tv ratio: {args.ts_tv_ratio}", file=sys.stderr)
    print(f"\nGap simulation parameters:", file=sys.stderr)
    print(f"  Target: {args.percent_remaining}% genome remaining", file=sys.stderr)
    print(f"  Gap removal: {100 - args.percent_remaining}% of genome", file=sys.stderr)
    print(f"  Beta distribution: Alpha={args.alpha}, Beta={args.beta}", file=sys.stderr)
    print(f"  Gap size range: [{args.min_gap_size}, {args.max_gap_size}] bp", file=sys.stderr)
    print(f"  Mode: Split sequences (no gap characters inserted)", file=sys.stderr)
    print(file=sys.stderr)
    
    # Read input FASTA file
    try:
        print("Reading input FASTA file...", file=sys.stderr)
        records = list(SeqIO.parse(args.input_fasta, "fasta"))
        print(f"Loaded {len(records)} sequence(s)", file=sys.stderr)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if len(records) == 0:
        print("Error: No sequences found in input file", file=sys.stderr)
        sys.exit(1)
    
    # Calculate total sequence length
    total_bases = sum(len(record.seq) for record in records)
    
    # Introduce point mutations (if mutation rate > 0)
    if args.mutation_rate > 0:
        print(f"\nIntroducing point mutations...", file=sys.stderr)
        records, total_mutations, transitions, transversions = introduce_mutations(
            records, args.mutation_rate, args.ts_tv_ratio
        )
        mutation_percentage = (total_mutations / total_bases * 100) if total_bases > 0 else 0
        observed_ts_tv = (transitions / transversions) if transversions > 0 else float('inf')
        
        print(f"  Total mutations: {total_mutations:,}", file=sys.stderr)
        print(f"  Mutation rate achieved: {mutation_percentage:.4f}%", file=sys.stderr)
        print(f"  Transitions: {transitions:,}", file=sys.stderr)
        print(f"  Transversions: {transversions:,}", file=sys.stderr)
        print(f"  Observed Ts/Tv ratio: {observed_ts_tv:.2f}", file=sys.stderr)
        print(f"  Total bases: {total_bases:,}", file=sys.stderr)
    else:
        print("\nSkipping mutations (rate = 0.0)", file=sys.stderr)
    
    # Introduce gaps by splitting sequences
    modified_records, mapping_data = introduce_gaps(
        records,
        args.percent_remaining,
        args.alpha,
        args.beta,
        args.min_gap_size,
        args.max_gap_size
    )
    
    # Write output FASTA file or to stdout
    try:
        if args.output:
            print(f"\nWriting output to {args.output}...", file=sys.stderr)
            SeqIO.write(modified_records, args.output, "fasta")
            print("Done!", file=sys.stderr)
            
            # Write mapping file
            mapping_file = args.output + ".mapping"
            print(f"Writing mapping file to {mapping_file}...", file=sys.stderr)
            with open(mapping_file, 'w') as f:
                f.write("# Position mapping: original -> gapped genome\n")
                f.write("# Format: original_seq_id\toriginal_start\toriginal_end\tnew_seq_id\tnew_start\tnew_end\n")
                for orig_id, orig_start, orig_end, new_id, new_start, new_end in mapping_data:
                    f.write(f"{orig_id}\t{orig_start}\t{orig_end}\t{new_id}\t{new_start}\t{new_end}\n")
            print("Done!", file=sys.stderr)
        else:
            print(f"\nWriting output to stdout...", file=sys.stderr)
            # Write to stdout without trailing newline
            from io import StringIO
            buffer = StringIO()
            SeqIO.write(modified_records, buffer, "fasta")
            output = buffer.getvalue()
            # Remove trailing newline if present
            if output.endswith('\n'):
                output = output[:-1]
            sys.stdout.write(output)
            
            # Write mapping to stderr when output goes to stdout
            print("\n# Position mapping (original -> gapped):", file=sys.stderr)
            print("# original_seq_id\toriginal_start\toriginal_end\tnew_seq_id\tnew_start\tnew_end", file=sys.stderr)
            for orig_id, orig_start, orig_end, new_id, new_start, new_end in mapping_data:
                print(f"{orig_id}\t{orig_start}\t{orig_end}\t{new_id}\t{new_start}\t{new_end}", file=sys.stderr)
    except Exception as e:
        print(f"Error writing output: {e}", file=sys.stderr)
        sys.exit(1)
    
    if args.output:
        print(f"\n{'='*60}\n", file=sys.stderr)


if __name__ == "__main__":
    main()

