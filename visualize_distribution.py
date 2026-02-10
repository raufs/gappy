#!/usr/bin/env python3
"""
Generate a visualization of the gap length distribution using gappy's default parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Default parameters from gappy
ALPHA = 1.09
BETA = 10.0
MIN_GAP_SIZE = 1
MAX_GAP_SIZE = 1000000

# Sample gap lengths (same logic as in gappy.py)
def sample_gap_lengths(num_samples=10000):
    """Sample gap lengths from beta distribution."""
    beta_samples = np.random.beta(ALPHA, BETA, num_samples)
    gap_range = MAX_GAP_SIZE - MIN_GAP_SIZE
    gap_lengths = MIN_GAP_SIZE + (beta_samples * gap_range)
    return gap_lengths.astype(int)

# Set random seed for reproducibility
np.random.seed(42)

# Generate samples
gap_lengths = sample_gap_lengths(10000)

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Format function for axis labels
def format_bp(x, pos):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    elif x >= 1e3:
        return f'{x/1e3:.0f}k'
    else:
        return f'{int(x)}'

# Left plot: Linear scale histogram
ax1.hist(gap_lengths, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
ax1.set_xlabel('Gap Length (bp)', fontsize=11)
ax1.set_ylabel('Frequency', fontsize=11)
ax1.set_title('Gap Length Distribution\n(Linear Scale)', fontsize=12, fontweight='bold')
ax1.xaxis.set_major_formatter(FuncFormatter(format_bp))
ax1.grid(True, alpha=0.3)

# Add statistics text box
stats_text = f'Parameters:\nα = {ALPHA}\nβ = {BETA}\n\nStatistics:\nMean: {np.mean(gap_lengths)/1e3:.1f}k bp\nMedian: {np.median(gap_lengths)/1e3:.1f}k bp\nMode ≈ 10k bp'
ax1.text(0.65, 0.97, stats_text, transform=ax1.transAxes,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
         fontsize=9, family='monospace')

# Right plot: Log scale histogram
ax2.hist(gap_lengths, bins=100, edgecolor='black', alpha=0.7, color='coral')
ax2.set_xlabel('Gap Length (bp, log scale)', fontsize=11)
ax2.set_ylabel('Frequency', fontsize=11)
ax2.set_title('Gap Length Distribution\n(Log Scale)', fontsize=12, fontweight='bold')
ax2.set_xscale('log')
ax2.xaxis.set_major_formatter(FuncFormatter(format_bp))
ax2.grid(True, alpha=0.3, which='both')

# Add range text box
range_text = f'Range:\n[{MIN_GAP_SIZE} bp, {MAX_GAP_SIZE/1e6:.0f}M bp]'
ax2.text(0.05, 0.97, range_text, transform=ax2.transAxes,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
         fontsize=9, family='monospace')

plt.suptitle('GAPPY Default Gap Length Distribution', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

# Save figure
output_file = 'gap_distribution.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Visualization saved to: {output_file}")

# Print statistics
print(f"\nDistribution Statistics:")
print(f"  Mean: {np.mean(gap_lengths):,.1f} bp")
print(f"  Median: {np.median(gap_lengths):,.1f} bp")
print(f"  Min: {np.min(gap_lengths):,} bp")
print(f"  Max: {np.max(gap_lengths):,} bp")
print(f"  Std Dev: {np.std(gap_lengths):,.1f} bp")
