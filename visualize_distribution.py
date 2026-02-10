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
fig, ax = plt.subplots(1, 1, figsize=(9, 6))

# Format function for axis labels
def format_bp(x, pos):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    elif x >= 1e3:
        return f'{x/1e3:.0f}k'
    else:
        return f'{int(x)}'

# Histogram
ax.hist(gap_lengths, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
ax.set_xlabel('Gap Length (bp)', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('GAPPY Default Gap Length Distribution', fontsize=14, fontweight='bold')
ax.xaxis.set_major_formatter(FuncFormatter(format_bp))
ax.grid(True, alpha=0.3)

# Add statistics text box
stats_text = f'Parameters:\nα = {ALPHA}\nβ = {BETA}\nRange: [{MIN_GAP_SIZE} bp, {MAX_GAP_SIZE/1e6:.0f}M bp]\n\nStatistics:\nMean: {np.mean(gap_lengths)/1e3:.1f}k bp\nMedian: {np.median(gap_lengths)/1e3:.1f}k bp\nMode ≈ 10k bp'
ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
         verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85),
         fontsize=10, family='monospace')

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
