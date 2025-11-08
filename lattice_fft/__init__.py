"""
Lattice FFT: A package for crystal lattice simulation and reciprocal space analysis.

This package provides tools for:
- Generating 2D crystal lattices with various symmetries
- Modeling atomic structures with Gaussian density distributions
- Computing 2D Fast Fourier Transforms for reciprocal space analysis
- Introducing and analyzing structural defects
- Visualizing real-space and reciprocal-space patterns
"""

__version__ = "0.1.0"
__author__ = "EigenJames"

from .lattice import make_lattice_positions, introduce_defects, add_jitter
from .atoms import gaussian_kernel, stamp_atoms, render_lattice
from .fft_utils import compute_fft, reciprocal_axes, get_fft_magnitude
from .visualization import (
    plot_2x2_comparison,
    plot_lattice,
    plot_fft,
    plot_unit_cell_comparison
)

__all__ = [
    # Lattice generation
    'make_lattice_positions',
    'introduce_defects',
    'add_jitter',
    # Atomic rendering
    'gaussian_kernel',
    'stamp_atoms',
    'render_lattice',
    # FFT utilities
    'compute_fft',
    'reciprocal_axes',
    'get_fft_magnitude',
    # Visualization
    'plot_2x2_comparison',
    'plot_lattice',
    'plot_fft',
    'plot_unit_cell_comparison',
]
