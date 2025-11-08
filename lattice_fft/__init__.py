"""Utilities for simulating 2D lattices and their FFT signatures."""
from .config import LatticeParameters, DefectParameters, JitterParameters
from .simulation import (
    SimulationResult,
    compute_fft_magnitude,
    generate_lattice_positions,
    render_lattice_image,
    simulate_lattice,
)
from .plotting import plot_real_and_fft

__all__ = [
    "LatticeParameters",
    "DefectParameters",
    "JitterParameters",
    "SimulationResult",
    "generate_lattice_positions",
    "render_lattice_image",
    "compute_fft_magnitude",
    "simulate_lattice",
    "plot_real_and_fft",
]
