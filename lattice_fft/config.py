"""Configuration dataclasses for lattice FFT simulations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class LatticeParameters:
    """Geometric description of a 2D lattice simulation domain."""

    ny: int = 512
    nx: int = 512
    px_size_angstrom: float = 0.5
    lattice_constant_angstrom: float = 5.43
    atom_sigma_angstrom: float = 0.8
    atom_amplitude: float = 1.0

    def __post_init__(self) -> None:
        if self.ny <= 0 or self.nx <= 0:
            raise ValueError("Simulation dimensions must be positive")
        if self.px_size_angstrom <= 0:
            raise ValueError("Pixel size must be positive")
        if self.lattice_constant_angstrom <= 0:
            raise ValueError("Lattice constant must be positive")
        if self.atom_sigma_angstrom <= 0:
            raise ValueError("Atomic Gaussian width must be positive")
        if self.atom_amplitude <= 0:
            raise ValueError("Atomic amplitude must be positive")

    @property
    def shape(self) -> Tuple[int, int]:
        return (self.ny, self.nx)

    @property
    def period_px(self) -> float:
        return self.lattice_constant_angstrom / self.px_size_angstrom

    @property
    def atom_sigma_px(self) -> float:
        return self.atom_sigma_angstrom / self.px_size_angstrom


@dataclass(frozen=True)
class DefectParameters:
    """Parameters controlling how defects are introduced."""

    defect_fraction: float = 0.05
    max_row_shift_px: float = 2.0
    remove_fraction: float = 0.03

    def __post_init__(self) -> None:
        if not 0.0 <= self.defect_fraction <= 1.0:
            raise ValueError("defect_fraction must lie in [0, 1]")
        if self.max_row_shift_px < 0:
            raise ValueError("max_row_shift_px must be non-negative")
        if not 0.0 <= self.remove_fraction <= 1.0:
            raise ValueError("remove_fraction must lie in [0, 1]")


@dataclass(frozen=True)
class JitterParameters:
    """Parameters describing random jitter/disorder."""

    sigma_px: float = 0.0

    def __post_init__(self) -> None:
        if self.sigma_px < 0:
            raise ValueError("sigma_px must be non-negative")

    @property
    def has_jitter(self) -> bool:
        return self.sigma_px > 0
