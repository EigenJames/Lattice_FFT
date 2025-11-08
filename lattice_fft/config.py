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


@dataclass(frozen=True)
class JitterParameters:
    """Parameters describing random jitter/disorder."""

    sigma_px: float = 0.0

    @property
    def has_jitter(self) -> bool:
        return self.sigma_px > 0
