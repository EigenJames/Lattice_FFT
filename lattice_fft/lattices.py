"""Lattice basis helpers for generating atomic coordinates."""
from __future__ import annotations

from typing import Literal, Tuple

import numpy as np

Basis = Literal["square", "hex"]


def square_lattice(period_px: float, shape: Tuple[int, int]) -> np.ndarray:
    """Return square lattice positions within the image bounds."""
    ny, nx = shape
    ys = np.arange(0, ny, period_px)
    xs = np.arange(0, nx, period_px)
    grid_y, grid_x = np.meshgrid(ys, xs, indexing="ij")
    return np.column_stack((grid_y.ravel(), grid_x.ravel()))


def hexagonal_lattice(period_px: float, shape: Tuple[int, int]) -> np.ndarray:
    """Return hexagonal lattice positions within the image bounds."""
    ny, nx = shape
    dy = period_px * np.sqrt(3.0) / 2.0
    ys = np.arange(0, ny + dy, dy)
    positions: list[Tuple[float, float]] = []
    for idx, y in enumerate(ys):
        if y >= ny:
            continue
        offset = (period_px / 2.0) if idx % 2 else 0.0
        xs = np.arange(offset, nx + period_px, period_px)
        for x in xs:
            if x >= nx:
                continue
            positions.append((y, x))
    return np.asarray(positions, dtype=float)


def generate_lattice(period_px: float, shape: Tuple[int, int], basis: Basis) -> np.ndarray:
    """Dispatch lattice generation according to the requested basis."""
    if basis == "square":
        return square_lattice(period_px, shape)
    if basis == "hex":
        return hexagonal_lattice(period_px, shape)
    raise ValueError(f"Unsupported basis '{basis}'")


__all__ = [
    "Basis",
    "generate_lattice",
    "hexagonal_lattice",
    "square_lattice",
]
