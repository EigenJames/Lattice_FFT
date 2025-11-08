"""Simulation primitives for lattice FFT demonstrations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Optional, Tuple

import numpy as np
from numpy.fft import fft2, fftfreq, fftshift

from .config import DefectParameters, JitterParameters, LatticeParameters

Basis = Literal["square", "hex"]


@dataclass
class SimulationResult:
    """Container for the outputs of a lattice simulation."""

    positions: np.ndarray
    image: np.ndarray
    fft_magnitude: np.ndarray
    kx: np.ndarray
    ky: np.ndarray


def gaussian_kernel(size_px: int, sigma_px: float) -> np.ndarray:
    """Create a normalized 2D Gaussian kernel."""
    if size_px <= 0:
        raise ValueError("size_px must be positive")
    half = size_px // 2
    y, x = np.mgrid[-half : half + 1, -half : half + 1]
    kernel = np.exp(-(x**2 + y**2) / (2.0 * sigma_px**2))
    kernel_sum = kernel.sum()
    if kernel_sum <= 0:
        raise ValueError("Gaussian kernel sum must be positive")
    return kernel / kernel_sum


def stamp_atoms(image: np.ndarray, positions: Iterable[Tuple[float, float]], kernel: np.ndarray) -> None:
    """Stamp Gaussian atoms into an image at sub-pixel locations."""
    ky, kx = kernel.shape
    hy, hx = ky // 2, kx // 2
    ny, nx = image.shape

    for (y0, x0) in positions:
        iy = int(np.floor(y0))
        ix = int(np.floor(x0))
        ys = max(0, iy - hy)
        ye = min(ny, iy + hy + 1)
        xs = max(0, ix - hx)
        xe = min(nx, ix + hx + 1)

        ky0 = ys - (iy - hy)
        ky1 = ky0 + (ye - ys)
        kx0 = xs - (ix - hx)
        kx1 = kx0 + (xe - xs)

        if ys < ye and xs < xe:
            image[ys:ye, xs:xe] += kernel[ky0:ky1, kx0:kx1]


def generate_lattice_positions(params: LatticeParameters, basis: Basis = "square") -> np.ndarray:
    """Generate lattice site coordinates in pixel units."""
    ny, nx = params.shape
    period_px = params.period_px
    positions = []

    if basis == "square":
        ys = np.arange(0, ny, period_px)
        xs = np.arange(0, nx, period_px)
        for y in ys:
            for x in xs:
                positions.append((y, x))
    elif basis == "hex":
        dy = period_px * np.sqrt(3.0) / 2.0
        ys = np.arange(0, ny, dy)
        for idx, y in enumerate(ys):
            offset = (period_px / 2.0) if (idx % 2 == 1) else 0.0
            xs = np.arange(offset, nx, period_px)
            for x in xs:
                positions.append((y, x))
    else:
        raise ValueError(f"Unsupported basis '{basis}'")

    return np.array(positions, dtype=float)


def apply_defects(
    positions: np.ndarray,
    params: LatticeParameters,
    defect_params: Optional[DefectParameters] = None,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Apply row shifts and random removals to lattice positions."""
    if defect_params is None:
        return positions.copy()

    rng = rng or np.random.default_rng()
    ny, nx = params.shape
    pos = positions.copy()

    rows = np.unique(np.round(pos[:, 0]).astype(int))
    n_shift_rows = max(1, int(defect_params.defect_fraction * len(rows)))
    shift_rows = rng.choice(rows, size=min(n_shift_rows, len(rows)), replace=False)
    for r in shift_rows:
        shift = rng.uniform(-defect_params.max_row_shift_px, defect_params.max_row_shift_px)
        mask = np.round(pos[:, 0]).astype(int) == r
        pos[mask, 1] += shift

    n_remove = int(defect_params.remove_fraction * len(pos))
    if n_remove > 0:
        remove_idx = rng.choice(len(pos), size=n_remove, replace=False)
        pos = np.delete(pos, remove_idx, axis=0)

    pos[:, 0] = np.clip(pos[:, 0], 0, ny - 1)
    pos[:, 1] = np.clip(pos[:, 1], 0, nx - 1)
    return pos


def apply_jitter(
    positions: np.ndarray,
    jitter: Optional[JitterParameters] = None,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Apply isotropic Gaussian jitter to lattice positions."""
    if jitter is None or not jitter.has_jitter:
        return positions.copy()

    rng = rng or np.random.default_rng()
    noise = rng.normal(0.0, jitter.sigma_px, size=positions.shape)
    return positions + noise


def render_lattice_image(params: LatticeParameters, positions: np.ndarray) -> np.ndarray:
    """Render a lattice as a float32 image."""
    image = np.zeros(params.shape, dtype=float)
    ksize = int(max(7, np.ceil(6.0 * params.atom_sigma_px)))
    if ksize % 2 == 0:
        ksize += 1
    kernel = gaussian_kernel(ksize, params.atom_sigma_px)
    stamp_atoms(image, positions, kernel)
    return image


def compute_fft_magnitude(image: np.ndarray, params: LatticeParameters) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute log-magnitude FFT and reciprocal-space axes."""
    F = fftshift(fft2(image))
    mag = np.abs(F)
    mag_log = np.log10(mag + 1e-6)
    ky = fftshift(fftfreq(image.shape[0], d=params.px_size_angstrom))
    kx = fftshift(fftfreq(image.shape[1], d=params.px_size_angstrom))
    return mag_log, kx, ky


def simulate_lattice(
    params: LatticeParameters,
    *,
    basis: Basis = "square",
    defect_params: Optional[DefectParameters] = None,
    jitter: Optional[JitterParameters] = None,
    rng: Optional[np.random.Generator] = None,
) -> SimulationResult:
    """Generate a lattice image and its FFT signature."""
    rng = rng or np.random.default_rng()

    base_positions = generate_lattice_positions(params, basis=basis)
    defected = apply_defects(base_positions, params, defect_params, rng=rng)
    jittered = apply_jitter(defected, jitter, rng=rng)

    image = render_lattice_image(params, jittered)
    fft_mag, kx, ky = compute_fft_magnitude(image, params)

    return SimulationResult(positions=jittered, image=image, fft_magnitude=fft_mag, kx=kx, ky=ky)
