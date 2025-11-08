"""
Lattice generation and defect introduction utilities.

This module provides functions to create periodic lattice structures
and introduce various types of crystal defects.
"""

import numpy as np


def make_lattice_positions(ny, nx, period_px, basis="square", rng=None):
    """
    Generate 2D lattice positions in pixels for a given basis within the image.
    
    Parameters
    ----------
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    period_px : float
        Lattice spacing in pixels
    basis : str, optional
        Lattice type: 'square' or 'hex' (hexagonal/triangular)
        Default is 'square'
    rng : numpy.random.Generator, optional
        Random number generator for reproducibility
        
    Returns
    -------
    positions : ndarray
        Array of shape (N, 2) with (y, x) positions in pixels
        
    Examples
    --------
    >>> positions = make_lattice_positions(512, 512, 10.0, basis='square')
    >>> positions.shape
    (2601, 2)
    """
    positions = []
    
    if basis == "square":
        ys = np.arange(0, ny, period_px)
        xs = np.arange(0, nx, period_px)
        for y in ys:
            for x in xs:
                positions.append((y, x))
                
    elif basis == "hex":
        # Hexagonal (triangular) lattice via row offsets
        dy = period_px * np.sqrt(3) / 2
        ys = np.arange(0, ny, dy)
        for i, y in enumerate(ys):
            offset = (period_px / 2) if (i % 2 == 1) else 0
            xs = np.arange(offset, nx, period_px)
            for x in xs:
                positions.append((y, x))
    else:
        raise ValueError(f"Unknown basis: {basis}. Choose 'square' or 'hex'.")
        
    return np.array(positions, dtype=float)


def introduce_defects(positions, ny, nx, defect_fraction=0.05, 
                      max_row_shift_px=2.0, remove_fraction=0.03, rng=None):
    """
    Introduce structural defects into lattice positions.
    
    This function simulates two types of crystal defects:
    - Dislocations: Random row shifts
    - Vacancies: Random atom removal
    
    Parameters
    ----------
    positions : ndarray
        Array of shape (N, 2) with (y, x) positions
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    defect_fraction : float, optional
        Fraction of rows to shift (0 to 1), default 0.05
    max_row_shift_px : float, optional
        Maximum shift magnitude in pixels, default 2.0
    remove_fraction : float, optional
        Fraction of atoms to remove (0 to 1), default 0.03
    rng : numpy.random.Generator, optional
        Random number generator for reproducibility
        
    Returns
    -------
    defected_positions : ndarray
        Modified positions array with defects introduced
        
    Examples
    --------
    >>> pos = make_lattice_positions(512, 512, 10.0)
    >>> rng = np.random.default_rng(42)
    >>> defected = introduce_defects(pos, 512, 512, rng=rng)
    """
    if rng is None:
        rng = np.random.default_rng()
        
    pos = positions.copy()
    
    # Row shifts: group by rounded y index (row)
    rows = np.unique(np.round(pos[:, 0]).astype(int))
    n_shift_rows = max(1, int(defect_fraction * len(rows)))
    shift_rows = rng.choice(rows, size=n_shift_rows, replace=False)
    
    for r in shift_rows:
        shift = rng.uniform(-max_row_shift_px, max_row_shift_px)
        mask = (np.round(pos[:, 0]).astype(int) == r)
        pos[mask, 1] += shift  # shift x positions within that row

    # Remove random atoms
    n_remove = int(remove_fraction * len(pos))
    if n_remove > 0:
        remove_idx = rng.choice(len(pos), size=n_remove, replace=False)
        pos = np.delete(pos, remove_idx, axis=0)

    # Keep within bounds
    pos[:, 0] = np.clip(pos[:, 0], 0, ny - 1)
    pos[:, 1] = np.clip(pos[:, 1], 0, nx - 1)
    
    return pos


def add_jitter(positions, sigma_jitter_px, rng=None):
    """
    Add random Gaussian jitter to atomic positions (thermal disorder).
    
    Parameters
    ----------
    positions : ndarray
        Array of shape (N, 2) with (y, x) positions
    sigma_jitter_px : float
        Standard deviation of Gaussian jitter in pixels
    rng : numpy.random.Generator, optional
        Random number generator for reproducibility
        
    Returns
    -------
    jittered_positions : ndarray
        Positions with added random displacement
        
    Examples
    --------
    >>> pos = make_lattice_positions(512, 512, 10.0)
    >>> rng = np.random.default_rng(42)
    >>> jittered = add_jitter(pos, sigma_jitter_px=0.5, rng=rng)
    """
    if rng is None:
        rng = np.random.default_rng()
        
    jittered = positions.copy()
    jittered += rng.normal(0, sigma_jitter_px, size=positions.shape)
    
    return jittered
