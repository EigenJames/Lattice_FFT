"""
Atomic structure modeling with Gaussian density distributions.

This module provides functions to create and render atomic structures
as Gaussian peaks in 2D images.
"""

import numpy as np


def gaussian_kernel(size_px, sigma_px):
    """
    Create a normalized 2D Gaussian kernel centered in a square.
    
    Parameters
    ----------
    size_px : int
        Size of the square kernel in pixels
    sigma_px : float
        Standard deviation of the Gaussian in pixels
        
    Returns
    -------
    kernel : ndarray
        2D normalized Gaussian kernel of shape (size_px, size_px)
        
    Examples
    --------
    >>> kern = gaussian_kernel(15, 2.0)
    >>> kern.shape
    (15, 15)
    >>> np.isclose(kern.sum(), 1.0)
    True
    """
    half = size_px // 2
    y, x = np.mgrid[-half:half+1, -half:half+1]
    g = np.exp(-(x**2 + y**2) / (2 * sigma_px**2))
    g /= g.sum() + 1e-12
    return g


def stamp_atoms(image, positions, kernel):
    """
    Add Gaussian atoms to an image at specified positions.
    
    Uses simple interpolation for subpixel positioning.
    
    Parameters
    ----------
    image : ndarray
        2D array to stamp atoms onto (modified in-place)
    positions : ndarray
        Array of shape (N, 2) with (y, x) positions in pixels
    kernel : ndarray
        2D Gaussian kernel to stamp at each position
        
    Returns
    -------
    image : ndarray
        Modified image with atoms stamped
        
    Examples
    --------
    >>> img = np.zeros((100, 100))
    >>> positions = np.array([[50.0, 50.0], [25.0, 75.0]])
    >>> kernel = gaussian_kernel(15, 2.0)
    >>> img = stamp_atoms(img, positions, kernel)
    """
    ky, kx = kernel.shape
    hy, hx = ky // 2, kx // 2
    ny, nx = image.shape

    for (y0, x0) in positions:
        iy = int(np.floor(y0))
        ix = int(np.floor(x0))
        
        # Bounds for pasting
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
            
    return image


def render_lattice(ny, nx, positions, sigma_px):
    """
    Render a complete lattice image with Gaussian atoms.
    
    Parameters
    ----------
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    positions : ndarray
        Array of shape (N, 2) with (y, x) positions in pixels
    sigma_px : float
        Standard deviation of Gaussian atoms in pixels
        
    Returns
    -------
    image : ndarray
        2D array with rendered lattice
        
    Examples
    --------
    >>> from lattice_fft.lattice import make_lattice_positions
    >>> positions = make_lattice_positions(512, 512, 10.0)
    >>> img = render_lattice(512, 512, positions, sigma_px=1.6)
    >>> img.shape
    (512, 512)
    """
    img = np.zeros((ny, nx), dtype=float)
    
    # Kernel size ~ 6 sigma to capture most mass
    ksize = int(max(7, np.ceil(6 * sigma_px)))
    if ksize % 2 == 0:
        ksize += 1
        
    kern = gaussian_kernel(ksize, sigma_px)
    stamp_atoms(img, positions, kern)
    
    return img
