"""
Fast Fourier Transform utilities for reciprocal space analysis.

This module provides functions to compute and analyze 2D FFTs
of crystal lattice images.
"""

import numpy as np
from numpy.fft import fft2, fftshift, fftfreq


def compute_fft(image, log_scale=True):
    """
    Compute 2D FFT of an image with optional log scaling.
    
    Parameters
    ----------
    image : ndarray
        2D input image
    log_scale : bool, optional
        If True, return log10 of magnitude, default True
        
    Returns
    -------
    fft_result : ndarray
        FFT magnitude (log-scaled if log_scale=True)
        
    Examples
    --------
    >>> img = np.random.rand(512, 512)
    >>> fft_mag = compute_fft(img)
    >>> fft_mag.shape
    (512, 512)
    """
    F = fftshift(fft2(image))
    mag = np.abs(F)
    
    if log_scale:
        return np.log10(mag + 1e-6)
    else:
        return mag


def get_fft_magnitude(image, log_scale=False):
    """
    Get raw FFT magnitude without log scaling (alias for compute_fft).
    
    Parameters
    ----------
    image : ndarray
        2D input image
    log_scale : bool, optional
        If True, return log10 of magnitude, default False
        
    Returns
    -------
    magnitude : ndarray
        FFT magnitude
    """
    return compute_fft(image, log_scale=log_scale)


def reciprocal_axes(ny, nx, px_size_A):
    """
    Generate reciprocal space axes in inverse Angstroms.
    
    Parameters
    ----------
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    px_size_A : float
        Pixel size in Angstroms per pixel
        
    Returns
    -------
    kx : ndarray
        Horizontal reciprocal space axis (1/Angstrom)
    ky : ndarray
        Vertical reciprocal space axis (1/Angstrom)
        
    Examples
    --------
    >>> kx, ky = reciprocal_axes(512, 512, px_size_A=0.5)
    >>> len(kx), len(ky)
    (512, 512)
    >>> kx.max()  # Nyquist frequency
    1.0
    """
    # Frequency in cycles per pixel -> convert to inverse Angstrom
    ky = fftshift(fftfreq(ny, d=px_size_A))  # 1/Angstrom
    kx = fftshift(fftfreq(nx, d=px_size_A))  # 1/Angstrom
    return kx, ky


def get_reciprocal_extent(ny, nx, px_size_A):
    """
    Get extent for matplotlib imshow in reciprocal space.
    
    Parameters
    ----------
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    px_size_A : float
        Pixel size in Angstroms per pixel
        
    Returns
    -------
    extent : list
        [left, right, bottom, top] for imshow extent parameter
        
    Examples
    --------
    >>> extent = get_reciprocal_extent(512, 512, 0.5)
    >>> len(extent)
    4
    """
    kx, ky = reciprocal_axes(ny, nx, px_size_A)
    return [kx[0], kx[-1], ky[0], ky[-1]]
