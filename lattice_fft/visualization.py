"""
Visualization utilities for lattice and reciprocal space analysis.

This module provides plotting functions for visualizing lattice structures
and their Fourier transforms.
"""

import numpy as np
import matplotlib.pyplot as plt
from .fft_utils import reciprocal_axes


def plot_lattice(image, title="Lattice (real space)", figsize=(6, 6), 
                 cmap="gray", show_colorbar=True):
    """
    Plot a single lattice image.
    
    Parameters
    ----------
    image : ndarray
        2D lattice image
    title : str, optional
        Plot title
    figsize : tuple, optional
        Figure size (width, height)
    cmap : str, optional
        Colormap name, default 'gray'
    show_colorbar : bool, optional
        Whether to show colorbar, default True
        
    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(image, cmap=cmap, origin="lower")
    ax.set_title(title)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    
    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Intensity (a.u.)")
        
    plt.tight_layout()
    return fig, ax


def plot_fft(fft_image, px_size_A, title="FFT magnitude (log)", 
             figsize=(6, 6), cmap="inferno", show_colorbar=True, 
             klim=None):
    """
    Plot FFT with reciprocal space axes.
    
    Parameters
    ----------
    fft_image : ndarray
        2D FFT magnitude (typically log-scaled)
    px_size_A : float
        Pixel size in Angstroms per pixel
    title : str, optional
        Plot title
    figsize : tuple, optional
        Figure size (width, height)
    cmap : str, optional
        Colormap name, default 'inferno'
    show_colorbar : bool, optional
        Whether to show colorbar, default True
    klim : float or None, optional
        Reciprocal space axis limit (±klim in 1/Å)
        
    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    ny, nx = fft_image.shape
    kx, ky = reciprocal_axes(ny, nx, px_size_A)
    extent = [kx[0], kx[-1], ky[0], ky[-1]]
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(fft_image, cmap=cmap, origin="lower", 
                   extent=extent, aspect="auto")
    ax.set_title(title)
    ax.set_xlabel("kx (1/Å)")
    ax.set_ylabel("ky (1/Å)")
    
    if klim is not None:
        ax.set_xlim(-klim, klim)
        ax.set_ylim(-klim, klim)
    
    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log10 |F|")
        
    plt.tight_layout()
    return fig, ax


def plot_2x2_comparison(perfect_img, perfect_fft, defect_img, defect_fft, 
                       px_size_A, title_suffix="", figsize=(10, 10)):
    """
    Plot 2x2 comparison of perfect vs defected lattice and their FFTs.
    
    Parameters
    ----------
    perfect_img : ndarray
        Perfect lattice real-space image
    perfect_fft : ndarray
        FFT of perfect lattice (log-scaled)
    defect_img : ndarray
        Defected lattice real-space image
    defect_fft : ndarray
        FFT of defected lattice (log-scaled)
    px_size_A : float
        Pixel size in Angstroms per pixel
    title_suffix : str, optional
        Additional text for main title
    figsize : tuple, optional
        Figure size (width, height)
        
    Returns
    -------
    fig, axes : matplotlib Figure and array of Axes
    """
    kx, ky = reciprocal_axes(*perfect_img.shape, px_size_A)
    extent_fft = [kx[0], kx[-1], ky[0], ky[-1]]

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # 1. Perfect lattice
    ax = axes[0, 0]
    im0 = ax.imshow(perfect_img, cmap="gray", origin="lower")
    ax.set_title("1) Perfect lattice (real space)")
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    fig.colorbar(im0, ax=ax, fraction=0.046, pad=0.04, label="Intensity (a.u.)")

    # 2. FFT of perfect lattice
    ax = axes[0, 1]
    im1 = ax.imshow(perfect_fft, cmap="inferno", origin="lower", 
                    extent=extent_fft, aspect="auto")
    ax.set_title("2) FFT: Bragg reflections (kx, ky)")
    ax.set_xlabel("kx (1/Å)")
    ax.set_ylabel("ky (1/Å)")
    fig.colorbar(im1, ax=ax, fraction=0.046, pad=0.04, label="log10 |F|")

    # 3. Defected lattice
    ax = axes[1, 0]
    im2 = ax.imshow(defect_img, cmap="gray", origin="lower")
    ax.set_title("3) Defected lattice (dislocations/vacancies)")
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label="Intensity (a.u.)")

    # 4. FFT of defected lattice
    ax = axes[1, 1]
    im3 = ax.imshow(defect_fft, cmap="inferno", origin="lower", 
                    extent=extent_fft, aspect="auto")
    ax.set_title("4) FFT: symmetry breaking / diffuse scattering")
    ax.set_xlabel("kx (1/Å)")
    ax.set_ylabel("ky (1/Å)")
    fig.colorbar(im3, ax=ax, fraction=0.046, pad=0.04, label="log10 |F|")

    # Add caption
    caption = ("Perfect periodicity → discrete symmetric peaks; "
               "defects → smeared/streaked intensity (order → disorder)")
    if title_suffix:
        caption += f"\n{title_suffix}"
        
    fig.suptitle(caption, y=0.98)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    return fig, axes


def plot_unit_cell_comparison(images, ffts, periods_px, labels, px_size_A,
                              view_px=None, kmax=1.0, figsize=(12, 13)):
    """
    Plot multiple lattices with different periodicities side-by-side.
    
    Parameters
    ----------
    images : list of ndarray
        List of real-space lattice images
    ffts : list of ndarray
        List of corresponding FFT images (log-scaled)
    periods_px : list of float
        List of lattice periods in pixels
    labels : list of str
        List of descriptive labels for each case
    px_size_A : float
        Pixel size in Angstroms per pixel
    view_px : float or None, optional
        Real-space window size for zoomed view
    kmax : float, optional
        Reciprocal space axis limit, default 1.0 (1/Å)
    figsize : tuple, optional
        Figure size (width, height)
        
    Returns
    -------
    fig, axes : matplotlib Figure and array of Axes
    """
    n = len(images)
    fig, axes = plt.subplots(n, 2, figsize=figsize)
    
    ny, nx = images[0].shape
    kx, ky = reciprocal_axes(ny, nx, px_size_A)
    extent_fft = [kx[0], kx[-1], ky[0], ky[-1]]
    
    for i, (img, fft_img, per, lab) in enumerate(zip(images, ffts, periods_px, labels)):
        # Real-space panel (left)
        axL = axes[i, 0] if n > 1 else axes[0]
        axL.imshow(img, cmap='gray', origin='lower')
        axL.set_title(f'Real space: {lab}')
        
        if view_px is not None:
            axL.set_xlim(0, view_px)
            axL.set_ylim(0, view_px)
            
        axL.set_xlabel('x (px)')
        axL.set_ylabel('y (px)')
        
        # Reciprocal-space panel (right)
        axR = axes[i, 1] if n > 1 else axes[1]
        axR.imshow(fft_img, cmap='inferno', origin='lower', 
                   extent=extent_fft, aspect='auto')
        axR.set_title('FFT magnitude (log)')
        axR.set_xlabel('kx (1/Å)')
        axR.set_ylabel('ky (1/Å)')
        axR.set_xlim(-kmax, kmax)
        axR.set_ylim(-kmax, kmax)
    
    fig.tight_layout()
    return fig, axes
