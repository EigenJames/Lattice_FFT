"""Plotting helpers for lattice FFT simulations."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_real_and_fft(
    image: np.ndarray,
    fft_magnitude: np.ndarray,
    *,
    kx: np.ndarray,
    ky: np.ndarray,
    title_suffix: str = "",
    cmap_real: str = "gray",
    cmap_fft: str = "inferno",
    show_colorbar: bool = True,
    figsize: tuple[int, int] = (10, 5),
    fft_is_log: bool = True,
    log_eps: float = 1e-6,
) -> plt.Figure:
    """Create a side-by-side plot of real-space image and FFT magnitude.

    If ``fft_is_log`` is ``False`` the FFT magnitude is log-scaled internally
    using ``log_eps`` for numerical stability before plotting.
    """
    if fft_is_log:
        fft_display = fft_magnitude
        colorbar_label = "log10 |F|"
    else:
        fft_display = np.log10(np.abs(fft_magnitude) + log_eps)
        colorbar_label = "log10 |F|" if show_colorbar else ""

    extent = [kx.min(), kx.max(), ky.min(), ky.max()]

    fig, (ax_real, ax_fft) = plt.subplots(1, 2, figsize=figsize)

    im_real = ax_real.imshow(image, cmap=cmap_real, origin="lower")
    ax_real.set_title("Real space lattice" + (f" — {title_suffix}" if title_suffix else ""))
    ax_real.set_xlabel("x (px)")
    ax_real.set_ylabel("y (px)")
    if show_colorbar:
        fig.colorbar(im_real, ax=ax_real, fraction=0.046, pad=0.04, label="Intensity (a.u.)")

    im_fft = ax_fft.imshow(
        fft_display,
        cmap=cmap_fft,
        origin="lower",
        extent=extent,
        aspect="auto",
    )
    ax_fft.set_title("FFT magnitude" + (f" — {title_suffix}" if title_suffix else ""))
    ax_fft.set_xlabel("kx (1/Å)")
    ax_fft.set_ylabel("ky (1/Å)")
    if show_colorbar:
        fig.colorbar(im_fft, ax=ax_fft, fraction=0.046, pad=0.04, label=colorbar_label)

    fig.tight_layout()
    return fig
