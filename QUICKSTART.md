# Lattice FFT Package - Quick Reference

## Installation

```bash
# Clone the repository
git clone https://github.com/EigenJames/Lattice_FFT.git
cd Lattice_FFT

# Install in editable mode
pip install -e .
```

## Quick Start

```python
from lattice_fft import (
    make_lattice_positions,
    introduce_defects,
    render_lattice,
    compute_fft,
    plot_2x2_comparison
)
import numpy as np

# Setup parameters
ny, nx = 512, 512
px_size_A = 0.5
period_px = 10.86
sigma_px = 1.6
rng = np.random.default_rng(42)

# Generate perfect lattice
pos_perfect = make_lattice_positions(ny, nx, period_px, basis='square')
img_perfect = render_lattice(ny, nx, pos_perfect, sigma_px)
fft_perfect = compute_fft(img_perfect)

# Generate defected lattice
pos_defect = introduce_defects(pos_perfect, ny, nx, 
                               defect_fraction=0.05,
                               remove_fraction=0.03, 
                               rng=rng)
img_defect = render_lattice(ny, nx, pos_defect, sigma_px)
fft_defect = compute_fft(img_defect)

# Visualize comparison
plot_2x2_comparison(img_perfect, fft_perfect, 
                    img_defect, fft_defect, 
                    px_size_A)
```

## Module Overview

### lattice.py - Lattice Generation
- `make_lattice_positions(ny, nx, period_px, basis='square', rng=None)`
  - Creates square or hexagonal lattice positions
- `introduce_defects(positions, ny, nx, defect_fraction=0.05, max_row_shift_px=2.0, remove_fraction=0.03, rng=None)`
  - Adds dislocations and vacancies
- `add_jitter(positions, sigma_jitter_px, rng=None)`
  - Simulates thermal disorder

### atoms.py - Atomic Rendering
- `gaussian_kernel(size_px, sigma_px)`
  - Creates 2D Gaussian kernel
- `stamp_atoms(image, positions, kernel)`
  - Places atoms at positions
- `render_lattice(ny, nx, positions, sigma_px)`
  - Generates complete lattice image

### fft_utils.py - Reciprocal Space
- `compute_fft(image, log_scale=True)`
  - Computes 2D FFT with optional log scaling
- `reciprocal_axes(ny, nx, px_size_A)`
  - Generates k-space axes in 1/Å
- `get_fft_magnitude(image, log_scale=False)`
  - Extracts FFT amplitude

### visualization.py - Plotting
- `plot_lattice(image, title, figsize, cmap, show_colorbar)`
  - Plots single lattice image
- `plot_fft(fft_image, px_size_A, title, figsize, cmap, show_colorbar, klim)`
  - Plots FFT with reciprocal axes
- `plot_2x2_comparison(perfect_img, perfect_fft, defect_img, defect_fft, px_size_A, title_suffix, figsize)`
  - 2x2 perfect vs defected comparison
- `plot_unit_cell_comparison(images, ffts, periods_px, labels, px_size_A, view_px, kmax, figsize)`
  - Multi-period comparison

## Common Workflows

### 1. Compare Different Lattice Types
```python
# Square lattice
pos_sq = make_lattice_positions(512, 512, 10.86, basis='square')
img_sq = render_lattice(512, 512, pos_sq, 1.6)

# Hexagonal lattice
pos_hex = make_lattice_positions(512, 512, 10.86, basis='hex')
img_hex = render_lattice(512, 512, pos_hex, 1.6)
```

### 2. Study Defect Effects
```python
defect_levels = [0.0, 0.05, 0.1, 0.2]
for frac in defect_levels:
    pos_def = introduce_defects(pos_perfect, ny, nx, 
                                defect_fraction=frac, 
                                rng=rng)
    img = render_lattice(ny, nx, pos_def, sigma_px)
    fft = compute_fft(img)
    # Analyze...
```

### 3. Thermal Disorder Analysis
```python
from lattice_fft import add_jitter

jitter_values = [0.0, 0.5, 1.0, 2.0]
for sigma in jitter_values:
    pos_jit = add_jitter(pos_perfect, sigma, rng=rng)
    img = render_lattice(ny, nx, pos_jit, sigma_px)
    fft = compute_fft(img)
    # Analyze peak broadening...
```

### 4. Unit Cell Scaling Study
```python
scales = [0.7, 1.0, 1.5]
for scale in scales:
    pos = make_lattice_positions(ny, nx, period_px * scale, basis='square')
    img = render_lattice(ny, nx, pos, sigma_px)
    fft = compute_fft(img)
    # Observe 1/d relationship...
```

## Notebooks

1. **`FFT and lattice.ipynb`**: Original comprehensive analysis
   - Mathematical derivations
   - Step-by-step explanations
   - Extensive visualizations

2. **`demo.ipynb`**: Package demonstration
   - Concise examples
   - Clean API usage
   - Quick reference

## Parameters Guide

### Physical Parameters
- `px_size_A`: Pixel size in Ångströms (typical: 0.1-1.0)
- `a_A`: Lattice constant in Ångströms (Si: 5.43)
- `atom_sigma_A`: Atomic width in Ångströms (typical: 0.5-1.5)

### Defect Parameters
- `defect_fraction`: Fraction of rows with dislocations (0.0-0.3)
- `max_row_shift_px`: Maximum dislocation shift in pixels (1-5)
- `remove_fraction`: Vacancy rate (0.0-0.1)
- `sigma_jitter_px`: Thermal jitter std dev in pixels (0.0-3.0)

### Image Parameters
- `ny, nx`: Image dimensions (256-1024 typical)
- Larger images → better frequency resolution
- More unit cells → sharper Bragg peaks

## Advanced: Custom Analysis

```python
def analyze_peak_quality(fft_magnitude, kx, ky):
    """Custom function to extract peak metrics."""
    cy, cx = len(ky) // 2, len(kx) // 2
    
    # Find peak positions
    from scipy.signal import find_peaks
    profile_x = fft_magnitude[cy, :]
    peaks_x, _ = find_peaks(profile_x, height=threshold)
    
    # Calculate peak widths, spacings, etc.
    # ...
    return metrics

# Use in analysis loop
for condition in conditions:
    img = generate_image(condition)
    fft_mag = get_fft_magnitude(img, log_scale=False)
    kx, ky = reciprocal_axes(ny, nx, px_size_A)
    metrics = analyze_peak_quality(fft_mag, kx, ky)
```

## Tips

1. **Reproducibility**: Always pass `rng=np.random.default_rng(seed)` for reproducible defects/jitter
2. **Memory**: Large images (>1024²) may require significant RAM for FFT
3. **Visualization**: Use `klim` parameter in `plot_fft()` to focus on central peaks
4. **Performance**: Pre-compute kernels for batch rendering

## Troubleshooting

**Import errors**: Ensure package is installed (`pip install -e .`)

**Memory issues**: Reduce image size or use smaller batches

**Weak FFT peaks**: Increase number of unit cells or reduce atom sigma

**Noisy FFT**: Ensure proper normalization and log scaling

## References

See main README.md for theoretical background and literature references.
