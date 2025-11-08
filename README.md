# Crystal Lattice and Diffraction: From Order to Disorder

## Abstract

This repository provides a comprehensive computational demonstration of Fourier-space analysis of periodic crystal lattices and the effects of structural defects on reciprocal-space symmetry. Using two-dimensional Fast Fourier Transform (FFT) techniques, we illustrate how periodic atomic arrangements produce discrete Bragg reflections in reciprocal space, and how disorder—in the form of dislocations, vacancies, and thermal motion—manifests as diffuse scattering and symmetry breaking.

## Theoretical Background

### Fourier Analysis and Reciprocal Space

The relationship between real-space crystal structure and reciprocal-space diffraction patterns is governed by the Fourier transform. For a discrete 2D lattice, the forward transform is:

$$F[u,v] = \sum_{y=0}^{N_y-1} \sum_{x=0}^{N_x-1} f[y,x]\, e^{-\,i\,2\pi \left( \tfrac{u}{N_x}x + \tfrac{v}{N_y}y \right)}$$

where $f[y,x]$ represents the real-space intensity distribution and $F[u,v]$ its reciprocal-space representation. The spatial frequencies are given by:

$$f_x(u) = \frac{u}{N_x\,\Delta x}, \quad f_y(v) = \frac{v}{N_y\,\Delta y}$$

with pixel size $\Delta x = \Delta y$ in Ångströms per pixel.

### Reciprocal Lattice and Bragg Reflections

A perfectly periodic lattice with real-space unit cell dimension $d$ produces discrete diffraction peaks at reciprocal lattice positions:

$$\mathbf{G}_{hk} = h\,\mathbf{b}_1 + k\,\mathbf{b}_2$$

where $|\mathbf{b}_i| = 2\pi/a_i$ for primitive lattice vectors $\mathbf{a}_i$. This inverse relationship is fundamental:

- **Smaller real-space periodicity** ($d \downarrow$) → **larger reciprocal-space separation** ($1/d \uparrow$)
- **Larger real-space periodicity** ($d \uparrow$) → **smaller reciprocal-space separation** ($1/d \downarrow$)

### Atomic Form Factor and Gaussian Approximation

Individual atoms are modeled as Gaussian electron density distributions:

$$g(\mathbf{r}) \propto \exp\!\left(-\frac{\|\mathbf{r}\|^2}{2\sigma^2}\right)$$

The Fourier transform of a Gaussian remains Gaussian:

$$G(\mathbf{f}) \propto \exp\!\left(-2\pi^2\,\sigma^2\,\|\mathbf{f}\|^2\right)$$

This produces a natural high-frequency envelope that attenuates Bragg peak intensities at larger scattering vectors, consistent with atomic form factor behavior.

### Disorder and Diffuse Scattering

Structural imperfections introduce deviations from perfect periodicity:

1. **Thermal Disorder (Debye-Waller Effect)**: Random atomic displacements with standard deviation $\sigma_j$ attenuate coherent Bragg intensity:
   $$I \propto \exp\left(-2\pi^2\sigma_j^2\,\|\mathbf{f}\|^2\right)$$

2. **Dislocations**: Systematic row shifts introduce phase ramps that manifest as streaking along specific crystallographic directions in reciprocal space.

3. **Vacancies and Point Defects**: Random removal of lattice sites creates diffuse background scattering, breaking the sharp Bragg condition.

## Implementation

### Simulation Parameters

The notebook employs the following physical and computational parameters:

- **Pixel Size**: 0.5 Å/pixel (subatomic resolution)
- **Lattice Spacing**: 5.43 Å (approximating silicon diamond cubic projected to 2D)
- **Atomic Width**: σ = 0.8 Å (Gaussian standard deviation)
- **Image Dimensions**: 512 × 512 pixels
- **Lattice Types**: Square (simple cubic projection) and hexagonal (close-packed projection)

### Key Functions

```python
gaussian_kernel(size_px, sigma_px)
```
Generates normalized 2D Gaussian representing atomic electron density.

```python
make_lattice_positions(ny, nx, period_px, basis)
```
Creates periodic lattice with specified symmetry (square or hexagonal).

```python
introduce_defects(positions, defect_fraction, max_row_shift_px, remove_fraction)
```
Simulates crystal imperfections: dislocations (row shifts) and vacancies (atom removal).

```python
compute_fft(image)
```
Performs 2D FFT with log-scaling for visualization of dynamic range in diffraction intensities.

## Demonstrations

### 1. One-Dimensional Fourier Analysis

**Objective**: Establish fundamental FFT concepts with simple periodic signals.

- Single sine wave → sharp symmetric peaks at ±f₀
- Composite signal (sum of multiple frequencies) → distinct peaks revealing each periodic component
- **Analogy**: Multiple crystal plane spacings produce multiple diffraction orders

### 2. Two-Dimensional Perfect Lattice

**Objective**: Visualize reciprocal space for ideal periodic structures.

**Results**:
- Square lattice → square reciprocal lattice with Bragg peaks at $(h/d, k/d)$
- Hexagonal lattice → hexagonal reciprocal pattern with six-fold symmetry
- Log-scale FFT reveals both fundamental and higher-order reflections

### 3. Unit Cell and Reciprocal Spacing Relationship

**Objective**: Quantitatively demonstrate inverse relationship between real and reciprocal space.

**Experimental Design**:
- Generate lattices with varied periodicity: 0.7d, d, 1.5d
- Fixed real-space field of view for direct comparison
- Standardized reciprocal-space axes (±1.0 Å⁻¹)

**Observations**:
- Smaller d → Bragg peaks farther apart in k-space
- Larger d → peaks closer together
- Peak positions measured along kₓ, k_y, and diagonal directions validate $\mathbf{k} = 2\pi/d$

### 4. Defect Analysis

**Objective**: Characterize how structural disorder affects diffraction patterns.

**Defect Types Introduced**:
- **Row Dislocations**: 5% of rows shifted by ±2 pixels
- **Vacancies**: 3% random atom removal
- **Thermal Jitter**: Gaussian position perturbations (σ_jitter controllable)

**Diffraction Consequences**:
- Smearing and broadening of Bragg peaks
- Anisotropic streaking along dislocation directions
- Elevated diffuse background intensity
- Loss of higher-order reflections due to disorder-induced coherence loss

### 5. Systematic Disorder: Jitter Parameter Sweep

**Objective**: Quantify relationship between atomic displacement magnitude and reciprocal-space peak degradation.

**Method**: Apply Gaussian jitter with increasing σ_jitter from 0 to 2 Å, observe:
- Progressive peak broadening (Debye-Waller-like attenuation)
- Reduction in peak-to-background ratio
- Transition from crystalline to amorphous-like diffraction patterns

## Physical Significance

### Semiconductor Metrology Applications

This FFT-based analysis methodology is directly applicable to:

1. **Wafer Crystal Quality Assessment**: Detection of strain, dislocations, and epitaxial mismatch through reciprocal-space mapping
2. **Device Fabrication Monitoring**: Identification of process-induced defects (ion implantation damage, etching artifacts)
3. **Transmission Electron Microscopy (TEM)**: Selected-area diffraction pattern interpretation
4. **X-ray Diffraction (XRD)**: Correlation between reciprocal-space peaks and real-space crystal structure

### Connection to Experimental Techniques

The simulated FFT patterns directly correspond to:
- **XRD**: Bragg peak positions, widths (crystallite size via Scherrer equation), and integrated intensities
- **Electron Diffraction**: Ring patterns (polycrystalline), spot patterns (single crystal)
- **High-Resolution TEM**: Fourier filtering for lattice imaging and strain mapping

## Dependencies

```python
numpy          # Array operations and FFT
matplotlib     # Visualization
scipy          # (Optional) Peak fitting and analysis
```

## Usage

Open and execute the Jupyter notebook sequentially:

```bash
jupyter notebook "FFT and lattice.ipynb"
```

Each cell is self-contained with explanatory markdown. Adjust parameters in the configuration cells:
- `px_size_A`: Spatial resolution (Å/pixel)
- `a_A`: Lattice constant
- `atom_sigma_A`: Atomic width
- `defect_fraction`, `remove_fraction`: Disorder magnitude

## Educational Value

This notebook serves as:
- **Pedagogical tool** for solid-state physics and crystallography courses
- **Visualization platform** for reciprocal space concepts often challenging to grasp from theory alone
- **Computational lab** for exploring disorder effects without requiring expensive experimental equipment
- **Foundation** for more advanced topics: structure factor calculations, kinematical diffraction theory, orientation relationships

## Future Extensions

Potential enhancements include:
- 3D lattice generation with volumetric FFT
- Dynamical diffraction theory (multiple scattering)
- Orientation texture and polycrystalline averaging
- Quantitative structure factor calculations
- Strain field simulation and geometric phase analysis
- Integration with experimental diffraction data (comparison workflows)

## References

- Ashcroft, N. W., & Mermin, N. D. (1976). *Solid State Physics*. Holt, Rinehart and Winston.
- Cullity, B. D., & Stock, S. R. (2001). *Elements of X-Ray Diffraction* (3rd ed.). Prentice Hall.
- Kittel, C. (2004). *Introduction to Solid State Physics* (8th ed.). Wiley.
- Williams, D. B., & Carter, C. B. (2009). *Transmission Electron Microscopy: A Textbook for Materials Science* (2nd ed.). Springer.

## License

MIT License - Free for educational and research purposes.

## Author

EigenJames  
Repository: https://github.com/EigenJames/Lattice_FFT

## Acknowledgments

This work demonstrates fundamental crystallographic principles applicable to semiconductor device characterization, materials science research, and solid-state physics education.
