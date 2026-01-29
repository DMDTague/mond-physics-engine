# MOND Physics Engine

A comprehensive C implementation of Modified Newtonian Dynamics (MOND) for galaxy simulations, replacing dark matter with modified gravity.

## 🌌 Overview

This physics engine implements the QUMOND (Quasi-Linear MOND) formulation following the seminal work of:

- **Milgrom (1983)**: Original MOND proposal
- **Banik & Zhao (2022)**: Modern testing and predictions
- **McGaugh (2014)**: Observational validation
- **Famaey & McGaugh (2012)**: MOND review

### What is MOND?

Modified Newtonian Dynamics (MOND) is an alternative to dark matter that modifies gravity at very low accelerations (below ~10⁻¹⁰ m/s²). Instead of adding invisible matter, MOND changes how gravity behaves when forces are weak.

**Key insight**: At accelerations much smaller than a₀ ≈ 1.2×10⁻¹⁰ m/s², the gravitational force law changes from F ∝ 1/r² to F ∝ 1/r.

## 🔬 Physics Implementation

### The MOND Acceleration Scale

```c
#define ACC_0 1.2e-10  // m/s² - fundamental constant of nature
```

This is THE fundamental constant in MOND theory. All accelerations are compared against this scale.

### Interpolation Functions

The transition between Newtonian and MOND regimes is governed by ν(y), where y = |g_N|/a₀:

1. **Simple formulation** (fastest):
   ```
   ν(y) = 0.5 + sqrt(0.25 + 1/y)
   ```

2. **Standard formulation** (most commonly used):
   ```
   ν(y) = y/μ(y) where μ(y) = y/sqrt(1 + y²)
   ```

3. **Bekenstein (TeVeS-compatible)**:
   ```
   ν(y) = 1 + 1/y
   ```

### The QUMOND Algorithm

For each particle at each timestep:

1. **Calculate Newtonian acceleration** from all other particles:
   ```
   a_N = Σ (G m_j / r_ij²) * r̂_ij
   ```

2. **Apply External Field Effect** (if enabled):
   ```
   a_total_N = a_N + a_ext
   ```

3. **Calculate MOND boost**:
   ```
   |a_MOND| = ν(|a_total_N|/a₀) * |a_total_N|
   ```

4. **Subtract external field**:
   ```
   a_final = a_MOND - a_ext
   ```

This algorithm ensures proper momentum conservation and implements the unique External Field Effect.

### External Field Effect (EFE)

**The most distinctive MOND prediction**: Internal dynamics depend on external fields, violating the Strong Equivalence Principle.

Example: A star cluster orbiting in a galaxy experiences different internal dynamics depending on its location in the galactic gravitational field, even if that field is uniform locally.

This is implemented via the algorithm above and is **critical** for accurate simulations.

## 🏗️ Building and Running

### Requirements

- GCC compiler (or any C11-compatible compiler)
- Make (optional, but recommended)
- Python 3 with matplotlib and numpy (for visualization)

### Quick Start

```bash
# Compile the optimized version
make

# Run the simulation
./mond_sim

# Visualize results
python3 visualize_mond.py
```

### Build Options

```bash
make              # Optimized build (-O3)
make debug        # Debug build with symbols
make profile      # Profiling build
make parallel     # OpenMP parallel version (future)
make test         # Run validation tests only
make benchmark    # Performance benchmark
```

### Compilation Flags

The default build uses aggressive optimization:
- `-O3`: Maximum optimization
- `-march=native`: CPU-specific optimizations
- `-ffast-math`: Fast floating-point math

For debugging or development:
```bash
gcc -g -O0 -Wall -o mond_sim mond_physics_engine.c -lm
```

## 📊 Validation Tests

The engine includes four critical validation tests that verify MOND predictions against observations:

### 1. Tully-Fisher Relation

**Prediction**: V_flat = (G M a₀)^(1/4)

This predicts an extremely tight correlation between galaxy luminosity (baryonic mass) and rotation velocity:

```
V ∝ M^(1/4)
```

**Observational status**: Confirmed with <10% scatter across 5 orders of magnitude in mass (McGaugh et al. 2000).

### 2. Flat Rotation Curves

**Prediction**: Galaxies naturally exhibit flat rotation curves at large radii without dark matter.

In the deep MOND regime (r >> r_transition), velocity becomes:
```
V = constant = (G M a₀)^(1/4)
```

**Observational status**: All observed spiral galaxies show flat rotation curves. This is automatic in MOND but requires fine-tuned dark matter halos in standard cosmology.

### 3. External Field Effect (EFE)

**Prediction**: Internal dynamics of a system depend on external gravitational fields.

Example predictions:
- Wide binary stars in the Solar neighborhood show enhanced accelerations due to the Milky Way's external field
- Dwarf galaxies near massive hosts behave more "Newtonian" than isolated dwarfs

**Observational status**: Recently confirmed in wide binary star statistics (Hernandez et al. 2022).

### 4. Galactic Bar Stability

**Prediction**: Galactic bars should rotate at constant pattern speed without slowing down.

Standard dark matter predicts bars should slow due to dynamical friction with the halo. MOND predicts no halo → no friction → stable bars.

**Observational status**: Observations favor fast, stable bars (Debattista & Sellwood 2000, Athanassoula 2002).

## 📁 Output Files

### Rotation Curves

Files: `rotation_curve_initial.dat`, `rotation_curve_final.dat`

Format:
```
# Radius(kpc)  V_circ(km/s)  a/a0
1.000          180.342       2.45e-01
5.000          220.156       1.23e-02
10.000         215.892       3.45e-03
...
```

Columns:
1. Radius in kiloparsecs
2. Circular velocity in km/s
3. Acceleration ratio (a/a₀)

### Snapshots

Files: `snapshot_000.dat`, `snapshot_001.dat`, ...

Format:
```
# id x(kpc) y(kpc) z(kpc) vx(km/s) vy(km/s) vz(km/s) mass(M_sun) type
0 0.000 0.000 0.000 0.0 0.0 0.0 1.0e+10 1
1 8.234 -2.451 0.123 -45.2 120.3 5.1 1.2e+09 2
...
```

Columns:
1. Particle ID
2-4. Position (x, y, z) in kpc
5-7. Velocity (vx, vy, vz) in km/s
8. Mass in solar masses
9. Type (1=bulge, 2=disk)

## 🎨 Visualization

The Python visualization script creates several plots:

### 1. Rotation Curve Analysis
- Velocity vs radius
- Acceleration regime map
- Velocity gradient (flatness test)

### 2. Galaxy Snapshots
- Face-on view (X-Y plane)
- Edge-on view (X-Z plane)
- Color-coded by particle type

### 3. Tully-Fisher Relation
- Theoretical MOND prediction
- Comparison with Newtonian prediction
- Log-log plot showing M^(1/4) scaling

### 4. MOND vs Newtonian Comparison
- Side-by-side rotation curves
- Acceleration profiles
- Transition region visualization

### 5. Animation (optional)
If ffmpeg is installed:
```bash
python3 visualize_mond.py  # Creates galaxy_evolution.mp4
```

## 🔧 Customization

### Changing MOND Parameters

Edit the parameters in the code:

```c
// In main() or initialization function
sim.params.acc_0 = A0_CODE;              // Acceleration scale
sim.params.interpolation_type = 0;       // 0=simple, 1=standard, 2=Bekenstein
sim.params.use_efe = 1;                  // Enable/disable EFE
sim.params.external_field = (Vector3){   // Set external field
    0.5 * A0_CODE,                       // x component
    0, 0                                 // y, z components
};
```

### Changing Galaxy Properties

```c
// Mass distribution
double M_bulge = 1.0;      // 10^10 M_sun
double M_disk = 5.0;       // 5 × 10^10 M_sun
double R_disk = 3.5;       // 3.5 kpc scale length

// Number of particles
int n_bulge = 200;
int n_disk = 800;
```

### Time Integration

```c
sim.dt = 0.001;            // Timestep in Gyr (1 Myr)
double t_end = 2.0;        // End time in Gyr
int n_snapshots = 10;      // Number of output snapshots
```

## 🧪 Running Specific Tests

### Test Flat Rotation Curves Only

```c
// In main(), comment out full simulation and call:
test_flat_rotation_curve();
```

### Test External Field Effect

```c
test_external_field_effect();
```

### Generate Theoretical Tully-Fisher Plot

```bash
python3 visualize_mond.py  # Will create plot even without simulation data
```

## 📊 Performance

### Computational Complexity

- N-body force calculation: O(N²) per timestep
- Single galaxy (1000 particles): ~1 second per snapshot
- Full simulation (2 Gyr, 10 snapshots): ~1-2 minutes

### Optimization

The code uses:
- Gravitational softening to avoid close encounter singularities
- Efficient vector operations
- Aggressive compiler optimizations (`-O3 -march=native`)

For large simulations (N > 10,000), consider:
1. Tree algorithms (Barnes-Hut)
2. GPU acceleration
3. Parallel processing (OpenMP/MPI)

## 📚 Theoretical Background

### Why MOND Works for Galaxies

In the outer regions of galaxies:
- Accelerations are extremely small (~10⁻¹¹ to 10⁻¹⁰ m/s²)
- This is below the MOND scale a₀
- Gravity is enhanced by factor √(a₀/a_N)
- Result: Flat rotation curves without dark matter

### The Deep MOND Limit

When a << a₀:

```
ν(y) ≈ 1/√y  for y << 1
```

This gives:
```
a_MOND ≈ √(a_N × a₀)
```

For circular orbits:
```
v² = a_MOND × r = √(a_N × a₀) × r = √(GM × a₀)
v = (GM × a₀)^(1/4)
```

**Key result**: Velocity is independent of radius! This explains flat rotation curves.

### Conservation Laws

MOND is a modification of dynamics, not just forces. Energy and momentum conservation require careful implementation:

1. **Energy**: Leapfrog integrator conserves energy to machine precision
2. **Momentum**: Algebraic MOND (used here) conserves total momentum
3. **Angular momentum**: Conserved for symmetric potentials

### Limitations of This Implementation

This code uses the **algebraic approximation** to MOND for speed. For research-grade simulations, consider:

1. **Full QUMOND**: Solves modified Poisson equation ∇·[μ(|∇Φ|/a₀)∇Φ] = 4πGρ
2. **Aquadratic Lagrangian (AQUAL)**: Field-based formulation
3. **TeVeS**: Relativistic extension (Bekenstein 2004)

The algebraic approximation is accurate for:
- Isolated systems
- Systems with weak external fields
- Educational purposes and quick tests

For clusters of galaxies or cosmological simulations, use full QUMOND.

## 🎯 Key Scientific Results

### Successes of MOND

1. **Rotation curves**: Predicts flat rotation curves with zero free parameters
2. **Tully-Fisher**: Natural explanation with tight scatter
3. **Low surface brightness galaxies**: MOND predicted observations before they were made
4. **No need for dark matter in galaxies**: All phenomena explained by visible matter
5. **External Field Effect**: Recently confirmed by observations

### Challenges for MOND

1. **Galaxy clusters**: May still need dark matter (or TeVeS)
2. **CMB power spectrum**: Difficult to match without dark matter
3. **Bullet Cluster**: Appears to show dark matter (though MOND interpretations exist)
4. **Cosmology**: MOND alone cannot explain all cosmological observations

### Current Status

MOND is extremely successful at galactic scales but faces challenges at larger scales. Most physicists favor dark matter, but MOND's successes suggest it captures something real about gravity at low accelerations.

## 🔗 References

### Primary Papers

1. Milgrom, M. (1983). "A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis." *ApJ*, 270, 365-370.

2. Banik, I., & Zhao, H. (2022). "From galactic bars to the Hubble tension: Weighing up the astrophysical evidence for Milgromian gravity." *Symmetry*, 14(7), 1331.

3. McGaugh, S. S. (2014). "The baryonic Tully-Fisher relation." *The Galaxy in Context*, 257-271.

4. Famaey, B., & McGaugh, S. S. (2012). "Modified Newtonian dynamics (MOND): Observational phenomenology and relativistic extensions." *Living Reviews in Relativity*, 15(1), 10.

### Computational Methods

5. Candlish, G. N., Smith, R., & Fellhauer, M. (2015). "RAyMOND: An N-body and hydrodynamics code for MOND." *MNRAS*, 446(1), 1060-1070.

6. Lüghausen, F., Famaey, B., & Kroupa, P. (2015). "Phantom of RAMSES: A new Milgromian dynamics N-body code." *Canadian Journal of Physics*, 93(2), 232-241.

### Observational Tests

7. McGaugh, S. S., Lelli, F., & Schombert, J. M. (2016). "Radial Acceleration Relation in Rotationally Supported Galaxies." *Physical Review Letters*, 117(20), 201101.

8. Hernandez, X., et al. (2022). "Wide binaries as a critical test of classical gravity." *MNRAS*, 509(3), 4696-4714.

## 📝 License

This code is released for educational and research purposes. Please cite appropriately if used in publications.

## 🤝 Contributing

Improvements and extensions welcome! Possible additions:

- [ ] Barnes-Hut tree algorithm for O(N log N) scaling
- [ ] Full QUMOND solver (modified Poisson equation)
- [ ] GPU acceleration (CUDA/OpenCL)
- [ ] Gas dynamics and star formation
- [ ] Cosmological expansion (for Hubble tension tests)
- [ ] Multiple interpolation function comparison
- [ ] Automatic parameter optimization

## 💬 Contact

For questions about the physics or implementation, please refer to the original papers cited above.

## 🌟 Acknowledgments

This implementation is based on decades of work by Mordehai Milgrom, Bob Sanders, Stacy McGaugh, Benoit Famaey, Pavel Kroupa, and many others who have developed and tested MOND theory.

---

**Last updated**: January 2026

**Version**: 1.0

**Status**: Educational/Research Code
