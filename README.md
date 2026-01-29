# 🌌 MOND 3D Physics Engine
> A high-performance C-based physics simulator implementing **Modified Newtonian Dynamics**.

## 📊 Simulation Results
The engine simulates galactic evolution without the need for Dark Matter by modifying the gravitational law at low accelerations ( < a_0$).

### Galaxy Evolution (2.0 Gyr)
| Initial State (T=0) | Final State (T=2.0 Gyr) |
| :---: | :---: |
| <img src="https://raw.githubusercontent.com/DMDTague/mond-physics-engine/main/snapshot_000.png" width="400"> | <img src="https://raw.githubusercontent.com/DMDTague/mond-physics-engine/main/snapshot_009.png" width="400"> |

### Physics Validation: Rotation Curves
<img src="https://raw.githubusercontent.com/DMDTague/mond-physics-engine/main/rotation_curves_analysis.png" width="800">

*The plot above demonstrates the "Flat Rotation Curve" phenomenon, where star velocity remains constant even at extreme radii.*

---

## ➗ Core Theory
We utilize the **QUMOND** (Quasi-Linear MOND) formulation:

45193 \mathbf{g} = \nu\left(\frac{|\mathbf{g}_N|}{a_0}\right) \mathbf{g}_N 45193

Using the simple interpolation function:
45193 \nu(y) = \frac{1}{2} + \sqrt{\frac{1}{4} + \frac{1}{y}} 45193

## 🛠 Tech Stack
- **Language:** C11 (Engine), Python 3 (Visualization)
- **Math:** Symplectic Leapfrog Integration
- **Platform:** Optimized for Apple Silicon (ARM64)

---
