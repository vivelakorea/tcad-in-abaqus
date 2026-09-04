# TCAD in Abaqus — Semiconductor Device Simulation with User Elements

**Abaqus 사용자 요소(UEL)로 구현한 반도체 소자 시뮬레이션(TCAD): pn 다이오드와 3D
MOSFET의 드리프트-확산(drift–diffusion) 해석.**

Semiconductor device simulation (drift–diffusion TCAD) implemented as **Abaqus user
elements (UEL)** — the structural FE solver becomes a device simulator. Nodal degrees of
freedom are electrostatic potential and carrier densities instead of displacements; the
discretization is the same edge-based Scharfetter–Gummel box method used by real TCAD
codes (Silvaco, Sentaurus). Everything is verified against classic papers and independent
Python reference implementations.

## Demo — MOSFET turn-on, live from Abaqus

Gate-voltage sweep 0 → 3 V (V_D = 50 mV) in a single Abaqus job: the inversion
channel switches on at the surface as V_G crosses threshold (~1.3 V), and the
UEL drain current traces the Pao–Sah (1966) exact transfer curve point by point.

![MOSFET turn-on demo](docs/demo_turnon.gif)

Reproduce with `python mosfet/make_demo.py`.

![3D NMOS results](docs/fig_tcad_uel.png)

## Contents

| | model | UEL | verification |
|---|---|---|---|
| `diode/` | 1D pn junction, drift–diffusion (ψ, n, p) | 2-node element, Scharfetter–Gummel (1969) exponential-fitting flux, analytic unsymmetric Jacobian | node-wise 5×10⁻⁸ agreement with an independent Python solver; built-in potential, mass-action law, current conservation |
| `mosfet/` | 3D long-channel NMOS (6×0.5×2 µm, n⁺ 10¹⁹ S/D, p-substrate 10¹⁷, 10 nm oxide, L = 4 µm) | 8-node hex box-method element with SG fluxes on all 12 edges, quasi-Fermi variables (ψ, φₙ, φₚ) | drain current vs. Pao–Sah (1966) exact double integral and Brews (1978) charge-sheet model |

### MOSFET drain current: 3D UEL vs. the papers

| V_G [V] | V_D [V] | I_D UEL [µA] | Pao–Sah 1966 | Brews 1978 | dev (UEL/PS) |
|---|---|---|---|---|---|
| 2.0 | 0.05 | 3.61 | 3.53 | 3.41 | 2.3 % |
| 2.0 | 0.50 | 19.1 | 18.0 | 17.2 | 6.3 % |
| 3.0 | 0.50 | 83.8 | 82.4 | 80.7 | 1.7 % |
| 3.0 | 0.05 | 10.2 | 10.1 | 9.9 | 1.0 % |

The largest deviation sits at the saturation knee (V_G − V_T0 ≈ 0.7 V), exactly where the
gradual-channel assumption behind Pao–Sah starts to strain — the 2D device solution does
not make that assumption, so diverging there is correct behavior, not error.

## How to run

Requires Abaqus (tested with 2024) with a linked Fortran compiler, plus Python 3 with
numpy (and matplotlib for the figure).

```
cd diode
python run_diode.py     # ~1 min: writes inp, runs abaqus job=... user=uel_dd.f, checks

cd mosfet
python run_mosfet.py    # ~3 min: full bias sweep in one job (6 steps), figure + checks
```

Each driver generates the mesh/inp, launches Abaqus, parses the `.dat` output, and
asserts the physics checks — if it prints `check passed`, everything reproduced.

## Implementation notes (the parts that bite)

- **No initial conditions in Abaqus statics** → starting from U = 0 with full doping
  charge makes Newton explode. Fix: ramp the doping with step time *inside the UEL*
  (the classic TCAD doping-continuation trick).
- **(ψ, n, p) unknowns diverge at real doping levels** (n/nᵢ ~ 10⁹). Fix: switch to
  quasi-Fermi potentials, n = nᵢ e^(ψ−φₙ), p = nᵢ e^(φₚ−ψ) — all unknowns become
  bounded voltages, and contact BCs become plain voltage ramps. Also standard practice
  in production device solvers.
- **`DEXP(T) - 1.D0` is not `expm1`**: catastrophic cancellation at |t| ~ 10⁻⁸ puts a
  10⁻⁴ floor under the Bernoulli-function fluxes. Use a series branch for small |t|.
- **Terminal current = reaction force.** At Dirichlet-constrained contact dofs, the
  Abaqus reaction (RF) of the continuity equations *is* the discrete carrier flux into
  the contact: I = q·nᵢ·ΣRF. No flux post-processing needed, and it is exactly
  conservative.
- The coupled Jacobian is unsymmetric — `UNSYMM` on the `*USER ELEMENT` line and
  `UNSYMM=YES` on every step, or convergence quietly degrades.

## References

- D. L. Scharfetter, H. K. Gummel, *Large-signal analysis of a silicon Read diode
  oscillator*, IEEE Trans. Electron Devices 16 (1969) 64.
- H. C. Pao, C. T. Sah, *Effects of diffusion current on characteristics of
  metal-oxide (insulator)-semiconductor transistors*, Solid-State Electron. 9 (1966) 927.
- J. R. Brews, *A charge-sheet model of the MOSFET*, Solid-State Electron. 21 (1978) 345.

## Author

**심규장 (Gyu-Jang Sim)** — 서울대학교 재료역학연구실 ([MAMEL](https://mamel.snu.ac.kr), Seoul National University)
[github.com/vivelakorea](https://github.com/vivelakorea) · gyujang95@gmail.com

Formulation, Fortran UELs, Python reference solvers, and verification in this
repository are by the author. Background: computational solid mechanics
(crystal-plasticity FEM, UMAT/VUMAT development for anisotropic yield and
distortional hardening models) — this project applies the same FE machinery,
Newton solvers, and consistent-linearization discipline to semiconductor devices.

## License

MIT
