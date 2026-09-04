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

![Abaqus UEL MOSFET simulation: inversion channel forming during gate voltage sweep, drain current following the Pao-Sah transfer curve](docs/demo_turnon.gif)

Reproduce with `python mosfet/make_demo.py`.

![3D NMOS TCAD simulation in Abaqus: doping map, electrostatic potential, electron inversion channel, drain current vs Pao-Sah and Brews models](docs/fig_tcad_uel.png)

## Contents

| | model | UEL | verification |
|---|---|---|---|
| `diode/` | 1D pn junction, drift–diffusion (ψ, n, p) | 2-node element, Scharfetter–Gummel (1969) exponential-fitting flux, analytic unsymmetric Jacobian | node-wise 5×10⁻⁸ agreement with an independent Python solver; built-in potential, mass-action law, current conservation |
| `mosfet/` | 3D long-channel NMOS (6×0.5×2 µm, n⁺ 10¹⁹ S/D, p-substrate 10¹⁷, 10 nm oxide, L = 4 µm) | 8-node hex box-method element with SG fluxes on all 12 edges, quasi-Fermi variables (ψ, φₙ, φₚ) | drain current vs. Pao–Sah (1966) exact double integral and Brews (1978) charge-sheet model |
| `mosfet/` (electro-thermal) | same NMOS with **self-heating**: lattice temperature as a 4th nodal dof, steady-state heat equation with edge-lumped Joule heating, V_T(T) and µ(T) feedback (a minimal Wachutka 1990 thermodynamic model) | **monolithic** (ψ, φₙ, φₚ, ΔT) Newton — one unsymmetric 4×4-block Jacobian, not the usual staggered TCAD↔thermal loop | energy balance (heatsink reaction heat = ΣI·V, Tellegen) to 0.01 %; isothermal limit reproduces Pao–Sah; I_D–V_D droop with negative output conductance |
| `mosfet/` (electro-thermo-mechanical) | + small-strain **thermoelasticity** (trilinear hex, 2×2×2 Gauss, σ = C:(ε − αΔT·I)) and **piezoresistive** mobility feedback from the element stress (Smith 1954 n-Si coefficients) | **7-dof monolithic** (ψ, φₙ, φₚ, ΔT, uₓ, u_y, u_z) — displacements ride Abaqus dof slots 5–7 (UR2/UR3/WARP) | uniform-ΔT free expansion vs closed form (0.08 %); uniaxial −100 MPa → ΔI_D/I_D = −10.1 % vs Smith's −π₁₁σ = −10.2 %; energy balance 0.01 %; thermal-stress droop on top of self-heating |

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

## Self-heating MOSFET: monolithic electro-thermal UEL

Lattice temperature rise ΔT is added as a **4th nodal degree of freedom** and solved in
the *same* Newton matrix as the device equations — (ψ, φₙ, φₚ, ΔT) fully coupled. This
is the differentiator: published electro-thermal device studies overwhelmingly couple a
TCAD solver to a thermal FE solver in a staggered loop; here the whole thing is one
unsymmetric Jacobian inside Abaqus. The heat equation is discretized with the same
box method (edge conductances κA/h), the Joule source J·E is lumped edge-wise from the
SG fluxes, and temperature feeds back through V_T(T) = k_BT/q in the SG exponent and
µ(T) ∝ (T/300)⁻³ᐟ².

![Monolithic electro-thermal MOSFET simulation in Abaqus: drain current droop from self-heating and lattice temperature hotspot at the drain end of the channel](docs/fig_selfheating.png)

V_G = 3 V, substrate bottom held at 300 K, all other surfaces adiabatic. The device is a
0.5 µm-wide toy, so the physical dissipation (~µW) gives mK heating; a heat multiplier
HSCALE = 200 emulates a multi-finger power layout (equivalent to scaling κ down) to make
the classic signatures visible:

| V_D [V] | I_D isothermal [µA] | I_D self-heating [µA] | ΔT_max [K] | energy balance err |
|---|---|---|---|---|
| 0.5 | 10.5 | 10.0 | 5.8 | 0.00 % |
| 1.5 | 17.9 | 15.0 | 37.7 | 0.01 % |
| 3.0 | 18.4 | 14.0 | 118.9 | 0.00 % |

Verification (no closed form exists for the coupled problem):
- **Energy balance / Tellegen**: total heat flowing out through the heatsink (reaction
  "moment" of the ΔT dof at the Dirichlet nodes) equals HSCALE·I_D·V_D to 0.01 % — the
  discrete Joule sum over edges telescopes exactly to the terminal power.
- **Weak-coupling limit**: with HSCALE = 0 the element reduces to the isothermal UEL
  equation-for-equation; I_D matches Pao–Sah within 3.7 % (V_D ≤ 1 V, where the
  gradual-channel assumption holds) and ΔT ≡ 0 to machine precision.
- **Self-heating signature**: 23 % I_D droop at V_D = 3 V with *negative* output
  conductance beyond V_D ≈ 1.5 V, and the temperature hotspot sits at the drain end of
  the channel — the textbook qualitative picture.

Reproduce with `python mosfet/run_selfheating.py` (two Abaqus jobs, heating on/off).

## Electro-thermo-mechanical: the full three-field loop

`mosfet/uel_mos_etm.f` extends the element to **seven nodal dofs**
(ψ, φₙ, φₚ, ΔT, uₓ, u_y, u_z), closing the loop
*current → Joule heat → temperature → thermal strain → stress → piezoresistance →
mobility → current* in a single monolithic Newton iteration. Displacements occupy
Abaqus's remaining structural dof slots (5, 6, 7 = UR2, UR3, warping — a `*STATIC`
step accepts all of them for user elements). Mechanics is standard small-strain
thermoelasticity on the same hexahedra (2×2×2 Gauss, consistent K_uu and
thermal-expansion K_uT blocks); the electron mobility carries a piezoresistance factor
1 − (π₁₁σₓₓ + π₁₂(σ_yy + σ_zz)) with Smith (1954) n-Si coefficients evaluated from the
element-center stress.

![Electro-thermo-mechanical MOSFET simulation in Abaqus: 7-dof monolithic user element with self-heating and piezoresistive feedback](docs/fig_etm.png)

Verification (`python mosfet/run_etm.py`, three Abaqus jobs):

| check | result | reference |
|---|---|---|
| uniform ΔT = 100 K free expansion (bottom rollers) | u_z(top) error 0.08 % | closed form u_z = −αΔT·H |
| uniaxial σₓₓ = −100 MPa via face displacement | ΔI_D/I_D = −10.13 % | Smith (1954): −π₁₁σₓₓ = −10.22 % |
| full loop, bottom clamped, HSCALE = 200 | energy balance ≤ 0.01 %, ΔT_max 118 K | Tellegen, as above |
| thermal-stress feedback | I_D(3 V) = 13.99 µA < 14.04 µA (thermal-only) | compressive hotspot lowers µₙ — extra droop on top of self-heating |

The piezoresistive coupling enters the residual only (its Jacobian columns are
omitted — a deliberate quasi-Newton shortcut for a ~3 % effect; the V_D continuation
ladder carries convergence, and all four checks still pass).

## How to run

Requires Abaqus (tested with 2024) with a linked Fortran compiler, plus Python 3 with
numpy (and matplotlib for the figure).

```
cd diode
python run_diode.py     # ~1 min: writes inp, runs abaqus job=... user=uel_dd.f, checks

cd mosfet
python run_mosfet.py         # ~3 min: full bias sweep in one job (6 steps), figure + checks
python run_selfheating.py    # ~4 min: electro-thermal, 2 jobs (HSCALE 0/200), figure + checks
python run_etm.py            # ~6 min: electro-thermo-mechanical, 3 jobs, figure + checks
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
- **NaN passes every convergence check.** Coupling the weakly-conducting temperature
  equation to exponentials is a trap: during Newton excursions the SG fluxes ride the
  exp-clamp (~e⁸⁰), an unguarded Joule source kicks ΔT to 10¹⁸ K, the residual overflows
  — and Abaqus reports the step *converged*, because every comparison against NaN is
  false. It then sails through all steps in one iteration each and writes a `.dat` full
  of NaN with zero error messages. Guards: clamp the per-edge Joule power and drop the
  thermal-coupling Jacobian entries while fluxes are unphysical (both inactive at the
  converged solution, so consistency is preserved where it matters).

## FAQ

**Can Abaqus simulate semiconductor devices?**
Yes — a `*USER ELEMENT` is just "residual + Jacobian per element", and nothing forces
the degrees of freedom to be displacements. Here they are electrostatic potential and
quasi-Fermi potentials, the residuals are Poisson + electron/hole continuity, and
Abaqus's Newton loop, unsymmetric solver, and load stepping do the rest.

**How is the drift–diffusion equation discretized?**
With the Scharfetter–Gummel (1969) exponential-fitting flux on every mesh edge — the
same box-integration scheme used by commercial TCAD tools. It reduces to central
differencing for pure diffusion, becomes upwind automatically under strong drift, and
gives exactly zero flux at equilibrium.

**Why quasi-Fermi variables instead of carrier densities?**
With (ψ, n, p) at realistic doping (n/nᵢ ~ 10⁹) the coupled Newton diverges from any
cold start. With (ψ, φₙ, φₚ) every unknown is a bounded voltage, contacts become plain
voltage boundary conditions, and the same problem converges with ordinary ramping.

**Is this a full TCAD replacement?**
No — no recombination models, mobility models, or impact ionization; constant mobility
and Boltzmann statistics only. It is a transparent, verified reference implementation
of the numerical core: discretization, linearization, and solver.

## 한국어 소개

반도체 소자 시뮬레이션(TCAD)의 수치 코어 — 포아송 방정식과 전자/정공 연속
방정식의 드리프트-확산(drift–diffusion) 모델 — 를 구조해석 코드 Abaqus의 사용자
요소(UEL)로 구현한 프로젝트입니다. 상용 TCAD(Sentaurus, Silvaco Atlas)와 같은
Scharfetter–Gummel box method 이산화, 준페르미 퍼텐셜 변수, 완전 연성 Newton
풀이를 사용하며, 1D pn 접합 다이오드(내장전위·질량작용법칙·전류보존 검증)와
3D 긴채널 NMOS MOSFET(문턱전압, 반전 채널, 전달특성 I_D–V_G)을 Pao–Sah(1966)
정확해 및 Brews(1978) charge-sheet 모델과 정량 대조해 1–6% 수준으로
재현합니다. 유한요소 정식화·일관 선형화·사용자 서브루틴(UMAT/UEL) 개발
경험을 소자 시뮬레이션으로 확장한 작업입니다.

## References

- D. L. Scharfetter, H. K. Gummel, *Large-signal analysis of a silicon Read diode
  oscillator*, IEEE Trans. Electron Devices 16 (1969) 64.
- H. C. Pao, C. T. Sah, *Effects of diffusion current on characteristics of
  metal-oxide (insulator)-semiconductor transistors*, Solid-State Electron. 9 (1966) 927.
- J. R. Brews, *A charge-sheet model of the MOSFET*, Solid-State Electron. 21 (1978) 345.
- G. K. Wachutka, *Rigorous thermodynamic treatment of heat generation and conduction in
  semiconductor device modeling*, IEEE Trans. CAD 9 (1990) 1141.
- C. S. Smith, *Piezoresistance effect in germanium and silicon*, Phys. Rev. 94
  (1954) 42.

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
