# TCAD in Abaqus — Semiconductor Device Simulation with User Elements

**Abaqus 사용자 요소(UEL)로 구현한 반도체 소자 시뮬레이션(TCAD).**

Semiconductor device physics implemented as Abaqus user elements — the structural FE
solver becomes a device simulator. Drift–diffusion (Scharfetter–Gummel box method,
quasi-Fermi variables) for a 1D pn diode and a 3D MOSFET, extended monolithically to
electro-thermal and electro-thermo-mechanical coupling, metal-line electromigration
(Korhonen model), and a density-gradient quantum correction — each physics added as an
extra nodal degree of freedom in one Newton matrix.

Every module is verified against published results (Pao–Sah 1966, Brews 1978, Smith
1954, Wachutka 1990, Korhonen 1993, Blech 1976, Lee 2009) or an independent Python
reference implementation. Each `run_*.py` driver generates the mesh and input deck,
launches Abaqus, parses the output, and asserts its checks — `check passed` means it
reproduced.

Requires Abaqus (tested with 2024) with a linked Fortran compiler and Python 3 with
numpy/scipy/matplotlib.

## Author

심규장 (Gyu-Jang Sim) — 서울대학교 재료역학연구실 ([MAMEL](https://mamel.snu.ac.kr), Seoul National University)
[github.com/vivelakorea](https://github.com/vivelakorea) · gyujang95@gmail.com

## License

MIT
