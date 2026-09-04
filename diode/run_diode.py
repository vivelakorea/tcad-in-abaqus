# -*- coding: utf-8 -*-
"""run_diode: pn 접합 드리프트-확산을 Abaqus UEL로 푼다.

Scharfetter-Gummel (1969) 이산화를 2절점 UEL(uel_dd.f, 절점 자유도
psi,n,p)로 구현하고 reference_dd1d.py(파이썬 참조해)와 같은 문제를 돈다:
  Step 1: 접촉 BC를 0에서 평형값까지 램프 (V=0)
  Step 2: p쪽 접촉에 순바이어스 V=10 V_T 램프
검증: (1) Vbi vs ln(Na*Nd), 질량법칙  (2) 전류 보존
      (3) 파이썬 coupled Newton 해와 절점별 비교 (같은 이산화 -> 같은 해)
사용: python run_diode.py   (Abaqus 필요)
"""
import io
import os
import re
import subprocess

import numpy as np

import reference_dd1d as ref

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
NA, ND, L, NN = 100.0, 100.0, 4.0, 201
NINC1, NINC2, VFWD = 20, 10, 10.0


def write_inp(job):
    xs = np.linspace(-L/2, L/2, NN)
    (psiL, nL, pL), (psiR, nR, pR) = ref.contacts(0.0)
    Lns = ['*HEADING', 'ex16 drift-diffusion UEL (Scharfetter-Gummel 1969)',
           '*USER ELEMENT, NODES=2, TYPE=U1, PROPERTIES=2, COORDINATES=1,'
           ' VARIABLES=1, UNSYMM',
           '1,2,3', '*NODE']
    for a, x in enumerate(xs, 1):
        Lns.append(f'{a}, {x:.10f}')
    Lns.append('*ELEMENT, TYPE=U1, ELSET=EALL')
    for e in range(1, NN):
        Lns.append(f'{e}, {e}, {e+1}')
    Lns += ['*UEL PROPERTY, ELSET=EALL', f'{NA}, {ND}',
            '*NSET, NSET=NALL, GENERATE', f'1, {NN}, 1',
            # Step 1: 평형 BC 램프
            '*STEP, INC=200, UNSYMM=YES', '*STATIC, DIRECT',
            f'{1.0/NINC1}, 1.0',
            '*CONTROLS, PARAMETERS=FIELD', '1e-6, 1e-6,,,,,,',
            '*BOUNDARY',
            f'1, 1, 1, {psiL}', f'1, 2, 2, {nL}', f'1, 3, 3, {pL}',
            f'{NN}, 1, 1, {psiR}', f'{NN}, 2, 2, {nR}', f'{NN}, 3, 3, {pR}',
            '*NODE PRINT, NSET=NALL, FREQUENCY=1', 'U', '*END STEP',
            # Step 2: 순바이어스 램프
            '*STEP, INC=200, UNSYMM=YES', '*STATIC, DIRECT',
            f'{1.0/NINC2}, 1.0',
            '*CONTROLS, PARAMETERS=FIELD', '1e-6, 1e-6,,,,,,',
            '*BOUNDARY', f'1, 1, 1, {psiL + VFWD}',
            '*NODE PRINT, NSET=NALL, FREQUENCY=1', 'U', '*END STEP']
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(Lns) + '\n')


def parse_blocks(job):
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    blocks = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        rows = re.findall(r'^\s+(\d+)\s+' + (num + r'\s+')*2 + num, blk, re.M)
        d = np.zeros((NN, 3))
        for r_ in rows:
            d[int(r_[0]) - 1] = [float(v) for v in r_[1:]]
        blocks.append(d)
    return blocks


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex16dd'
    write_inp(job)
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = f'abaqus job={job} user={os.path.join(HERE, "uel_dd.f")} interactive'
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True, text=True)
    assert 'COMPLETED' in r.stdout, r.stdout[-500:]
    blocks = parse_blocks(job)
    assert len(blocks) == NINC1 + NINC2, f'{len(blocks)} blocks'
    eq, fw = blocks[NINC1 - 1], blocks[-1]

    # (1) 평형
    Vbi = eq[-1, 0] - eq[0, 0]
    mass = np.abs(eq[:, 1] * eq[:, 2] - 1).max()
    print(f'equilibrium: Vbi = {Vbi:.6f}  vs  ln(Na*Nd) = {np.log(NA*ND):.6f}')
    print(f'mass law  max|np-1| = {mass:.2e}')

    # (2) 순바이어스 전류 보존
    psi, n, p = fw[:, 0], fw[:, 1], fw[:, 2]
    t = np.diff(psi)
    h = L / (NN - 1)
    Jn = (ref.bern(t)*n[1:] - ref.bern(-t)*n[:-1]) / h
    Jp = (ref.bern(-t)*p[1:] - ref.bern(t)*p[:-1]) / h
    Jtot = Jn - Jp
    dev = (Jtot.max() - Jtot.min()) / np.abs(Jtot).max()
    print(f'forward V=10VT: current constancy dev = {dev:.2e}')

    # (3) 파이썬 참조해와 절점별 비교
    d0, _ = ref.gummel(ref.init_guess(), 0.0)
    dn, _ = ref.newton_ramped(d0.copy(), VFWD)
    dpsi = np.abs(psi - dn[:NN]).max()
    dn_n = np.abs(n - dn[NN:2*NN]).max() / dn[NN:2*NN].max()
    print(f'|psi_UEL - psi_py|max = {dpsi:.2e},  rel|n_UEL - n_py|max = {dn_n:.2e}')
    assert abs(Vbi - np.log(NA*ND)) < 1e-3 and mass < 1e-6
    assert dev < 1e-4 and dpsi < 1e-6 and dn_n < 1e-6
    print('check passed: 같은 SG 이산화 -> Abaqus UEL과 파이썬이 같은 해.')


if __name__ == '__main__':
    main()
