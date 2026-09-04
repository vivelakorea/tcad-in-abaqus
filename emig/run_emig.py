# -*- coding: utf-8 -*-
"""run_emig: 1D electromigration UEL (uel_em.f, Korhonen 모델) 검증.

Al 배선 (L=200um, A=1x0.5um^2, j=1e5 A/cm^2, ~500K 가속시험 물성),
양끝 원자 플럭스 차단(blocked), (V, sigma) monolithic + backward Euler.

문헌 수치 검증:
1. 과도: sigma(끝단, t) = G*sqrt(4 kappa t / pi), G = e Z* rho j / Omega
   -- Korhonen et al., J. Appl. Phys. 73 (1993) 3790 닫힌형. sqrt(t) 법칙.
2. 정상상태: 선형 back-stress 프로파일, delta_sigma = e Z* rho j L / Omega
   (Blech-Herring back stress). 캐소드 인장 / 애노드 압축.
3. Blech (1976) 임계곱: (jL)_c ~ 1260 A/cm 실측(Al). 본 파라미터로 환산한
   임계 back-stress를 출력해 대조 (jL 이하면 EM 정지 = Blech effect).

사용: python run_emig.py   (Abaqus job 1개, ~1 min)
"""
import io
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
# uel_em.f 과 동일 물성 (Al, ~500K)
RHO, XKAP = 4.9e-6, 1.8e-9                   # Ohm*cm, cm^2/s
ZST, OMEG, QE = 4.0, 1.66e-23, 1.602e-19
CEZ = QE*ZST/OMEG                            # C/cm^3; CEZ*E = MPa/cm
LX = 200e-4                                  # cm (200 um)
AR = 1.0e-4*0.5e-4                           # cm^2
J = 1.0e5                                    # A/cm^2
VAPP = J*RHO*LX                              # V
NEL = 80
T1, T2 = 2500.0, 5000.0                      # 과도 체크 시각 [s]


def write_inp(job):
    xs = np.linspace(0, LX, NEL + 1)
    L = ['*HEADING', 'ex21 Al line electromigration (Korhonen model UEL)',
         '*USER ELEMENT, NODES=2, TYPE=U1, PROPERTIES=1, COORDINATES=1,'
         ' VARIABLES=1, UNSYMM', '1,2', '*NODE']
    L += [f'{i+1}, {xs[i]:.10e}' for i in range(NEL + 1)]
    L.append('*ELEMENT, TYPE=U1, ELSET=ELINE')
    L += [f'{i+1}, {i+1}, {i+2}' for i in range(NEL)]
    L += ['*UEL PROPERTY, ELSET=ELINE', f'{AR:.6e}',
          '*NSET, NSET=NA', '1', f'*NSET, NSET=NC', f'{NEL+1}',
          '*NSET, NSET=NALL, GENERATE', f'1, {NEL+1}, 1',
          '*AMPLITUDE, NAME=INST, TIME=TOTAL TIME',
          '0., 1., 1e9, 1.']
    bc = ['*BOUNDARY, AMPLITUDE=INST',
          f'NA, 1, 1, {VAPP:.8e}', 'NC, 1, 1, 0.']
    ninc1 = int(T2/50)
    L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
          f'50., {T2:.1f}, 1., 50.'] + bc
    L += [f'*NODE PRINT, NSET=NALL, FREQUENCY={ninc1//2}', 'U', '*END STEP']
    L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
          '1e5, 2e6, 1e3, 1e5'] + bc + ['*END STEP']
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs


def parse(job, nn):
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    out = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        for pp in blk.split('NODE FOOT-')[1:]:
            if 'U1' in pp.split('\n', 1)[0]:
                rows = re.findall(r'^\s+(\d+)\s+' + num + r'\s+' + num,
                                  pp, re.M)
                d = np.zeros((nn, 2))
                for r_ in rows:
                    d[int(r_[0]) - 1] = [float(r_[1]), float(r_[2])]
                out.append(d)
    return out


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex21em'
    xs = write_inp(job)
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = (f'abaqus job={job} user={os.path.join(HERE, "uel_em.f")}'
           f' interactive cpus=1')
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                       text=True, errors='replace')
    assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    ev = parse(job, NEL + 1)
    assert len(ev) == 3, f'{len(ev)} print blocks'

    G = CEZ*VAPP/LX                                  # MPa/cm
    # 1. 과도 sqrt(t): 캐소드 끝 인장
    s1, s2 = ev[0][-1, 1], ev[1][-1, 1]
    th1 = G*np.sqrt(4*XKAP*T1/np.pi)
    th2 = G*np.sqrt(4*XKAP*T2/np.pi)
    print(f'[sqrt(t)] sigma_c({T1:.0f}s) = {s1:.1f} MPa (닫힌형 {th1:.1f}), '
          f'sigma_c({T2:.0f}s) = {s2:.1f} MPa (닫힌형 {th2:.1f})')
    print(f'          비율 {s2/s1:.4f} vs sqrt(2) = {np.sqrt(2):.4f}')
    # 2. 정상상태 선형 back-stress
    prof = ev[2][:, 1]
    dsig = prof[-1] - prof[0]
    dsig_th = CEZ*VAPP                               # = e Z* rho j L / Omega
    lin = np.polyfit(xs, prof, 1)
    resid = np.abs(prof - np.polyval(lin, xs)).max()/abs(dsig)
    print(f'[정상상태] delta_sigma = {dsig:.1f} MPa, eZ*rho*j*L/Omega = '
          f'{dsig_th:.1f} MPa; 선형성 잔차 {resid*100:.2f}%')
    print(f'          캐소드 {prof[-1]:+.1f} MPa (인장), '
          f'애노드 {prof[0]:+.1f} MPa (압축)')
    # 3. Blech (1976)
    jl = J*LX
    jlc = 1260.0                                     # A/cm, Blech 실측(Al)
    sig_c = CEZ*RHO*jlc/2*1                          # (jL)_c 의 최대 back-stress
    print(f'[Blech] 본 시험 jL = {jl:.0f} A/cm (> (jL)_c ~ {jlc:.0f}: EM 진행).'
          f' (jL)_c 환산 임계 back-stress = {sig_c:.0f} MPa'
          f' (void 핵생성 문턱 O(100 MPa) 문헌과 정합)')

    # 그림
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    xu = xs*1e4
    for d, t, c in ((ev[0], T1, 'tab:blue'), (ev[1], T2, 'tab:orange')):
        ax.plot(xu, d[:, 1], '-', color=c, label=f't = {t:.0f} s')
    ax.plot(xu, prof, '-', color='crimson', label='steady state')
    ax.plot(xu, prof[0] + (dsig_th/LX)*xs, 'k--', lw=1,
            label='$eZ^*\\rho jL/\\Omega$ (Blech-Herring)')
    ax.set_xlabel('x [$\\mu$m] (anode $\\to$ cathode)')
    ax.set_ylabel('$\\sigma$ [MPa]')
    ax.set_title('Electromigration back-stress (Korhonen model UEL)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_emig')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    print('figure ->', os.path.abspath(out) + '.png')

    assert abs(s2/th2 - 1) < 0.10                    # Korhonen sqrt(t)
    assert abs(s2/s1 - np.sqrt(2)) < 0.07
    assert abs(dsig/dsig_th - 1) < 0.01              # Blech-Herring 선형
    assert resid < 0.01
    print('check passed: EM UEL — Korhonen sqrt(t) / 정상 back-stress / Blech.')


if __name__ == '__main__':
    main()
