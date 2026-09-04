# -*- coding: utf-8 -*-
"""run_emig: 1D electromigration UEL (uel_em.f, Korhonen 모델) 검증.

Al 배선 (L=200um, A=1x0.5um^2, j=1e5 A/cm^2, ~500K 가속시험 물성),
양끝 원자 플럭스 차단(blocked), (V, sigma) monolithic + backward Euler.

검증 = 논문 결과 재현:
1. Korhonen et al., J. Appl. Phys. 73 (1993) 3790 이 발표한 유한 배선
   응력 전개 해(그 논문 그림들의 원천, reference_korhonen.py 독립 구현)를
   kappa*t/L^2 = 0.01 ~ 0.5 및 정상상태에서 절점별 프로파일 대조.
2. 정상상태 Blech-Herring back-stress: delta_sigma = e Z* rho j L / Omega,
   선형 프로파일 (캐소드 인장 / 애노드 압축).
3. Blech, J. Appl. Phys. 47 (1976) 1203 실측 임계곱 (jL)_c ~ 1260 A/cm 를
   본 파라미터의 back-stress 스케일로 환산 대조.

사용: python run_emig.py   (Abaqus job 1개, ~1 min)
"""
import io
import os
import re
import subprocess

import numpy as np

import reference_korhonen as ref

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
TAU = LX*LX/XKAP                             # 확산 시간 L^2/kappa [s]


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

    def step(dt, period, freq):
        return (['*STEP, INC=400, UNSYMM=YES', '*STATIC',
                 f'{dt:.6e}, {period:.6e}, {dt*0.99:.6e}, {dt:.6e}'] + bc
                + [f'*NODE PRINT, NSET=NALL, FREQUENCY={freq}', 'U',
                   '*END STEP'])

    # 논문 그림의 시각축: kappa*t/L^2 = 0.01 / 0.055 / 0.1 / 0.3 / 0.5 / 정상
    L += step(5e-4*TAU, 0.01*TAU, 999)               # -> 0.01 tau
    L += step(5e-3*TAU, 0.09*TAU, 9)                 # -> 0.055, 0.1 tau
    L += step(2.5e-2*TAU, 0.4*TAU, 8)                # -> 0.3, 0.5 tau
    L += step(0.45*TAU, 9.0*TAU, 999)                # -> ~9.5 tau (정상)
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
    tlist = [0.01, 0.055, 0.1, 0.3, 0.5]             # [tau]
    assert len(ev) == len(tlist) + 1, f'{len(ev)} print blocks'

    G = CEZ*VAPP/LX                                  # MPa/cm
    sig_ss = G*LX/2                                  # 캐소드 정상 back-stress
    # 1. Korhonen (1993) 유한 배선 해와 프로파일 대조
    print(' kappa*t/L^2   max|UEL - Korhonen1993| / sig_ss')
    worst = 0.0
    for d, tn in zip(ev[:-1], tlist):
        th = ref.sigma(xs, tn*TAU, LX, XKAP, G)
        dev = np.abs(d[:, 1] - th).max()/sig_ss
        worst = max(worst, dev)
        print(f'   {tn:5.3f}        {dev*100:.2f}%')
    print(f'worst profile dev vs Korhonen(1993) = {worst*100:.2f}%')
    # 2. 정상상태 Blech-Herring
    prof = ev[-1][:, 1]
    dsig = prof[-1] - prof[0]
    dsig_th = CEZ*VAPP                               # = e Z* rho j L / Omega
    lin = np.polyfit(xs, prof, 1)
    resid = np.abs(prof - np.polyval(lin, xs)).max()/abs(dsig)
    print(f'[정상상태] delta_sigma = {dsig:.1f} MPa vs eZ*rho*j*L/Omega = '
          f'{dsig_th:.1f} MPa; 선형성 잔차 {resid*100:.2f}%; '
          f'캐소드 {prof[-1]:+.1f} / 애노드 {prof[0]:+.1f} MPa')
    # 3. Blech (1976)
    jl, jlc = J*LX, 1260.0
    sig_c = CEZ*RHO*jlc/2
    print(f'[Blech] jL = {jl:.0f} A/cm (> (jL)_c ~ {jlc:.0f}: EM 진행); '
          f'(jL)_c 환산 임계 back-stress = {sig_c:.0f} MPa')

    # 그림: 논문 해(점선) 위에 UEL(실선) 겹치기
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    xu = xs*1e4
    cmap = plt.get_cmap('viridis')
    for k, (d, tn) in enumerate(zip(ev[:-1], tlist)):
        c = cmap(k/(len(tlist)))
        ax.plot(xu, d[:, 1], '-', color=c, lw=1.8,
                label=f'$\\kappa t/L^2$ = {tn:g}')
        ax.plot(xu, ref.sigma(xs, tn*TAU, LX, XKAP, G), 'k--', lw=0.9)
    ax.plot(xu, prof, '-', color='crimson', lw=1.8, label='steady state')
    ax.plot(xu, G*(xs - LX/2), 'k--', lw=0.9,
            label='Korhonen (1993) solution')
    ax.set_xlabel('x [$\\mu$m] (anode $\\to$ cathode)')
    ax.set_ylabel('$\\sigma$ [MPa]')
    ax.set_title('EM stress evolution: UEL vs Korhonen (1993)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_emig')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    print('figure ->', os.path.abspath(out) + '.png')

    assert worst < 0.02                              # Korhonen 1993 전 구간
    assert abs(dsig/dsig_th - 1) < 0.01              # Blech-Herring
    assert resid < 0.01
    print('check passed: EM UEL — Korhonen(1993) 프로파일 전개 / '
          'Blech-Herring / Blech(1976).')


if __name__ == '__main__':
    main()
