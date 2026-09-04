# -*- coding: utf-8 -*-
"""run_etm: monolithic 전기-열-기계 NMOS UEL (uel_mos_etm.f).

7-dof (psi, phi_n, phi_p, dT, ux, uy, uz) 를 한 Newton 행렬로 푼다.
- 열탄성: sigma = C:(eps - alpha*dT*I), 2x2x2 Gauss, K_uu/K_uT 일관 선형화.
- 압저항: 요소 중심 응력 -> 전자 이동도 FPZ = 1 - (pi11 sxx + pi12(syy+szz)),
  Smith (1954) n-Si <100> 계수. (잔차만 반영, ponytail 주석 참조.)

검증:
1. 균일 dT=100K 자유팽창(바닥 롤러): u_z(top) = -alpha*dT*H 닫힌형.
2. 단축 압축 sxx=-100MPa (x-면 변위 규정): 선형영역 I_D 변화
   dI/I = FPZ-1 = -pi11*sxx = -10.2% (Smith 1954 직접 대조).
3. full loop (HSCALE=200, 바닥 고정): 에너지수지 + I_D-V_D 음의 기울기
   (열+압저항 droop; ET-only 결과와 비교 출력).

사용: python run_etm.py   (Abaqus job 3개, ~수 분)
"""
import io
import os
import re
import subprocess

import numpy as np

from run_mosfet import (grids, contacts, Q, VT, NI, NA_SUB, ND_PLUS,
                        UM, XG0, XG1, XJ, WM)

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
VG = 3.0
HSCALE = 200.0
ALSI, EYSI = 2.6e-6, 1.30e11
P11, P12 = -102.2e-11, 53.4e-11              # Smith (1954) n-Si [1/Pa]
DT0 = 100.0                                  # 팽창 테스트 dT [K]
EPS0 = -1.0e8 / EYSI                         # 압저항 테스트: sxx=-100MPa


def write_inp(job, hscale, steps, mech_bc):
    """steps: (vg, vd, ninc, extra_bc_lines) 목록. mech_bc: 공통 역학 BC."""
    xs, zs, zox = grids()
    gate_ix = np.where((xs >= XG0 - 1e-9) & (xs <= XG1 + 1e-9))[0]
    nid = {}
    lines = []
    a = 0
    for iz, z in enumerate(np.r_[zox, zs]):
        for iy in (0, 1):
            for ix, x in enumerate(xs):
                zlev = iz - len(zox)
                if zlev < 0 and ix not in gate_ix:
                    continue
                a += 1
                nid[(ix, zlev, iy)] = a
                lines.append(f'{a}, {x*UM:.10e}, {iy*WM*UM:.10e}, '
                             f'{np.r_[zox, zs][iz]*UM:.10e}')
    els = {'EOX': [], 'ENP': [], 'EPSUB': []}
    enum = 0
    conn_lines = {k: [] for k in els}
    for zlev in range(-len(zox), len(zs) - 1):
        for ix in range(len(xs) - 1):
            def n8():
                return [nid.get((ix, zlev, 0)), nid.get((ix+1, zlev, 0)),
                        nid.get((ix+1, zlev, 1)), nid.get((ix, zlev, 1)),
                        nid.get((ix, zlev+1, 0)), nid.get((ix+1, zlev+1, 0)),
                        nid.get((ix+1, zlev+1, 1)), nid.get((ix, zlev+1, 1))]
            nn = n8()
            if None in nn:
                continue
            xc = 0.5*(xs[ix] + xs[ix+1])
            if zlev < 0:
                key = 'EOX'
            else:
                zc = 0.5*(zs[zlev] + zs[zlev+1])
                np_reg = (xc < XG0 or xc > XG1) and zc < XJ
                key = 'ENP' if np_reg else 'EPSUB'
            enum += 1
            els[key].append(enum)
            conn_lines[key].append(f'{enum}, ' + ', '.join(map(str, nn)))
    psn, psp = contacts()
    L = ['*HEADING', 'ex19 3D NMOS electro-thermo-mechanical UEL',
         '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=3, COORDINATES=3,'
         ' VARIABLES=1, UNSYMM', '1,2,3,4,5,6,7', '*NODE'] + lines
    for key in ('EOX', 'ENP', 'EPSUB'):
        if conn_lines[key]:
            L.append(f'*ELEMENT, TYPE=U1, ELSET={key}')
            L += conn_lines[key]
    L += ['*UEL PROPERTY, ELSET=EOX', f'0, 0., {hscale:.6e}',
          f'*UEL PROPERTY, ELSET=ENP', f'1, {ND_PLUS/NI:.6e}, {hscale:.6e}',
          f'*UEL PROPERTY, ELSET=EPSUB',
          f'1, {-NA_SUB/NI:.6e}, {hscale:.6e}']
    sets = {'SRC': [], 'DRN': [], 'BLK': [], 'GATE': [], 'OXI': [],
            'XL': [], 'XR': [], 'YL': [], 'YR': []}
    for (ix, zlev, iy), aa in nid.items():
        x = xs[ix]
        if zlev == 0 and x <= 0.9 + 1e-9:
            sets['SRC'].append(aa)
        if zlev == 0 and x >= 5.1 - 1e-9:
            sets['DRN'].append(aa)
        if zlev == len(zs) - 1:
            sets['BLK'].append(aa)
        if zlev == -len(zox):
            sets['GATE'].append(aa)
        if zlev < 0:
            sets['OXI'].append(aa)
        if ix == 0:
            sets['XL'].append(aa)
        if ix == len(xs) - 1:
            sets['XR'].append(aa)
        sets['YL' if iy == 0 else 'YR'].append(aa)
    for nm, ids in sets.items():
        L.append(f'*NSET, NSET={nm}')
        ids = sorted(ids)
        for i in range(0, len(ids), 12):
            L.append(', '.join(map(str, ids[i:i+12])))
    L.append('*NSET, NSET=NALL, GENERATE')
    L.append(f'1, {a}, 1')
    nzb = len(zs) - 1
    L += ['*NSET, NSET=PINA', str(nid[(0, nzb, 0)]),
          '*NSET, NSET=PINB', str(nid[(0, nzb, 1)])]
    psi_bulk = psp

    first = True
    for vg, vd, ninc, extra in steps:
        s = ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
             f'{1.0/ninc}, 1.0, 1e-9, {1.0/ninc}',
             '*CONTROLS, PARAMETERS=FIELD', '1e-6,,,,,,,',
             '*BOUNDARY',
             f'GATE, 1, 1, {vg/VT + psi_bulk:.8e}',
             'GATE, 2, 3, 0.', 'OXI, 2, 3, 0.',
             f'SRC, 1, 1, {psn:.8e}', 'SRC, 2, 3, 0.',
             f'DRN, 1, 1, {psn + vd/VT:.8e}',
             f'DRN, 2, 2, {vd/VT:.8e}', f'DRN, 3, 3, {vd/VT:.8e}',
             f'BLK, 1, 1, {psp:.8e}', 'BLK, 2, 3, 0.'] + mech_bc + extra
        if first:
            s += ['*NODE PRINT, NSET=NALL, FREQUENCY=999', 'U',
                  '*NODE PRINT, NSET=DRN, FREQUENCY=999, TOTALS=YES', 'RF',
                  '*NODE PRINT, NSET=BLK, FREQUENCY=999, TOTALS=YES', 'RF']
            first = False
        s.append('*END STEP')
        L += s
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs, zs, nid


def parse(job, nn):
    """스텝별 [U[nn,7], [RF_DRN 총합(7), RF_BLK 총합(7)]]."""
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    tot = r'(-?\d+\.\d*(?:E[+-]\d+)?|NaN)'
    events = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        for pp in blk.split('NODE FOOT-')[1:]:
            hdr = pp.split('\n', 1)[0]
            if 'U1' in hdr:
                rows = re.findall(r'^\s+(\d+)\s+' + (num + r'\s+')*6 + num,
                                  pp, re.M)
                d = np.zeros((nn, 7))
                for r_ in rows:
                    d[int(r_[0]) - 1] = [float(v) for v in r_[1:]]
                events.append([d, []])
            elif 'RF1' in hdr and events:
                m = re.search(r'^ TOTAL\s+' + (tot + r'\s+')*6 + tot,
                              pp, re.M)
                if m:
                    events[-1][1].append([float(m.group(k))
                                          for k in range(1, 8)])
    return events


def run_job(job, hscale, steps, mech_bc):
    xs, zs, nid = write_inp(job, hscale, steps, mech_bc)
    nn = max(nid.values())
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = (f'abaqus job={job} user={os.path.join(HERE, "uel_mos_etm.f")}'
           f' interactive cpus=1')
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                       text=True, errors='replace')
    assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    ev = parse(job, nn)
    assert len(ev) == len(steps), f'{len(ev)} print blocks != {len(steps)}'
    return xs, zs, nid, ev


def main():
    os.makedirs(RUN, exist_ok=True)
    xs, zs, _ = grids()
    H = zs[-1] * UM                                  # 기판 두께 [cm]
    LX = xs[-1] * UM

    # ---- 1. 자유 열팽창: 바닥 롤러, 균일 dT -> u_z(top) = -alpha*dT*H ----
    rollers = ['BLK, 7, 7, 0.', 'PINA, 5, 6, 0.', 'PINB, 5, 5, 0.',
               f'NALL, 4, 4, {DT0:.4e}']
    steps = [(0., 0., 20, [])]
    xs, zs, nid, ev = run_job('ex19exp', 0.0, steps, rollers)
    top = nid[(len(xs)//2, 0, 0)] - 1
    uz = ev[0][0][top, 6]
    uz_exact = -ALSI*DT0*H
    err_exp = abs(uz/uz_exact - 1)
    print(f'[팽창] u_z(top) = {uz:.4e} cm, 닫힌형 {uz_exact:.4e} cm, '
          f'오차 {err_exp*100:.2f}%')

    # ---- 2. 압저항: x-면 변위로 sxx = E*eps0, dI/I = -pi11*sxx ----
    pins = ['PINA, 6, 7, 0.', 'PINB, 7, 7, 0.', 'BLK, 4, 4, 0.']
    steps = [(0., 0., 20, ['XL, 5, 5, 0.', 'XR, 5, 5, 0.']),
             (VG, 0.05, 10, ['XL, 5, 5, 0.', 'XR, 5, 5, 0.']),
             (VG, 0.05, 4, ['XL, 5, 5, 0.',
                            f'XR, 5, 5, {EPS0*LX:.6e}'])]
    _, _, _, ev = run_job('ex19pz', 0.0, steps, pins)
    i0 = Q*NI*ev[1][1][0][1]
    i1 = Q*NI*ev[2][1][0][1]
    di = i1/i0 - 1
    di_exact = -P11*EYSI*EPS0                        # dI/I = -pi11*sxx
    print(f'[압저항 종] dI/I = {di*100:.2f}%, Smith(1954) pi11 예측 '
          f'{di_exact*100:.2f}% (sxx = {EYSI*EPS0/1e6:.0f} MPa)')

    # ---- 2b. 횡방향 압저항: y-면 변위로 syy = E*eps0, dI/I = -pi12*syy ----
    pins_t = ['XL, 5, 5, 0.', 'PINA, 7, 7, 0.', 'PINB, 7, 7, 0.',
              'BLK, 4, 4, 0.']
    y0 = ['YL, 6, 6, 0.', 'YR, 6, 6, 0.']
    yc = ['YL, 6, 6, 0.', f'YR, 6, 6, {EPS0*WM*UM:.6e}']
    steps = [(0., 0., 20, y0), (VG, 0.05, 10, y0), (VG, 0.05, 4, yc)]
    _, _, _, ev = run_job('ex19pzt', 0.0, steps, pins_t)
    i0t = Q*NI*ev[1][1][0][1]
    i2 = Q*NI*ev[2][1][0][1]
    dit = i2/i0t - 1
    dit_exact = -P12*EYSI*EPS0                       # dI/I = -pi12*syy
    print(f'[압저항 횡] dI/I = {dit*100:.2f}%, Smith(1954) pi12 예측 '
          f'{dit_exact*100:.2f}% (syy = {EYSI*EPS0/1e6:.0f} MPa)')

    # ---- 3. full loop: 바닥 고정, HSCALE=200, VD 소인 ----
    clamp = ['BLK, 5, 7, 0.', 'BLK, 4, 4, 0.']
    steps = [(0., 0., 20, []), (VG, 0., 10, [])]
    vds = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    steps += [(VG, vd, 8, []) for vd in vds]
    _, _, nid3, ev = run_job('ex19hot', HSCALE, steps, clamp)
    print(' VD    I_ETM         dTmax    balance   (ET-only I: ex18hot)')
    ietm, bal_w = [], 0.0
    for k, vd in enumerate(vds):
        d, rf = ev[2+k]
        i_d = abs(Q*NI*rf[0][1])
        q_out = abs(rf[1][3])
        bal = abs(q_out - HSCALE*i_d*vd)/(HSCALE*i_d*vd)
        bal_w = max(bal_w, bal)
        ietm.append(i_d)
        print(f'{vd:4.1f}  {i_d:.4e}   {d[:, 3].max():6.1f}   {bal*100:5.2f}%')

    # ---- 그림: I_D-V_D + dT 맵 (ETM) ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax0.plot(vds, np.array(ietm)*1e6, 's-', color='crimson',
             label='electro-thermo-mechanical')
    ax0.set_xlabel('$V_D$ [V]')
    ax0.set_ylabel('$I_D$ [$\\mu$A]')
    ax0.set_title(f'$I_D$-$V_D$ (VG={VG:g}V, 7-dof monolithic)')
    ax0.legend()
    ax0.grid(alpha=0.3)
    ax0.text(0.98, 0.05,
             f'expansion check: {err_exp*100:.2f}% err\n'
             f'piezo $\\pi_{{11}}$: {di*100:.1f}% vs {di_exact*100:.1f}%, '
             f'$\\pi_{{12}}$: {dit*100:.1f}% vs {dit_exact*100:.1f}% '
             f'(Smith 1954)',
             transform=ax0.transAxes, ha='right', fontsize=8)
    dh = ev[-1][0]
    y0 = [nid3[(ix, iz, 0)] - 1 for iz in range(len(zs))
          for ix in range(len(xs))]
    dT = dh[y0, 3].reshape(len(zs), len(xs))
    Xf, Zf = np.meshgrid(xs, zs)
    pc = ax1.pcolormesh(Xf, -Zf, dT, cmap='hot', shading='gouraud')
    fig.colorbar(pc, ax=ax1, label='$\\Delta T$ [K]')
    ax1.set_title(f'$\\Delta T$ (VD={vds[-1]:g}V, bottom clamped)')
    ax1.set_xlabel('x [$\\mu$m]')
    ax1.set_ylabel('z [$\\mu$m]')
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_etm')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    print('figure ->', os.path.abspath(out) + '.png')

    assert err_exp < 0.02                            # 1. 열팽창 닫힌형
    assert abs(di - di_exact) < 0.02                 # 2. Smith 1954 pi11
    assert abs(dit - dit_exact) < 0.015              # 2b. Smith 1954 pi12
    assert bal_w < 0.05                              # 3. 에너지 수지
    assert ietm[-1] < ietm[-2]                       # 3. 음의 출력 컨덕턴스
    print('check passed: monolithic 전기-열-기계 UEL — '
          '열팽창/압저항(Smith)/에너지수지/droop.')


if __name__ == '__main__':
    main()
