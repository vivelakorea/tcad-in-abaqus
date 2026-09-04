# -*- coding: utf-8 -*-
"""run_selfheating: monolithic 전기-열 NMOS UEL (uel_mos_et.f).

기존 3D NMOS(run_mosfet)와 같은 소자/격자에 격자온도 상승 dT(자유도 4)를
얹어 (psi, phi_n, phi_p, dT) 를 한 Newton 행렬로 푼다 (staggered 아님).
- 열: 정상상태 -k lap(T) = q_J (box method), Joule = 모서리별 (Jn+Jp)*dpsi.
- 피드백: SG의 t = dpsi/V_T(T), mu(T) ∝ (T/300)^-1.5.
- 열 BC: 기판 바닥(BLK) dT=0 (히트싱크), 나머지 단열.
- HSCALE(발열 배율): 이 토이 소자는 W=0.5um 라 실제 발열이 mK 수준 ->
  다중 핑거 전력소자처럼 발열을 HSCALE배 (κ 축소와 등가) 해서 droop 재현.

검증 (닫힌형 없음 -> CLAUDE 계획 1,2,5):
1. 약결합 극한: HSCALE=0 이면 등온 UEL과 같은 방정식 -> I_D vs Pao-Sah 수 %.
2. 에너지 수지: 히트싱크 반력(RF dof4) 총합 = HSCALE * I_D*V_D (Tellegen).
3. self-heating on/off I_D-V_D 비교: 고 V_D 에서 I_D droop.

사용: python run_selfheating.py   (Abaqus 필요, job 2개, ~수 분)
"""
import io
import os
import re
import subprocess

import numpy as np

import reference_mosfet as ref
from run_mosfet import (grids, contacts, Q, VT, NI, NA_SUB, ND_PLUS,
                        UM, XG0, XG1, XJ, WM, LG, TOXU)

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
VG = 3.0
VDS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
HSCALE = 200.0        # 발열 배율 (핑거 수 유사); 0 = 등온 기준


def write_inp(job, hscale):
    xs, zs, zox = grids()
    gate_ix = np.where((xs >= XG0 - 1e-9) & (xs <= XG1 + 1e-9))[0]
    nid = {}
    lines = []
    a = 0
    for iz, z in enumerate(np.r_[zox, zs]):
        for iy in (0, 1):
            for ix, x in enumerate(xs):
                zlev = iz - len(zox)                 # 음수 = 산화막층
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
    L = ['*HEADING', 'ex18 3D NMOS electro-thermal UEL (monolithic)',
         '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=3, COORDINATES=3,'
         ' VARIABLES=1, UNSYMM', '1,2,3,4', '*NODE'] + lines
    for key in ('EOX', 'ENP', 'EPSUB'):
        if conn_lines[key]:
            L.append(f'*ELEMENT, TYPE=U1, ELSET={key}')
            L += conn_lines[key]
    L += ['*UEL PROPERTY, ELSET=EOX', f'0, 0., {hscale:.6e}',
          f'*UEL PROPERTY, ELSET=ENP', f'1, {ND_PLUS/NI:.6e}, {hscale:.6e}',
          f'*UEL PROPERTY, ELSET=EPSUB',
          f'1, {-NA_SUB/NI:.6e}, {hscale:.6e}']
    sets = {'SRC': [], 'DRN': [], 'BLK': [], 'GATE': [], 'OXI': []}
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
    for nm, ids in sets.items():
        L.append(f'*NSET, NSET={nm}')
        ids = sorted(ids)
        for i in range(0, len(ids), 12):
            L.append(', '.join(map(str, ids[i:i+12])))
    L.append('*NSET, NSET=NALL, GENERATE')
    L.append(f'1, {a}, 1')
    psi_bulk = psp

    def bstep(vg, vd, ninc, prt):
        s = ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
             f'{1.0/ninc}, 1.0, 1e-9, {1.0/ninc}',
             '*CONTROLS, PARAMETERS=FIELD', '1e-6,,,,,,,',
             '*BOUNDARY',
             f'GATE, 1, 1, {vg/VT + psi_bulk:.8e}',
             'GATE, 2, 3, 0.', 'OXI, 2, 3, 0.',
             f'SRC, 1, 1, {psn:.8e}', 'SRC, 2, 3, 0.',
             f'DRN, 1, 1, {psn + vd/VT:.8e}',
             f'DRN, 2, 2, {vd/VT:.8e}', f'DRN, 3, 3, {vd/VT:.8e}',
             f'BLK, 1, 1, {psp:.8e}', 'BLK, 2, 3, 0.',
             'BLK, 4, 4, 0.']                        # 히트싱크 dT=0
        if prt:
            s += ['*NODE PRINT, NSET=NALL, FREQUENCY=999', 'U',
                  '*NODE PRINT, NSET=DRN, FREQUENCY=999, TOTALS=YES', 'RF',
                  '*NODE PRINT, NSET=BLK, FREQUENCY=999, TOTALS=YES', 'RF']
        s.append('*END STEP')
        return s

    L += bstep(0.0, 0.0, 20, True)                   # 1: 평형(도핑 램프)
    L += bstep(VG, 0.0, 10, False)                   # 2: VG 인가 (출력 상속)
    for vd in VDS:                                   # 3...: VD 램프
        L += bstep(VG, vd, 8, False)
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs, zs, nid


def parse(job, nn):
    """스텝별 (U[nn,4], [RF_DRN 총합, RF_BLK 총합]) 목록. dof4 반력은 RM1."""
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    tot = r'(-?\d+\.\d*(?:E[+-]\d+)?)'
    events = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        for pp in blk.split('NODE FOOT-')[1:]:
            hdr = pp.split('\n', 1)[0]
            if 'U1' in hdr:
                rows = re.findall(r'^\s+(\d+)\s+' + (num + r'\s+')*3 + num,
                                  pp, re.M)
                d = np.zeros((nn, 4))
                for r_ in rows:
                    d[int(r_[0]) - 1] = [float(v) for v in r_[1:]]
                events.append([d, []])
            elif 'RF1' in hdr and events:
                m = re.search(r'^ TOTAL\s+' + (tot + r'\s+')*3 + tot,
                              pp, re.M)
                if m:
                    events[-1][1].append([float(m.group(k))
                                          for k in range(1, 5)])
    return events


def run_job(job, hscale):
    xs, zs, nid = write_inp(job, hscale)
    nn = max(nid.values())
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = (f'abaqus job={job} user={os.path.join(HERE, "uel_mos_et.f")}'
           f' interactive cpus=1')
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                       text=True, errors='replace')
    assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    ev = parse(job, nn)
    assert len(ev) == 2 + len(VDS), f'{len(ev)} print blocks'
    return xs, zs, nid, ev[2:]                       # VD 스텝만


def main():
    os.makedirs(RUN, exist_ok=True)
    res = {}
    for job, hs in (('ex18iso', 0.0), ('ex18hot', HSCALE)):
        res[job] = run_job(job, hs)

    print(f'\n VG={VG:g}V, HSCALE={HSCALE:g}')
    print(' VD     I_iso         I_hot         dTmax    I_PaoSah      '
          'balance')
    worst_ps, worst_bal = 0.0, 0.0
    iso_i, hot_i = [], []
    for k, vd in enumerate(VDS):
        _, _, _, evi = res['ex18iso']
        _, _, _, evh = res['ex18hot']
        di, rfi = evi[k][0], evi[k][1]
        dh, rfh = evh[k][0], evh[k][1]
        i_iso = abs(Q*NI*rfi[0][1])                  # 드레인 전자 반력 -> A
        i_hot = abs(Q*NI*rfh[0][1])
        iso_i.append(i_iso)
        hot_i.append(i_hot)
        dtmax = dh[:, 3].max()
        assert abs(di[:, 3]).max() < 1e-8            # 등온 job 은 dT=0
        ip = ref.id_paosah(VG, vd) * WM / LG         # 소자 전류로 환산
        if vd <= 1.0:      # GCA 성립 영역만 (포화영역은 CLM으로 갈라짐이 정상)
            worst_ps = max(worst_ps, abs(i_iso - ip)/ip)
        # 에너지 수지: 히트싱크로 나간 열 = HSCALE * I_D*V_D
        q_out = abs(rfh[1][3])
        p_in = HSCALE * i_hot * vd
        bal = abs(q_out - p_in) / p_in
        worst_bal = max(worst_bal, bal)
        print(f'{vd:4.1f}  {i_iso:.4e}    {i_hot:.4e}   {dtmax:6.1f}   '
              f'{ip:.4e}   {bal*100:5.2f}%')
    droop = 1.0 - hot_i[-1]/iso_i[-1]
    print(f'worst dev vs Pao-Sah (등온) = {worst_ps*100:.2f}%')
    print(f'worst energy-balance error  = {worst_bal*100:.2f}%')
    print(f'I_D droop @ VD={VDS[-1]:g}V        = {droop*100:.2f}%')

    # ---- 그림: I_D-V_D droop + dT 맵 ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    xs, zs, nid, evh = res['ex18hot']
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax0.plot(VDS, np.array(iso_i)*1e6, 'o-', label='isothermal (HSCALE=0)')
    ax0.plot(VDS, np.array(hot_i)*1e6, 's-',
             label=f'self-heating (HSCALE={HSCALE:g})')
    ax0.set_xlabel('$V_D$ [V]')
    ax0.set_ylabel('$I_D$ [$\\mu$A]')
    ax0.set_title(f'$I_D$-$V_D$ droop (VG={VG:g}V, monolithic)')
    ax0.legend()
    ax0.grid(alpha=0.3)
    dh = evh[-1][0]
    y0 = [nid[(ix, iz, 0)] - 1 for iz in range(len(zs))
          for ix in range(len(xs))]
    dT = dh[y0, 3].reshape(len(zs), len(xs))
    Xf, Zf = np.meshgrid(xs, zs)
    pc = ax1.pcolormesh(Xf, -Zf, dT, cmap='hot', shading='gouraud')
    fig.colorbar(pc, ax=ax1, label='$\\Delta T$ [K]')
    ax1.set_title(f'lattice temperature rise (VD={VDS[-1]:g}V)')
    ax1.set_xlabel('x [$\\mu$m]')
    ax1.set_ylabel('z [$\\mu$m]')
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_selfheating')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    print('figure ->', os.path.abspath(out) + '.png')

    assert worst_ps < 0.10                           # 1. 약결합 극한
    assert worst_bal < 0.05                          # 2. 에너지 수지
    assert droop > 0.03                              # 3. droop 재현
    print('check passed: monolithic 전기-열 UEL — 약결합/에너지수지/droop.')


if __name__ == '__main__':
    main()
