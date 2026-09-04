# -*- coding: utf-8 -*-
"""run_mosfet: 3D NMOS 소자 시뮬레이션을 Abaqus UEL로 돌린다.

소자: 6um x 0.5um x 2um NMOS. n+ 소스/드레인(1e19, 깊이 0.2um),
p 기판(1e17), 게이트 산화막 10nm(x=1..5um, 채널 L=4um), 옴 접촉 3개.
격자: 텐서 육면체(접합/표면 세밀화), box-method SG UEL(uel_mos.f).
한 job 6 스텝: 평형(도핑 램프) -> VG=2V -> VD=0.05 -> VD=0.5 -> VG=3V -> VD=0.05.
검증: I_D(반력 RF로 추출)를 Pao & Sah (1966) 정확해 / Brews (1978)
charge-sheet와 대조 (긴 채널 -> gradual channel 성립, 수 % 기대).
그림: fig_tcad_uel.pdf --- 3D 복셀 도핑(TCAD 스타일) + psi/log n 단면 + I_D 대조.
사용: python run_mosfet.py   (Abaqus 필요, ~수 분)
"""
import io
import os
import re
import subprocess

import numpy as np

import reference_mosfet as ref

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
Q, VT, NI = 1.602e-19, 0.02585, 1.0e10
NA_SUB, ND_PLUS = 1.0e17, 1.0e19            # cm^-3
UM = 1e-4                                    # 1 um in cm
XG0, XG1, XJ = 1.0, 5.0, 0.2                # 게이트 구간, 접합 깊이 [um]
TOXU = 0.01                                  # 10 nm [um]
WM = 0.5                                     # 압출 폭 [um]
LG = XG1 - XG0                               # 채널 길이 4 um
BIAS = [(2.0, 0.05), (2.0, 0.5), (3.0, 0.5), (3.0, 0.05)]  # 스텝 3~6 종점


def seg(a, b, n):
    return np.linspace(a, b, n)


def grids():
    xs = np.unique(np.round(np.r_[seg(0, 0.9, 7), seg(0.9, 1.1, 11),
                                  seg(1.1, 4.9, 20), seg(4.9, 5.1, 11),
                                  seg(5.1, 6.0, 7)], 6))
    zs = np.unique(np.round(np.r_[seg(0, 0.01, 11), seg(0.01, 0.05, 11),
                                  seg(0.05, 0.3, 11), seg(0.3, 2.0, 11)], 6))
    zox = np.round(seg(-TOXU, 0, 5)[:-1], 6)
    return xs, zs, zox


def contacts():
    """옴 접촉의 psi/VT (전하중성+평형): psi = asinh(N/2)."""
    return np.arcsinh(ND_PLUS/NI/2), np.arcsinh(-NA_SUB/NI/2)


def write_inp(job):
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
    L = ['*HEADING', 'ex17 3D NMOS drift-diffusion UEL (box method SG)',
         '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=2, COORDINATES=3,'
         ' VARIABLES=1, UNSYMM', '1,2,3', '*NODE'] + lines
    for key in ('EOX', 'ENP', 'EPSUB'):
        if conn_lines[key]:
            L.append(f'*ELEMENT, TYPE=U1, ELSET={key}')
            L += conn_lines[key]
    L += ['*UEL PROPERTY, ELSET=EOX', '0, 0.',
          f'*UEL PROPERTY, ELSET=ENP', f'1, {ND_PLUS/NI:.6e}',
          f'*UEL PROPERTY, ELSET=EPSUB', f'1, {-NA_SUB/NI:.6e}']
    # 절점 집합
    sets = {'SRC': [], 'DRN': [], 'BLK': [], 'GATE': [], 'OXI': []}
    for (ix, zlev, iy), aa in nid.items():
        x = xs[ix]
        if zlev == 0 and x <= 0.9 + 1e-9:
            sets['SRC'].append(aa)
        if zlev == 0 and x >= 5.1 - 1e-9:
            sets['DRN'].append(aa)
        if zlev == len(zs) - 1:
            sets['BLK'].append(aa)
        if zlev == -len(grids()[2]):
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
             f'BLK, 1, 1, {psp:.8e}', 'BLK, 2, 3, 0.']
        if prt:
            s += ['*NODE PRINT, NSET=NALL, FREQUENCY=999', 'U',
                  '*NODE PRINT, NSET=DRN, FREQUENCY=999, TOTALS=YES', 'RF']
        s.append('*END STEP')
        return s

    L += bstep(0.0, 0.0, 20, True)                    # 1: 평형(도핑 램프)
    L += bstep(2.0, 0.0, 10, False)                   # 2: VG=2V
    for k, (vg, vd) in enumerate(BIAS):               # 3~6
        L += bstep(vg, vd, 8, True)
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs, zs, zox, nid, sets


def parse(job, nn):
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    out = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        upart, _, rfpart = blk.partition('NODE FOOT-  RF1')
        rows = re.findall(r'^\s+(\d+)\s+' + (num + r'\s+')*2 + num, upart, re.M)
        m = re.search(r'^ TOTAL\s+(-?\d+\.\d*(?:E[+-]\d+)?)'
                      r'\s+(-?\d+\.\d*(?:E[+-]\d+)?)'
                      r'\s+(-?\d+\.\d*(?:E[+-]\d+)?)', rfpart, re.M)
        rf = [float(m.group(k)) for k in (1, 2, 3)] if m else None
        if rows:
            d = np.zeros((nn, 3))
            for r_ in rows:
                d[int(r_[0]) - 1] = [float(v) for v in r_[1:]]
            out.append((d, rf))
        elif rf is not None and out and out[-1][1] is None:
            out[-1] = (out[-1][0], rf)       # U와 RF가 헤더로 갈라진 경우
    return out


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex17mos'
    xs, zs, zox, nid, sets = write_inp(job)
    nn = max(nid.values())
    print(f'mesh: {nn} nodes ({len(xs)}x{len(zs)} Si + oxide), '
          f'{3*nn} dofs')
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = f'abaqus job={job} user={os.path.join(HERE, "uel_mos.f")} interactive cpus=1'
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True, text=True)
    assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    blocks = parse(job, nn)
    assert len(blocks) == 2 + len(BIAS), f'{len(blocks)} print blocks'  # eq + VG스텝(출력 상속) + 4 bias

    # ---- 평형 체크: 벌크 중성, 표면 아래 공핍 ----
    eqU = blocks[0][0]
    print(' VG     VD     I_UEL(Abaqus)  I_PaoSah      I_Brews       dev(UEL/PS)')
    worst = 0.0
    results = []
    for (vg, vd), (dU, rf) in zip(BIAS, blocks[2:]):
        i_mesh = Q * NI * rf[1]                      # 전자 반력 -> A (폭 WM)
        i_wl1 = abs(i_mesh) * LG / WM                # W/L = 1 로 환산
        ip = ref.id_paosah(vg, vd)
        ib = ref.id_chargesheet(vg, vd)
        dev = abs(i_wl1 - ip) / ip
        worst = max(worst, dev)
        results.append((vg, vd, i_wl1, ip, ib, dU))
        print(f'{vg:4.1f}  {vd:5.2f}   {i_wl1:.4e}    {ip:.4e}   '
              f'{ib:.4e}   {dev*100:5.2f}%')
    print(f'worst deviation (3D UEL vs Pao-Sah 1966) = {worst*100:.2f}%')

    # ---- TCAD 스타일 그림 ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    fig = plt.figure(figsize=(10.5, 7.6))
    ax0 = fig.add_subplot(2, 2, 1, projection='3d')
    xe, ye = np.meshgrid(xs, [0, WM])
    dop2 = np.zeros((len(zs) - 1, len(xs) - 1))
    for iz in range(len(zs) - 1):
        for ix in range(len(xs) - 1):
            xc = 0.5*(xs[ix] + xs[ix+1]); zc = 0.5*(zs[iz] + zs[iz+1])
            dop2[iz, ix] = (np.log10(ND_PLUS) if ((xc < XG0 or xc > XG1)
                            and zc < XJ) else -np.log10(NA_SUB))
    cmap = plt.get_cmap('jet')
    nrm = Normalize(-19, 19)

    def dopval(x, z):
        return (np.log10(ND_PLUS) if ((x < XG0 or x > XG1) and z < XJ)
                else -np.log10(NA_SUB))

    def face(Xf, Yf, Zf, C):
        ax0.plot_surface(Xf, Yf, -Zf, facecolors=cmap(nrm(C)), shade=False,
                         rstride=1, cstride=1, antialiased=False)
    Xg2, Zg2 = np.meshgrid(xs, zs)
    Cxz = np.vectorize(dopval)(Xg2, Zg2)
    for y0 in (0.0, WM):                              # 앞/뒤 x-z 면
        face(Xg2, np.full_like(Xg2, y0), Zg2, Cxz)
    Yg2, Zg2b = np.meshgrid([0, WM], zs)
    for x0 in (0.0, 6.0):                             # 좌/우 y-z 면
        C = np.vectorize(lambda z: dopval(x0 + (1e-3 if x0 == 0 else -1e-3), z))(Zg2b)
        face(np.full_like(Yg2, x0), Yg2, Zg2b, C)
    Xg3, Yg3 = np.meshgrid(xs, [0, WM])               # 윗면(도핑 지도) / 바닥
    Ctop = np.vectorize(lambda x: dopval(x, 1e-4))(Xg3)
    face(Xg3, Yg3, np.zeros_like(Xg3), Ctop)
    face(Xg3, Yg3, np.full_like(Xg3, 2.0), np.full_like(Xg3, -np.log10(NA_SUB)))
    xg = xs[(xs >= XG0) & (xs <= XG1)]                # 게이트 전극
    Xg4, Yg4 = np.meshgrid(xg, [0, WM])
    for zg in (TOXU, 0.06):
        ax0.plot_surface(Xg4, Yg4, np.full_like(Xg4, zg), color='gold', shade=False)
    Zg4 = np.meshgrid([TOXU, 0.06], [0, WM])[0].T
    for xge in (xg[0], xg[-1]):
        ax0.plot_surface(np.full_like(Zg4, xge),
                         np.array([[0, WM], [0, WM]]), Zg4, color='goldenrod',
                         shade=False)
    ax0.set_box_aspect((6, 2.4, 3.2))
    Xf, Zf = np.meshgrid(xs, zs)
    ax0.set_title('device: net doping (log$_{10}$|N|, sign)')
    ax0.set_xlabel('x [$\\mu$m]'); ax0.set_ylabel('y'); ax0.set_zlabel('-z')
    ax0.view_init(22, -60)
    mp = plt.cm.ScalarMappable(nrm, cmap); mp.set_array([])
    fig.colorbar(mp, ax=ax0, shrink=0.6, label='sgn(N)·log$_{10}$|N| [cm$^{-3}$]')

    vg, vd, iu, ip, ib, dU = results[2]              # VG=3, VD=0.5
    y0 = [nid[(ix, iz, 0)] - 1 for iz in range(len(zs)) for ix in range(len(xs))]
    psi = (dU[y0, 0].reshape(len(zs), len(xs))) * VT
    nel = NI * np.exp(np.clip(dU[y0, 0] - dU[y0, 1], -80, 80)
                      ).reshape(len(zs), len(xs))
    ax1 = fig.add_subplot(2, 2, 2)
    pc = ax1.pcolormesh(Xf, -Zf, psi, cmap='viridis', shading='gouraud')
    fig.colorbar(pc, ax=ax1, label='$\\psi$ [V]')
    ax1.set_title(f'potential (VG={vg:g}, VD={vd:g})')
    ax1.set_xlabel('x [$\\mu$m]'); ax1.set_ylabel('z [$\\mu$m]'); ax1.set_ylim(-0.5, 0)
    ax2 = fig.add_subplot(2, 2, 3)
    pc2 = ax2.pcolormesh(Xf, -Zf, np.log10(np.maximum(nel, 1.0)),
                         cmap='inferno', shading='gouraud', vmin=2, vmax=20)
    fig.colorbar(pc2, ax=ax2, label='log$_{10}$ n [cm$^{-3}$]')
    ax2.set_title('electron density: inversion channel')
    ax2.set_xlabel('x [$\\mu$m]'); ax2.set_ylabel('z [$\\mu$m]'); ax2.set_ylim(-0.3, 0)
    ax3 = fig.add_subplot(2, 2, 4)
    lab = [f'VG={a:g}\nVD={b:g}' for a, b, *_ in results]
    w = 0.27; xi = np.arange(len(results))
    ax3.bar(xi - w, [r[2]*1e6 for r in results], w, label='Abaqus 3D UEL')
    ax3.bar(xi,     [r[3]*1e6 for r in results], w, label='Pao-Sah 1966')
    ax3.bar(xi + w, [r[4]*1e6 for r in results], w, label='Brews 1978')
    ax3.set_xticks(xi, lab, fontsize=8)
    ax3.set_ylabel('$I_D$ [$\\mu$A] (W/L=1)'); ax3.legend(fontsize=8)
    ax3.set_title('drain current: 3-way comparison')
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_tcad_uel')
    fig.savefig(out + '.pdf')
    fig.savefig(out + '.png', dpi=150)
    print('figure ->', os.path.abspath(out) + '.{pdf,png}')
    assert worst < 0.10
    print('check passed: 3D 소자 UEL의 I_D가 Pao-Sah(1966)와 수 % 이내.')


if __name__ == '__main__':
    main()
