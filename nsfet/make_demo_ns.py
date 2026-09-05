# -*- coding: utf-8 -*-
"""make_demo_ns: 나노시트 GAAFET 턴온 데모 GIF.

한 job 으로 VG 를 -0.4 -> 0.8 V (VDS=0.7V) 30 증분 램프하며 매 증분의
절점해/드레인 반력을 찍고 두 패널 GIF 를 만든다:
  (좌) A-A 단면 (y=0) 전자밀도 log10 n --- 시트 3장이 켜지는 순간
  (우) I_D-V_G' (LP 정렬) --- UEL 점들이 쌓이고, 옆에 Wang et al. (2023)
       논문 전달곡선 디지타이즈 점 (CC BY, 출처 표기) 을 나란히.
ODB 에도 매 증분 U/UVARM 필드가 실려 Abaqus Viewer 애니메이션도 가능.
사용: python make_demo_ns.py   (~10 min)  출력: ../docs/demo_nsfet.gif
"""
import io
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'jlfet'))
import run_nsfet as ns
import paper_wang2023 as paper

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
VG0, VG1, NINC = -0.4, 0.8, 30
VDD = ns.VDD


def write_inp(job):
    nn, nid, xs, ys, zs = ns.write_inp(job, 0.0)
    txt = io.open(os.path.join(RUN, job + '.inp')).read()
    head = txt[:txt.index('*STEP')]

    def bc(vg, vd):
        return ['*BOUNDARY',
                f'SRC, 1, 1, {ns.PSN:.8e}', 'SRC, 2, 3, 0.',
                f'DRN, 1, 1, {ns.PSN + vd/ns.VT:.8e}',
                f'DRN, 2, 2, {vd/ns.VT:.8e}', f'DRN, 3, 3, {vd/ns.VT:.8e}',
                f'GATE, 1, 1, {vg/ns.VT:.8e}', 'OXI, 2, 3, 0.']

    L = [head.rstrip()]
    for vg, vd, ninc in [(0., 0., 20), (VG0, VDD, 10)]:
        L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
              f'{1.0/ninc}, 1.0, 1e-9, {1.0/ninc}',
              '*CONTROLS, PARAMETERS=FIELD', '1e-4,,,,,,,1e-4'] \
             + bc(vg, vd) + ['*END STEP']
    L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
          f'{1.0/NINC}, 1.0, 1e-9, {1.0/NINC}',
          '*CONTROLS, PARAMETERS=FIELD', '1e-4,,,,,,,1e-4'] \
         + bc(VG1, VDD) + [
          '*NODE PRINT, NSET=NALL, FREQUENCY=1', 'U',
          '*NODE PRINT, NSET=DRN, FREQUENCY=1, TOTALS=YES', 'RF',
          '*OUTPUT, FIELD, FREQUENCY=1',
          '*NODE OUTPUT', 'U',
          '*ELEMENT OUTPUT, ELSET=EVIS', 'UVARM',
          '*END STEP']
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return nn, nid, xs, ys, zs


def parse3(job, nn):
    """증분별 [U[nn,3], RF_DRN 총합]."""
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    num = r'(-?\d\.\d+E[+-]\d+)'
    tot = r'(-?\d+\.\d*(?:E[+-]\d+)?)'
    ev = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        for pp in blk.split('NODE FOOT-')[1:]:
            hdr = pp.split('\n', 1)[0]
            if 'U1' in hdr:
                rows = re.findall(r'^\s+(\d+)\s+' + (num + r'\s+')*2 + num,
                                  pp, re.M)
                d = np.zeros((nn, 3))
                for r_ in rows:
                    d[int(r_[0]) - 1] = [float(v) for v in r_[1:]]
                ev.append([d, None])
            elif 'RF1' in hdr and ev:
                m = re.search(r'^ TOTAL\s+' + (tot + r'\s+')*2 + tot,
                              pp, re.M)
                if m:
                    ev[-1][1] = float(m.group(2))
    return [e for e in ev if e[1] is not None]


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex25demo'
    nn, nid, xs, ys, zs = write_inp(job)
    if not os.path.exists(os.path.join(RUN, job + '.dat')):
        lck = os.path.join(RUN, job + '.lck')
        if os.path.exists(lck):
            os.remove(lck)
        cmd = (f'abaqus job={job} user={os.path.abspath(ns.UELF)}'
               f' interactive cpus=1')
        print('>', cmd)
        r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                           text=True, errors='replace')
        assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    ev = parse3(job, nn)
    assert len(ev) >= NINC - 2, f'{len(ev)} frames'
    vgs = VG0 + (VG1 - VG0)*np.arange(1, len(ev) + 1)/len(ev)
    ids = 4.0*np.abs(ns.Q*ns.NI*np.array([e[1] for e in ev]))
    lg10 = np.log10(np.maximum(ids, 1e-30))
    vg_lp = np.interp(np.log10(ns.IOFF_LP), lg10, vgs)   # LP 정렬 (WF 역할)
    # A-A 단면 (y=0) 격자 (구멍은 NaN)
    zmap = np.full((len(zs), len(xs)), np.nan)
    idx = {}
    for iz in range(len(zs)):
        for ix in range(len(xs)):
            if (ix, 0, iz) in nid:
                idx[(iz, ix)] = nid[(ix, 0, iz)] - 1
    Xf, Zf = np.meshgrid(xs, zs)
    zm2 = 2*(ns.TNS + ns.TSP + ns.TNS/2)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.7))
    pv, pi = zip(*paper.FIG3B_IDEAL)

    def frame(k):
        ax1.clear()
        ax2.clear()
        dU = ev[k][0]
        nel = zmap.copy()
        for (iz, ix), j in idx.items():
            nel[iz, ix] = np.log10(max(
                ns.NI*np.exp(min(max(dU[j, 0] - dU[j, 1], -80), 80)), 1.0))
        # 거울 복원 (z)
        Zfull = np.r_[zs, 2*(ns.TNS + ns.TSP + ns.TNS/2) - zs[::-1][1:]]
        Nfull = np.vstack([nel, nel[::-1][1:]])
        Xg, Zg = np.meshgrid(xs, Zfull)
        ax1.pcolormesh(Xg, -Zg, np.ma.masked_invalid(Nfull), cmap='inferno',
                       shading='gouraud', vmin=4, vmax=21)
        ax1.set_facecolor('0.85')
        ax1.set_xlabel('x [nm]')
        ax1.set_ylabel('-z [nm]')
        ax1.set_title(f'electron density (A-A),  $V_G\'$ = '
                      f'{vgs[k]-vg_lp:+.2f} V')
        ax2.semilogy(np.array(pv), np.array(pi), 'ks', ms=4, mfc='none',
                     label='Wang et al. 2023 (digitized)')
        ax2.semilogy(vgs[:k+1]-vg_lp, ids[:k+1], 'ro', ms=4,
                     label='Abaqus UEL')
        ax2.semilogy(vgs[k]-vg_lp, ids[k], 'ro', ms=9, mfc='none')
        ax2.set_xlim(VG0-vg_lp, VG1-vg_lp)
        ax2.set_ylim(1e-14, 1e-3)
        ax2.set_xlabel("$V_G'$ [V] (LP-aligned)")
        ax2.set_ylabel('$I_D$ [A]')
        ax2.set_title(f'transfer @ $V_{{DS}}$ = {VDD} V')
        ax2.legend(fontsize=8, loc='lower right')
        ax2.grid(alpha=0.3, which='both')
        fig.tight_layout()

    anim = FuncAnimation(fig, frame, frames=len(ev))
    out = os.path.join(HERE, '..', 'docs', 'demo_nsfet.gif')
    anim.save(out, writer=PillowWriter(fps=5))
    print('demo ->', os.path.abspath(out))
    span = lg10.max() - lg10.min()
    assert span > 8.0, f'{span:.1f} decades'
    print(f'demo check passed: {span:.1f} decades swing.')


if __name__ == '__main__':
    main()
