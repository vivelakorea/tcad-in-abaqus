# -*- coding: utf-8 -*-
"""make_demo: 게이트 전압 소인 데모 GIF 생성.

한 job으로 VG를 0 -> 3V로 30 증분 램프(VD=50mV)하며 매 증분의 절점해와
드레인 반력을 찍고, 두 패널 GIF로 만든다:
  (좌) 전자 밀도 log10 n --- 문턱(~1.3V)을 넘으며 반전 채널이 켜지는 순간
  (우) I_D-V_G 전달곡선 --- UEL 점들이 Pao-Sah(1966) 곡선을 따라간다
사용: python make_demo.py   (Abaqus + matplotlib/Pillow 필요, ~수 분)
출력: ../docs/demo_turnon.gif
"""
import io
import os
import subprocess

import numpy as np

import reference_mosfet as ref
import run_mosfet as m

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
VD, VGMAX, NINC = 0.05, 3.0, 30


def write_inp(job):
    xs, zs, zox = m.grids()
    # 격자/절점/요소/집합 부분은 run_mosfet의 inp를 재사용하고 스텝만 교체
    m.write_inp(job)
    txt = io.open(os.path.join(RUN, job + '.inp')).read()
    head = txt[:txt.index('*STEP')]
    psn, psp = m.contacts()

    def bc(vg, vd):
        return ['*BOUNDARY',
                f'GATE, 1, 1, {vg/m.VT + psp:.8e}',
                'GATE, 2, 3, 0.', 'OXI, 2, 3, 0.',
                f'SRC, 1, 1, {psn:.8e}', 'SRC, 2, 3, 0.',
                f'DRN, 1, 1, {psn + vd/m.VT:.8e}',
                f'DRN, 2, 2, {vd/m.VT:.8e}', f'DRN, 3, 3, {vd/m.VT:.8e}',
                f'BLK, 1, 1, {psp:.8e}', 'BLK, 2, 3, 0.']
    L = [head.rstrip()]
    # 1: 평형(도핑 램프), 2: VD=50mV --- 출력 없음
    for vg, vd, ninc in [(0.0, 0.0, 20), (0.0, VD, 2)]:
        L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
              f'{1.0/ninc}, 1.0, 1e-9, {1.0/ninc}',
              '*CONTROLS, PARAMETERS=FIELD', '1e-6,,,,,,,'] + bc(vg, vd) + ['*END STEP']
    # 3: VG 램프, 매 증분 출력
    L += ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
          f'{1.0/NINC}, 1.0, 1e-9, {1.0/NINC}',
          '*CONTROLS, PARAMETERS=FIELD', '1e-6,,,,,,,'] + bc(VGMAX, VD) + [
          '*NODE PRINT, NSET=NALL, FREQUENCY=1', 'U',
          '*NODE PRINT, NSET=DRN, FREQUENCY=1, TOTALS=YES', 'RF', '*END STEP']
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs, zs, zox


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex17demo'
    xs, zs, zox = write_inp(job)
    nid = {}
    a = 0
    gate_ix = np.where((xs >= m.XG0 - 1e-9) & (xs <= m.XG1 + 1e-9))[0]
    for iz in range(-len(zox), len(zs)):
        for iy in (0, 1):
            for ix in range(len(xs)):
                if iz < 0 and ix not in gate_ix:
                    continue
                a += 1
                nid[(ix, iz, iy)] = a
    lck = os.path.join(RUN, job + '.lck')
    if os.path.exists(lck):
        os.remove(lck)
    cmd = f'abaqus job={job} user={os.path.join(HERE, "uel_mos.f")} interactive cpus=1'
    print('>', cmd)
    r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True, text=True)
    assert 'COMPLETED' in r.stdout, r.stdout[-500:]
    blocks = m.parse(job, a)
    assert len(blocks) == NINC, f'{len(blocks)} frames'

    # Pao-Sah 기준 곡선
    vg_ref = np.linspace(0.6, VGMAX, 25)
    id_ref = np.array([ref.id_paosah(v, VD) for v in vg_ref]) * 1e6
    y0 = [nid[(ix, iz, 0)] - 1 for iz in range(len(zs)) for ix in range(len(xs))]
    Xf, Zf = np.meshgrid(xs, zs)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.5))
    vgs, ids = [], []
    for k, (dU, rf) in enumerate(blocks, 1):
        vgs.append(VGMAX * k / NINC)
        ids.append(abs(m.Q * m.NI * rf[1]) * m.LG / m.WM * 1e6)

    def frame(k):
        ax1.clear(); ax2.clear()
        dU, rf = blocks[k]
        nel = m.NI * np.exp(np.clip(dU[y0, 0] - dU[y0, 1], -80, 80)
                            ).reshape(len(zs), len(xs))
        ax1.pcolormesh(Xf, -Zf, np.log10(np.maximum(nel, 1.0)),
                       cmap='inferno', shading='gouraud', vmin=2, vmax=20)
        ax1.set_ylim(-0.3, 0); ax1.set_xlabel('x [$\\mu$m]')
        ax1.set_ylabel('z [$\\mu$m]')
        ax1.set_title(f'electron density,  $V_G$ = {vgs[k]:.1f} V')
        ax2.plot(vg_ref, id_ref, 'k-', lw=1.2, label='Pao-Sah 1966 (exact)')
        ax2.plot(vgs[:k+1], ids[:k+1], 'ro', ms=4, label='Abaqus 3D UEL')
        ax2.plot(vgs[k], ids[k], 'ro', ms=9, mfc='none')
        ax2.set_xlim(0, VGMAX); ax2.set_ylim(-0.3, max(ids)*1.08)
        ax2.set_xlabel('$V_G$ [V]'); ax2.set_ylabel('$I_D$ [$\\mu$A] (W/L=1)')
        ax2.set_title(f'transfer curve, $V_D$ = {VD*1000:.0f} mV')
        ax2.legend(fontsize=8, loc='upper left')
        fig.tight_layout()

    anim = FuncAnimation(fig, frame, frames=len(blocks))
    out = os.path.join(HERE, '..', 'docs', 'demo_turnon.gif')
    anim.save(out, writer=PillowWriter(fps=5))
    print('demo ->', os.path.abspath(out))
    # 정합성: 마지막 프레임(VG=3, VD=0.05)을 Pao-Sah와 대조
    dev = abs(ids[-1]*1e-6 - ref.id_paosah(VGMAX, VD)) / ref.id_paosah(VGMAX, VD)
    print(f'final point vs Pao-Sah: {dev*100:.2f}%')
    assert dev < 0.05
    print('demo check passed.')


if __name__ == '__main__':
    main()
