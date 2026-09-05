# -*- coding: utf-8 -*-
"""run_nsfet: Wang et al., Electronics 12, 770 (2023) 3단 적층 나노시트
GAAFET 를 Table 1 입력으로 재현 (jlfet/paper_wang2023.py).

기하 (nm, Table 1): Lg=18, Lsp=5, NS 36x6 x3장(간격 12), S/D 2e20,
채널 1e15(p), Vdd=0.7. EOT=1 가정(미명시), WF는 논문 절차대로
Ioff=0.1nA 타깃팅 -- 상수 mu 모델에선 강체 이동이라 후처리로 충분.

모델링:
- 이중 거울 대칭: z=가운데 시트 중앙, y=시트 폭 중앙 -> 전류 x4.
- 금속 게이트는 메쉬 구멍: EOT 셸 외면이 게이트 Dirichlet.
- 게이트 창 = 시트+셸(너머 금속), 스페이서 창 = 갭 전체 유전체,
  S/D 창 = 갭까지 통짜 n+ 에피. Lrcs = 시트 도핑 경계의 게이트쪽 이동.
- 바닥 서브핀/PTS 생략 (논문: Hrcs=0 이상 구조에선 Isub>97% 지배).

비교 목표 (논문 본문 명시): Lrcs=5nm 에서 Ioff 가 이상 구조의 10배 이상,
SS/on-off 열화, Ion 증가 트렌드. (Ion 절대값은 vsat/ballistic 보정
스택이 없어 비교 불가 -- 문서화.)

사용: python run_nsfet.py [lrcs_nm]   (기본 0; job 1개, ~수 분)
"""
import io
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
UELF = os.path.join(HERE, 'uel_ns.f')       # jlfet/uel_jl.f + UVARM 브리지
Q, VT, NI = 1.602e-19, 0.02585, 1.0e10
NM = 1e-7
LG, LSP, LSD = 18.0, 5.0, 10.0               # nm
EOT = 1.0                                    # 가정 (미명시)
WNS2 = 18.0                                  # 시트 반폭
TNS, TSP = 6.0, 12.0
X0 = LSD + LSP + LG/2                        # 게이트 중심 = 24
XT = 2*X0                                    # 전장 48
NSD, NCH = 2e20/NI, -1e15/NI                 # ni 단위 (채널 p형)
PSN = np.arcsinh(NSD/2)
VDD = 0.7
VGS = list(np.round(np.arange(-0.4, 0.81, 0.1), 3))
IOFF_LP = 1e-10                              # 0.1 nA (논문 LP 기준)


def grids():
    seg = np.linspace
    e0, e1 = X0 - LG/2, X0 + LG/2            # 게이트 모서리 15, 33
    xs = np.unique(np.round(np.r_[
        seg(0, LSD, 6), seg(LSD, e0, 11), seg(e0, e1, 19),
        seg(e1, XT - LSD, 11), seg(XT - LSD, XT, 6)], 4))
    ys = np.unique(np.round(np.r_[seg(0, 14, 8), seg(14, WNS2, 9),
                                  seg(WNS2, WNS2 + EOT, 3)], 4))
    zs = np.unique(np.round(np.r_[
        seg(-EOT, 0, 3), seg(0, TNS, 13), seg(TNS, TNS + EOT, 3),
        seg(TNS + TSP - EOT, TNS + TSP, 3),
        seg(TNS + TSP, TNS + TSP + TNS/2, 7)], 4))
    return xs, ys, zs


SH = [(0.0, TNS), (TNS + TSP, TNS + TSP + TNS)]  # 시트 z 범위 (반쪽: 2번째 절반)


def in_sheet(z):
    return any(a - 1e-6 <= z <= b + 1e-6 for a, b in SH)


def in_shell_z(z):
    return any((a - EOT - 1e-6 <= z <= a + 1e-6)
               or (b - 1e-6 <= z <= b + EOT + 1e-6) for a, b in SH)


def region(xc, yc, zc, lrcs):
    """요소 중심 -> 'NSD'/'NCH'/'EOX'/None(금속·외부)."""
    dxl, dxr = xc, XT - xc                   # 좌우 대칭 거리
    zone = ('sd' if min(dxl, dxr) < LSD
            else 's' if min(dxl, dxr) < LSD + LSP else 'g')
    sheet = in_sheet(zc) and yc < WNS2
    if zone == 'sd':
        return 'NSD' if (yc < WNS2 and 0 < zc) else None
    if sheet:                                # 채널 시트 (Lrcs 만큼 n+ 침식)
        return 'NSD' if min(dxl, dxr) < LSD + lrcs else 'NCH'
    if zone == 's':                          # 스페이서: 갭·셸 전부 유전체
        return 'EOX' if (yc < WNS2 + EOT and zc > -EOT) else None
    # 게이트 창: 시트 주변 EOT 셸만 유전체, 나머지는 금속(구멍)
    shell = ((yc < WNS2 and in_shell_z(zc))
             or (WNS2 < yc < WNS2 + EOT and (in_sheet(zc) or in_shell_z(zc))))
    return 'EOX' if shell else None


def write_inp(job, lrcs):
    xs, ys, zs = grids()
    nid = {}
    conn = {'ENSD': [], 'ENCH': [], 'EOX': []}
    touch_si = set()
    a = [0]

    def node(ix, iy, iz):
        key = (ix, iy, iz)
        if key not in nid:
            a[0] += 1
            nid[key] = a[0]
        return nid[key]

    e = 0
    for iz in range(len(zs)-1):
        for iy in range(len(ys)-1):
            for ix in range(len(xs)-1):
                xc = 0.5*(xs[ix] + xs[ix+1])
                yc = 0.5*(ys[iy] + ys[iy+1])
                zc = 0.5*(zs[iz] + zs[iz+1])
                r = region(xc, yc, zc, lrcs)
                if r is None:
                    continue
                nn = [node(ix, iy, iz), node(ix+1, iy, iz),
                      node(ix+1, iy+1, iz), node(ix, iy+1, iz),
                      node(ix, iy, iz+1), node(ix+1, iy, iz+1),
                      node(ix+1, iy+1, iz+1), node(ix, iy+1, iz+1)]
                e += 1
                key = {'NSD': 'ENSD', 'NCH': 'ENCH', 'EOX': 'EOX'}[r]
                conn[key].append(f'{e}, ' + ', '.join(map(str, nn)))
                if key != 'EOX':
                    touch_si.update(nn)
    lines = [f'{v}, {xs[ix]*NM:.8e}, {ys[iy]*NM:.8e}, {zs[iz]*NM:.8e}'
             for (ix, iy, iz), v in sorted(nid.items(), key=lambda kv: kv[1])]
    L = ['*HEADING', f'ex25 stacked-NS GAAFET (Wang 2023) Lrcs={lrcs}nm',
         '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=2, COORDINATES=3,'
         ' VARIABLES=1, UNSYMM', '1,2,3', '*NODE'] + lines
    for key in ('ENSD', 'ENCH', 'EOX'):
        if conn[key]:
            L.append(f'*ELEMENT, TYPE=U1, ELSET={key}')
            L += conn[key]
    L += ['*UEL PROPERTY, ELSET=ENSD', f'1, {NSD:.6e}',
          '*UEL PROPERTY, ELSET=ENCH', f'1, {NCH:.6e}',
          '*UEL PROPERTY, ELSET=EOX', '0, 0.']
    # ODB 가시화용 더미 오버레이: 같은 연결성의 C3D8 을 초저강성으로 겹침
    # -> Viewer 에 메쉬가 렌더되고 U1/U2/U3 = psi/phi_n/phi_p [VT] 컨투어.
    # 재료별 elset (VNSD/VNCH/VOX) 이라 색상 구분도 가능.
    for key in ('ENSD', 'ENCH', 'EOX'):
        if conn[key]:
            L.append(f'*ELEMENT, TYPE=C3D8, ELSET=V{key[1:]}')
            L += [f'{int(c.split(",")[0]) + 500000},'
                  + c.split(',', 1)[1] for c in conn[key]]
    L += ['*ELSET, ELSET=EVIS', 'VNSD, VNCH, VOX',
          '*MATERIAL, NAME=MDUM', '*ELASTIC', '1e-20, 0.',
          '*USER OUTPUT VARIABLES', '2',     # UVARM1/2 = log10 n / log10 p
          '*SOLID SECTION, ELSET=EVIS, MATERIAL=MDUM']
    e0, e1 = X0 - LG/2, X0 + LG/2
    sets = {'SRC': [], 'DRN': [], 'GATE': [], 'OXI': []}
    for (ix, iy, iz), v in nid.items():
        x, y, z = xs[ix], ys[iy], zs[iz]
        if ix == 0:
            sets['SRC'].append(v)
        if ix == len(xs)-1:
            sets['DRN'].append(v)
        if (e0 - 1e-6 <= x <= e1 + 1e-6
            and (abs(z - (-EOT)) < 1e-6 or abs(y - (WNS2+EOT)) < 1e-6
                 or (abs(z - (TNS+EOT)) < 1e-6)
                 or (abs(z - (TNS+TSP-EOT)) < 1e-6))):
            sets['GATE'].append(v)
        if v not in touch_si:
            sets['OXI'].append(v)
    for nm_, ids in sets.items():
        L.append(f'*NSET, NSET={nm_}')
        ids = sorted(set(ids))
        for i in range(0, len(ids), 12):
            L.append(', '.join(map(str, ids[i:i+12])))
    L.append('*NSET, NSET=NALL, GENERATE')
    L.append(f'1, {a[0]}, 1')

    def bstep(vg, vd, ninc, prt):
        s = ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
             f'{1.0/ninc}, 1.0, 1e-9, {1.0/ninc}',
             '*CONTROLS, PARAMETERS=FIELD', '1e-4,,,,,,,1e-4',
             '*BOUNDARY',
             f'SRC, 1, 1, {PSN:.8e}', 'SRC, 2, 3, 0.',
             f'DRN, 1, 1, {PSN + vd/VT:.8e}',
             f'DRN, 2, 2, {vd/VT:.8e}', f'DRN, 3, 3, {vd/VT:.8e}',
             f'GATE, 1, 1, {vg/VT:.8e}',
             'OXI, 2, 3, 0.']
        if prt:
            s += ['*NODE PRINT, NSET=DRN, FREQUENCY=999, TOTALS=YES', 'RF',
                  '*OUTPUT, FIELD, FREQUENCY=999',      # ODB: 스텝 말 필드
                  '*NODE OUTPUT', 'U',
                  '*ELEMENT OUTPUT, ELSET=EVIS', 'UVARM']
        s.append('*END STEP')
        return s

    L += bstep(0., 0., 20, False)                    # 1: 평형 (도핑 램프)
    L += bstep(VGS[0], VDD, 10, True)                # 2: VD=0.7 + VG 시작
    for vg in VGS[1:]:
        L += bstep(vg, VDD, 4, False)
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return a[0], nid, xs, ys, zs


def parse_rf(job):
    txt = io.open(os.path.join(RUN, job + '.dat'), errors='ignore').read()
    tot = r'(-?\d+\.\d*(?:E[+-]\d+)?|NaN)'
    out = []
    for blk in txt.split('N O D E   O U T P U T')[1:]:
        for pp in blk.split('NODE FOOT-')[1:]:
            if 'RF1' in pp.split('\n', 1)[0]:
                m = re.search(r'^ TOTAL\s+' + (tot + r'\s+')*2 + tot,
                              pp, re.M)
                if m:
                    out.append([float(m.group(k)) for k in (1, 2, 3)])
    return out


def run_job(job, lrcs):
    dat = os.path.join(RUN, job + '.dat')
    if os.path.exists(dat):
        print(f'> {job}: cached')
    else:
        nn = write_inp(job, lrcs)[0]
        print(f'> {job}: Lrcs={lrcs}nm, {nn} nodes')
        lck = os.path.join(RUN, job + '.lck')
        if os.path.exists(lck):
            os.remove(lck)
        cmd = (f'abaqus job={job} user={os.path.abspath(UELF)}'
               f' interactive cpus=1')
        r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                           text=True, errors='replace')
        assert 'COMPLETED' in r.stdout, r.stdout[-600:]
    rf = parse_rf(job)
    assert len(rf) == len(VGS), f'{len(rf)} blocks != {len(VGS)}'
    return 4.0*np.abs(Q*NI*np.array([r[1] for r in rf]))  # 이중 대칭 x4


def lp_metrics(i_d):
    """논문 LP 절차: VG축을 Ioff=0.1nA 에 정렬 -> Ioff/Ion/SS."""
    lg10 = np.log10(np.maximum(i_d, 1e-30))
    vg0 = np.interp(np.log10(IOFF_LP), lg10, VGS)    # I=0.1nA 인 raw VG
    ion = 10**np.interp(vg0 + VDD, VGS, lg10)        # VG' = Vdd
    m = (i_d > 1e-13) & (i_d < 1e-9)
    ss = 1000/np.polyfit(np.array(VGS)[m], lg10[m], 1)[0] if m.sum() >= 3 \
        else float('nan')
    return vg0, ion, ss


def compare():
    """논문 본문 명시 결과와 대조 (Lrcs = 0/3/5 nm)."""
    cur = {l: run_job(f'ex25r{l}', float(l)) for l in (0, 3, 5)}
    met = {l: lp_metrics(i) for l, i in cur.items()}
    vg0_ideal = met[0][0]                            # 이상 소자 WF 정렬
    print('\n Lrcs  SS[mV/dec]  Ion_LP[uA]  Ioff(WF고정)[A]  /이상')
    ioffs = {}
    for l in (0, 3, 5):
        lg10 = np.log10(np.maximum(cur[l], 1e-30))
        ioffs[l] = 10**np.interp(vg0_ideal, VGS, lg10)
        print(f'  {l}    {met[l][2]:6.1f}     {met[l][1]*1e6:8.2f}'
              f'    {ioffs[l]:.3e}   {ioffs[l]/ioffs[0]:6.1f}x')
    r5 = ioffs[5]/ioffs[0]
    ss = [met[l][2] for l in (0, 3, 5)]
    figure(cur, met, ioffs, vg0_ideal)
    # 논문 명시: Lrcs=5 에서 Ioff 10배 이상, SS/on-off 열화, Ion 증가
    assert r5 > 10.0, f'Ioff(5nm)/Ioff(0) = {r5:.1f}'
    assert ss[0] < ss[1] < ss[2]
    assert met[0][1] < met[3][1] < met[5][1]
    print(f'check passed: Wang et al. (2023) 명시 결과 재현 — '
          f'Ioff(Lrcs=5) = {r5:.1f}x (>10x), SS 열화 {ss[0]:.1f}->{ss[2]:.1f}.'
          f' (Ion 절대값은 vsat/ballistic 미구현으로 비교 제외)')


def _regmap(xx, yy, zz, lrcs=0.0):
    """가시화용: 거울 대칭 펼친 좌표 -> 0빈/1금속/2산화막/3채널/4 n+."""
    code = {'EOX': 2, 'NCH': 3, 'NSD': 4}
    zm = TNS + TSP + TNS/2                   # z 거울면 21
    out = np.zeros(np.broadcast(xx, yy, zz).shape)
    it = np.nditer([xx, yy, zz, out],
                   op_flags=[['readonly']]*3 + [['writeonly']])
    for x, y, z, o in it:
        ze = float(z) if z <= zm else 2*zm - float(z)
        r = region(float(x), abs(float(y)), ze, lrcs)
        if r is not None:
            o[...] = code[r]
        else:
            dx = min(float(x), XT - float(x))
            ingate = dx >= LSD + LSP
            inbox = (abs(float(y)) < WNS2 + EOT + 1e-6
                     and -EOT - 1e-6 < ze)
            o[...] = 1 if (ingate and inbox) else 0
    return out


def figure(cur, met, ioffs, vg0_ideal):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(['white', '#909090', '#9ecae1',
                           '#a1d99b', '#c1392b'])
    fig = plt.figure(figsize=(14.5, 7.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1])
    # (0) 3D 전체 구조: 1nm 복셀, y=0 컷어웨이 + 반투명 게이트 금속
    ax3 = fig.add_subplot(gs[:, 0], projection='3d')
    zm2 = 2*(TNS + TSP + TNS/2)              # 42
    xe = np.arange(0., XT + 0.5, 1.0)
    ye = np.arange(0., WNS2 + EOT + 0.5, 1.0)    # 컷어웨이: y >= 0 만
    ze = np.arange(-EOT, zm2 + EOT + 0.5, 1.0)
    xc = 0.5*(xe[:-1] + xe[1:])
    yc = 0.5*(ye[:-1] + ye[1:])
    zc = 0.5*(ze[:-1] + ze[1:])
    XC, YC, ZC = np.meshgrid(xc, yc, zc, indexing='ij')
    RM = _regmap(XC, YC, ZC)
    filled = RM > 0
    col = np.zeros(RM.shape + (4,))
    col[RM == 1] = (0.42, 0.42, 0.45, 0.60)  # 금속 (약반투명 진회색)
    col[RM == 2] = (0.62, 0.79, 0.88, 0.45)  # 산화막
    col[RM == 3] = (0.42, 0.72, 0.40, 1.0)   # 채널
    col[RM == 4] = (0.76, 0.22, 0.17, 1.0)   # S/D 에피
    X3, Y3, Z3 = np.meshgrid(xe, ye, -ze, indexing='ij')
    ax3.voxels(X3, Y3, Z3, filled, facecolors=col, edgecolor=None)
    ax3.set_box_aspect((XT, 2*(WNS2+EOT), zm2 + 2*EOT))
    ax3.set_xlabel('x [nm]')
    ax3.set_ylabel('y [nm]')
    ax3.set_zlabel('-z [nm]')
    ax3.set_title('3-stacked NS GAAFET (cutaway at y=0)')
    ax3.view_init(16, -78)
    # (a) x-z 단면 (y=0, 시트 중앙): 적층 시트 + S/D 에피 + 게이트/스페이서
    ax = fig.add_subplot(gs[0, 1])
    x = np.linspace(0, XT, 241)
    z = np.linspace(-EOT, 2*(TNS+TSP+TNS/2)+EOT, 221)
    X, Z = np.meshgrid(x, z)
    ax.pcolormesh(X, -Z, _regmap(X, np.zeros_like(X), Z), cmap=cmap,
                  vmin=0, vmax=4, shading='nearest')
    ax.set_title('cut A-A (y=0): 3-stacked nanosheets')
    ax.set_xlabel('x [nm] (source $\\to$ drain)')
    ax.set_ylabel('-z [nm]')
    # (b) y-z 단면 (게이트 중앙): gate-all-around 랩
    ax = fig.add_subplot(gs[0, 2])
    y = np.linspace(-(WNS2+EOT), WNS2+EOT, 241)
    Y, Z2 = np.meshgrid(y, z)
    ax.pcolormesh(Y, -Z2, _regmap(np.full_like(Y, X0), Y, Z2), cmap=cmap,
                  vmin=0, vmax=4, shading='nearest')
    ax.set_title('cut B-B (gate center): gate-all-around')
    ax.set_xlabel('y [nm]')
    ax.set_ylabel('-z [nm]')
    for lab, c in (('metal gate', '#909090'), ('oxide', '#9ecae1'),
                   ('channel 10$^{15}$', '#a1d99b'),
                   ('S/D 2$\\times$10$^{20}$', '#c1392b')):
        ax.plot([], [], 's', color=c, label=lab)
    ax.legend(fontsize=7, loc='center')
    # (c) 전달특성 (WF = 이상 소자 정렬로 고정)
    ax = fig.add_subplot(gs[1, 1])
    for l, c in ((0, 'tab:blue'), (3, 'tab:orange'), (5, 'crimson')):
        ax.semilogy(np.array(VGS) - vg0_ideal, cur[l], 'o-', ms=3,
                    color=c, label=f'$L_{{rcs}}$ = {l} nm')
    ax.axvline(0, color='k', lw=0.7, ls=':')
    ax.axhline(IOFF_LP, color='k', lw=0.7, ls=':')
    ax.annotate(f'$I_{{off}}$ @ fixed WF:\n1.0x / {ioffs[3]/ioffs[0]:.1f}x /'
                f' {ioffs[5]/ioffs[0]:.1f}x\n(paper: >10x at 5 nm)',
                xy=(0, ioffs[5]), xytext=(0.25, 1e-14), fontsize=8,
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.set_xlabel("$V_G'$ [V] (LP-aligned)")
    ax.set_ylabel('$I_D$ [A]')
    ax.set_title(f'transfer @ $V_{{DS}}$ = {VDD} V')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which='both')
    # (d) SS 열화
    ax = fig.add_subplot(gs[1, 2])
    ls = [0, 3, 5]
    ax.bar([str(v) for v in ls], [met[v][2] for v in ls],
           color=['tab:blue', 'tab:orange', 'crimson'], width=0.5)
    ax.axhline(64.0, color='k', ls='--', lw=1,
               label='paper calibration baseline ~64')
    for i, v in enumerate(ls):
        ax.text(i, met[v][2] + 0.3, f'{met[v][2]:.1f}', ha='center',
                fontsize=9)
    ax.set_ylim(55, 75)
    ax.set_xlabel('$L_{rcs}$ [nm]')
    ax.set_ylabel('SS [mV/dec]')
    ax.set_title('subthreshold swing degradation')
    ax.legend(fontsize=8)
    fig.suptitle('Stacked-nanosheet GAAFET rebuilt from Wang et al., '
                 'Electronics 12, 770 (2023), Table 1', y=0.99)
    fig.tight_layout()
    out = os.path.join(HERE, '..', 'docs', 'fig_nsfet')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    print('figure ->', os.path.abspath(out) + '.png')


def main():
    os.makedirs(RUN, exist_ok=True)
    arg = sys.argv[1] if len(sys.argv) > 1 else '0'
    if arg == 'all':
        compare()
        return
    lrcs = float(arg)
    i_d = run_job(f'ex25r{int(lrcs)}', lrcs)
    print(f'\nLrcs = {lrcs} nm (VDS = {VDD} V):')
    for vg, i in zip(VGS, i_d):
        print(f'  VGraw={vg:5.2f}  I_D={i:.3e} A')
    vg0, ion, ss = lp_metrics(i_d)
    print(f'LP 정렬: VG0(raw@0.1nA) = {vg0:.3f} V -> '
          f'Ion(VG\'=0.7) = {ion*1e6:.2f} uA, SS = {ss:.1f} mV/dec')


if __name__ == '__main__':
    main()
