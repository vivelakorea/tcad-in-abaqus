# -*- coding: utf-8 -*-
"""run_staggered: monolithic vs staggered형(분할반복) 수렴 비교 (검증계획 3).

같은 잔차(=같은 물리·같은 해)를 두 Jacobian 으로 푼다:
- monolithic: uel_mos_et.f (전기<->열 교차 블록 포함, 일관 Newton)
- staggered형: uel_mos_et_stag.f (교차 블록 제거 -> block-Jacobi 분할반복;
  문헌의 TCAD<->열 FE staggered 루프가 매 반복 하는 일과 동일)

결합 세기(HSCALE)를 올리며 반복수/증분수/cutback/완주 여부를 비교.
기대: 약결합에선 비슷, 강결합에선 staggered 반복수 급증 또는 발산 --
monolithic 을 쓰는 이유의 정량 증거. 두 모드가 완주하면 해(I_D)는
동일해야 함(잔차 동일) -> 상호 검증.

사용: python run_staggered.py   (Abaqus job 6개, ~20 min)
"""
import os
import re
import subprocess

from run_mosfet import Q, NI
from run_selfheating import write_inp, parse, VDS, VG, RUN, HERE

HSCALES = [5.0, 50.0, 200.0, 800.0]
MODES = (('mono', 'uel_mos_et.f'), ('stag', 'uel_mos_et_stag.f'))


def run_case(job, hscale, uel):
    xs, zs, nid = write_inp(job, hscale)
    nn = max(nid.values())
    msgp = os.path.join(RUN, job + '.msg')
    if os.path.exists(msgp) and os.path.exists(os.path.join(RUN, job + '.dat')):
        print(f'> {job}: cached')                    # 배치 실행용 (STAG_HS 참조)
        ok = 'THE ANALYSIS HAS BEEN COMPLETED' in open(msgp,
                                                       errors='ignore').read()
    else:
        lck = os.path.join(RUN, job + '.lck')
        if os.path.exists(lck):
            os.remove(lck)
        cmd = (f'abaqus job={job} user={os.path.join(HERE, uel)}'
               f' interactive cpus=1')
        print('>', cmd)
        r = subprocess.run(f'cmd /c "{cmd}"', cwd=RUN, capture_output=True,
                           text=True, errors='replace')
        ok = 'COMPLETED' in (r.stdout or '')
    msg = open(os.path.join(RUN, job + '.msg'), errors='ignore').read()
    incs = int(re.search(r'TOTAL OF\s+(\d+)\s+INCREMENTS', msg).group(1))
    iters = int(re.search(r'(\d+)\s+ITERATIONS INCLUDING', msg).group(1))
    cuts = int(re.search(r'(\d+)\s+CUTBACKS', msg).group(1))
    steps = len(re.findall(r'FRACTION OF STEP COMPLETED\s+1\.00', msg))
    i_d = None
    if ok:
        ev = parse(job, nn)
        i_d = abs(Q*NI*ev[-1][1][0][1])              # 마지막 VD 스텝 드레인 전류
    return dict(ok=ok, incs=incs, iters=iters, cuts=cuts, steps=steps,
                i_d=i_d)


def main():
    os.makedirs(RUN, exist_ok=True)
    hss = HSCALES
    if os.environ.get('STAG_HS'):                    # 배치 실행: 부분만 돌리기
        hss = [float(v) for v in os.environ['STAG_HS'].split(',')]
        for hs in hss:
            for mode, uel in MODES:
                run_case(f'ex22{mode}{int(hs)}', hs, uel)
        print('batch done:', hss)
        return
    res = {}
    for hs in HSCALES:
        for mode, uel in MODES:
            job = f'ex22{mode}{int(hs)}'
            res[(hs, mode)] = run_case(job, hs, uel)

    nstep = 2 + len(VDS)
    print(f'\n VG={VG:g}, VD 사다리 {VDS} ({nstep} steps)')
    print(' HSCALE  mode  완주  steps  incs  iters  cutbacks  I_D(VD=3)')
    for hs in HSCALES:
        for mode, _ in MODES:
            r = res[(hs, mode)]
            i_s = f'{r["i_d"]:.4e}' if r['i_d'] else '   -  (발산)'
            print(f'{hs:6.0f}  {mode}  {"O" if r["ok"] else "X":>3}  '
                  f'{r["steps"]:5d}  {r["incs"]:4d}  {r["iters"]:5d}  '
                  f'{r["cuts"]:8d}  {i_s}')

    # monolithic 은 전 구간 완주
    for hs in HSCALES:
        assert res[(hs, 'mono')]['ok'], f'monolithic 발산 @HSCALE={hs}'
    # 두 모드가 다 완주한 경우 같은 해로 수렴 (잔차 동일 -> 상호 검증)
    both = [hs for hs in HSCALES if res[(hs, 'stag')]['ok']]
    for hs in both:
        im, s = res[(hs, 'mono')]['i_d'], res[(hs, 'stag')]['i_d']
        dev = abs(im - s)/im
        print(f'HSCALE={hs:g}: I_D(mono/stag) 상대차 {dev:.2e}')
        assert dev < 1e-3
    # 강결합에서 staggered 비용 증가 또는 발산
    top = HSCALES[-1]
    if res[(top, 'stag')]['ok']:
        assert (res[(top, 'stag')]['iters']
                > res[(top, 'mono')]['iters']), '강결합 staggered 우위?'
        print(f'강결합(HSCALE={top:g}): staggered 반복수 '
              f'{res[(top, "stag")]["iters"]} > monolithic '
              f'{res[(top, "mono")]["iters"]}')
    else:
        print(f'강결합(HSCALE={top:g}): staggered 발산, monolithic 완주.')
    print('check passed: staggered vs monolithic 비교 (검증계획 3).')


if __name__ == '__main__':
    main()
