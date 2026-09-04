# -*- coding: utf-8 -*-
"""run_resistor: 전기-열-기계 UEL(uel_mos_etm.f)의 요소 수준 정량 검증.

균일 n+ 도핑 막대 저항기 (4 x 0.5 x 0.5 um, N_D = 1e19 cm^-3):
1. 저항: R = L/(q n mu_n A) 표준 반도체 저항 공식 (드리프트-확산의
   옴 극한, 균일 막대에서 정확).
2. 자기발열 온도 프로파일: 균일 Joule + 양끝 T 고정 -> 포물선
   dT(x) = (P/(2 kappa A)) x (L-x)/L,  dT_max = P L / (8 kappa A)
   (Carslaw & Jaeger 1959, 정상상태 열전도 닫힌형).
3. 에너지 수지: 양끝 히트싱크 반력열 합 = I*V.

MOSFET 검증(run_etm)과 달리 여기는 전 항목이 닫힌형 -> 요소 단위 정답 대조.
사용: python run_resistor.py   (Abaqus job 1개, ~1 min)
"""
import io
import os
import subprocess

import numpy as np

from run_mosfet import Q, VT, NI, UM
from run_etm import parse

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, 'abq_run')
ND = 1.0e19                                  # 도핑 [cm^-3]
XMUN = 400.0                                 # UEL 상수와 동일
XKSI = 1.5                                   # W/cm/K
LX, WY, WZ = 4.0, 0.5, 0.5                   # um
VAPP = 0.1                                   # V (옴 영역)


def write_inp(job):
    # 접촉 경계층(Debye ~1.3nm @1e19) 해상: 양끝 nm급 세밀화 (MOSFET과 동일 관례)
    seg = np.linspace
    xs = np.unique(np.round(np.r_[seg(0, 0.02, 11), seg(0.02, 0.2, 10),
                                  seg(0.2, LX - 0.2, 40),
                                  seg(LX - 0.2, LX - 0.02, 10),
                                  seg(LX - 0.02, LX, 11)], 7))
    nx = len(xs) - 1
    nid = {}
    lines = []
    a = 0
    for iz in (0, 1):
        for iy in (0, 1):
            for ix in range(nx + 1):
                a += 1
                nid[(ix, iy, iz)] = a
                lines.append(f'{a}, {xs[ix]*UM:.10e}, {iy*WY*UM:.10e}, '
                             f'{iz*WZ*UM:.10e}')
    conn = []
    for ix in range(nx):
        nn = [nid[(ix, 0, 0)], nid[(ix+1, 0, 0)], nid[(ix+1, 1, 0)],
              nid[(ix, 1, 0)], nid[(ix, 0, 1)], nid[(ix+1, 0, 1)],
              nid[(ix+1, 1, 1)], nid[(ix, 1, 1)]]
        conn.append(f'{ix+1}, ' + ', '.join(map(str, nn)))
    psn = np.arcsinh(ND/NI/2)
    L = ['*HEADING', 'ex20 uniform n+ resistor bar (ETM UEL element test)',
         '*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES=3, COORDINATES=3,'
         ' VARIABLES=1, UNSYMM', '1,2,3,4,5,6,7', '*NODE'] + lines
    L += ['*ELEMENT, TYPE=U1, ELSET=EBAR'] + conn
    L += ['*UEL PROPERTY, ELSET=EBAR', f'1, {ND/NI:.6e}, 1.0']
    xl = sorted(v for (ix, iy, iz), v in nid.items() if ix == 0)
    xr = sorted(v for (ix, iy, iz), v in nid.items() if ix == nx)
    for nm, ids in (('XL', xl), ('XR', xr)):
        L.append(f'*NSET, NSET={nm}')
        L.append(', '.join(map(str, ids)))
    L.append('*NSET, NSET=NALL, GENERATE')
    L.append(f'1, {a}, 1')

    # 접촉 psi BC 를 UEL 도핑 램프와 일관되게: psi_face(t) = asinh(t*N/2ni)
    # (선형 램프 BC 는 log 로 크는 벌크 중성 psi 와 어긋나 콜드스타트 발산)
    # log 간격 테이블: 어느 t 에서도 보간값 ~ asinh(t*N/2ni) (선형 구간 불일치 제거)
    L.append('*AMPLITUDE, NAME=DOPRAMP, TIME=TOTAL TIME')
    tpts = [0.] + list(np.logspace(-4, 0, 25))
    for t in tpts:
        L.append(f'{t:.6e}, {np.arcsinh(t*ND/NI/2)/psn:.8e}')

    def bstep(v, ninc, amp):
        dt0 = 1e-4 if amp else 1.0/ninc
        # 평균력 스케일 ~4 로 작아 1e-6 비율은 잔차 노이즈 바닥 아래 (함정 7 계열)
        s = ['*STEP, INC=400, UNSYMM=YES', '*STATIC',
             f'{dt0}, 1.0, 1e-9, {1.0/ninc}',
             '*CONTROLS, PARAMETERS=FIELD', '1e-4,,,,,,,']
        if amp:
            s += ['*BOUNDARY, AMPLITUDE=DOPRAMP',
                  f'XL, 1, 1, {psn:.8e}', f'XR, 1, 1, {psn:.8e}']
        else:
            s += ['*BOUNDARY',
                  f'XL, 1, 1, {psn:.8e}', f'XR, 1, 1, {psn + v/VT:.8e}']
        s += ['*BOUNDARY', 'XL, 2, 2, 0.', f'XR, 2, 2, {v/VT:.8e}',
              'NALL, 3, 3, 0.',      # n+ 막대: 정공 죽은 장 -> 전역 고정
              'XL, 4, 4, 0.', 'XR, 4, 4, 0.',
              'XL, 5, 7, 0.',
              '*END STEP']
        return s

    s1 = bstep(0.0, 20, True)
    s1 = s1[:-1] + ['*NODE PRINT, NSET=NALL, FREQUENCY=999', 'U',
                    '*NODE PRINT, NSET=XR, FREQUENCY=999, TOTALS=YES', 'RF',
                    '*NODE PRINT, NSET=XL, FREQUENCY=999, TOTALS=YES', 'RF',
                    '*END STEP']
    L += s1                                          # 1: 평형(도핑 램프)
    L += bstep(VAPP, 8, False)                       # 2: V 인가 (HS=0 게이트)
    L += bstep(VAPP, 4, False)                       # 3: 발열 on (KSTEP>=3)
    io.open(os.path.join(RUN, job + '.inp'), 'w').write('\n'.join(L) + '\n')
    return xs, nid


def main():
    os.makedirs(RUN, exist_ok=True)
    job = 'ex20bar'
    xs, nid = write_inp(job)
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
    assert len(ev) == 3, f'{len(ev)} print blocks'

    A = WY*WZ*UM*UM                                  # cm^2
    Lc = LX*UM                                       # cm
    d, rf = ev[2]
    i_uel = abs(Q*NI*rf[0][1])                       # 전자 반력 -> A
    r_uel = VAPP/i_uel
    r_th = Lc/(Q*ND*XMUN*A)                          # R = L/(q n mu A)
    err_r = abs(r_uel/r_th - 1)
    print(f'[저항] R_UEL = {r_uel:.2f} Ohm, R = L/(q n mu A) = {r_th:.2f} '
          f'Ohm, 오차 {err_r*100:.2f}%')

    p_in = i_uel*VAPP
    dt_max_th = p_in*Lc/(8*XKSI*A)                   # Carslaw-Jaeger 포물선
    nx = len(xs) - 1
    prof = np.array([d[nid[(ix, 0, 0)]-1, 3] for ix in range(nx+1)])
    dt_mid = np.interp(LX/2, xs, prof)
    dt_q = np.interp(LX/4, xs, prof)
    th = 4*dt_max_th*(xs/LX)*(1 - xs/LX)             # dT(x) 닫힌형
    err_t = abs(dt_mid/dt_max_th - 1)
    err_shape = np.abs(prof - th).max()/dt_max_th    # 전 절점 프로파일 오차
    print(f'[포물선] dT_max = {dt_mid:.4f} K, 닫힌형 PL/8kA = '
          f'{dt_max_th:.4f} K, 오차 {err_t*100:.2f}%; '
          f'dT(L/4)/dT(L/2) = {dt_q/dt_mid:.4f} (이론 0.75), '
          f'프로파일 최대오차 {err_shape*100:.2f}%')

    q_out = abs(rf[0][3]) + abs(rf[1][3])            # 양끝 RM1
    err_e = abs(q_out - p_in)/p_in
    print(f'[에너지] 히트싱크 반력열 {q_out:.3e} W vs I*V = {p_in:.3e} W, '
          f'오차 {err_e*100:.2f}%')

    assert err_r < 0.02
    assert err_t < 0.03 and err_shape < 0.03
    assert err_e < 0.02
    print('check passed: 저항기 요소검증 — R 공식/포물선 dT/에너지수지.')


if __name__ == '__main__':
    main()
