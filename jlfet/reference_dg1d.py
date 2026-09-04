# -*- coding: utf-8 -*-
"""reference_dg1d: density-gradient 평형 슬랩의 독립 파이썬 참조해.

Ancona & Tiersten (1987) / Ancona & Iafrate (1989) DG 모델, 평형(phi_n=0):
  Poisson:  psi'' = -(q ni / eps VT) (dop - S^2)      [psi 는 VT 단위, y cm]
  DG:       (2b/VT) S'' = (2 sigma - psi) S,   S = e^sigma = sqrt(n/ni)
  b = gamma hbar^2/(12 q m*),  m* = 0.32 m0
BC: 하드월 sigma(0)=sigma(T)=-10 (양자 구속), psi'=0 (게이트 없음).
풀이: 1D 유한차분 + gamma-continuation 감쇠 Newton (0 -> gamma).
UEL(uel_jl_dg.f) 슬랩 단면과 절점별 대조용.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

Q, VT, NI = 1.602e-19, 0.02585, 1.0e10
EPSI = 11.7*8.854e-14
BN0 = 1.985e-16                              # hbar^2/(12 q m*) [V cm^2]
SIGW = -10.0


def solve(T_cm, dop, gamma, ny=201):
    """절점 (psi_i, sigma_i) Newton. 반환: y, psi, sigma."""
    y = np.linspace(0, T_cm, ny)
    h = y[1] - y[0]
    cp = Q*NI/(EPSI*VT)

    psi = np.full(ny, np.arcsinh(dop/2))
    sig = np.full(ny, 0.5*np.log(dop))

    def newton(gam, psi, sig, wall, itmax=100):
        cb = 2.0*gam*BN0/VT

        def resid(psi, sig):
            S = np.exp(np.clip(sig, -60, 60))
            n = S*S
            rp = np.zeros(ny)
            rs = np.zeros(ny)
            rp[1:-1] = (psi[:-2] - 2*psi[1:-1] + psi[2:])/h**2 \
                + cp*(dop - n[1:-1])
            rp[0] = 2*(psi[1] - psi[0])/h**2 + cp*(dop - n[0])
            rp[-1] = 2*(psi[-2] - psi[-1])/h**2 + cp*(dop - n[-1])
            rs[1:-1] = cb*(S[:-2] - 2*S[1:-1] + S[2:])/h**2 \
                - (2*sig[1:-1] - psi[1:-1])*S[1:-1]
            rs[0] = sig[0] - wall
            rs[-1] = sig[-1] - wall
            return rp, rs, S, n

        import os
        dbg = os.environ.get('DGDBG')
        for it in range(itmax):
            rp, rs, S, n = resid(psi, sig)
            res = max(np.abs(rp).max()/cp/dop, np.abs(rs).max()/max(S.max(), 1))
            if dbg and it % 5 == 0:
                print(f'  it{it} res={res:.3e} sig[1]={sig[1]:.2f}')
            if res < 1e-9:                   # 1e-12 는 반올림 바닥 아래(함정 15)
                return psi, sig, True
            # Jacobian (2x2 블록 삼중대각)
            rows, cols, vals = [], [], []

            def add(r, c, v):
                rows.append(r)
                cols.append(c)
                vals.append(v)
            for i in range(ny):
                ip, isg = 2*i, 2*i + 1
                if i in (0, ny-1):
                    j = 1 if i == 0 else ny-2
                    add(ip, ip, -2/h**2)
                    add(ip, 2*j, 2/h**2)
                    add(ip, isg, -cp*2*n[i])
                    add(isg, isg, 1.0)
                    continue
                add(ip, ip, -2/h**2)
                add(ip, 2*(i-1), 1/h**2)
                add(ip, 2*(i+1), 1/h**2)
                add(ip, isg, -cp*2*n[i])
                add(isg, isg, cb*(-2*S[i])/h**2
                    - (2 + 2*sig[i] - psi[i])*S[i])
                add(isg, 2*(i-1)+1, cb*S[i-1]/h**2)
                add(isg, 2*(i+1)+1, cb*S[i+1]/h**2)
                add(isg, ip, S[i])
            J = sp.csr_matrix((vals, (rows, cols)), shape=(2*ny, 2*ny))
            r = np.empty(2*ny)
            r[0::2] = rp
            r[1::2] = rs
            # 행 평형화: 1/h^2 급 행과 O(1) 행의 스케일 격차로 인한
            # 직접해 오염 방지 (조건수 문제)
            rowmax = np.maximum(np.abs(J).max(axis=1).toarray().ravel(),
                                1e-300)
            Dinv = sp.diags(1.0/rowmax)
            d = spla.spsolve((Dinv @ J).tocsr(), Dinv @ (-r))
            dpsi, dsig = d[0::2], d[1::2]
            # 감쇠 순수 Newton. 단조 백트래킹은 금지: 수렴 상태에서 벽값을
            # 옮기면 어떤 스텝도 잔차를 일시 증가시키므로(지수 2차 오차)
            # 단조 요구가 정상 Newton 경로를 막아버린다.
            damp = min(1.0, 2.0/max(np.abs(dpsi).max(), np.abs(dsig).max(),
                                    1e-30))
            psi = psi + damp*dpsi
            sig = sig + damp*dsig
        return psi, sig, False

    # 벽값 continuation: sigma_wall 을 벌크값 -> SIGW 로 램프 (gamma 는 전량).
    # 지수 비선형이라 스텝 |d sigma| ~ 0.3 이내로 잘게 (한 번에 내리면
    # cb/h^2 * dS 가 이웃 잔차를 폭발시켜 Newton 선형화가 깨짐).
    sigb = 0.5*np.log(dop)
    ok = True
    for wall in np.linspace(sigb, SIGW, 81)[1:]:
        psi, sig, ok = newton(gamma, psi, sig, wall)
        assert ok, f'Newton 발산 @wall={wall:.2f}'
    return y, psi, sig


if __name__ == '__main__':
    # 자가 점검: 전역 중성(무게이트 -> 총전하 0), 하드월, 대칭, 격자수렴
    T = 5e-7
    y, psi, sig = solve(T, 8e9, 3.6)
    n = np.exp(2*sig)
    net = np.trapezoid(8e9 - n, y)/(8e9*T)
    assert abs(net) < 1e-4                   # 전역 중성
    assert n[0]/n.max() < 1e-8               # 하드월
    imax = np.argmax(n)
    assert abs(y[imax]/T - 0.5) < 0.01       # 중앙 최대(대칭)
    y2, _, sig2 = solve(T, 8e9, 3.6, ny=401)  # 2배 세분 (그 이상은 1/h^2
    d = np.abs(np.interp(y, y2, np.exp(2*sig2)) - n).max()/n.max()  # 조건수 벽)
    assert d < 5e-3                          # 격자 수렴 (2차: 오차비 4)
    print(f'check passed: DG 참조해 — n_max/dop = {n.max()/8e9:.3f}, '
          f'중앙 n = {n[len(y)//2]/8e9:.3f} dop, 격자수렴 {d:.1e}')
