# -*- coding: utf-8 -*-
"""reference_dd1d: 1D pn 접합 드리프트-확산 순수 파이썬 참조해 (SG + Gummel vs coupled Newton).

무차원화(V_T=1, n_i=1, q=eps=1). 미지수: 절점마다 (psi, n, p).
    Poisson:  -psi'' = p - n + N(x)
    연속:     J_n' = 0,  J_p' = 0   (재결합 무시, 정상상태)
SG 플럭스:  J_n = [B(t) n_R - B(-t) n_L]/h,  t = psi_R - psi_L
검증 3종:
  (1) 평형(V=0): built-in potential = ln(Na*Nd) 해석값과 비교, np=1 질량법칙
  (2) 순바이어스에서 J_n + J_p 가 소자 전체에서 상수인가 (전류 보존)
  (3) Gummel(staggered) vs coupled Newton(monolithic, 3x3 블록) 반복 수 비교
"""
import numpy as np

Nd = Na = 1.0e2
nnode, L = 201, 4.0
x = np.linspace(-L/2, L/2, nnode)
h = x[1] - x[0]
Ndop = np.where(x > 0, Nd, -Na)


def bern(t):
    """Bernoulli B(t)=t/(e^t-1), 안정 구현."""
    t = np.asarray(t, float)
    out = np.empty_like(t)
    t = np.clip(t, -500, 500)
    small = np.abs(t) < 1e-10
    out[small] = 1.0 - t[small]/2
    ts = t[~small]
    out[~small] = ts / np.expm1(ts)
    return out


def dbern(t):
    """B'(t)."""
    t = np.asarray(t, float)
    out = np.empty_like(t)
    small = np.abs(t) < 1e-6
    out[small] = -0.5 + t[small]/6
    ts = t[~small]
    e = np.expm1(ts)
    out[~small] = (e - ts*np.exp(ts)) / e**2
    return out


def contacts(V):
    """옴 접촉: np=1, n-p=N, psi = 인가전압 + ln(n)."""
    nL = (-Na/2) + np.sqrt((Na/2)**2 + 1)   # p쪽: n 작음
    pL = 1.0/nL
    nR = (Nd/2) + np.sqrt((Nd/2)**2 + 1)
    pR = 1.0/nR
    psiL = np.log(nL) + V                    # p쪽 접촉에 +V (순바이어스)
    psiR = np.log(nR)
    return (psiL, nL, pL), (psiR, nR, pR)


def residual(d, V):
    """d = [psi(0..), n(..), p(..)] 순서로 블록 배열."""
    psi, n, p = d[:nnode], d[nnode:2*nnode], d[2*nnode:]
    R = np.zeros(3*nnode)
    t = np.diff(psi)                          # 요소별 psi_R - psi_L
    # SG 플럭스 (요소별)
    Jn = (bern(t)*n[1:] - bern(-t)*n[:-1]) / h
    Jp = (bern(-t)*p[1:] - bern(t)*p[:-1]) / h
    # Poisson (FD와 동일한 1D FEM): -psi'' = p-n+N
    R[1:nnode-1] = (-(psi[2:] - 2*psi[1:-1] + psi[:-2]) / h
                    - h*(p[1:-1] - n[1:-1] + Ndop[1:-1]))
    # 연속: J' = 0  (절점 잔차 = 유입-유출)
    R[nnode+1:2*nnode-1] = Jn[1:] - Jn[:-1]
    R[2*nnode+1:3*nnode-1] = -(Jp[1:] - Jp[:-1])
    # 접촉 BC
    (psiL, nL, pL), (psiR, nR, pR) = contacts(V)
    R[0] = psi[0] - psiL; R[nnode-1] = psi[nnode-1] - psiR
    R[nnode] = n[0] - nL; R[2*nnode-1] = n[nnode-1] - nR
    R[2*nnode] = p[0] - pL; R[3*nnode-1] = p[nnode-1] - pR
    return R


def jacobian_fd(d, V):
    """검증용이 아니라 교육용 명료함을 위해 해석 Jacobian 대신 희소 FD.
    (열 그룹핑: 3점 스텐실이라 6개 색이면 충분하지만, 여기선 단순 FD로 충분히 빠름)"""
    m = 3*nnode
    J = np.zeros((m, m))
    R0 = residual(d, V)
    for j in range(m):
        dp = d.copy()
        step = 1e-7 * max(1.0, abs(d[j]))
        dp[j] += step
        J[:, j] = (residual(dp, V) - R0) / step
    return J, R0


def newton(d, V, tol=1e-9, maxit=40):
    """감쇠 Newton(Bank-Rose식): 잔차가 줄 때까지 스텝을 반토막."""
    hist = []
    for it in range(maxit):
        J, R = jacobian_fd(d, V)
        nr = np.abs(R).max()
        hist.append(nr)
        if nr < tol:
            return d, hist
        dd = np.linalg.solve(J, -R)
        s = 1.0
        for _ in range(30):
            dt_ = d + s*dd
            dt_[nnode:] = np.maximum(dt_[nnode:], 1e-30)
            if np.abs(residual(dt_, V)).max() < nr:
                break
            s *= 0.5
        d = dt_
    return d, hist


def newton_ramped(d, V, nramp=5):
    """바이어스를 nramp 단계로 램프(4장 하중 스텝). 단계별 반복 수 목록 반환."""
    counts = []
    for k in range(1, nramp+1):
        d, hist = newton(d, V*k/nramp)
        counts.append(len(hist)-1)
    return d, counts


def gummel(d, V, tol=1e-9, maxit=400):
    """staggered: psi(비선형 Poisson, 준페르미 고정) -> n -> p 순환."""
    psi, n, p = d[:nnode].copy(), d[nnode:2*nnode].copy(), d[2*nnode:].copy()
    (psiL, nL, pL), (psiR, nR, pR) = contacts(V)
    phin = np.log(n) - psi                    # 준페르미(스케일) 고정용
    phip = np.log(p) + psi
    for it in range(maxit):
        psi_old = psi.copy()
        # (a) 비선형 Poisson: n=exp(psi+phin), p=exp(-psi+phip)로 대입해 psi만 Newton
        for k in range(60):
            nn_ = np.exp(psi + phin); pp_ = np.exp(-psi + phip)
            Rp = np.zeros(nnode)
            Rp[1:-1] = (-(psi[2:] - 2*psi[1:-1] + psi[:-2]) / h
                        - h*(pp_[1:-1] - nn_[1:-1] + Ndop[1:-1]))
            Rp[0] = psi[0]-psiL; Rp[-1] = psi[-1]-psiR
            if np.abs(Rp).max() < 1e-11:
                break
            Jp_ = np.zeros((nnode, nnode))
            for i in range(1, nnode-1):
                Jp_[i, i-1] = -1/h; Jp_[i, i+1] = -1/h
                Jp_[i, i] = 2/h + h*(pp_[i] + nn_[i])
            Jp_[0, 0] = Jp_[-1, -1] = 1.0
            dpsi = np.linalg.solve(Jp_, -Rp)
            dpsi = np.clip(dpsi, -2, 2)
            psi += dpsi
        # (b) n 연속 (psi 고정, 선형 3중대각)
        A = np.zeros((nnode, nnode)); b = np.zeros(nnode)
        t = np.diff(psi)
        for i in range(1, nnode-1):
            A[i, i-1] = -bern(-t[i-1])/h
            A[i, i]   = (bern(t[i-1]) + bern(-t[i]))/h
            A[i, i+1] = -bern(t[i])/h
        A[0, 0] = A[-1, -1] = 1.0; b[0] = nL; b[-1] = nR
        n = np.linalg.solve(A, b)
        # (c) p 연속
        A = np.zeros((nnode, nnode)); b = np.zeros(nnode)
        for i in range(1, nnode-1):
            A[i, i-1] = -bern(t[i-1])/h
            A[i, i]   = (bern(-t[i-1]) + bern(t[i]))/h
            A[i, i+1] = -bern(-t[i])/h
        A[0, 0] = A[-1, -1] = 1.0; b[0] = pL; b[-1] = pR
        p = np.linalg.solve(A, b)
        n = np.maximum(n, 1e-30); p = np.maximum(p, 1e-30)
        phin = np.log(n) - psi; phip = np.log(p) + psi
        d = np.r_[psi, n, p]
        if np.abs(psi - psi_old).max() < tol:
            return d, it + 1
    return d, maxit


def init_guess(V=0.0):
    n0 = np.where(x > 0, Nd, 1.0/Na)
    p0 = 1.0/n0
    psi0 = np.log(n0)
    return np.r_[psi0, n0, p0]


def demo():
    # ---- (1) 평형: Vbi 해석값 ----
    d0, it_g = gummel(init_guess(), 0.0)
    psi = d0[:nnode]; n = d0[nnode:2*nnode]; p = d0[2*nnode:]
    Vbi = psi[-1] - psi[0]
    Vbi_exact = np.log(Na*Nd)
    print(f"equilibrium: Vbi = {Vbi:.6f}  vs  ln(Na*Nd) = {Vbi_exact:.6f}")
    mass = np.abs(n*p - 1).max()
    print(f"mass law  max|np-1| = {mass:.2e}")

    # ---- (2)+(3) 순바이어스 V=10 V_T: 전류 보존 + 반복 수 비교 ----
    V = 10.0
    d_g, itg = gummel(d0.copy(), V)
    d_n, hist = newton(d0.copy(), V)
    t = np.diff(d_n[:nnode])
    n_ = d_n[nnode:2*nnode]; p_ = d_n[2*nnode:]
    Jn = (bern(t)*n_[1:] - bern(-t)*n_[:-1])/h
    Jp = (bern(-t)*p_[1:] - bern(t)*p_[:-1])/h
    Jtot = Jn - Jp
    dev = (Jtot.max()-Jtot.min())/np.abs(Jtot).max()
    print(f"forward bias V=10VT: total current constancy dev = {dev:.2e}")
    print(f"Gummel iterations   : {itg}")
    print(f"coupled Newton |R|  : " + "  ".join(f"{v:.1e}" for v in hist))
    dpsi = np.abs(d_g[:nnode]-d_n[:nnode]).max()
    print(f"|psi_gummel - psi_newton|max = {dpsi:.2e}  (같은 해)")
    assert abs(Vbi-Vbi_exact) < 1e-3 and mass < 1e-6 and dev < 1e-8
    print("check passed.")


if __name__ == "__main__":
    demo()
