# -*- coding: utf-8 -*-
"""reference_korhonen: Korhonen et al., J. Appl. Phys. 73 (1993) 3790.

양끝 플럭스 차단 유한 배선의 EM 응력 전개 — 논문이 그림으로 발표한 해:
  dsig/dt = kappa * d2sig/dx2,  플럭스 kappa*(dsig/dx - G) = 0 at x=0,L
  sig(x,0) = 0,  G = e Z* rho j / Omega
정상해 G*(x-L/2) 를 빼면 순수 확산 + 제로플럭스 -> 코사인 고유전개:
  sig(x,t) = G(x-L/2) + sum_{n odd} (4GL/(n pi)^2) cos(n pi x/L)
             * exp(-(n pi)^2 kappa t / L^2)
t=0 에서 0, t->inf 에서 Blech-Herring 선형 프로파일로 수렴.
(UEL 부호 관례: 캐소드(x=L) 인장 양수 — 위 식 그대로.)
"""
import numpy as np


def sigma(x, t, L, kappa, G, nmax=399):
    """논문 급수해 sig(x,t). x: array [cm], t: scalar [s], G: [MPa/cm]."""
    x = np.asarray(x, dtype=float)
    s = G*(x - L/2)
    n = np.arange(1, nmax + 1, 2, dtype=float)       # 홀수 항
    lam = (n*np.pi/L)**2 * kappa
    coef = 4*G*L/(n*np.pi)**2 * np.exp(-lam*t)
    s = s + (coef[None, :] * np.cos(np.outer(x, n)*np.pi/L)).sum(axis=1)
    return s


if __name__ == '__main__':
    # 자가 점검: t=0 -> 0, t=inf -> 선형, 초기 성장 sqrt(t) 극한
    L, kap, G = 2e-2, 1.8e-9, 1.891e4
    x = np.linspace(0, L, 81)
    assert np.abs(sigma(x, 0.0, L, kap, G)).max() < G*L*2e-3  # 1/n^2 절단
    assert np.abs(sigma(x, 1e9, L, kap, G) - G*(x - L/2)).max() < 1e-6
    t = 2000.0                                       # kappa*t/L^2 ~ 0.009
    early = G*np.sqrt(4*kap*t/np.pi)
    assert abs(sigma(np.array([L]), t, L, kap, G)[0]/early - 1) < 0.01
    print('check passed: Korhonen(1993) 급수해 자가 점검.')
