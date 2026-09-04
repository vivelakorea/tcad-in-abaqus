# -*- coding: utf-8 -*-
"""reference_mosfet: 긴채널 NMOS I-V 논문 재현 (Pao-Sah 1966 vs Brews 1978).

두 고전 논문을 똑같이 구현해 서로 대조한다:
  [1] Pao & Sah, Solid-State Electronics 9 (1966) 927
      --- 수직 Poisson의 정확 적분 + 이중 적분 드레인 전류 (기준해)
  [2] Brews, Solid-State Electronics 21 (1978) 345
      --- charge-sheet 모델 (닫힌형). 논문의 주장: Pao-Sah와 ~1% 이내 일치.
검증: 같은 소자 파라미터로 두 모델의 I_D(V_G, V_D)를 계산, 강반전에서 수 % 일치하는지.
게이트 방정식의 표면전위는 '비선형 Poisson + Newton'으로 푼다.
"""
import numpy as np

# ---- 소자/물리 상수 (실단위) ----
q, kT = 1.602e-19, 0.02585 * 1.602e-19
VT = 0.02585
eps_si, eps_ox = 11.7 * 8.854e-14, 3.9 * 8.854e-14   # F/cm
ni = 1.0e10                                            # cm^-3
Na = 1.0e17
tox = 10e-7                                            # 10 nm in cm
Cox = eps_ox / tox
mu = 400.0                                             # cm^2/Vs
WL = 1.0                                               # W/L
phiF = VT * np.log(Na / ni)
gam = np.sqrt(2 * eps_si * q * Na) / Cox               # body factor
VFB = 0.0                                              # 이상 MOS


def F_field(psis, Vch):
    """Kingston-Neustadter 정규화 전기장 함수 (Pao-Sah 식 그대로).
    Vch: 채널 준페르미 분리(드레인 방향 위치의 V(y))."""
    a = np.exp(-psis / VT) + psis / VT - 1.0
    b = (ni / Na) ** 2 * (np.exp((psis - Vch) / VT) - np.exp(-Vch / VT) + psis / VT * 0)
    # 전자항: (ni/Na)^2 [exp((psi-V)/VT) - exp(-V/VT)] - psi/VT 항은 홀쪽에 이미 포함
    b = (ni / Na) ** 2 * (np.exp((psis - Vch) / VT) - np.exp(-Vch / VT))
    val = a + b
    return np.sqrt(np.maximum(val, 1e-300))


def psis_of(VG, Vch):
    """게이트 방정식 VG = VFB + psis + gam*sqrt(VT)*F 를 Newton으로 (ex10과 같은 요령)."""
    ps = phiF + Vch  # 초기 추정
    for it in range(100):
        Fv = F_field(ps, Vch)
        g = VFB + ps + gam * np.sqrt(VT) * Fv - VG
        h = 1e-7
        Fp = F_field(ps + h, Vch)
        dg = 1 + gam * np.sqrt(VT) * (Fp - Fv) / h
        dps = -g / dg
        ps += np.clip(dps, -0.2, 0.2)
        if abs(dps) < 1e-12:
            break
    return ps


def id_paosah(VG, VD, nV=600, npts=800):
    """Pao-Sah 이중 적분: I_D = mu W/L q ∫0^VD ∫ n(psi) / F dpsi dV  (표준형)."""
    LD = np.sqrt(eps_si * VT / (q * Na))
    Vs = np.linspace(0, VD, nV)
    inner = np.zeros_like(Vs)
    for k, V in enumerate(Vs):
        ps = psis_of(VG, V)
        # 표면 근처에 지수적으로 몰린 적분: u=(ps-psi)/VT 치환 후 균일 u 격자
        u = np.linspace(0.0, min(40.0, (ps - 1e-6) / VT), npts)
        psi = ps - VT * u
        n_over_F = (ni**2 / Na) * np.exp((psi - V) / VT) / (
            np.sqrt(2) * VT / LD * F_field(psi, V))
        inner[k] = np.trapezoid(n_over_F * VT, u)
    return mu * WL * q * np.trapezoid(inner, Vs)


def id_chargesheet(VG, VD):
    """Brews charge-sheet 닫힌형 (표면전위 판): 논문 (또는 Tsividis 표준형).
    I = mu Cox W/L [ (VG-VFB-psi)dpsi 항 - 2/3 gam (psi^{3/2}) 항 + VT(공핍보정) ]"""
    ps0 = psis_of(VG, 0.0)
    psL = psis_of(VG, VD)
    t1 = (VG - VFB + VT) * (psL - ps0) - 0.5 * (psL**2 - ps0**2)
    t2 = -2.0 / 3.0 * gam * (psL**1.5 - ps0**1.5) + gam * VT * (psL**0.5 - ps0**0.5)
    return mu * Cox * WL * (t1 + t2)


def demo():
    print(f"NMOS: Na=1e17, tox=10nm, VT0(approx) = "
          f"{VFB + 2*phiF + gam*np.sqrt(2*phiF):.3f} V")
    print(" VG     VD     I_PaoSah      I_ChargeSheet   dev")
    worst_strong = 0.0
    for VG in (1.5, 2.0, 3.0):
        for VD in (0.05, 0.5, 1.5):
            i1 = id_paosah(VG, VD)
            i2 = id_chargesheet(VG, VD)
            dev = abs(i1 - i2) / i1
            if VG >= 2.0:
                worst_strong = max(worst_strong, dev)
            print(f"{VG:4.1f}  {VD:5.2f}   {i1:.4e}   {i2:.4e}   {dev*100:5.2f}%")
    print(f"strong-inversion worst deviation = {worst_strong*100:.2f}%")
    print("(문헌의 알려진 경향 그대로: 강반전 수 %, 중반전(VG=1.5)에서 커짐)")
    assert worst_strong < 0.05
    print("check passed: 두 논문의 모델이 수 % 이내로 서로를 재현한다.")


if __name__ == "__main__":
    demo()
