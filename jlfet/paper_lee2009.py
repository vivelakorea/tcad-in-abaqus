# -*- coding: utf-8 -*-
"""paper_lee2009: C.-W. Lee et al., Appl. Phys. Lett. 94, 053511 (2009)
"Junctionless multigate field-effect transistor" — 재현 대상 논문 데이터.

Table I (본문 그대로, junctionless 열):
  채널 도핑 8e19 cm^-3 (N형, S/D 와 동일 -> 접합 없음)
  게이트 산화막 2 nm, 게이트 일함수 5.5 eV
  T_si = W_fin = 5 nm (단면 5x5 nm^2), L_gate = 10-30 nm, 트라이게이트
  시뮬레이터: Atlas 3D (고전 드리프트-확산; 이동도 모델 미명시)

Fig. 3 (V_DS 50 mV) 디지타이즈 값 — junctionless 곡선:
  원본 그림 292x211 px 해상도 한계로 SS +-1.5 mV/dec, DIBL +-4 mV.
  DIBL = Vth(0.05V) - Vth(1V).
"""

# L_gate [nm] -> (SS [mV/dec], DIBL [mV])
FIG3_JL = {
    10: (64.0, 73.0),
    15: (61.0, 39.0),
    20: (60.0, 14.0),
    30: (60.0, 8.0),
}
SS_TOL, DIBL_TOL = 1.5, 4.0                  # 디지타이즈 불확도

# Table I
ND = 8.0e19                                  # cm^-3
TOX = 2.0                                    # nm
WF_GATE = 5.5                                # eV
TSI = WSI = 5.0                              # nm
PHI_I_SI = 4.05 + 0.56                       # Si 진성 준위 일함수 [eV]
