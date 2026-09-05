# -*- coding: utf-8 -*-
"""paper_wang2023: D. Wang et al., Electronics 12, 770 (2023), MDPI 오픈액세스
"Investigation of Source/Drain Recess Engineering and Its Impacts on
FinFET and GAA Nanosheet FET at 5 nm Node" (CC BY).

복잡 형상 재현 대상 2호: 5 nm 노드 3단 적층 나노시트 GAAFET vs 2-fin FinFET.
Table 1 (본문 그대로):
"""

COMMON = dict(CPP=51., LG=18., LSP0=5., LCNT=14.,            # nm
              NSD=2e20, NCH=1e15, NPTS=2e18)                 # cm^-3
FIN = dict(NFIN=2, WFIN=6., HFIN=56., PFIN=30.)              # nm
NS = dict(NNS=3, WNS=36., TNS=6., TSP=12.)                   # nm
VDD = 0.7                                                    # V

# 본문 명시 결과값 (재현 비교 목표):
# - 유효폭: Weff(Fin)=236 nm, Weff(NS)=252 nm  (식 (1),(2))
# - 보정 기준(TSMC 7nm FinFET): SS ~64 mV/dec, DIBL ~30 mV/V
# - LP 설계: WF 튜닝으로 Ioff = 0.1 nA 고정 (WF는 파생량 -> 우리도 동일 절차)
# - Lrcs=5 nm 에서 Ioff 가 이상 구조 대비 10배 이상 증가
# - NSFET Ion 이 FinFET 대비 ~10% 우위 (Weff + (100)면 이동도)
# - Lrcs+Hrcs(2/10 nm) 종합: NSFET Ioff 가 FinFET 대비 37% 작음
#
# Fig. 3(b) NSFET 이상(Lrcs=0) 전달곡선 디지타이즈 (CC BY, 캔버스 픽셀 판독).
# Vds=0.7V. (0,1e-10)은 본문 LP 스펙(외삽과 일치). 0.15-0.45V 구간은
# 레전드와 겹쳐 제외, 0.8V 점은 프레임 경계 의심으로 제외.
# 서브문턱 fit -> 논문 NSFET SS ~ 69.6 mV/dec (DG 양자보정 포함 값).
FIG3B_IDEAL = [
    (0.000, 1.00e-10), (0.025, 1.74e-10), (0.050, 4.01e-10),
    (0.075, 9.23e-10), (0.100, 2.02e-09), (0.125, 4.72e-09),
    (0.500, 1.41e-05), (0.550, 1.93e-05), (0.600, 2.46e-05),
    (0.650, 2.99e-05), (0.700, 3.63e-05), (0.750, 4.21e-05),
]

# 미명시 -> 가정 필요(재현 시 명시할 것):
# - EOT/tox (5nm 노드 관례 ~1 nm 가정), 게이트 WF (Ioff 타깃팅으로 결정)
# - 이동도/양자보정 모델 스택 (Sentaurus IAL+DG+ballistic+vsat, TSMC 7nm 보정)
#   -> 고전 상수-mu DD 로는 Ion 절대값 재현 불가. SS/전기정역학 비교가 목표.
