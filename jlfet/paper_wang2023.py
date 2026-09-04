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
# 미명시 -> 가정 필요(재현 시 명시할 것):
# - EOT/tox (5nm 노드 관례 ~1 nm 가정), 게이트 WF (Ioff 타깃팅으로 결정)
# - 이동도/양자보정 모델 스택 (Sentaurus IAL+DG+ballistic+vsat, TSMC 7nm 보정)
#   -> 고전 상수-mu DD 로는 Ion 절대값 재현 불가. SS/전기정역학 비교가 목표.
