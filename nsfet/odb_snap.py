# -*- coding: utf-8 -*-
"""odb_snap: ODB 를 Abaqus/CAE 로 렌더해 PNG 저장 (noGUI).
사용: abaqus cae noGUI=odb_snap.py -- <odb> <out_png_base>
UVARM1(log10 n) 컨투어, 실리콘 elset 만 표시, 마지막 프레임.
"""
import sys
from abaqus import session
from abaqusConstants import (CONTOURS_ON_UNDEF, INTEGRATION_POINT, PNG, OFF,
                             ON, FEATURE, PARALLEL)
import visualization                          # session.openOdb 표시형 등록
try:
    import displayGroupOdbToolset as dgo
except ImportError:
    dgo = None

odbp = sys.argv[-2]
outp = sys.argv[-1]
o = session.openOdb(name=odbp)
vp = session.viewports['Viewport: 1']
vp.setValues(displayedObject=o)
vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF,))
vp.odbDisplay.setPrimaryVariable(variableLabel='UVARM1',
                                 outputPosition=INTEGRATION_POINT)
nstep = len(o.steps)
try:
    ifr = int(sys.argv[-3])                  # 선택: 마지막 스텝의 프레임 번호
except (ValueError, IndexError):
    ifr = -1
vp.odbDisplay.setFrame(step=nstep - 1, frame=ifr)
# 실리콘(더미 오버레이) 만 표시 (모듈 없으면 전체 표시로 폴백)
if dgo is not None:
    try:
        inst = o.rootAssembly.instances.keys()[0]
        leaf = dgo.LeafFromElementSets(elementSets=(inst + '.VNSD',
                                                    inst + '.VNCH'))
        vp.odbDisplay.displayGroup.replace(leaf=leaf)
    except Exception as e:
        print('displaygroup fallback:', e)
vp.odbDisplay.contourOptions.setValues(minAutoCompute=OFF, minValue=4.0,
                                       maxAutoCompute=OFF, maxValue=21.0)
vp.odbDisplay.commonOptions.setValues(visibleEdges=FEATURE)
# 파트 외 전부 제거 + 뷰포트 최대화 + fit (사용자 제안)
vp.viewportAnnotationOptions.setValues(triad=OFF, legend=OFF, title=OFF,
                                       state=OFF, compass=OFF,
                                       annotations=OFF)
vp.maximize()
vp.view.setValues(session.views['Iso'])
# noGUI 는 원근 투영 기본 카메라가 단위 스케일이라 nm 소자가 점이 됨
# -> 평행 투영 + 수치 카메라 (소자 중심, 폭 8e-6 cm)
vp.view.setValues(projection=PARALLEL)
# 좌표 nm 단위 (Viewer 줌 바닥 회피) -> 소자 중심 (24, 9.5, 21), 폭 80
vp.view.setValues(cameraTarget=(24.0, 9.5, 21.0), width=80.0)
session.pngOptions.setValues(imageSize=(1600, 1100))
session.printOptions.setValues(vpDecorations=OFF, vpBackground=OFF,
                               reduceColors=False)
session.printToFile(fileName=outp, format=PNG, canvasObjects=(vp,))
print('saved', outp, 'view width =', vp.view.width)
