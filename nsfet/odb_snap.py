# -*- coding: utf-8 -*-
"""odb_snap: ODB 를 Abaqus/CAE 로 렌더해 PNG 저장 (noGUI).
사용: abaqus cae noGUI=odb_snap.py -- <odb> <out_png_base>
UVARM1(log10 n) 컨투어, 실리콘 elset 만 표시, 마지막 프레임.
"""
import sys
from abaqus import session
from abaqusConstants import (CONTOURS_ON_UNDEF, INTEGRATION_POINT, PNG, OFF,
                             ON, FEATURE)
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
vp.odbDisplay.setFrame(step=nstep - 1, frame=-1)
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
# noGUI 에서 fitView 가 무시됨 -> 카메라 수치 지정 (cm 좌표, 소자 중심)
vp.view.setValues(session.views['Iso'])
vp.view.setValues(cameraTarget=(2.4e-06, 9.5e-07, 2.1e-06), width=8.0e-06)
vp.view.zoom(20.0)
session.printOptions.setValues(vpDecorations=ON, reduceColors=False)
session.printToFile(fileName=outp, format=PNG, canvasObjects=(vp,))
print('saved', outp)
