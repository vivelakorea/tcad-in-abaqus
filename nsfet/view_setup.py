# -*- coding: utf-8 -*-
"""view_setup: Abaqus Viewer(GUI) 시작 스크립트 — 나노시트 GAAFET ODB 를
전자밀도(UVARM1 = log10 n) 컨투어로 열어준다.
사용: abaqus viewer script=view_setup.py  (nsfet/abq_run 에서)
보기: F 키(fit), Result->Step/Frame 로 VG 스텝 이동, 애니메이션 재생 가능.
U1/U2/U3 = psi/phi_n/phi_p [VT], UVARM2 = log10 p.
"""
from abaqus import session
from abaqusConstants import CONTOURS_ON_UNDEF, INTEGRATION_POINT, OFF, FEATURE
import visualization
import displayGroupOdbToolset as dgo

o = session.openOdb(name='ex25demo.odb')
vp = session.viewports['Viewport: 1']
vp.setValues(displayedObject=o)
vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF,))
vp.odbDisplay.setPrimaryVariable(variableLabel='UVARM1',
                                 outputPosition=INTEGRATION_POINT)
vp.odbDisplay.setFrame(step=len(o.steps) - 1, frame=-1)
try:
    inst = o.rootAssembly.instances.keys()[0]
    leaf = dgo.LeafFromElementSets(elementSets=(inst + '.VNSD',
                                                inst + '.VNCH'))
    vp.odbDisplay.displayGroup.replace(leaf=leaf)
except Exception:
    pass
vp.odbDisplay.contourOptions.setValues(minAutoCompute=OFF, minValue=4.0,
                                       maxAutoCompute=OFF, maxValue=21.0)
vp.odbDisplay.commonOptions.setValues(visibleEdges=FEATURE)
vp.view.setValues(session.views['Iso'])
vp.view.fitView()
