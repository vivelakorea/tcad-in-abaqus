# -*- coding: utf-8 -*-
"""odb_movie: 나노시트 GAAFET ODB 의 시간이력 애니메이션을 동영상(AVI)으로.
사용: abaqus cae noGUI=odb_movie.py -- <odb> <out_avi_base>
odb_snap 과 같은 화면 구성(파트만, 전체화면, 평행투영 nm 카메라)으로
전 스텝(도핑 램프 -> VD -> VG 램프)의 UVARM1(log10 n) 컨투어를 기록.
"""
import sys
from abaqus import session
from abaqusConstants import (CONTOURS_ON_UNDEF, INTEGRATION_POINT, OFF, ON,
                             FEATURE, PARALLEL, PNG)
import visualization
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
vp.viewportAnnotationOptions.setValues(triad=OFF, legend=OFF, title=OFF,
                                       state=OFF, compass=OFF,
                                       annotations=OFF)
vp.maximize()
vp.view.setValues(session.views['Iso'])
vp.view.setValues(projection=PARALLEL)
vp.view.setValues(cameraTarget=(24.0, 9.5, 21.0), width=80.0)
# noGUI 에는 animationController 가 없음 -> 마지막 스텝(VG 램프)의 프레임을
# 개별 PNG 로 찍고, 파이썬(Pillow)에서 동영상으로 조립한다.
session.pngOptions.setValues(imageSize=(960, 660))
session.printOptions.setValues(vpDecorations=OFF, vpBackground=OFF,
                               reduceColors=False)
nstep = len(o.steps)
nfr = len(o.steps[o.steps.keys()[nstep - 1]].frames)
for i in range(1, nfr):
    vp.odbDisplay.setFrame(step=nstep - 1, frame=i)
    session.printToFile(fileName='%s_f%03d' % (outp, i), format=PNG,
                        canvasObjects=(vp,))
print('frames saved:', nfr - 1)
