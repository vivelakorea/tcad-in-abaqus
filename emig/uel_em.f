C uel_em.f --- 1D electromigration UEL (Korhonen model, TCAD in Abaqus).
C 금속 배선의 EM 응력 축적을 (V, sigma) 2-자유도 monolithic으로 푼다:
C   전기:  d/dx( (1/rho) dV/dx ) = 0   (옴 전도)
C   응력:  dsig/dt = d/dx[ kappa ( dsig/dx + (e Z* rho j)/Omega ) ]
C          (Korhonen et al., J. Appl. Phys. 73 (1993) 3790)
C   전자풍 항은 rho j = -dV/dx 로 환원 -> G = (e Z*/Omega) * E [MPa/cm].
C 끝단 원자 플럭스 차단(자연 BC), 과도해석은 backward Euler:
C   sigma_old = U - DU (SVARS 불필요). 정상상태: dsig/dx = -G (선형).
C 자유도: 1 = V [V], 2 = sigma [MPa].  단위: cm, s, MPa.
C 물성 (Al, ~500K 가속시험): rho=4.9e-6 Ohm*cm, kappa=D_a*B*Omega/kT
C   = 1.8e-9 cm^2/s (D_a=1e-11, B=7.5e10 Pa), Z*=4, Omega=1.66e-23 cm^3.
C PROPS: (1) 단면적 A [cm^2]
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),ENERGY(8),
     1 PROPS(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
      PARAMETER (RHO=4.9D-6, XKAP=1.8D-9)
      PARAMETER (ZST=4.D0, OMEG=1.66D-23, QE=1.602D-19)
C     CEZ = e Z*/Omega [C/cm^3]; CEZ*E [J/cm^4] = [MPa/cm]
      PARAMETER (CEZ=QE*ZST/OMEG)
C
      AR = PROPS(1)
      H  = DABS(COORDS(1,2) - COORDS(1,1))
      DT = DTIME
      IF (DT .LT. 1.D-30) DT = 1.D-30
      VL = AR*H/2.D0
      DO I = 1, NDOFEL
        RHS(I,1) = 0.D0
        DO J = 1, NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
      END DO
C     ---- 전기 (선형 옴 전도) ----
      GE = AR/(RHO*H)
      RHS(1,1) = RHS(1,1) - GE*(U(1) - U(3))
      RHS(3,1) = RHS(3,1) - GE*(U(3) - U(1))
      AMATRX(1,1) = AMATRX(1,1) + GE
      AMATRX(1,3) = AMATRX(1,3) - GE
      AMATRX(3,3) = AMATRX(3,3) + GE
      AMATRX(3,1) = AMATRX(3,1) - GE
C     ---- 응력: backward Euler + 확산/전자풍 플럭스 ----
C     sigma = 인장 양수. 정상상태 dsig/dx = +CEZ*E -> 캐소드(저전위) 인장,
C     애노드 압축 (전자풍이 원자를 애노드로 밀어냄; void는 캐소드).
      SAOLD = U(2) - DU(2,1)
      SBOLD = U(4) - DU(4,1)
      EF = (U(1) - U(3))/H
      GW = CEZ*EF
      PHI = XKAP*AR*((U(4) - U(2))/H - GW)
      RHS(2,1) = RHS(2,1) - VL*(U(2) - SAOLD)/DT + PHI
      RHS(4,1) = RHS(4,1) - VL*(U(4) - SBOLD)/DT - PHI
C     Jacobian (AMATRX = -dRHS/dU)
      CK = XKAP*AR/H
      CV = XKAP*AR*CEZ/H
      AMATRX(2,2) = AMATRX(2,2) + VL/DT + CK
      AMATRX(2,4) = AMATRX(2,4) - CK
      AMATRX(2,1) = AMATRX(2,1) + CV
      AMATRX(2,3) = AMATRX(2,3) - CV
      AMATRX(4,4) = AMATRX(4,4) + VL/DT + CK
      AMATRX(4,2) = AMATRX(4,2) - CK
      AMATRX(4,1) = AMATRX(4,1) - CV
      AMATRX(4,3) = AMATRX(4,3) + CV
      RETURN
      END
