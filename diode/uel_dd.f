C uel_dd.f --- 1D pn-junction drift-diffusion UEL (TCAD in Abaqus).
C Scharfetter & Gummel (1969, IEEE Trans. Electron Devices 16:64)의
C 지수 맞춤 플럭스를 그대로 2절점 요소로 구현. 무차원화는 reference_dd1d.py와 동일
C (V_T=1, n_i=1). 절점 자유도: psi(1), n(2), p(3). 접선 비대칭 -> UNSYMM.
C   Poisson:  -psi'' = p - n + N(x),  N = +Nd (x>0), -Na (x<=0)
C   연속:     J_n' = 0, J_p' = 0
C   SG:       J_n = [B(t) n_R - B(-t) n_L]/h,  t = psi_R - psi_L
C PROPS: (1)=Na, (2)=Nd
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),ENERGY(8),
     1 PROPS(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      DIMENSION R(6),AK(6,6)
C
      XNA = PROPS(1)
      XND = PROPS(2)
      H   = COORDS(1,2) - COORDS(1,1)
C     절점 도핑 (ex10과 동일: x>0 -> +Nd, 아니면 -Na)
      DOPL = -XNA
      IF (COORDS(1,1) .GT. 0.D0) DOPL = XND
      DOPR = -XNA
      IF (COORDS(1,2) .GT. 0.D0) DOPR = XND
C     0 초기값에서 출발하므로 1단계에서는 도핑을 시간 램프
C     (TCAD의 doping continuation; 접촉 BC 램프와 보조를 맞춘다)
      IF (KSTEP .EQ. 1) THEN
        SCALE = TIME(1) + DTIME
        IF (SCALE .GT. 1.D0) SCALE = 1.D0
        DOPL = DOPL*SCALE
        DOPR = DOPR*SCALE
      END IF
C
      PSIL = U(1)
      ENL  = U(2)
      PL   = U(3)
      PSIR = U(4)
      ENR  = U(5)
      PR   = U(6)
      T = PSIR - PSIL
      CALL EX16B(T,  BT,  BPT)
      CALL EX16B(-T, BMT, BPMT)
C     SG 플럭스와 t-미분
      XJN  = (BT*ENR - BMT*ENL)/H
      XJP  = (BMT*PR - BT*PL)/H
      DJNT = (BPT*ENR + BPMT*ENL)/H
      DJPT = (-BPMT*PR - BPT*PL)/H
C     잔차 (ex10의 절점 조립과 동일한 부호)
      R(1) = (PSIL-PSIR)/H - 0.5D0*H*(PL-ENL+DOPL)
      R(4) = (PSIR-PSIL)/H - 0.5D0*H*(PR-ENR+DOPR)
      R(2) = -XJN
      R(5) =  XJN
      R(3) =  XJP
      R(6) = -XJP
C     Jacobian
      DO I = 1, 6
        DO J = 1, 6
          AK(I,J) = 0.D0
        END DO
      END DO
C     Poisson 행
      AK(1,1) =  1.D0/H
      AK(1,4) = -1.D0/H
      AK(1,2) =  0.5D0*H
      AK(1,3) = -0.5D0*H
      AK(4,4) =  1.D0/H
      AK(4,1) = -1.D0/H
      AK(4,5) =  0.5D0*H
      AK(4,6) = -0.5D0*H
C     n 연속 행 (row2 = -Jn, row5 = +Jn)
      AK(2,1) =  DJNT
      AK(2,4) = -DJNT
      AK(2,2) =  BMT/H
      AK(2,5) = -BT/H
      AK(5,1) = -DJNT
      AK(5,4) =  DJNT
      AK(5,2) = -BMT/H
      AK(5,5) =  BT/H
C     p 연속 행 (row3 = +Jp, row6 = -Jp)
      AK(3,1) = -DJPT
      AK(3,4) =  DJPT
      AK(3,3) = -BT/H
      AK(3,6) =  BMT/H
      AK(6,1) =  DJPT
      AK(6,4) = -DJPT
      AK(6,3) =  BT/H
      AK(6,6) = -BMT/H
C
      DO I = 1, 6
        RHS(I,1) = -R(I)
        DO J = 1, 6
          AMATRX(I,J) = AK(I,J)
        END DO
      END DO
      RETURN
      END
C
      SUBROUTINE EX16B(T,B,BP)
C     Bernoulli B(t)=t/(e^t-1) 와 B'(t), 안정 구현
      INCLUDE 'ABA_PARAM.INC'
      TC = T
      IF (TC .GT.  100.D0) TC =  100.D0
      IF (TC .LT. -100.D0) TC = -100.D0
C     |t| 작으면 급수로 (DEXP(t)-1 의 상쇄 오차 회피; expm1 대용)
      IF (DABS(TC) .LT. 1.D-4) THEN
        B  = 1.D0 - TC/2.D0 + TC*TC/12.D0
        BP = -0.5D0 + TC/6.D0
      ELSE
        E  = DEXP(TC) - 1.D0
        B  = TC/E
        BP = (E - TC*DEXP(TC))/(E*E)
      END IF
      RETURN
      END
