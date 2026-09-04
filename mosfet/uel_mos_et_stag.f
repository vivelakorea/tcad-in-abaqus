C uel_mos_et_stag.f --- staggered형(분할반복) 비교용 변형판.
C uel_mos_et.f 와 잔차(물리)는 완전히 동일하고, Jacobian 의 전기<->열
C 교차 블록만 제거: DJNTH/DJPTH(전기행의 T열), DPJ1~7(T행의 전기열 +
C Joule 의 T-대각). 즉 각 Newton 반복이 "전기 풀고 열 풀기"를 한 번에
C 하는 block-Jacobi 분할반복이 된다 -- staggered TCAD<->열 루프의
C 수렴 특성(선형 수렴, 강결합 발산)을 같은 코드베이스에서 재현.
C run_staggered.py 가 HSCALE 을 올리며 monolithic 과 반복수/수렴 비교.
C PROPS: (1) 0=산화막/1=실리콘, (2) net doping / ni, (3) 발열 배율 HSCALE
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
      PARAMETER (Q=1.602D-19, VT=0.02585D0, XNI=1.0D10, T0=300.D0)
      PARAMETER (EPSI=11.7D0*8.854D-14, EPOX=3.9D0*8.854D-14)
      PARAMETER (XMUN=400.D0, XMUP=150.D0)
      PARAMETER (XKSI=1.5D0, XKOX=0.014D0)
      DIMENSION IEA(12),IEB(12),IDIR(12),EN(8),P(8)
      DATA IEA /1,4,5,8, 1,2,5,6, 1,2,3,4/
      DATA IEB /2,3,6,7, 4,3,8,7, 5,6,7,8/
      DATA IDIR/1,1,1,1, 2,2,2,2, 3,3,3,3/
C
      ISI  = NINT(PROPS(1))
      DOP  = PROPS(2)
C     발열 배율 * q*ni*VT -> 모서리 Joule [W] 환산 계수.
C     스텝 1(평형)·2(VG, VD=0)는 발열 0 (물리적으로도 0) -> 열 결합을 끄고
C     등온 해에서 웜스타트. cold start 의 e80 플럭스 폭주가 T장에 못 들어감.
      HS   = PROPS(3)*Q*XNI*VT
      IF (KSTEP .LE. 2) HS = 0.D0
C     0 초기값 출발 -> 1단계에서 도핑 램프 (continuation)
      IF (KSTEP .EQ. 1) THEN
        SCALE = TIME(1) + DTIME
        IF (SCALE .GT. 1.D0) SCALE = 1.D0
        DOP = DOP*SCALE
      END IF
      DX = DABS(COORDS(1,2) - COORDS(1,1))
      DY = DABS(COORDS(2,4) - COORDS(2,1))
      DZ = DABS(COORDS(3,5) - COORDS(3,1))
      VOL8 = DX*DY*DZ/8.D0
      EPS = EPOX
      XK  = XKOX
      IF (ISI .EQ. 1) THEN
        EPS = EPSI
        XK  = XKSI
      END IF
      CEP = EPS*VT/Q
C
      DO I = 1, NDOFEL
        RHS(I,1) = 0.D0
        DO J = 1, NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
      END DO
C     ---- 미친 시도해 감지 -> 강제 cutback (오버플로->NaN 좀비수렴 차단)
      DO IA = 1, 8
        KP = 4*(IA-1) + 1
        IF (DABS(U(KP)) .GT. 1.D5 .OR. DABS(U(KP+3)) .GT. 1.D4)
     1    PNEWDT = 0.25D0
      END DO
C     ---- 절점 농도 (준페르미 변수, 300K 스케일 유지) ----
      DO IA = 1, 8
        KP = 4*(IA-1) + 1
        ARG1 = U(KP) - U(KP+1)
        ARG2 = U(KP+2) - U(KP)
        EN(IA) = DEXP(DMIN1(DMAX1(ARG1, -80.D0), 80.D0))
        P(IA)  = DEXP(DMIN1(DMAX1(ARG2, -80.D0), 80.D0))
      END DO
C     ---- 전하 (실리콘): R_psi 에 -ni(p-n+N)vol/8 ----
      IF (ISI .EQ. 1) THEN
        DO IA = 1, 8
          KP = 4*(IA-1) + 1
          RHS(KP,1) = RHS(KP,1) + XNI*(P(IA) - EN(IA) + DOP)*VOL8
          AMATRX(KP,KP)   = AMATRX(KP,KP)   + XNI*(P(IA)+EN(IA))*VOL8
          AMATRX(KP,KP+1) = AMATRX(KP,KP+1) - XNI*EN(IA)*VOL8
          AMATRX(KP,KP+2) = AMATRX(KP,KP+2) - XNI*P(IA)*VOL8
        END DO
      END IF
C     ---- 12 모서리 순회 ----
      DO 100 IE = 1, 12
        IA = IEA(IE)
        IB = IEB(IE)
        IF (IDIR(IE) .EQ. 1) THEN
          H = DX
          AF = DY*DZ/4.D0
        ELSE IF (IDIR(IE) .EQ. 2) THEN
          H = DY
          AF = DX*DZ/4.D0
        ELSE
          H = DZ
          AF = DX*DY/4.D0
        END IF
        KA = 4*(IA-1) + 1
        KB = 4*(IB-1) + 1
C       Poisson 모서리 전도
        G = CEP*AF/H
        RHS(KA,1) = RHS(KA,1) - G*(U(KA) - U(KB))
        RHS(KB,1) = RHS(KB,1) - G*(U(KB) - U(KA))
        AMATRX(KA,KA) = AMATRX(KA,KA) + G
        AMATRX(KA,KB) = AMATRX(KA,KB) - G
        AMATRX(KB,KB) = AMATRX(KB,KB) + G
        AMATRX(KB,KA) = AMATRX(KB,KA) - G
C       열전도 (산화막 포함): g_T = k*A/h, 선형
        GT = XK*AF/H
        RHS(KA+3,1) = RHS(KA+3,1) - GT*(U(KA+3) - U(KB+3))
        RHS(KB+3,1) = RHS(KB+3,1) - GT*(U(KB+3) - U(KA+3))
        AMATRX(KA+3,KA+3) = AMATRX(KA+3,KA+3) + GT
        AMATRX(KA+3,KB+3) = AMATRX(KA+3,KB+3) - GT
        AMATRX(KB+3,KB+3) = AMATRX(KB+3,KB+3) + GT
        AMATRX(KB+3,KA+3) = AMATRX(KB+3,KA+3) - GT
        IF (ISI .NE. 1) GO TO 100
C       모서리 온도, V_T(T)/mu(T) 인자
        TK = T0 + 0.5D0*(U(KA+3) + U(KB+3))
        IF (TK .LT. 200.D0)  TK = 200.D0
        IF (TK .GT. 1500.D0) TK = 1500.D0
        RT = T0/TK
        FMU = RT*DSQRT(RT)
C       SG 플럭스 (a -> b), t = dpsi * T0/T  (dpsi 는 VT(300K) 스케일)
        DPS = U(KB) - U(KA)
        TT = RT*DPS
        CALL EX18B(TT,  BT,  BPT)
        CALL EX18B(-TT, BMT, BPMT)
        ENA = EN(IA)
        ENB = EN(IB)
        PA  = P(IA)
        PB  = P(IB)
        CN = XMUN*VT*AF/H*FMU
        CP = XMUP*VT*AF/H*FMU
        XJN  = CN*(BT*ENB - BMT*ENA)
        XJP  = CP*(BMT*PB - BT*PA)
        DJNT = CN*(BPT*ENB + BPMT*ENA)
        DJPT = CP*(-BPMT*PB - BPT*PA)
        RHS(KA+1,1) = RHS(KA+1,1) + XJN
        RHS(KB+1,1) = RHS(KB+1,1) - XJN
        RHS(KA+2,1) = RHS(KA+2,1) - XJP
        RHS(KB+2,1) = RHS(KB+2,1) + XJP
C       dXJN: psi (t-사슬 x RT + 농도-사슬), phi_n, dT
        DJNPA = -RT*DJNT - CN*BMT*ENA
        DJNPB =  RT*DJNT + CN*BT*ENB
        DJNFA =  CN*BMT*ENA
        DJNFB = -CN*BT*ENB
C       dXJP
        DJPPA = -RT*DJPT + CP*BT*PA
        DJPPB =  RT*DJPT - CP*BMT*PB
        DJPFA = -CP*BT*PA
        DJPFB =  CP*BMT*PB
C       dT 사슬: dt/dT_node = -0.5 t/TK, dmu-인자 -> -0.75 XJ/TK.
C       Newton 과도상태(플럭스 폭주, 수렴값 ~1e3)에서는 T열 결합을 끊어
C       T피벗(g_T~1e-3) 소거 폭발/NaN 방지. 해 근방 비활성 -> 일관성 유지.
C       [staggered] 전기행의 T열 결합 제거
        DJNTH = 0.D0
        DJPTH = 0.D0
C       n 행 (row a = -XJN, row b = +XJN)
        AMATRX(KA+1,KA)   = AMATRX(KA+1,KA)   - DJNPA
        AMATRX(KA+1,KB)   = AMATRX(KA+1,KB)   - DJNPB
        AMATRX(KA+1,KA+1) = AMATRX(KA+1,KA+1) - DJNFA
        AMATRX(KA+1,KB+1) = AMATRX(KA+1,KB+1) - DJNFB
        AMATRX(KA+1,KA+3) = AMATRX(KA+1,KA+3) - DJNTH
        AMATRX(KA+1,KB+3) = AMATRX(KA+1,KB+3) - DJNTH
        AMATRX(KB+1,KA)   = AMATRX(KB+1,KA)   + DJNPA
        AMATRX(KB+1,KB)   = AMATRX(KB+1,KB)   + DJNPB
        AMATRX(KB+1,KA+1) = AMATRX(KB+1,KA+1) + DJNFA
        AMATRX(KB+1,KB+1) = AMATRX(KB+1,KB+1) + DJNFB
        AMATRX(KB+1,KA+3) = AMATRX(KB+1,KA+3) + DJNTH
        AMATRX(KB+1,KB+3) = AMATRX(KB+1,KB+3) + DJNTH
C       p 행 (row a = +XJP, row b = -XJP)
        AMATRX(KA+2,KA)   = AMATRX(KA+2,KA)   + DJPPA
        AMATRX(KA+2,KB)   = AMATRX(KA+2,KB)   + DJPPB
        AMATRX(KA+2,KA+2) = AMATRX(KA+2,KA+2) + DJPFA
        AMATRX(KA+2,KB+2) = AMATRX(KA+2,KB+2) + DJPFB
        AMATRX(KA+2,KA+3) = AMATRX(KA+2,KA+3) + DJPTH
        AMATRX(KA+2,KB+3) = AMATRX(KA+2,KB+3) + DJPTH
        AMATRX(KB+2,KA)   = AMATRX(KB+2,KA)   - DJPPA
        AMATRX(KB+2,KB)   = AMATRX(KB+2,KB)   - DJPPB
        AMATRX(KB+2,KA+2) = AMATRX(KB+2,KA+2) - DJPFA
        AMATRX(KB+2,KB+2) = AMATRX(KB+2,KB+2) - DJPFB
        AMATRX(KB+2,KA+3) = AMATRX(KB+2,KA+3) - DJPTH
        AMATRX(KB+2,KB+3) = AMATRX(KB+2,KB+3) - DJPTH
C       모서리 Joule: PJ = HS*(Jp - Jn)*dpsi [W], 절점 반절 lumping.
C       Newton 과도상태에서 T장 폭주 방지: |PJ| <= 1e-2 W 클램프
C       (물리 최대 모서리 Joule ~1e-3 W -> 해 근방 비활성, 일관성 유지)
        PJ = HS*(XJP - XJN)*DPS
        FCLMP = 1.D0
        IF (DABS(PJ) .GT. 1.D-2) THEN
          PJ = DSIGN(1.D-2, PJ)
          FCLMP = 0.D0
        END IF
        RHS(KA+3,1) = RHS(KA+3,1) + 0.5D0*PJ
        RHS(KB+3,1) = RHS(KB+3,1) + 0.5D0*PJ
C       [staggered] T행의 전기열/Joule-대각 결합 제거 (열은 소스 동결)
        DPJ1 = 0.D0
        DPJ2 = 0.D0
        DPJ3 = 0.D0
        DPJ4 = 0.D0
        DPJ5 = 0.D0
        DPJ6 = 0.D0
        DPJ7 = 0.D0
        AMATRX(KA+3,KA)   = AMATRX(KA+3,KA)   - 0.5D0*DPJ1
        AMATRX(KA+3,KB)   = AMATRX(KA+3,KB)   - 0.5D0*DPJ2
        AMATRX(KA+3,KA+1) = AMATRX(KA+3,KA+1) - 0.5D0*DPJ3
        AMATRX(KA+3,KB+1) = AMATRX(KA+3,KB+1) - 0.5D0*DPJ4
        AMATRX(KA+3,KA+2) = AMATRX(KA+3,KA+2) - 0.5D0*DPJ5
        AMATRX(KA+3,KB+2) = AMATRX(KA+3,KB+2) - 0.5D0*DPJ6
        AMATRX(KA+3,KA+3) = AMATRX(KA+3,KA+3) - 0.5D0*DPJ7
        AMATRX(KA+3,KB+3) = AMATRX(KA+3,KB+3) - 0.5D0*DPJ7
        AMATRX(KB+3,KA)   = AMATRX(KB+3,KA)   - 0.5D0*DPJ1
        AMATRX(KB+3,KB)   = AMATRX(KB+3,KB)   - 0.5D0*DPJ2
        AMATRX(KB+3,KA+1) = AMATRX(KB+3,KA+1) - 0.5D0*DPJ3
        AMATRX(KB+3,KB+1) = AMATRX(KB+3,KB+1) - 0.5D0*DPJ4
        AMATRX(KB+3,KA+2) = AMATRX(KB+3,KA+2) - 0.5D0*DPJ5
        AMATRX(KB+3,KB+2) = AMATRX(KB+3,KB+2) - 0.5D0*DPJ6
        AMATRX(KB+3,KA+3) = AMATRX(KB+3,KA+3) - 0.5D0*DPJ7
        AMATRX(KB+3,KB+3) = AMATRX(KB+3,KB+3) - 0.5D0*DPJ7
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE EX18B(T,B,BP)
C     Bernoulli B(t)=t/(e^t-1), B'(t). 급수 분기로 상쇄 오차 회피.
      INCLUDE 'ABA_PARAM.INC'
      TC = T
      IF (TC .GT.  100.D0) TC =  100.D0
      IF (TC .LT. -100.D0) TC = -100.D0
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
