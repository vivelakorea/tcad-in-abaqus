C uel_mos_etm.f --- 3D NMOS electro-thermo-mechanical monolithic UEL.
C uel_mos_et.f 에 소변형 열탄성 + 압저항 피드백을 얹은 7-dof monolithic:
C   1=psi/VT, 2=phi_n/VT, 3=phi_p/VT, 4=dT[K], 5=ux, 6=uy, 7=uz [cm]
C (Abaqus 슬롯: U1-U3, UR1(dT), UR2/UR3/WARP(ux,uy,uz). *STATIC 그대로.)
C 역학: 트라이리니어 육면체 2x2x2 Gauss, sigma = C:(eps - alpha*dT*I).
C   K_uu, K_uT(열팽창) 일관 선형화. 단위: 응력 Pa, 길이 cm.
C 압저항: 요소 중심 응력으로 전자 이동도 배율
C   FPZ = 1 - (pi11*sxx + pi12*(syy+szz)),  Smith, Phys. Rev. 94 (1954) 42
C   (n-Si <100> 채널. 정공 생략 -- 채널 전류가 전자.)
C   ponytail: FPZ 는 잔차에만 반영(전기행의 u/dT-사슬 열 생략, 약결합 ~수 %.
C   수렴 문제가 생기면 그 열을 추가할 것.)
C 열/전기 부분과 함정 대책(HS 게이트, Joule 클램프, FTH)은 uel_mos_et.f 동일.
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
C     역학 물성: Si <100> / SiO2 (Pa, 1/K)
      PARAMETER (EYSI=1.30D11, PNSI=0.278D0, ALSI=2.6D-6)
      PARAMETER (EYOX=7.0D10,  PNOX=0.17D0,  ALOX=5.0D-7)
C     Smith(1954) n-Si 압저항 [1/Pa]
      PARAMETER (P11=-102.2D-11, P12=53.4D-11)
      DIMENSION IEA(12),IEB(12),IDIR(12),EN(8),P(8)
      DIMENSION XIL(8),ETL(8),ZEL(8),GP(2)
      DIMENSION DNX(8),DNY(8),DNZ(8),SN(8)
      DATA IEA /1,4,5,8, 1,2,5,6, 1,2,3,4/
      DATA IEB /2,3,6,7, 4,3,8,7, 5,6,7,8/
      DATA IDIR/1,1,1,1, 2,2,2,2, 3,3,3,3/
      DATA XIL /-1.D0, 1.D0, 1.D0,-1.D0, -1.D0, 1.D0, 1.D0,-1.D0/
      DATA ETL /-1.D0,-1.D0, 1.D0, 1.D0, -1.D0,-1.D0, 1.D0, 1.D0/
      DATA ZEL /-1.D0,-1.D0,-1.D0,-1.D0,  1.D0, 1.D0, 1.D0, 1.D0/
      DATA GP  /-0.577350269189626D0, 0.577350269189626D0/
C
      ISI  = NINT(PROPS(1))
      DOP  = PROPS(2)
C     스텝 1(평형)·2(VD=0)는 발열 0 -> 등온 해에서 웜스타트
      HS   = PROPS(3)*Q*XNI*VT
      IF (KSTEP .LE. 2) HS = 0.D0
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
      EYM = EYOX
      PRAT = PNOX
      ALP = ALOX
      IF (ISI .EQ. 1) THEN
        EPS = EPSI
        XK  = XKSI
        EYM = EYSI
        PRAT = PNSI
        ALP = ALSI
      END IF
      CEP = EPS*VT/Q
      XLAM = EYM*PRAT/((1.D0+PRAT)*(1.D0-2.D0*PRAT))
      XMUG = EYM/(2.D0*(1.D0+PRAT))
      BKA  = (3.D0*XLAM + 2.D0*XMUG)*ALP
C
      DO I = 1, NDOFEL
        RHS(I,1) = 0.D0
        DO J = 1, NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
      END DO
C     ---- 미친 시도해 감지 -> 강제 cutback (오버플로->NaN 좀비수렴 차단)
C     수렴 스케일: |psi|~1e2, |dT|~1e2. 그 1e3 배를 넘으면 포기가 답.
      DO IA = 1, 8
        KP = 7*(IA-1) + 1
        IF (DABS(U(KP)) .GT. 1.D5 .OR. DABS(U(KP+3)) .GT. 1.D4)
     1    PNEWDT = 0.25D0
      END DO
C     ---- 압저항 배율: 요소 중심 응력 (실리콘만) ----
      FPZ = 1.D0
      IF (ISI .EQ. 1) THEN
        EXX = 0.D0
        EYY = 0.D0
        EZZ = 0.D0
        TAV = 0.D0
        DO IA = 1, 8
          KU = 7*(IA-1) + 5
          EXX = EXX + XIL(IA)/(4.D0*DX)*U(KU)
          EYY = EYY + ETL(IA)/(4.D0*DY)*U(KU+1)
          EZZ = EZZ + ZEL(IA)/(4.D0*DZ)*U(KU+2)
          TAV = TAV + 0.125D0*U(7*(IA-1)+4)
        END DO
        TRE = EXX + EYY + EZZ
        SXX = XLAM*TRE + 2.D0*XMUG*EXX - BKA*TAV
        SYY = XLAM*TRE + 2.D0*XMUG*EYY - BKA*TAV
        SZZ = XLAM*TRE + 2.D0*XMUG*EZZ - BKA*TAV
        FPZ = 1.D0 - (P11*SXX + P12*(SYY + SZZ))
        IF (FPZ .LT. 0.3D0) FPZ = 0.3D0
        IF (FPZ .GT. 2.0D0) FPZ = 2.0D0
      END IF
C     ---- 절점 농도 (준페르미 변수, 300K 스케일 유지) ----
      DO IA = 1, 8
        KP = 7*(IA-1) + 1
        ARG1 = U(KP) - U(KP+1)
        ARG2 = U(KP+2) - U(KP)
        EN(IA) = DEXP(DMIN1(DMAX1(ARG1, -80.D0), 80.D0))
        P(IA)  = DEXP(DMIN1(DMAX1(ARG2, -80.D0), 80.D0))
      END DO
C     ---- 전하 (실리콘): R_psi 에 -ni(p-n+N)vol/8 ----
      IF (ISI .EQ. 1) THEN
        DO IA = 1, 8
          KP = 7*(IA-1) + 1
          RHS(KP,1) = RHS(KP,1) + XNI*(P(IA) - EN(IA) + DOP)*VOL8
          AMATRX(KP,KP)   = AMATRX(KP,KP)   + XNI*(P(IA)+EN(IA))*VOL8
          AMATRX(KP,KP+1) = AMATRX(KP,KP+1) - XNI*EN(IA)*VOL8
          AMATRX(KP,KP+2) = AMATRX(KP,KP+2) - XNI*P(IA)*VOL8
        END DO
      END IF
C     ---- 역학: 2x2x2 Gauss, K_uu + K_uT (양 재료) ----
      WV = VOL8
      DO 50 IG = 1, 2
      DO 50 JG = 1, 2
      DO 50 KG = 1, 2
        XI = GP(IG)
        ET = GP(JG)
        ZE = GP(KG)
        DO IA = 1, 8
          SN(IA) = (1.D0+XI*XIL(IA))*(1.D0+ET*ETL(IA))
     1             *(1.D0+ZE*ZEL(IA))/8.D0
          DNX(IA) = XIL(IA)*(1.D0+ET*ETL(IA))*(1.D0+ZE*ZEL(IA))
     1              /(4.D0*DX)
          DNY(IA) = ETL(IA)*(1.D0+XI*XIL(IA))*(1.D0+ZE*ZEL(IA))
     1              /(4.D0*DY)
          DNZ(IA) = ZEL(IA)*(1.D0+XI*XIL(IA))*(1.D0+ET*ETL(IA))
     1              /(4.D0*DZ)
        END DO
        EXX = 0.D0
        EYY = 0.D0
        EZZ = 0.D0
        GXY = 0.D0
        GYZ = 0.D0
        GZX = 0.D0
        TGP = 0.D0
        DO IA = 1, 8
          KU = 7*(IA-1) + 5
          EXX = EXX + DNX(IA)*U(KU)
          EYY = EYY + DNY(IA)*U(KU+1)
          EZZ = EZZ + DNZ(IA)*U(KU+2)
          GXY = GXY + DNY(IA)*U(KU) + DNX(IA)*U(KU+1)
          GYZ = GYZ + DNZ(IA)*U(KU+1) + DNY(IA)*U(KU+2)
          GZX = GZX + DNX(IA)*U(KU+2) + DNZ(IA)*U(KU)
          TGP = TGP + SN(IA)*U(7*(IA-1)+4)
        END DO
        TRE = EXX + EYY + EZZ
        SXX = XLAM*TRE + 2.D0*XMUG*EXX - BKA*TGP
        SYY = XLAM*TRE + 2.D0*XMUG*EYY - BKA*TGP
        SZZ = XLAM*TRE + 2.D0*XMUG*EZZ - BKA*TGP
        SXY = XMUG*GXY
        SYZ = XMUG*GYZ
        SZX = XMUG*GZX
        DO IA = 1, 8
          KU = 7*(IA-1) + 5
          RHS(KU,1)   = RHS(KU,1)
     1      - (SXX*DNX(IA) + SXY*DNY(IA) + SZX*DNZ(IA))*WV
          RHS(KU+1,1) = RHS(KU+1,1)
     1      - (SXY*DNX(IA) + SYY*DNY(IA) + SYZ*DNZ(IA))*WV
          RHS(KU+2,1) = RHS(KU+2,1)
     1      - (SZX*DNX(IA) + SYZ*DNY(IA) + SZZ*DNZ(IA))*WV
          DO IB = 1, 8
            KV = 7*(IB-1) + 5
            KT = 7*(IB-1) + 4
            DD = DNX(IA)*DNX(IB) + DNY(IA)*DNY(IB) + DNZ(IA)*DNZ(IB)
C           K_uu 3x3 블록: lam*dNi_a dNj_b + mu*(dNj_a dNi_b + dij*DD)
            AMATRX(KU,KV) = AMATRX(KU,KV) + WV*(XLAM*DNX(IA)*DNX(IB)
     1        + XMUG*(DNX(IA)*DNX(IB) + DD))
            AMATRX(KU,KV+1) = AMATRX(KU,KV+1) + WV*(XLAM*DNX(IA)
     1        *DNY(IB) + XMUG*DNY(IA)*DNX(IB))
            AMATRX(KU,KV+2) = AMATRX(KU,KV+2) + WV*(XLAM*DNX(IA)
     1        *DNZ(IB) + XMUG*DNZ(IA)*DNX(IB))
            AMATRX(KU+1,KV) = AMATRX(KU+1,KV) + WV*(XLAM*DNY(IA)
     1        *DNX(IB) + XMUG*DNX(IA)*DNY(IB))
            AMATRX(KU+1,KV+1) = AMATRX(KU+1,KV+1) + WV*(XLAM*DNY(IA)
     1        *DNY(IB) + XMUG*(DNY(IA)*DNY(IB) + DD))
            AMATRX(KU+1,KV+2) = AMATRX(KU+1,KV+2) + WV*(XLAM*DNY(IA)
     1        *DNZ(IB) + XMUG*DNZ(IA)*DNY(IB))
            AMATRX(KU+2,KV) = AMATRX(KU+2,KV) + WV*(XLAM*DNZ(IA)
     1        *DNX(IB) + XMUG*DNX(IA)*DNZ(IB))
            AMATRX(KU+2,KV+1) = AMATRX(KU+2,KV+1) + WV*(XLAM*DNZ(IA)
     1        *DNY(IB) + XMUG*DNY(IA)*DNZ(IB))
            AMATRX(KU+2,KV+2) = AMATRX(KU+2,KV+2) + WV*(XLAM*DNZ(IA)
     1        *DNZ(IB) + XMUG*(DNZ(IA)*DNZ(IB) + DD))
C           K_uT: d f_u / d dT_b = -BKA * dNi_a * N_b
            AMATRX(KU,KT)   = AMATRX(KU,KT)   - BKA*DNX(IA)*SN(IB)*WV
            AMATRX(KU+1,KT) = AMATRX(KU+1,KT) - BKA*DNY(IA)*SN(IB)*WV
            AMATRX(KU+2,KT) = AMATRX(KU+2,KT) - BKA*DNZ(IA)*SN(IB)*WV
          END DO
        END DO
   50 CONTINUE
C     ---- 12 모서리 순회 (Poisson / 열전도 / SG / Joule) ----
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
        KA = 7*(IA-1) + 1
        KB = 7*(IB-1) + 1
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
        CALL EX19B(TT,  BT,  BPT)
        CALL EX19B(-TT, BMT, BPMT)
        ENA = EN(IA)
        ENB = EN(IB)
        PA  = P(IA)
        PB  = P(IB)
C       전자 이동도에 압저항 배율 FPZ
        CN = XMUN*VT*AF/H*FMU*FPZ
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
C       dT 사슬. Newton 과도상태(플럭스 폭주, 수렴값 ~1e3)에서는 T열
C       결합을 끊어 T피벗(g_T 소거) 폭발/NaN 방지. 해 근방 비활성.
        FTH = 1.D0
        IF (DMAX1(DABS(XJN), DABS(XJP)) .GT. 1.D7) FTH = 0.D0
        DJNTH = FTH*(-0.5D0*TT/TK*DJNT - 0.75D0*XJN/TK)
        DJPTH = FTH*(-0.5D0*TT/TK*DJPT - 0.75D0*XJP/TK)
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
C       모서리 Joule: |PJ| <= 1e-2 W 클램프 (물리 최대 ~1e-3 W)
        PJ = HS*(XJP - XJN)*DPS
        FCLMP = 1.D0
        IF (DABS(PJ) .GT. 1.D-2) THEN
          PJ = DSIGN(1.D-2, PJ)
          FCLMP = 0.D0
        END IF
        RHS(KA+3,1) = RHS(KA+3,1) + 0.5D0*PJ
        RHS(KB+3,1) = RHS(KB+3,1) + 0.5D0*PJ
C       dPJ (T행 8열 선형화)
        DPJ1 = FCLMP*HS*((DJPPA - DJNPA)*DPS - (XJP - XJN))
        DPJ2 = FCLMP*HS*((DJPPB - DJNPB)*DPS + (XJP - XJN))
        DPJ3 = -FCLMP*HS*DJNFA*DPS
        DPJ4 = -FCLMP*HS*DJNFB*DPS
        DPJ5 = FCLMP*HS*DJPFA*DPS
        DPJ6 = FCLMP*HS*DJPFB*DPS
        DPJ7 = FCLMP*HS*(DJPTH - DJNTH)*DPS
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
      SUBROUTINE EX19B(T,B,BP)
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
