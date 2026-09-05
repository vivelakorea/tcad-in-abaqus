C uel_ns.f --- 나노시트 GAAFET 용 3-dof DD UEL + ODB 가시화 브리지.
C 물리는 jlfet/uel_jl.f 와 동일 (SG box method, 준페르미, mu_n=100).
C 가시화: UEL 이 요소 중심 log10(n), log10(p) 를 COMMON /VISCOM/ 에
C 적재하고, 같은 절점을 공유하는 더미 C3D8(요소번호 +500000)의 UVARM 이
C 읽어 ODB 필드(UVARM1=log10 n, UVARM2=log10 p [cm^-3])로 출력한다.
C -> Abaqus Viewer 에서 U1/U2/U3(psi/phi_n/phi_p [VT]) + 캐리어 밀도
C    컨투어/스텝 애니메이션 가능. (UEL 자체는 뷰어에 형상이 없음)
C PROPS: (1) 0=산화막/1=실리콘, (2) net doping / ni
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
      PARAMETER (Q=1.602D-19, VT=0.02585D0, XNI=1.0D10)
      PARAMETER (EPSI=11.7D0*8.854D-14, EPOX=3.9D0*8.854D-14)
      PARAMETER (XMUN=100.D0, XMUP=40.D0)
      COMMON /VISCOM/ VLOG(600000,2)
      DIMENSION IEA(12),IEB(12),IDIR(12),EN(8),P(8)
      DATA IEA /1,4,5,8, 1,2,5,6, 1,2,3,4/
      DATA IEB /2,3,6,7, 4,3,8,7, 5,6,7,8/
      DATA IDIR/1,1,1,1, 2,2,2,2, 3,3,3,3/
C
      ISI  = NINT(PROPS(1))
      DOP  = PROPS(2)
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
      IF (ISI .EQ. 1) EPS = EPSI
      CEP = EPS*VT/Q
      DN = XMUN*VT
      DP = XMUP*VT
C
      DO I = 1, NDOFEL
        RHS(I,1) = 0.D0
        DO J = 1, NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
      END DO
C     ---- 절점 농도 (준페르미 변수) ----
      DO IA = 1, 8
        KP = 3*(IA-1) + 1
        ARG1 = U(KP) - U(KP+1)
        ARG2 = U(KP+2) - U(KP)
        EN(IA) = DEXP(DMIN1(DMAX1(ARG1, -80.D0), 80.D0))
        P(IA)  = DEXP(DMIN1(DMAX1(ARG2, -80.D0), 80.D0))
      END DO
C     ---- 가시화: 요소 중심 log10 농도 [cm^-3] -> UVARM 브리지 ----
C     (실리콘만 의미 있음; 산화막은 센티널 -1 -> Viewer 에서 하한 컷)
      IF (JELEM .LE. 600000) THEN
        IF (ISI .EQ. 1) THEN
          ENC = 0.D0
          PC  = 0.D0
          DO IA = 1, 8
            ENC = ENC + EN(IA)/8.D0
            PC  = PC  + P(IA)/8.D0
          END DO
          VLOG(JELEM,1) = DLOG10(DMAX1(XNI*ENC, 1.D-30))
          VLOG(JELEM,2) = DLOG10(DMAX1(XNI*PC,  1.D-30))
        ELSE
          VLOG(JELEM,1) = -1.D0
          VLOG(JELEM,2) = -1.D0
        END IF
      END IF
C     ---- 전하 (실리콘) ----
      IF (ISI .EQ. 1) THEN
        DO IA = 1, 8
          KP = 3*(IA-1) + 1
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
        KA = 3*(IA-1) + 1
        KB = 3*(IB-1) + 1
        G = CEP*AF/H
        RHS(KA,1) = RHS(KA,1) - G*(U(KA) - U(KB))
        RHS(KB,1) = RHS(KB,1) - G*(U(KB) - U(KA))
        AMATRX(KA,KA) = AMATRX(KA,KA) + G
        AMATRX(KA,KB) = AMATRX(KA,KB) - G
        AMATRX(KB,KB) = AMATRX(KB,KB) + G
        AMATRX(KB,KA) = AMATRX(KB,KA) - G
        IF (ISI .NE. 1) GO TO 100
        T = U(KB) - U(KA)
        CALL EXNSB(T,  BT,  BPT)
        CALL EXNSB(-T, BMT, BPMT)
        ENA = EN(IA)
        ENB = EN(IB)
        PA  = P(IA)
        PB  = P(IB)
        CN = DN*AF/H
        CP = DP*AF/H
        XJN  = CN*(BT*ENB - BMT*ENA)
        XJP  = CP*(BMT*PB - BT*PA)
        DJNT = CN*(BPT*ENB + BPMT*ENA)
        DJPT = CP*(-BPMT*PB - BPT*PA)
        RHS(KA+1,1) = RHS(KA+1,1) + XJN
        RHS(KB+1,1) = RHS(KB+1,1) - XJN
        RHS(KA+2,1) = RHS(KA+2,1) - XJP
        RHS(KB+2,1) = RHS(KB+2,1) + XJP
        DJNPA = -DJNT - CN*BMT*ENA
        DJNPB =  DJNT + CN*BT*ENB
        DJNFA =  CN*BMT*ENA
        DJNFB = -CN*BT*ENB
        DJPPA = -DJPT + CP*BT*PA
        DJPPB =  DJPT - CP*BMT*PB
        DJPFA = -CP*BT*PA
        DJPFB =  CP*BMT*PB
        AMATRX(KA+1,KA)   = AMATRX(KA+1,KA)   - DJNPA
        AMATRX(KA+1,KB)   = AMATRX(KA+1,KB)   - DJNPB
        AMATRX(KA+1,KA+1) = AMATRX(KA+1,KA+1) - DJNFA
        AMATRX(KA+1,KB+1) = AMATRX(KA+1,KB+1) - DJNFB
        AMATRX(KB+1,KA)   = AMATRX(KB+1,KA)   + DJNPA
        AMATRX(KB+1,KB)   = AMATRX(KB+1,KB)   + DJNPB
        AMATRX(KB+1,KA+1) = AMATRX(KB+1,KA+1) + DJNFA
        AMATRX(KB+1,KB+1) = AMATRX(KB+1,KB+1) + DJNFB
        AMATRX(KA+2,KA)   = AMATRX(KA+2,KA)   + DJPPA
        AMATRX(KA+2,KB)   = AMATRX(KA+2,KB)   + DJPPB
        AMATRX(KA+2,KA+2) = AMATRX(KA+2,KA+2) + DJPFA
        AMATRX(KA+2,KB+2) = AMATRX(KA+2,KB+2) + DJPFB
        AMATRX(KB+2,KA)   = AMATRX(KB+2,KA)   - DJPPA
        AMATRX(KB+2,KB)   = AMATRX(KB+2,KB)   - DJPPB
        AMATRX(KB+2,KA+2) = AMATRX(KB+2,KA+2) - DJPFA
        AMATRX(KB+2,KB+2) = AMATRX(KB+2,KB+2) - DJPFB
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE EXNSB(T,B,BP)
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
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
C     더미 C3D8 (요소번호 = UEL + 500000) -> UEL 이 적재한 농도 읽기
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2),COORD(*)
      DIMENSION JMAC(*),JMATYP(*)
      COMMON /VISCOM/ VLOG(600000,2)
      JE = NOEL - 500000
      IF (JE .GE. 1 .AND. JE .LE. 600000) THEN
        UVAR(1) = VLOG(JE,1)
        UVAR(2) = VLOG(JE,2)
      END IF
      RETURN
      END
