C uel_jl_dg.f --- density-gradient 양자보정 드리프트-확산 UEL (4-dof).
C uel_jl.f 에 DG 양자보정(Ancona & Tiersten 1987, Ancona & Iafrate 1989;
C Sentaurus 의 "quantum correction")을 4번째 절점 자유도로 얹은 monolithic:
C   1=psi/VT, 2=phi_n/VT, 3=phi_p/VT, 4=sigma = ln sqrt(n/ni)
C n = exp(2*sigma) 이고 Lambda = 2*sigma - psi + phi_n 이 DG 퍼텐셜(VT단위).
C DG 방정식 (box method): (2b/VT) lap(S) = Lambda*S,  S = sqrt(n/ni) = e^sigma
C   b = gamma * hbar^2/(12 q m*),  m* = 0.32 m0, gamma = 보정계수(PROPS(3)).
C 전자 SG 플럭스: 유효 드리프트 t = d(2*sigma + phi_n)  (n = e^{2sigma} 라
C   균일 phi_n 에서 B(t)e^t = B(-t) 항등식으로 플럭스 정확히 0 = 평형 보존).
C gamma=0 극한: DG식이 sigma=(psi-phi_n)/2 를 강제 -> 고전 DD 로 환원.
C 경계(드라이버): Si/SiO2 계면 = hard wall sigma=-10 (n~0, dark space),
C   옴 접촉 = 고전 sigma = psn/2. 정공은 고전 유지(전자 소자 v1).
C PROPS: (1) 0=산화막/1=실리콘, (2) net doping / ni, (3) gamma (DG 세기)
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
C     BN0 = hbar^2/(12 q m*), m*=0.32 m0  [V cm^2]
      PARAMETER (BN0=1.985D-16)
      DIMENSION IEA(12),IEB(12),IDIR(12),EN(8),P(8),S(8)
      DATA IEA /1,4,5,8, 1,2,5,6, 1,2,3,4/
      DATA IEB /2,3,6,7, 4,3,8,7, 5,6,7,8/
      DATA IDIR/1,1,1,1, 2,2,2,2, 3,3,3,3/
C
      ISI  = NINT(PROPS(1))
      DOP  = PROPS(2)
      GAM  = PROPS(3)
      CB   = 2.D0*GAM*BN0/VT
C     sigma-행 전역 스케일: 나노 격자의 vol~1e-22 행이 솔버 0-피벗 문턱
C     아래로 떨어지는 것 방지 (해 불변, 함정 15 계열)
      DGSCL = 1.D18
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
C     ---- 미친 시도해 감지 -> 강제 cutback ----
      DO IA = 1, 8
        KP = 4*(IA-1) + 1
        IF (DABS(U(KP)) .GT. 1.D5) PNEWDT = 0.25D0
      END DO
C     ---- 절점 농도: n = e^{2 sigma}, S = e^sigma, p 고전 ----
      DO IA = 1, 8
        KP = 4*(IA-1) + 1
        EN(IA) = DEXP(DMIN1(DMAX1(2.D0*U(KP+3), -80.D0), 80.D0))
        S(IA)  = DEXP(DMIN1(DMAX1(U(KP+3), -40.D0), 40.D0))
        P(IA)  = DEXP(DMIN1(DMAX1(U(KP+2) - U(KP), -80.D0), 80.D0))
      END DO
C     ---- 전하 + DG 절점항 (실리콘) ----
      IF (ISI .EQ. 1) THEN
        DO IA = 1, 8
          KP = 4*(IA-1) + 1
          RHS(KP,1) = RHS(KP,1) + XNI*(P(IA) - EN(IA) + DOP)*VOL8
          AMATRX(KP,KP)   = AMATRX(KP,KP)   + XNI*P(IA)*VOL8
          AMATRX(KP,KP+2) = AMATRX(KP,KP+2) - XNI*P(IA)*VOL8
          AMATRX(KP,KP+3) = AMATRX(KP,KP+3) + 2.D0*XNI*EN(IA)*VOL8
C         DG: R_sig -= vol*(2sig - psi + phi_n)*S   (Lambda*S 항).
C         gamma=0 (CB=0) 은 순수 대수식 -> S-가중을 떼서 선형화
C         (S-가중 대각 vol*S*(2+Lambda) 는 Lambda=-2 에서 0 = 특이)
          ARGL = 2.D0*U(KP+3) - U(KP) + U(KP+1)
          IF (GAM .GT. 0.D0) THEN
            RHS(KP+3,1) = RHS(KP+3,1) - VOL8*DGSCL*ARGL*S(IA)
            AMATRX(KP+3,KP+3) = AMATRX(KP+3,KP+3)
     1                          + VOL8*DGSCL*S(IA)*(2.D0 + ARGL)
            AMATRX(KP+3,KP)   = AMATRX(KP+3,KP)   - VOL8*DGSCL*S(IA)
            AMATRX(KP+3,KP+1) = AMATRX(KP+3,KP+1) + VOL8*DGSCL*S(IA)
          ELSE
            RHS(KP+3,1) = RHS(KP+3,1) - VOL8*DGSCL*ARGL
            AMATRX(KP+3,KP+3) = AMATRX(KP+3,KP+3) + 2.D0*VOL8*DGSCL
            AMATRX(KP+3,KP)   = AMATRX(KP+3,KP)   - VOL8*DGSCL
            AMATRX(KP+3,KP+1) = AMATRX(KP+3,KP+1) + VOL8*DGSCL
          END IF
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
        IF (ISI .NE. 1) GO TO 100
C       DG 모서리 확산: R_sig += CB*A/h*(S_b - S_a)
        GD = CB*AF/H*DGSCL
        RHS(KA+3,1) = RHS(KA+3,1) + GD*(S(IB) - S(IA))
        RHS(KB+3,1) = RHS(KB+3,1) + GD*(S(IA) - S(IB))
        AMATRX(KA+3,KA+3) = AMATRX(KA+3,KA+3) + GD*S(IA)
        AMATRX(KA+3,KB+3) = AMATRX(KA+3,KB+3) - GD*S(IB)
        AMATRX(KB+3,KB+3) = AMATRX(KB+3,KB+3) + GD*S(IB)
        AMATRX(KB+3,KA+3) = AMATRX(KB+3,KA+3) - GD*S(IA)
C       전자 SG: 유효 드리프트 t = d(2*sigma + phi_n)
        TT = (2.D0*U(KB+3) + U(KB+1)) - (2.D0*U(KA+3) + U(KA+1))
        CALL EXDGB(TT,  BT,  BPT)
        CALL EXDGB(-TT, BMT, BPMT)
        ENA = EN(IA)
        ENB = EN(IB)
        CN = DN*AF/H
        XJN  = CN*(BT*ENB - BMT*ENA)
        DJNT = CN*(BPT*ENB + BPMT*ENA)
        RHS(KA+1,1) = RHS(KA+1,1) + XJN
        RHS(KB+1,1) = RHS(KB+1,1) - XJN
C       dXJN: sigma (t-사슬 x2 + 농도-사슬 x2), phi_n (t-사슬)
        DJNSA = -2.D0*DJNT - 2.D0*CN*BMT*ENA
        DJNSB =  2.D0*DJNT + 2.D0*CN*BT*ENB
        DJNFA = -DJNT
        DJNFB =  DJNT
        AMATRX(KA+1,KA+3) = AMATRX(KA+1,KA+3) - DJNSA
        AMATRX(KA+1,KB+3) = AMATRX(KA+1,KB+3) - DJNSB
        AMATRX(KA+1,KA+1) = AMATRX(KA+1,KA+1) - DJNFA
        AMATRX(KA+1,KB+1) = AMATRX(KA+1,KB+1) - DJNFB
        AMATRX(KB+1,KA+3) = AMATRX(KB+1,KA+3) + DJNSA
        AMATRX(KB+1,KB+3) = AMATRX(KB+1,KB+3) + DJNSB
        AMATRX(KB+1,KA+1) = AMATRX(KB+1,KA+1) + DJNFA
        AMATRX(KB+1,KB+1) = AMATRX(KB+1,KB+1) + DJNFB
C       정공 SG (고전, t = d psi)
        TP = U(KB) - U(KA)
        CALL EXDGB(TP,  BTP,  BPTP)
        CALL EXDGB(-TP, BMTP, BPMTP)
        PA  = P(IA)
        PB  = P(IB)
        CP = DP*AF/H
        XJP  = CP*(BMTP*PB - BTP*PA)
        DJPT = CP*(-BPMTP*PB - BPTP*PA)
        RHS(KA+2,1) = RHS(KA+2,1) - XJP
        RHS(KB+2,1) = RHS(KB+2,1) + XJP
        DJPPA = -DJPT + CP*BTP*PA
        DJPPB =  DJPT - CP*BMTP*PB
        DJPFA = -CP*BTP*PA
        DJPFB =  CP*BMTP*PB
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
      SUBROUTINE EXDGB(T,B,BP)
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
