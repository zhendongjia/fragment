      SUBROUTINE FPOLY(I)
C
C
C          POLYNOMIAL INITIALIZATION.
C          --------------------------
C
      INCLUDE 'COMMONP.FOR'
      DIMENSION  A(6),P(3),PDOT(3)
C
C
      NBP1 = NBPERT + 1
      JMIN = 0
      RMIN2 = 100.0
C
      DO 1 K = 1,3
      P(K) = 0.0
      PDOT(K) = 0.0
    1 CONTINUE
      IBIN = ILIST(I)
      IBINR = LISTR(I)
      IBINJ = IBIN - NBPERT
      IF (IBINJ.LT.1)  IBINJ = IBINJ + NBTOT
C
C          FIRST OBTAIN THE PERTURBING FORCE AND ITS DERIVATIVE.
      DO 10 IDUM = 1,NBMAX
      NNB = LIST(1,IBINJ)
      IF (NNB.EQ.0)  GO TO 9
      IDIS = IABS (IDUM - NBP1)
      L = 1
    2 L = L + 1
      J = LIST(L,IBINJ)
      IF (J.EQ.I)  GO TO 8
C          INCLUDE ALL EMBRYOS.
      IF (BODY(J).GT.EMBRYO)  GO TO 4
      IF (IABS(IBINR - LISTR(J)).GT.NRPERT.OR.IDIS.GT.MPERT(J))  GO TO 8
C
    4 RIJ2 = 0.0
      RDOT = 0.0
C
      DO 5 K = 1,3
      A(K) = X(K,J) - X(K,I)
      A(K+3) = XDOT(K,J) - XDOT(K,I)
      RIJ2 = RIJ2 + A(K)**2
      RDOT = RDOT + A(K)*A(K+3)
    5 CONTINUE
C
      RIJ3 = RIJ2*DSQRT (RIJ2)
C          FIND PLANETESIMAL CLOSEST TO BODY #I.       
      IF (RIJ2.LT.RMIN2)  THEN
          JMIN = J
          RMIN2 = RIJ2
          VRMIN = RDOT
      END IF
      RDOT = 3.0D0*RDOT/RIJ2
      FIJ = BODY(J)/RIJ3
C
      DO 7 K = 1,3
      P(K) = P(K) + A(K)*FIJ
      PDOT(K) = PDOT(K) + (A(K+3) - A(K)*RDOT)*FIJ
    7 CONTINUE
C
    8 IF (L.LE.NNB)  GO TO 2
    9 IBINJ = IBINJ + 1
      IF (IBINJ.GT.NBTOT)  IBINJ = IBINJ - NBTOT
   10 CONTINUE
C
C          ADD SOLAR CONTRIBUTION TO FORCE & THREE TAYLOR SERIES DERIVATIVES.
      SUNPL = 1.0 + BODY(I)
      RIJ2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
      FS = -SUNPL/(RIJ2*DSQRT (RIJ2))
      RDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
C
      A1 = 0.0
      A2 = 0.0
      DO 15 K = 1,3
      FI(K,I) = X(K,I)*FS + P(K)
      FDOT(K,I) = FS*(XDOT(K,I) - 3.0D0*RDOT*X(K,I)/RIJ2) + PDOT(K)
      A1 = A1 + XDOT(K,I)**2 + X(K,I)*FI(K,I)
      A2 = A2 + X(K,I)*FDOT(K,I) + 3.0D0*XDOT(K,I)*FI(K,I)
C          PLANETARY PERTURBATIONS AND NON-DOMINANT SOLAR TERMS ARE OMITTED.
   15 CONTINUE
C
      DO 20 K = 1,3
      D2(K,I) = FS*(FI(K,I) - 3.0D0*X(K,I)*A1/RIJ2)
      D3(K,I) = FS*(FDOT(K,I) - 9.0D0*XDOT(K,I)*A1/RIJ2 -
     &                                             3.0D0*X(K,I)*A2/RIJ2)
   20 CONTINUE
C
C          SPECIFY A CONSERVATIVE TIME-STEP AND INITIALIZE VARIABLES.
      FI2 = FI(1,I)**2 + FI(2,I)**2 + FI(3,I)**2
      F2DOT2 = D2(1,I)**2 + D2(2,I)**2 + D2(3,I)**2
      F3DOT2 = D3(1,I)**2 + D3(2,I)**2 + D3(3,I)**2
      A1 = ETA*DSQRT (FI2/F2DOT2)
      A2 = 0.1*ETA*F2DOT2/F3DOT2
      STEP(I) = 0.1*DSQRT (MIN (A1,A2))
C          REDUCE STEP OF MERGED BODY TO COMPENSATE FOR DERIVATIVE ERRORS.
      IF (TIME.GT.0.0)  STEP(I) = 0.1*STEP(I)
C          REDUCE EVEN FURTHER IF NEIGHBOUR IS CLOSE.
      CNT = 1.0D-03
CC      IF (BODY(I).LT.0.1*ZMOON)  CNT = 1.0D-07
      IF (JMIN.GT.0)  STEP(I) = MIN (CNT*RMIN2/DABS(VRMIN),STEP(I))
      T0(I) = TIME
      T1(I) = TIME - STEP(I)
      T2(I) = TIME - 2.0*STEP(I)
      T3(I) = TIME - 3.0*STEP(I)
      DT1 = STEP(I)
C
C          CONVERT FROM TAYLOR SERIES DERIVATIVES TO DIVIDED DIFFERENCES.
      DO 30 K = 1,3
      D1(K,I) = (ONE6*D3(K,I)*DT1 - 0.5D0*D2(K,I))*DT1 + FDOT(K,I)
      D2(K,I) = 0.5D0*D2(K,I) - 0.5D0*D3(K,I)*DT1
      D3(K,I) = ONE6*D3(K,I)
C          INITIALIZE THE PRIMARY COORDINATES AND VELOCITIES.
      X0(K,I) = X(K,I)
      X0DOT(K,I) = XDOT(K,I)
      F(K,I) = 0.5D0*FI(K,I)
      FDOT(K,I) = ONE6*FDOT(K,I)
C          HALF THE FORCE AND SIXTH THE DERIVATIVE IS STORED AS USUAL.
   30 CONTINUE
C
C          SET THE CURRENT BINDING ENERGY AND DOMINANT TERMS OF DERIVATIVES.
      H(I) = 0.5D0*(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2) -
     &                   SUNPL/DSQRT (X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      HDOT(I) = 0.0
      D1HDOT(I) = 0.0
      D2HDOT(I) = 0.0
      D3HDOT(I) = 0.0
C
      DO 40 K = 1,3
      HDOT(I) = HDOT(I) + XDOT(K,I)*P(K)
      D1HDOT(I) = D1HDOT(I) + FI(K,I)*P(K) + XDOT(K,I)*PDOT(K)
      D2HDOT(I) = D2HDOT(I) + 6.0*FDOT(K,I)*P(K) + 2.0*FI(K,I)*PDOT(K)
      D3HDOT(I) = D3HDOT(I) + D2(K,I)*P(K) + 18.0*FDOT(K,I)*PDOT(K)
   40 CONTINUE
C
C          INITIALIZE SEMI-MAJOR AXIS TO TRACE FEEDING ZONES.
      A0(I) = -0.5*SUNPL/H(I)
      IF (TIME.EQ.0.0)  GO TO 50
C
C          SEE WHETHER STABILIZATION PROCEDURE CAN BE ACTIVATED AGAIN.
      GAMMA2 = (P(1)**2 + P(2)**2 + P(3)**2)*RIJ2**2
      IF (GAMMA2.LT.DFMIN2.AND.KZ(1).GT.0)  THEN
          ISTAB(I) = 1
      ELSE
          ISTAB(I) = 0
      END IF
C
   50 RETURN
C
      END
