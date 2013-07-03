      SUBROUTINE FPOLY(I)
C
C
C          POLYNOMIAL INITIALIZATION.
C          --------------------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/FPOLYC/ P(3),PDOT(3),RMIN2,VRMIN,JMIN
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
      IDIS = IABS (IDUM - NBP1)

      DO 8 L = 2, NNB+1
      J = LIST(L,IBINJ)
      IF (J.EQ.I) CYCLE
C     AVOID DUPLICATE COMPUTION FOR PERTURBING PLANET
      IF (KZ(3).GT.0 .AND. KZ(17).GT.0 .AND. J.EQ.JUPITER) CYCLE
C          INCLUDE ALL EMBRYOS.
      IF (BODY(J).LE.EMBRYO) THEN
         IF (IABS(IBINR - LISTR(J)).GT.NRPERT.OR.IDIS.GT.MPERT(J)) CYCLE
      END IF
C
      CALL CALC_FPOLY(I, J)
C
    8 CONTINUE

      IBINJ = IBINJ + 1
      IF (IBINJ.GT.NBTOT)  IBINJ = IBINJ - NBTOT
   10 CONTINUE

C     CALCULATE THE IMPACT OF THE PERTURBING PLANET.
      IF (KZ(3).GT.0 .AND. KZ(17).GT.0 .AND. I.NE.JUPITER) THEN
         CALL CALC_FPOLY(I, JUPITER)
      END IF
C     CALCULATE THE IMPACT OF GAS POTENTIAL.
      IF (KZ(18).GT.0) THEN
         CALL GAS_POTENTIAL(I)
       END IF
C
C     CALCULATE THE IMPACT OF GAS DAMPING.
       IF (KZ(19).GT.0.AND.
     &       SQRT(X(1,I)**2+X(2,I)**2+X(3,I)**2).LE.R_IN) THEN
          CALL GAS_DAMPING(I)
        END IF
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
      DT01(I) = STEP(I)
      DT02(I) = 2 * STEP(I)
      DT03(I) = 3 * STEP(I)
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


      SUBROUTINE CALC_FPOLY(I,J)
      INCLUDE 'COMMONP.FOR'
      COMMON/FPOLYC/ P(3),PDOT(3),RMIN2,VRMIN,JMIN
      REAL*8 A(6)
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
      END
C
C
C
      SUBROUTINE GAS_POTENTIAL(I)
      INCLUDE 'COMMONP.FOR'
      COMMON/FPOLYC/ P(3), PDOT(3),RMIN2,VRMIN,JMIN
      REAL *8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), DENS, DENS_DOT, 
     &     FIN_ABS, FIN_ABS_DOT
C
      R12 = 0.0
      R12_DOT = 0.0
      DO 100 K = 1, 3 
      R12 = R12 + X(K,I)**2
      R12_DOT = R_DOT12 + X(K,I)*XDOT(K,I)
 100  CONTINUE
      R12 = SQRT(R12)
      R12_DOT = R12_DOT/R12
C
      DO 101 K =1, 3
         THE(K) = X(K,I)/R12
         THE_DOT(K) = XDOT(K,I)/R12 - X(K,I)*R12_DOT/R12**2
 101     CONTINUE
C         
      YEAR = TIME/TWOPI
      DENS = DENS0*EXP(-YEAR/T_DEP)*R12**(-1.5)
      DENS_DOT = DENS0*(-1/T_DEP)*EXP(-YEAR/T_DEP)*R12**(-1.5) 
     &           + DENS0*(-1.5)*R12**(-2.5)*R12_DOT*EXP(-YEAR/T_DEP)
C      
      FIN_ABS = TWOPI*DENS*(0.2*(R12/R_EDGE)**(2.5)
     &          + 32*(R12/R_EDGE)**(4.5)) 
C
      FIN_ABS_DOT = TWOPI*DENS_DOT*(0.2*(R12/R_EDGE)**(2.5)
     &                              + 32*(R12/R_EDGE)**(4.5))
     &              + TWOPI*DENS*(0.5*(R12/R_EDGE)**(1.5)
     &                              + 144*(R12/R_EDGE)**(3.5))
     &              *(R12_DOT/R_EDGE)
C
      DO 102 K = 1, 3
         P(K) = P(K) + FIN_ABS*THE(K)
         PDOT(K) = PDOT(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 102     CONTINUE
C
      END  
C
C
C
      SUBROUTINE GAS_DAMPING(I)     
      INCLUDE 'COMMONP.FOR'
      COMMON/FPOLYC/ P(3), PDOT(3),RMIN2,VRMIN,JMIN
      REAL*8 R12, R12_DOT, THE(3), THE_DOT(3), VCI(3), VCI_DOT(3),
     &       DELT_V(3),DELT_V_DOT(3), ABS_DELT_V, ABS_DELT_V_DOT,
     &       HI, HI_DOT, DENS_S, DENS_S_DOT, DENS_G, DENS_G_DOT, DENS_P,
     &       T_TIDAL1_DOT, T_TIDAL2_DOT
C
      R12 = 0.0
      R12_DOT = 0.0
      V12 = 0.0
      V12_DOT = 0.0
      DO 200 K = 1, 3 
      R12 = R12 + X(K,I)**2
      R12_DOT = R_DOT12 + X(K,I)*XDOT(K,I)
      V12 = V12 + XDOT(K,I)**2
      V12_DOT = V12_DOT + XDOT(K,I)*P(K)
 200  CONTINUE
      R12 = SQRT(R12)
      R12_DOT = R12_DOT/R12
      V12 = SQRT(V12)
      V12_DOT = V12_DOT/V12
C
      DO 201 K =1, 3
         THE(K) = X(K,I)/R12
         THE_DOT(K) = XDOT(K,I)/R12 - X(K,I)*R12_DOT/R12**2
 201     CONTINUE
C
      VCI(1) = -R12**(-0.5)*THE(2)
      VCI(2) = R12**(-0.5)*THE(1)
      VCI(3) = R12**(-0.5)*THE(3)
C
      VCI_DOT(1) = 0.5*R12**(-1.5)*THE(2) - R12**(-0.5)*THE_DOT(2)
      VCI_DOT(2) = (-0.5)*R12**(-1.5)*THE(1) + R12**(-0.5)*THE_DOT(1)
      VCI_DOT(3) = (-0.5)*R12**(-1.5)*THE(3) + R12**(-0.5)*THE_DOT(3)
C
      DO 202 K = 1, 3
         DELT_V(K) = XDOT(K,I) - VCI(K)
         DELT_V_DOT(K) = P(K) - VCI_DOT(K)
 202     CONTINUE
C
      ABS_DELT_V = 0.0
      ABS_DELT_V_DOT = 0.0
      DO 203 K = 1, 3
         ABS_DELT_V = ABS_DELT_V + DELT_V(K)**2
         ABS_DELT_V_DOT = ABS_DELT_V_DOT + DELT_V(K)*DELT_V_DOT(K)
 203     CONTINUE
      ABS_DELT_V = SQRT(ABS_DELT_V)
      ABS_DELT_V_DOT = ABS_DELT_V_DOT/ABS_DELT_V
C
C
      HI = 0.05*R12**1.25
      HI_DOT = HI*1.25*R12_DOT/R12
C
      YEAR = TIME/TWOPI
      DENS_S = DENS0*EXP(-YEAR/T_DEP)
      DENS_S_DOT = DENS_S*(-1/T_DEP)
      DENS_G = DENS_S/HI
      DENS_G_DOT = DENS_G*DENS_S_DOT/DENS_S - DENS_G*HI_DOT/HI 
C
      DENS_P = BODY(I)/((2*TWOPI/3)*R(I)**3)
C
      T_TIDAL1(I) = (DENS_P/DENS_G)*R(I)*(8/3)/ABS_DELT_V
      T_TIDAL1(I) = T_TIDAL1(I)/TWOPI
      T_TIDAL1_DOT = T_TIDAL1(I)*(-DENS_G_DOT)/DENS_G
     &         + T_TIDAL1(I)*(-ABS_DELT_V_DOT)/ABS_DELT_V
C
      T_TIDAL2(I) = (1/BODY(I))*(1/(DENS_S*R12**2))*(HI/R12)**4
     &                   *(R12/V12)
      T_TIDAL2(I) = T_TIDAL2(I)/TWOPI
      T_TIDAL2_DOT = T_TIDAL2(I)*(-DENS_S_DOT)/DENS_S
     &           + T_TIDAL2(I)*(-2*R12_DOT)/R12
     &           + T_TIDAL2(I)*(4*HI_DOT)/HI
     &           + T_TIDAL2(I)*(-4*R12_DOT)/R12
     &           + T_TIDAL2(I)*R12_DOT/R12
     &           + T_TIDAL2(I)*(-V12_DOT)/V12
C
      DO 204 K = 1, 3
         P(K) = P(K) + (-DELT_V(K))*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
         PDOT(K) = PDOT(K) 
     &        + (-DELT_V_DOT(K))*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
     &        + (-DELT_V(K))*(-T_TIDAL1_DOT/T_TIDAL1(I)**2
     &           -T_TIDAL2_DOT/T_TIDAL2(I)**2)
 204     CONTINUE
     
      END

