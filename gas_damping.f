C     
C     
C     
      SUBROUTINE GAS_DAMPING(I, P, PDOT)     
      INCLUDE 'COMMONP.FOR'
      REAL*8 R12, R12_DOT, THE(3), THE_DOT(3), VCI(3), VCI_DOT(3),
     &     DELT_V(3),DELT_V_DOT(3), ABS_DELT_V, ABS_DELT_V_DOT,
     &     HI, HI_DOT, DENS_S, DENS_S_DOT, DENS_G, DENS_G_DOT, DENS_P,
     &     T_TIDAL1_DOT, T_TIDAL2_DOT, P(3), PDOT(3)
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
 201  CONTINUE
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
 202  CONTINUE
C     
      ABS_DELT_V = 0.0
      ABS_DELT_V_DOT = 0.0
      DO 203 K = 1, 3
         ABS_DELT_V = ABS_DELT_V + DELT_V(K)**2
         ABS_DELT_V_DOT = ABS_DELT_V_DOT + DELT_V(K)*DELT_V_DOT(K)
 203  CONTINUE
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
     &     + T_TIDAL1(I)*(-ABS_DELT_V_DOT)/ABS_DELT_V
C     
      T_TIDAL2(I) = (1/BODY(I))*(1/(DENS_S*R12**2))*(HI/R12)**4
     &     *(R12/V12)
      T_TIDAL2(I) = T_TIDAL2(I)/TWOPI
      T_TIDAL2_DOT = T_TIDAL2(I)*(-DENS_S_DOT)/DENS_S
     &     + T_TIDAL2(I)*(-2*R12_DOT)/R12
     &     + T_TIDAL2(I)*(4*HI_DOT)/HI
     &     + T_TIDAL2(I)*(-4*R12_DOT)/R12
     &     + T_TIDAL2(I)*R12_DOT/R12
     &     + T_TIDAL2(I)*(-V12_DOT)/V12
C     
      DO 204 K = 1, 3
         P(K) = P(K) + (-DELT_V(K))*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
         PDOT(K) = PDOT(K) 
     &        + (-DELT_V_DOT(K))*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
     &        + (-DELT_V(K))*(-T_TIDAL1_DOT/T_TIDAL1(I)**2
     &        -T_TIDAL2_DOT/T_TIDAL2(I)**2)
 204  CONTINUE
      
      END

