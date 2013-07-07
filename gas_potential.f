C     
C     
C     
      SUBROUTINE GAS_POTENTIAL(I, P, PDOT)
      INCLUDE 'COMMONP.FOR'
      REAL *8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), DENS, DENS_DOT,
     &     FIN_ABS, FIN_ABS_DOT, P(3), PDOT(3)
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
 101  CONTINUE
C     
      YEAR = TIME/TWOPI
      DENS = DENS0*EXP(-YEAR/T_DEP)*R12**(-1.5)
      DENS_DOT = DENS0*(-1/T_DEP)*EXP(-YEAR/T_DEP)*R12**(-1.5)/TWOPI
     &     + DENS0*(-1.5)*R12**(-2.5)*R12_DOT*EXP(-YEAR/T_DEP)
C     
      FIN_ABS = TWOPI*DENS*(0.2*(R12/R_EDGE)**(2.5)
     &     + 32*(R12/R_EDGE)**(4.5)) 
C     
      FIN_ABS_DOT = TWOPI*DENS_DOT*(0.2*(R12/R_EDGE)**(2.5)
     &     + 32*(R12/R_EDGE)**(4.5))
     &     + TWOPI*DENS*(0.5*(R12/R_EDGE)**(1.5)
     &     + 144*(R12/R_EDGE)**(3.5))
     &     *(R12_DOT/R_EDGE)
C     
      DO 102 K = 1, 3
         P(K) = P(K) + FIN_ABS*THE(K)
         PDOT(K) = PDOT(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 102  CONTINUE
C     
      END  
