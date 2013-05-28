      SUBROUTINE UPDATE_ORBIT(I)
C
C
C          UPDATE THE SEMI AND ECC FOR A PARTICLE.
C          ---------------------------------------
C
      INCLUDE 'COMMONP.FOR'

      REAL*8 RI, VI2, RDOT, INVSEMI

      RI = DSQRT (X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      INVSEMI = 2.0/RI - VI2/(1.0 + BODY(I))
      SEMI(I) = 1.0/INVSEMI
      RDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
      ECC(I) = DSQRT ((1.0 - RI/SEMI(I))**2 + RDOT**2/SEMI(I))

      END
