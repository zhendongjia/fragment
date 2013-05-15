      SUBROUTINE FRAGMT(I)
C
C
C          GENERATION OF FRAGMENTS.
C          ------------------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(8),VF(8),BF(8),AI,
     & 		  NF,IF(10),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
      REAL*4  RAN2,WK(2)
C
C
      BMAX = 0.5*ZMEJ*AI**1.24
      BETA = 4.0*BMAX/ZMEJ
C
C          SPECIFY FRAGMENT MASS FUNCTION FOR BODY #I.
      IF (BETA.GE.1.0)  THEN
          BF(1) = BMAX
          BF(2) = BMAX/BETA
          BF(3) = 0.5*(3.0 - BETA)*BMAX/BETA
          BF(4) = BF(3)
      ELSE
          DO 10 L = 1,4
   10     BF(L) = 0.25*ZMEJ
      END IF
C
C          REDUCE MASS & RADIUS OF BODY #I AND UPDATE GRAVITATIONAL ENERGY.
      RI = R(I)
      BODYI = BODY(I)
      BODY(I) = BODY(I) - ZMEJ
      R(I) = RI*(BODY(I)/BODYI)**0.3333
      EGRAV = EGRAV + 0.6*BODY(I)**2/R(I)
C
C          ASSIGN NEW LOCATIONS FOR THE FRAGMENTS AND INCREASE N.
      DO 20 L = 1,4
      NF = NF + 1
      J = N + 1
      IF(NF) = J			
      N = N + 1
   20 CONTINUE
C
C          GENERATE COORDINATES & VELOCITIES WITH RESPECT TO CENTRE OF MASS.
      RESC = 10.0*RI*(BCM/BODYI)**0.3333
      VREL = SQRT (2.0*BCM/RESC)
      KKK = 1
C
      DO 30 L = 1,4
      PHI = TWOPI*(FLOAT(L) - 0.5 + (FLOAT(NF) - 4.0)/8.0)/4.0	
      LK = 2*(L - 1)
      XF(LK+1) = RESC*DCOS (PHI)*(1.0 + 0.1*RAN2(KKK))
      XF(LK+2) = RESC*DSIN (PHI)*(1.0 + 0.1*RAN2(KKK))
      VF(LK+1) = VREL*XF(LK+1)/RESC
      VF(LK+2) = VREL*XF(LK+2)/RESC
   30 CONTINUE
C
C          SPECIFY GLOBAL VARIABLES FOR THE FRAGMENTS (NO Z-DISPERSION).
      DO 50 L = 1,4
      J = IF(L+NF-4)
C
      DO 45 K = 1,3
      LK = 2*(L - 1) + K
      IF (K.LE.2)  THEN
          X(K,J) = XCM(K) + XF(LK)
          XDOT(K,J) = VCM(K) + VF(LK)
      ELSE
          X(K,J) = XCM(K)
          XDOT(K,J) = VCM(K)
      END IF
   45 CONTINUE
C
      BODY(J) = BF(L)
      ISTAB(J) = 0
C          LINK FRAGMENTS TO THEIR ORIGIN.
      NAME(J) = NAME(I)
C          ALLOCATE FRACTIONAL SPIN.
      SPIN(J) = SPIN(I)*BODY(J)/BODYI
      R(J) = RI*(BODY(J)/BODYI)**0.3333
      EGRAV = EGRAV + 0.6*BODY(J)**2/R(J)
      LISTR(J) = ROUT2/(X(1,J)**2 + X(2,J)**2)
      WK(1) = X(1,J)
      WK(2) = X(2,J)
      THETA = ATAN2 (WK(2),WK(1))
      IF (THETA.LT.0.0)  THETA = THETA + TWOPI
      IBIN = 1 + THETA/ALPHA
      ILIST(J) = IBIN
      NNB = LIST(1,IBIN) + 1
      LIST(NNB+1,IBIN) = J
      LIST(1,IBIN) = NNB
      MPERT(J) = 1 + SQRT (BODY(J)/EMBRYO)*NBPERT
   50 CONTINUE
C
C          RETAIN FRACTIONAL PART OF INTRINSIC REMNANT SPIN.
      SPIN(I) = SPIN(I)*BODY(I)/BODYI
      NSTEPN(9) = NSTEPN(9) + 1
C
      RETURN
C
      END
