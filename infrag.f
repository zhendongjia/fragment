      SUBROUTINE INFRAG
C
C
C          INITIALIZATION OF FRAGMENTS.
C          ----------------------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(NFMAX),VF(NFMAX),BF(NFMAX),
     & 		  AI,NF,IF(NFMAX),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
C
C
C          ADD THE ORBITAL SPIN OF CORES BEFORE COALESCING TO C.M.
      SPIN(ICOMP) = SPIN(ICOMP) + BODY(ICOMP)*BODY(JCOMP)*
     &       ((X(1,ICOMP)-X(1,JCOMP))*(XDOT(1,ICOMP)-XDOT(1,JCOMP)) +
     &        (X(2,ICOMP)-X(2,JCOMP))*(XDOT(2,ICOMP)-XDOT(2,JCOMP)))/BCM
C
C          MOVE REMNANTS TO C.M. FOR ENERGY CALCULATION (MERGED LATER).
      DO 10 K = 1,3
      X(K,ICOMP) = XCM(K)
      X(K,JCOMP) = XCM(K)
      XDOT(K,ICOMP) = VCM(K)
      XDOT(K,JCOMP) = VCM(K)
   10 CONTINUE
      CALL UPDATE_ORBIT(ICOMP)
      CALL UPDATE_ORBIT(JCOMP)
C
C          OBTAIN TOTAL KINETIC & POTENTIAL ENERGY (INCLUDING CENTRAL CORES).
      ZKIN = 0.0
      POT = 0.0
      IF(NF+1) = ICOMP
      IF(NF+2) = JCOMP
      NF2 = NF + 2
C
      DO 20 L = 1,NF2
      J = IF(L)
      DO 12 K = 1,3
 12      ZKIN = ZKIN + BODY(J)*(XDOT(K,J) - VCM(K))**2
      POTK = 0.0
      DO 15 LK = 1,NF2
      K = IF(LK)
C          SKIP MUTUAL INTERACTION OF ICOMP & JCOMP.
      IF (J.EQ.K.OR.J + K.LT.2*JCOMP)  GO TO 15
      RIJ2 = (X(1,J)-X(1,K))**2 + (X(2,J)-X(2,K))**2 +(X(3,J)-X(3,K))**2
      POTK = POTK + BODY(K)/DSQRT (RIJ2)
   15 CONTINUE
      POT = POT + BODY(J)*POTK
   20 CONTINUE
C
C          SCALE THE VELOCITIES BY AVAILABLE ENERGY AND UPDATE TOTAL ENERGY.
      ZKIN = 0.5*ZKIN
      POT = 0.5*POT
      IF (EGRAV.LT.0.0)  EGRAV = 0.0
      FACTOR = SQRT ((EGRAV + POT)/ZKIN)
      ETOT = ETOT + EGRAV
C
C          INTRODUCE EXPANDING FRAGMENT VELOCITIES WITH RESPECT TO C.M.
      DO 50 L = 1,NF
      J = IF(L)
      DO 40 K = 1,3
   40 XDOT(K,J) = VCM(K) + (XDOT(K,J) - VCM(K))*FACTOR
C
C          SET ORBITAL ELEMENTS OF THE FRAGMENTS.
      RI = DSQRT (X(1,J)**2 + X(2,J)**2 + X(3,J)**2)
      CALL UPDATE_ORBIT(J)
      IF (L.EQ.1) WRITE (6,60)  FACTOR, EGRAV
      WRITE (6,70)  J,BODY(J)/ZMOON,SEMI(J),ECC(J),RI,(X(K,J),K=1,3),
     &                    (XDOT(K,J),K=1,3)
   50 CONTINUE
C
   60 FORMAT (3X,'  FRAGM. FACTOR ENERGY:', F7.2,1PE10.1)   
   70 FORMAT (3X,'  FRAGM. :',I5,F9.4,3F8.3,6F9.4)
C
C          INITIALIZE FORCE POLYNOMIALS AND UPDATE DIAGNOSTIC VARIABLES.
      DO 80 L = 1,NF
      J = IF(L)
      CALL FPOLY(J)
      IF (T0(J) + STEP(J).GT.TLIST)  GO TO 80
      NNB = NLIST(1) + 1
      NLIST(NNB+1) = J
      NLIST(1) = NNB
   80 CONTINUE
C
      RETURN
C
      END
