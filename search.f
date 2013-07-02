       SUBROUTINE SEARCH(I,ICASE,IBINT,RCL2)
C
C
C          COLLISION SEARCH.
C          -----------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(NFMAX),VF(NFMAX),BF(NFMAX),
     &		  AI,NF,IF(NFMAX),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
C
C
      IBIN = IBINT
      ITRY = 0
C          FIRST SEARCH CURRENT BIN OF BODY #I.
   10 NNB = LIST(1,IBINT)
      JCOMP = 0
      IF (NNB.EQ.0)  GO TO 30
      L = 1
   15 L = L + 1
      J = LIST(L,IBINT)
      IF (J.EQ.I)  GO TO 25
      S = TIME - T0(J)
      RIJ2 = 0.0
C
      DO 20 K = 1,3
      XIJ = ((FDOT(K,J)*S + F(K,J))*S + X0DOT(K,J))*S + X0(K,J) - X(K,I)
      RIJ2 = RIJ2 + XIJ**2
   20 CONTINUE 
C
      IF (RIJ2.GT.RCL2)  GO TO 25
      RCL2 = RIJ2
      JCOMP = J
   25 IF (L.LE.NNB)  GO TO 15
C
   30 IF (JCOMP.NE.0.OR.ITRY.EQ.2)  GO TO 40
C          TRY THE NEIGHBOURING BINS AS WELL.
      ITRY = ITRY + 1
      IBINT = IBIN + 2*ITRY - 3
      IF (IBINT.LT.1)  IBINT = IBINT + NBTOT
      IF (IBINT.GT.NBTOT)  IBINT = IBINT - NBTOT
      GO TO 10
C
   40 IF (JCOMP.EQ.0)  THEN
          NSTEPN(7) = NSTEPN(7) + 1
          GO TO 100
      END IF
C
C          SKIP PERICENTRE CALCULATION IF APPROXIMATE RDOT > 0.
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
      IF (RDOT.GT.0.0.AND.ICASE.EQ.1)  GO TO 100
      NSTEPN(3) = NSTEPN(3) + 1
C
C          PREDICT CURRENT COORDINATES AND VELOCITIES FOR BODY JCOMP.
      DT = TIME - T0(JCOMP)
      A1 = 0.05D0*DT
      A2 = 0.25D0*DT
      A3 = DT01(JCOMP) + DT02(JCOMP)
C
      RIJ2 = 0.0
      VREL2 = 0.0
      RDOT = 0.0
      DO 45 K = 1,3
      F2DOTK = D3(K,JCOMP)*A3 + D2(K,JCOMP)
      X(K,JCOMP) = ((((D3(K,JCOMP)*A1 + ONE12*F2DOTK)*DT + 
     &		                  FDOT(K,JCOMP))*DT + F(K,JCOMP))*DT + 
     &                                  X0DOT(K,JCOMP))*DT + X0(K,JCOMP)
      XDOT(K,JCOMP) =  (((D3(K,JCOMP)*A2 + ONE3*F2DOTK)*DT +
     &                  3.0D0*FDOT(K,JCOMP))*DT + 2.0D0*F(K,JCOMP))*DT +
     &                                                    X0DOT(K,JCOMP)
      RIJ2 = RIJ2 + (X(K,I) - X(K,JCOMP))**2
      VREL2 = VREL2 + (XDOT(K,I) - XDOT(K,JCOMP))**2
      RDOT = RDOT + (X(K,I) - X(K,JCOMP))*(XDOT(K,I) - XDOT(K,JCOMP))
   45 CONTINUE
      CALL UPDATE_ORBIT(JCOMP)
C
C          EVALUATE OSCULATING SEMI-MAJOR AXIS AND ECCENTRICITY.
      RIJ = SQRT (RIJ2)
      BCM = BODY(I) + BODY(JCOMP)
C          ECCENTRICITY EXPRESSION HOLDS FOR BOTH ELLIPTIC & HYPERBOLIC CASE.
      PERI = SEMI(I)*(1.0 - ECC(I))
      IF (KZ(5).LT.3.OR.I.GT.JCOMP)  GO TO 60
C
      YEARS = TIME/TWOPI
CC    GAMMA = DSQRT (GAMMA2)
      GAMMA = 0.0
      RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
      GAMMAS = 2.0*RIJ**3/(BCM*RI2*DSQRT (RI2))
      WRITE (6,50)  NAME(I),NAME(JCOMP),YEARS,RIJ,SEMI(I),ECC(I),PERI,
     &     GAMMA, GAMMAS,STEP(I)
   50 FORMAT (5X,'ENCOUNTER',2I4,F9.2,F12.6,F10.5,F10.4,F12.6,1P5E10.1)
C
   60 IF (ICASE.EQ.1)  GO TO 70
C
C          CHECK MERGER CRITERION FOR STABLE BINARY (ICASE = 2).
      IF (SEMI(I).LT.0.0)  GO TO 100
      PERIOD = SEMI(I)*DSQRT (SEMI(I)/BCM)
C          DEFINE CLOSE BINARY IF PERIOD < 0.05 YEARS.
      IF (PERIOD.GT.0.05)  GO TO 100
      WRITE (6,65)  NAME(I),NAME(JCOMP),SEMI(I),PERI,ECC(I),PERIOD
   65 FORMAT (/,5X,'CLOSE BINARY','   NAME = ',2I5,'  SEMI =',
     &      1PE9.1,'  PERI =',E9.1,'  ECC =',0PF6.2,'  PERIOD =',F7.3)
      NSTEPN(8) = NSTEPN(8) + 1
      IF (NSTEPN(8).LT.2)  GO TO 100
C          RESET COUNTER AND MERGE BINARY AFTER SECOND TRY.
      NSTEPN(8) = 0
      GO TO 75
C
C          CHECK COLLISION CONDITION.
   70 IF (RIJ.GT.R(I) + R(JCOMP))  GO TO 100
 75   IF (BODY(I).GT.BODY(JCOMP)) THEN
         ICOMP = I
      ELSE
         ICOMP = JCOMP
         JCOMP = I
      ENDIF
C          ROUTINE MERGE ASSUMES ICOMP < JCOMP.
C
      EPAIR = -BODY(ICOMP)*BODY(JCOMP)/(2.0*SEMI(I))
C          SUBTRACT BINDING ENERGY TO MAINTAIN CONSERVATION.
      ETOT = ETOT - EPAIR
C          TIDAL CORRECTION DUE TO THE SUN IS PERFORMED IN ROUTINE MERGE.
      IF (EPAIR.LT.0.0)  NSTEPN(5) = NSTEPN(5) + 1
C
C          FORM SQUARE RADIAL VELOCITY FOR FRAGMENTATION CRITERION (6/5/93).
      VR2 = RDOT**2/RIJ2
C
C          SET ESCAPE VELOCITY & REBOUND VELOCITY.
      VESC = DSQRT (2.0*BCM/RIJ)
      VREB = CI*DSQRT (VR2)
C
C          CHECK OPTION FOR DIAGNOSTIC COLLISION OUTPUT.
      IF (KZ(5).LT.2)  GO TO 90
C
      YEARS = TIME/TWOPI
CC    GAMMA = DSQRT (GAMMA2)
      GAMMA = 0.0
      VCOLL = BCM*(2.0/(R(ICOMP) + R(JCOMP)) - 1.0/SEMI(I))
      IF (VCOLL.LT.0.0)  VCOLL = 0.0
      VCOLL = DSQRT (VCOLL)
      POT1 = 1.2*(BODY(ICOMP)/(BODY(JCOMP)*R(ICOMP)) +
     &                               BODY(JCOMP)/(BODY(ICOMP)*R(JCOMP)))
      POT2 = 2.0/(R(ICOMP) + R(JCOMP))
      VCRIT = DSQRT (BCM*(POT1 + POT2))
      RATIO = VCOLL/VCRIT
      AMASS1 = 1.001*BODY(ICOMP)/ZMOON
      AMASS2 = 1.001*BODY(JCOMP)/ZMOON
      AMASS3 = AMASS1 + AMASS2
      PERRI = PERI/(R(ICOMP) + R(JCOMP))
      ACRAT = VREB/VESC
      WRITE (6,85) NAME(ICOMP), NAME(JCOMP), YEARS
      WRITE (6,86) NAME(ICOMP), BODY(ICOMP), SEMI(ICOMP), ECC(ICOMP),
     &     (X(K,ICOMP),K=1,3), (XDOT(K,ICOMP),K=1,3)
      WRITE (6,86) NAME(JCOMP), BODY(JCOMP), SEMI(JCOMP), ECC(JCOMP),
     &     (X(K,JCOMP),K=1,3), (XDOT(K,JCOMP),K=1,3)
   85 FORMAT (3X, 'COLLISION', 2I10, F15.1)
 86   FORMAT (5X, 'SOURCE :', I5, 9E12.4)
C
C          DECIDE ON ACCRETION, CRATERING OR FRAGMENTATION.
   90 CALL EVENT
C
  100 RETURN
C
      END
