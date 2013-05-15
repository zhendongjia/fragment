      SUBROUTINE EVENT
C
C
C          COLLISION EVENT.
C          ----------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(8),VF(8),BF(8),AI,
     &		  NF,IF(10),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
      REAL*4  VREL(3)
C
C
C          FIRST CHECK ACCRETION CRITERION.
      IF (VREB.LT.VESC)  GO TO 90
C
C          ACCRETE IF SMALLEST BODY <= 0.01 LUNAR MASS.
      BMIN = MIN (BODY(ICOMP),BODY(JCOMP))
      IF (BMIN.LE.0.01*ZMOON)  GO TO 90
C
C          FORM C.M. COORDINATES & VELOCITIES AND SAVE THE RELATIVE VELOCITY.
      DO 10 K = 1,3
      XCM(K) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/BCM
      VCM(K) = (BODY(ICOMP)*XDOT(K,ICOMP)+BODY(JCOMP)*XDOT(K,JCOMP))/BCM
      VREL(K) = XDOT(K,JCOMP) - XDOT(K,ICOMP)
   10 CONTINUE
C
C          INITIALIZE SCALARS AND SAVE ORIGINAL MASSES.
      NF = 0
      NC = 0
      DMC = 0.0
      BODYI = BODY(ICOMP)
      BODYJ = BODY(JCOMP)
C
C          SET PRE-COLLISION KINETIC ENERGY & POST-COLLISION TOTAL ENERGY.
      ZMU = BODYI*BODYJ/BCM
      ER = 0.5*ZMU*VR2
      EGRAV = 0.5*ZMU*(VREB**2 - VESC**2) - 
     &                       0.6*(BODYI**2/R(ICOMP) + BODYJ**2/R(JCOMP))
C
C          TREAT EACH COMPONENT IN TURN.
      I = ICOMP
C          DEFINE RATIO OF IMPACT STRENGTH OVER IMPACT ENERGY.
   20 AI = BSTAR*SI*R(I)**3/ER
      IF (I.EQ.ICOMP)  THEN
          WRITE (6,25)  NAME(ICOMP),NAME(JCOMP),AI,VREB,VESC,ER,EGRAV
   25     FORMAT (3X'DIAG  ',2I4,1P5E9.1)
      END IF
C
C          DECIDE BETWEEN CRATERING (AI > 1) OR FRAGMENTATION (AI < 1).
   30 IF (AI.GE.1.0)  THEN
C          SET NET MASS TRANSFER DURING CRATER FORMATION.
          DMC = 0.05*CEJ*ZKM*(VR2/VESC**2.25)*(BODYI - BODYJ)
C          MERGE IF ONE OF THE NEW MASSES < 0.01 LUNAR MASSES.
          IF (DABS (DMC).GE.BMIN - 0.01*ZMOON)  GO TO 50
C          MERGE IF MODIFIED REBOUND VELOCITY LESS THAN ESCAPE VELOCITY.
          IF (CII*VREB.LT.CI*VESC)  GO TO 50
          GO TO 70
      END IF
C
C          REVERT TO CRATERING IF NOT SUFFICIENT ENERGY FOR FRAGMENTATION.
CC      IF (EGRAV.LT.0.0)  THEN
CC          AI = 1.0/AI
CC          GO TO 30
CC      END IF
C
C          CALCULATE EFFECTIVE MASS LOSS FOR FRAGMENTATION.
      ZMEJ = CEJ*ZKM*BODYJ*(VR2/VESC**2.25)
      IF (ZMEJ.GT.0.9*BODY(I))  ZMEJ = 0.9*BODY(I)
C
C          FORM FRAGMENTS IF MASS LOSS > 0.05 LUNAR MASSES.
      IF (ZMEJ.GE.0.05*ZMOON)  THEN
          CALL FRAGMT(I)
      ELSE
C          REPLACE CASE OF SMALL EJECTED MASS BY CRATERING.
          AI = 1.0/AI
          GO TO 30
      END IF
C
C          TREAT THE SECOND BODY SIMILARLY.
   40 IF (I.EQ.ICOMP)  THEN
          I = JCOMP
          BODYJ = BODYI
          GO TO 20 
      END IF

C          SET ORBITAL CHARACTERISTICS OF FRAGMENTS AND MERGE CORES.
   50 IF (NF.GT.0)  THEN
          CALL INFRAG
          GO TO 90
      END IF
C
C          MERGE CORES IF NO CRATERING.
      IF (NC.EQ.0)  GO TO 90
C
C          MODIFY THE MASSES CAUSED BY CRATERING.
      BODY(ICOMP) = BODY(ICOMP) + DMC
      BODY(JCOMP) = BODY(JCOMP) - DMC
C
C          MERGE BODIES IF SMALLEST NEW MASS <= 0.01*ZMOON.
      IF (MIN (BODY(ICOMP),BODY(JCOMP)).LE.0.01*ZMOON)  GO TO 90
C      
C          ENFORCE MERGER IF SEMI-MAJOR AXIS IS INSIDE SUM OF RADII.
      VREL2 = VREL(1)**2 + VREL(2)**2 + VREL(3)**2
      SEMI = (VESC**2 - CII**2*VREL2)/BCM
      SEMI = 1.0/SEMI
      IF (SEMI.GT.0.0.AND.SEMI.LT.R(ICOMP) + R(JCOMP))  THEN
          WRITE (6,60)  ICOMP,JCOMP,SEMI,R(ICOMP)+R(JCOMP)
  60      FORMAT (3X,'  ENFORCED MERGER ',2I5,1P2E10.1)
          GO TO 90
      END IF
C
C          PERFORM APPROXIMATE ENERGY CORRECTION (IGNORE MASS EXCHANGE).
      ETOT = ETOT + 0.5*ZMU*VREL2*(CII**2 - 1.0)
C
C          FORM NEW FORCE POLYNOMIAL AND TIME-STEP OF CRATERED BODIES.
      CALL FPOLY(ICOMP)
      CALL FPOLY(JCOMP)
C
      NSTEPN(10) = NSTEPN(10) + 1
      GO TO 100
C
C	   CRATER FORMATION ON BODY #I.
   70 NC = NC + 1
C          MODIFY THE VELOCITY FOR REBOUND-TYPE COLLISION WITH ENERGY LOSS.
      DO 75 K = 1,3
      IF (I.EQ.JCOMP)  VREL(K) = -VREL(K)
      XDOT(K,I) = VCM(K) + CII*BODYJ*VREL(K)/BCM
   75 CONTINUE
C
C          PRINT ORBITAL ELEMENTS OF THE CRATERED BODY.
      RI = DSQRT (X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      SEMI = 2.0/RI - VI2/(1.0 + BODY(I))
      SEMI = 1.0/SEMI
      RDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
      ECC = DSQRT ((1.0 - RI/SEMI)**2 + RDOT**2/SEMI)
      WRITE (6,80)  NAME(I),DABS (DMC)/ZMOON,SEMI,ECC,RI
   80 FORMAT (3X,'  CRATER :',I5,F9.4,3F8.3)
      GO TO 40
C   
C          MERGE BODIES (INELASTIC, MINOR SECONDARY OR TWO CORE FRAGMENTS). 
   90 CALL MERGE(ICOMP,JCOMP)
C
  100 RETURN
C
      END
