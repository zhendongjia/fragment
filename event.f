      SUBROUTINE EVENT
C
C
C          COLLISION EVENT.
C          ----------------
C
C      IMPLICIT NONE
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(NFMAX),VF(NFMAX),BF(NFMAX),
     &                   AI,NF,NC,IF(NFMAX),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
      REAL*8  VREL(3), ZMU, ER
      INTEGER  I, TAR, PRO
      REAL*8 B_CRIT, B, XR2, XDR2, V_STAR, V_PRO, V_IMPACT, QR, L_LAP, L_ALPHA, MIU_ALPHA, M_TOT1, V_ESC1, M_TOT, V_ESC, ITA, C1, C2, C3, C4, V_HR
C
C          FORM C.M. COORDINATES & VELOCITIES AND SAVE THE RELATIVE VELOCITY.
      DO 10 K = 1,3
      XCM(K) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/BCM
      VCM(K) = (BODY(ICOMP)*XDOT(K,ICOMP)+BODY(JCOMP)*XDOT(K,JCOMP))/BCM
   10 CONTINUE
C
C     MAKE SURE THE TARGET AND PROJECTILE
      IF(BODY(ICOMP).GE.BODY(JCOMP)) THEN 
         TAR = ICOMP
         PRO = JCOMP
      ELSE 
         TAR = JCOMP
         PRO = ICOMP
       END IF
C
C     CALCULATE THE IMPACT PARAMETERS
      B_CRIT = R(TAR)/(R(TAR) + R(PRO))
      B = (X(1,ICOMP) - X(1,JCOMP))*(XDOT(1,ICOMP) - XDOT(1,JCOMP))
     &     + (X(2,ICOMP) - X(2,JCOMP))*(XDOT(2,ICOMP) - XDOT(2,JCOMP))
     &     + (X(3,ICOMP) - X(3,JCOMP))*(XDOT(3,ICOMP) - XDOT(3,JCOMP))
      XR2 = (X(1,ICOMP) - X(1,JCOMP))**2 + (X(2,ICOMP) - X(2,JCOMP))**2
     &     + (X(3,ICOMP) - X(3,JCOMP))**2
      XDR2 = (XDOT(1,ICOMP) - XDOT(1,JCOMP))**2 
     &     +   (XDOT(2,ICOMP) - XDOT(2,JCOMP))**2   
     &     +   (XDOT(3,ICOMP) - XDOT(3,JCOMP))**2
      B = B/(SQRT(XR2)*SQRT(XDR2))
      V_TAR = (XDOT(1,TAR)-VCM(1))**2 + (XDOT(2,TAR)-VCM(2))**2
     &     + (XDOT(3,TAR)-VCM(3))**2
      V_PRO = (XDOT(1,PRO)-VCM(1))**2 + (XDOT(2,PRO)-VCM(2))**2 
     &     + (XDOT(3,PRO)-VCM(3))**2
      V_IMPACT = SQRT(XDR2)
      QR = (BODY(TAR)*V_TAR + BODY(PRO)*V_PRO)/(2*(BODY(TAR)+BODY(PRO)))
      L_LAP = (R(TAR)+R(PRO))*(1-B)
      IF(R(TAR).GT.L_LAP+R(PRO)) THEN
         L_ALPHA = (3*R(PRO)*L_LAP**2 - L_LAP**3)/(4*R(PRO)**3)
      ELSE 
         L_ALPHA = 1
      ENDIF
      MIU_ALPHA = L_ALPHA*BODY(PRO)*BODY(TAR)/(L_ALPHA*BODY(PRO) + BODY(TAR))
C
      M_TOT1 = L_ALPHA*BODY(PRO) + BODY(TAR)
      V_ESC1 = SQRT(2*M_TOT1/(R(TAR)+R(PRO)))
      M_TOT = BODY(PRO) + BODY(TAR)
      V_ESC = SQRT(2*M_TOT/(R(TAR)+R(PRO)))
      ITA = (BODY(TAR) - BODY(PRO))/M_TOT
      C1 = 2.43
      C2 = -0.0408
      C3 = 1.86
      C4 = 1.08
      V_HR = V_ESC*(C1*ITA**2*(1-B)**2.5 + C2*ITA**2 + C3*(1-B)**2.5 + C4)
C
C          SET PRE-COLLISION KINETIC ENERGY & POST-COLLISION TOTAL ENERGY.
      ZMU = BODY(TAR)*BODY(PRO)/(BODY(TAR)+BODY(PRO))
      ER = 0.5*ZMU*VR2
      EGRAV = 0.5*ZMU*(VREB**2 - VESC**2) - 
     &                       0.6*(BODYI**2/R(ICOMP) + BODYJ**2/R(JCOMP))
C
      WRITE (6,25)  NAME(ICOMP),NAME(JCOMP),V_IMPACT,V_ESC1,V_HR 
   25     FORMAT (3X'DIAG  ',2I4,3F10.4)
C
C     FIRST CHECK THE PERFECT MERGE.
      IF (V_IMPACT.LT.V_ESC1)  THEN 
          CALL MERGE(TAR,PRO)
          GO TO 100
      END IF
C
C          DECIDE BETWEEN CRATERING OR FRAGMENTATION.
      NF = 0
      NC = 0
       IF (B.LT.B_CRIT) THEN
           CALL FRAGMT(TAR, PRO, B, L_LAP, MIU_ALPHA, QR)
           CALL REMOVE(TAR)
           IF (PRO.GT.TAR) PRO = PRO - 1
           CALL REMOVE(PRO)
       ELSE
            IF (V_IMPACT.LT.V_HR) THEN
               CALL MERGE(TAR, PRO)
               GO TO 100
            END IF
            NC = 1
            CALL CRATER(TAR, PRO, B, L_LAP, MIU_ALPHA, QR)
            CALL REMOVE(PRO)
       END IF
C
C
  100 RETURN
C
      END
