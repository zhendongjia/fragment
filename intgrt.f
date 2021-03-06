      SUBROUTINE INTGRT
C
C
C          N-BODY INTEGRATOR.
C          ------------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,BMAX,XF(NFMAX),VF(NFMAX),BF(NFMAX),
     &		  AI,NF,NC,IF(NFMAX),ICOMP,JCOMP
      REAL*8  BODYJ,FRAGM
      REAL*8  F2DOT(4),WK(10)
      REAL*8  F1(4),F1DOT(4),PDOT(3)
      COMMON/INTGC/ XI, YI, ZI, FRAGM, CHECK, FP(3), JMIN
C
C
C          DEFINE FRAGMENT MASS & CENTRAL PERTURBER BIN FOR FORCE LOOP.
      FRAGM = 0.5*ZMOON
      NBP1 = NBPERT + 1
C
C          FIND NEXT BODY TO BE TREATED AND SET NEW TIME.
    1 NNB = NLIST(1) + 1
      A1 = 1.0E+10
      DO 2 L = 2,NNB
      J = NLIST(L)
      A2 = T0(J) + STEP(J)
      IF (A2.GE.A1)  GO TO 2
      I = J
      A1 = A2
    2 CONTINUE
C
      TIME = A1
      IF (TIME.LT.TLIST)  GO TO 7
C
C          SET NEW TIME-STEP LIST AND MODIFY LIST INTERVAL.
    4 NNB = 1
      TLIST = TLIST + DTLIST
C
      DO 5 J = 1,N
      IF (T0(J) + STEP(J).GE.TLIST)  GO TO 5
      NNB = NNB + 1
      NLIST(NNB) = J
    5 CONTINUE
C
      IF (NNB.EQ.1)  GO TO 4
      NLIST(1) = NNB - 1
C
C          STABILIZE NLIST INTERVAL ON A MEMBERSHIP IN (0.5*N**0.5, N**0.5).
      IF (N.GT.10)  THEN
          IF (NLIST(1).LT.0.5*NSTABL)  DTLIST = 1.25*DTLIST
          IF (NLIST(1).GT.NSTABL)      DTLIST = 0.75*DTLIST
      ELSE
          DTLIST = 10.0*STEP(1)
      END IF
      GO TO 1

 7    DISTANCE2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
      IF (DISTANCE2.LT.CRIT_DISTANCE2 .OR. DISTANCE2.GT.ESCAPE_DISTANCE2) THEN
         CALL WRITE_OBJECT(I, 6, 'ESCAPE :')
         CALL REMOVE(I)
         GO TO 1
      END IF

      IBIN0 = ILIST(I)
      DT = STEP(I)
      T1PR = DT01(I)
      T2PR = DT02(I)
      T12PR = T1PR + T2PR
C          SET THE COMBINED SUN - PLANET MASS.
      SUNPL = 1.0 + BODY(I)
C
C         CHECK INDICATOR FOR BAUMGARTE-STIEFEL STABILIZATION.
      IF (ISTAB(I).EQ.0)  GO TO 9
C
C          NOTE THAT THE STABILIZER IS SWITCHED OFF IF GAMMA2 > DFMIN2.
      VELOC2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
      ENERGY = 0.5D0*VELOC2 - SUNPL/DSQRT (X0(1,I)**2 + X0(2,I)**2
     &                                                + X0(3,I)**2)
C          EXPLICIT ORBITAL ENERGY AT BEGINNING OF STEP.
      CORR = 0.2*(H(I) - ENERGY)/(DT*VELOC2)
      DO 8 K = 1,3
      F(K,I) = F(K,I) + CORR*X0DOT(K,I)
    8 CONTINUE
C
    9 DT06 = 0.6D0*DT
      DT19 = ONE9*DT
      DT12 = ONE12*DT
      DT34 = 0.75D0*DT
      DT32 = 1.5D0*DT
      DT20 = 2.0D0*DT  
C
C          OBTAIN CURRENT COORDINATES & VELOCITIES FOR BODY I TO THIRD ORDER.
      DO 10 K = 1,3
      F2DOTK = D3(K,I)*T12PR + D2(K,I)
      X(K,I) = ((((D3(K,I)*DT06 + F2DOTK)*DT12 + FDOT(K,I))*DT + 
     &     F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
      X0DOT(K,I) = (((D3(K,I)*DT34 + F2DOTK)*DT19 +
     &     FDOT(K,I))*DT32 + F(K,I))*DT20 + X0DOT(K,I) 
      FP(K) = 0.0
   10 CONTINUE
      CALL UPDATE_ORBIT(I)
C
      XI = X(1,I)
      YI = X(2,I)
      ZI = X(3,I)
      JMIN = 0
      CHECK = -1000000.0
      IF (KZ(11).EQ.0.AND.BODY(I).LT.EMBRYO)  GO TO 15
C
C          INCLUDE ALL PERTURBATIONS (DIRECT AND INDIRECT TERMS).
      DO J = 1,N
         IF (J.NE.I) CALL CALC_INTGRT(I,J)
      END DO
      GO TO 40
C
C          OBTAIN THE DOMINANT PLANETARY PERTURBATIONS (DIRECT TERMS ONLY).
   15 IBINJ = IBIN0 - NBPERT
C          CHECK BIN NUMBER NEAR 360 DEGREES.
      IF (IBINJ.LT.1)  IBINJ = NBTOT + IBINJ
      IBINR = LISTR(I)
C
C          LOOP OVER NEIGHBOURING PERTURBER BINS (NBMAX = 2*NBPERT + 1).
      DO 36 IDUM = 1,NBMAX
      NNB = LIST(1,IBINJ)
      IDIS = IABS (IDUM - NBP1)

      DO 34 L=2, NNB+1
      J = LIST(L,IBINJ)
      IF (J.EQ.I) CYCLE
C     AVOID DUPLICATE COMPUTION FOR PERTURBING PLANET
      IF (KZ(3).GT.0 .AND. KZ(17).GT.0 .AND. IS_JUPITER(J).EQ.1) CYCLE

C          INCLUDE ALL EMBRYOS INSIDE THE SPECIFIED PERTURBER BINS.
      BODYJ = BODY(J)
      IF (BODYJ.GT.EMBRYO) THEN
         IF (IDIS.GT.MPERT(J).OR.IABS (IBINR-LISTR(J)).GT.NRPERT) CYCLE
      ENDIF

      CALL CALC_INTGRT(I, J)

   34 CONTINUE
C
      IBINJ = IBINJ + 1
      IF (IBINJ.GT.NBTOT)  IBINJ = IBINJ - NBTOT
   36 CONTINUE

C     CALCULATE THE IMPACT OF THE PERTURBING PLANET.
      IF (KZ(3).GT.0 .AND. KZ(17).GT.0) THEN
         DO K=1,NJUPITER
            J = JUPITER(K)
            IF (I.EQ.J) CYCLE
            CALL CALC_INTGRT(I, J)
         END Do
      END IF

C     CALCULATE THE IMPACT OF GAS POTENTIAL.
 40   IF (KZ(18).GT.0) THEN
         CALL GAS_POTENTIAL(I, FP, PDOT)
      END IF
C
C     CALCULATE THE IMPACT OF GAS DAMPING.
      IF (KZ(19).GT.0.AND.
     &     SQRT(X(1,I)**2+X(2,I)**2+X(3,I)**2).LE.R_IN) THEN
         CALL GAS_DAMPING(I, FP, PDOT)
      END IF      
C
C          FORM THE PLANETARY PERTURBATIONS.
      PERT2 = FP(1)**2 + FP(2)**2 + FP(3)**2
      RI2 = XI**2 + YI**2 + ZI**2
      GAMMA2 = PERT2*RI2**2
C
C          ADD THE SUN - PLANET FORCE COMPONENT.
      FS = -SUNPL/(RI2*DSQRT (RI2))
      F1(1) = FP(1) + XI*FS
      F1(2) = FP(2) + YI*FS
      F1(3) = FP(3) + ZI*FS

C
C          SET TIME INTERVALS FOR CORRECTOR AND UPDATE THE BACKWARDS TIMES.
      DT1 = STEP(I) + DT01(I)
      DT2 = STEP(I) + DT02(I)
      DT3 = STEP(I) + DT03(I)
      T3PR = DT03(I)
      S1 = DT + DT1
      S2 = T1PR*T2PR
      S3 = S2*T3PR
      S4 = S2 + T3PR*T12PR
      S5 = T12PR + T3PR
      S6 = (((0.6666666666667D0*DT + S5)*DT06 + S4)*DT12 + ONE6*S3)*DT
      S7 = ((0.2D0*DT + 0.25D0*S5)*DT + ONE3*S4)*DT + 0.5D0*S3
      DT03(I) = DT02(I) + STEP(I)
      DT02(I) = DT01(I) + STEP(I)
      DT01(I) = STEP(I)
      T0(I) = TIME
      A1 = 1.0/DT
      A2 = 1.0/DT1
      A3 = 1.0/DT2
      A4 = DT*DT/DT3
C
C          FORM NEW DIFFERENCES AND INCLUDE FOURTH-ORDER CORRECTOR.
      DO 50 K = 1,3
      AK4 = (F1(K) - FI(K,I))*A1
      AK7 = (AK4 - D1(K,I))*A2
      AK10 = (AK7 - D2(K,I))*A3
      F4DOTK = (AK10 - D3(K,I))*A4
      FI(K,I) = F1(K)
      D1(K,I) = AK4
      D2(K,I) = AK7
      D3(K,I) = AK10
      X(K,I) = F4DOTK*S6 + X(K,I)
      X0(K,I) = X(K,I)
      X0DOT(K,I) = F4DOTK*S7 + X0DOT(K,I)
      XDOT(K,I) = X0DOT(K,I)
      F(K,I) = 0.5D0*F1(K)
      F1DOT(K) = (AK10*DT1 + AK7)*DT + AK4
      FDOT(K,I) = ONE6*F1DOT(K)
      F2DOT(K) = AK10*S1 + AK7
   50 CONTINUE
      CALL UPDATE_ORBIT(I)
C
C          CARRY ON NORMALLY FOR SMALL PERTURBATIONS.
      IF (GAMMA2.LT.DFMIN2.AND.ISTAB(I).GT.0)  GO TO 70
C
C          SKIP BINDING ENERGY INTEGRATION DURING LARGE PERTURBATIONS.
      IF (GAMMA2.GT.DFMIN2.AND.ISTAB(I).EQ.0)  GO TO 80
C
      IF (KZ(1).EQ.0)  GO TO 80
C
C          SEE WHETHER PERTURBATION PERMITS STABILIZATION PROCEDURE AGAIN.
      IF (GAMMA2.GT.DFMIN2)  GO TO 60
C
      ISTAB(I) = 1
      H(I) = 0.5D0*(X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2) -
     &                SUNPL/DSQRT (X0(1,I)**2 + X0(2,I)**2 + X0(3,I)**2)
      HDOT(I) = X0DOT(1,I)*FP(1) + X0DOT(2,I)*FP(2) + X0DOT(3,I)*FP(3)
      D1HDOT(I) = 0.0
      D2HDOT(I) = 0.0
      D3HDOT(I) = 0.0
      STEP(I) = 0.35*STEP(I)
C          THE NEXT TIME-STEP WILL BE LIMITED TO HALF THE CURRENT VALUE.
      GO TO 80
C
C          SWITCH OFF STABILIZATION FOR STRONG PERTURBATIONS.
   60 ISTAB(I) = 0
      HDOT(I) = 0.0
      D1HDOT(I) = 0.0
      D2HDOT(I) = 0.0
      D3HDOT(I) = 0.0
      NSTEPN(2) = NSTEPN(2) + 1
      GO TO 80
C
C          PREDICT THE CURRENT BINDING ENERGY (OMIT FACTORIALS).
   70 HDOT2 = (D3HDOT(I)*T2PR + D2HDOT(I))*T1PR + D1HDOT(I)
      HDOT3 = D3HDOT(I)*T12PR + D2HDOT(I)
      HDOT4 = D3HDOT(I)
      H(I) = (((0.25D0*HDOT4*DT + ONE3*HDOT3)*DT + 0.5D0*HDOT2)*DT +
     &                                                HDOT(I))*DT + H(I)
C
C          UPDATE ENERGY DIFFERENCES AND INCLUDE FOURTH-ORDER CORRECTOR.
      HDOT1 = X0DOT(1,I)*FP(1) + X0DOT(2,I)*FP(2) + X0DOT(3,I)*FP(3)
      HDOT2 = (HDOT1 - HDOT(I))*A1
      HDOT3 = (HDOT2 - D1HDOT(I))*A2
      HDOT4 = (HDOT3 - D2HDOT(I))*A3
      HDOT5 = (HDOT4 - D3HDOT(I))*A4
      HDOT(I) = HDOT1
      D1HDOT(I) = HDOT2
      D2HDOT(I) = HDOT3
      D3HDOT(I) = HDOT4
      H(I) = HDOT5*S7 + H(I)
C
C          DETERMINE NEW TIME-STEP.
   80 WK(1) = F1(1)**2 + F1(2)**2 + F1(3)**2
      WK(2) = F1DOT(1)**2 + F1DOT(2)**2 + F1DOT(3)**2
      WK(3) = 4.0*(F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2)
C          LEAVE OUT FACTORIAL AND SCALE BY 0.0000001 TO AVOID OVERFLOW.
      WK(6) = 0.0000001*D3(1,I)
      WK(7) = 0.0000001*D3(2,I)
      WK(8) = 0.0000001*D3(3,I)
      WK(4) = 60000000.0*SQRT (WK(6)**2 + WK(7)**2 + WK(8)**2)
C
C          STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2)**0.5.
      WK(1) = ETA*(SQRT(WK(1)*WK(3)) + WK(2))/(SQRT(WK(2))*WK(4) +WK(3))
      WK(1) = SQRT (WK(1))
      IF (ISTAB(I).GT.0)  WK(1) = WK(1)/(1.0 + 4.0*GAMMA2/DFMIN2)
C
C          REDUCE TIME-STEP BY FACTOR 4 DURING NON-STABILIZED INTEGRATION.
      IF (GAMMA2.GT.DFMIN2.AND.ISTAB(I).EQ.0)  WK(1) = 0.25*WK(1)
C
C          RESTRICT INCREASE BY INERTIAL FACTOR 1.4.
      IF (WK(1).GT.1.4*STEP(I))  THEN
          STEP(I) = 1.4*STEP(I)
      ELSE
          STEP(I) = WK(1)
      END IF
C
C          APPLY KINEMATIC TIME-STEP CRITERION FOR APPROACHING LOW-MASS BODY.
      IF (JMIN.GT.0)  THEN
          STEP(I) = MIN (0.25D0*DABS (CHECK),STEP(I))
      END IF
      IF (STEP(I).LT.CRIT_STEP) STEP(I) = CRIT_STEP 
C
C          OBTAIN CURRENT PHASE ANGLE AND PERTURBER BINS.
      WK(1) = X0(1,I)
      WK(2) = X0(2,I)
      THETA = ATAN2 (WK(2),WK(1))
      IF (THETA.LT.0.0)  THETA = THETA + TWOPI
      IBIN = 1 + THETA/ALPHA
      LISTR(I) = ROUT2/RI2
      IF (IBIN.EQ.IBIN0)  GO TO 100
C
C          UPDATE THE OLD AND NEW PERTURBER LIST.
      ILIST(I) = IBIN
C          FIRST REMOVE BODY #I FROM THE OLD LIST.
      NNB0 = LIST(1,IBIN0)
      L = 1
   90 L = L + 1
      IF (L.GT.NNB0 + 1)  GO TO 95
      IF (LIST(L,IBIN0).NE.I)  GO TO 90
C
C          MOVE UP ANY REMAINING LIST MEMBERS AND UPDATE THE MEMBERSHIP.
      DO 91 K = L,NNB0
   91 LIST(K,IBIN0) = LIST(K+1,IBIN0)
   95 LIST(1,IBIN0) = NNB0 - 1
C          INCLUDE BODY #I IN THE NEW LIST.
      NNB = LIST(1,IBIN)
      LIST(NNB+2,IBIN) = I
      LIST(1,IBIN) = NNB + 1
C
C          CHECK THE PERTURBATION FOR POSSIBLE CLOSE ENCOUNTER.
  100 IF (GAMMA2.LT.DFMIN2)  GO TO 200
C
C          COMPARE DOMINANT FORCE WITH THE SCALED COLLISION FORCE.
      IF (PERT2*BODY(I)**2.LT.DFMAX2)  GO TO 200
C
C          ASSUME MAXIMUM RADIUS TO ENSURE COLLISION DETECTION.
      RCL2 = 4.0*(R(I) + RMAX)**2
C
C          PERFORM COLLISION SEARCH AND IMPLEMENT MERGER OR FRAGMENTATION.
      CALL SEARCH(I,1,IBIN,RCL2)
C
  200 NSTEPN(1) = NSTEPN(1) + 1	
      NSTEPS = NSTEPS + 1
      NTIMER = NTIMER + 1
C
C          CHECK NEXT OUTPUT TIME.
      IF (TIME.GT.TNEXT)  RETURN
C
      IF (NTIMER.LT.1000)  GO TO 215
      NTIMER = 0
C
      OPEN(99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          WRITE (6,210)
  210     FORMAT (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
C
      CALL CPUTIM(TCOMP)
C          ELAPSED COMPUTING TIME IN MINUTES.
      IF (TCOMP.LT.CPU)  GO TO 1
C
      IF (KZ(2).NE.0)  THEN
          CPUTOT = CPUTOT + TCOMP
          CALL COMMON(1,1)
          YEARS = TIME/TWOPI
          WRITE (6,212)  TIME/TWOPI,TCOMP,CPUTOT/60.0
  212     FORMAT (//,9X,'COMMON SAVED AT YEARS =',F10.1,5X,
     &                          'TCOMP =',F8.1,5X,'CPUTOT =',F7.1)
      END IF
C
      STOP
C
  215 IF (NSTEPS.LT.100000)  GO TO 1
      NSTEPS = 0
C
C          PROCEDURE FOR MERGING STABLE CLOSE BINARY (M > ZMOON).
      I = 1
  220 IF (STEP(I).GT.0.0001*TWOPI.OR.BODY(I).LT.ZMOON)  GO TO 225
      IBIN = ILIST(I)
      RCL2 = 1600.0*R(I)**2
C
C          PERFORM CLOSE BINARY SEARCH INSIDE 40*R.
      CALL SEARCH(I,2,IBIN,RCL2)
C
  225 I = I + 1
      IF (I.LT.N)  GO TO 220
C
      GO TO 1
C
      END


      SUBROUTINE CALC_INTGRT(I, J)
      INCLUDE 'COMMONP.FOR'
      COMMON/INTGC/ XI, YI, ZI, FRAGM, CHECK, FP(3), JMIN
C
C          PREDICT RELATIVE COORDINATES OF BODY #J WITH RESPECT TO BODY #I.
   32 S = TIME - T0(J)
      XIJ = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J) - XI
      YIJ = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J) - YI
      ZIJ = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J) - ZI
C
C          SUM THE PLANETARY FORCE CONTRIBUTIONS. 
      RIJ2 = XIJ**2 + YIJ**2 + ZIJ**2
      FIJ = BODY(J)/(RIJ2*SQRT (RIJ2))
C
      FP(1) = FP(1) + XIJ*FIJ
      FP(2) = FP(2) + YIJ*FIJ
      FP(3) = FP(3) + ZIJ*FIJ
C
C          CHECK APPROACHING FRAGMENT (IGNORE VERTICAL MOTION).
      IF (BODY(J).GT.FRAGM)  RETURN
      RDOT = XIJ*(X0DOT(1,J)-X0DOT(1,I)) + YIJ*(X0DOT(2,J)-X0DOT(2,I)) 
      IF (RDOT.GT.0.0)  RETURN
      DELT = RIJ2/RDOT
      IF (DELT.LT.CHECK) RETURN
      CHECK = DELT
      JMIN = J

      END SUBROUTINE
