      SUBROUTINE OUTPUT
C
C
C          DATA ANALYSIS AND OUTPUT.
C          -------------------------
C
      INCLUDE 'COMMONP.FOR'
      REAL*8  A(12)
C
C
C          OBTAIN CURRENT COORDINATES AND VELOCITIES TO ORDER F3DOT.
      IF (TIME.EQ.0.0)  GO TO 3
      DO 2 I = 1,N
      DT = TIME - T0(I)
      A(1) = 0.05D0*DT
      A(2) = 0.25D0*DT
      A(3) = (T0(I) - T1(I)) + (T0(I) - T2(I))
      DO 1 K = 1,3
      F2DOTK = D3(K,I)*A(3) + D2(K,I)
      X(K,I) = ((((D3(K,I)*A(1) + ONE12*F2DOTK)*DT + FDOT(K,I))*DT +
     &                             F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
      XDOT(K,I) =  (((D3(K,I)*A(2) + ONE3*F2DOTK)*DT +
     &               3.0D0*FDOT(K,I))*DT + 2.0D0*F(K,I))*DT + X0DOT(K,I)
    1 CONTINUE
      CALL UPDATE_ORBIT(I)
    2 CONTINUE
    3 A(1) = 0.0
C          FIRST OBTAIN THE POTENTIAL ENERGY OF THE PROTO-PLANETS.
      I = 1
    4 JMIN = I + 1
      A(2) = 0.0
      DO 5 J = JMIN,N
      A(5) = X(1,I) - X(1,J)
      A(6) = X(2,I) - X(2,J)
      A(7) = X(3,I) - X(3,J)
      A(2) = A(2) + BODY(J)/DSQRT (A(5)**2 + A(6)**2 + A(7)**2)
    5 CONTINUE
      A(1) = A(1) + BODY(I)*A(2)
      I = I + 1
      IF (I.LT.N)  GO TO 4
    6 DO 7 K = 3,10
    7 A(K) = 0.0
      AZ = 0.0
C          NEXT LOOP IS KINETIC ENERGY, INTERACTION TERMS & ANGULAR MOMENTUM.
      DO 9 I = 1,N
      A(3) = A(3) + 0.5D0*BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                            XDOT(3,I)**2)
      A(4) = A(4) + BODY(I)/DSQRT (X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      A(6) = A(6) + BODY(I)
      AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
      DO 8 K = 1,2
      A(K+6) = A(K+6) + BODY(I)*X(K,I)
      A(K+8) = A(K+8) + BODY(I)*XDOT(K,I)
    8 CONTINUE
    9 CONTINUE
      AZ = AZ - (A(7)*A(10) - A(8)*A(9))/(A(6) + 1.0)
C          ANGULAR MOMENTUM IN THE INERTIAL FRAME.
      A(5) = A(3) - A(4) - A(1)
      A(5) = A(5) - 0.5*(A(9)**2 + A(10)**2)/(A(6) + 1.0)
C          TOTAL ENERGY CORRECTED FOR INDIRECT TERMS (ADDED 6 FEBRUARY 1981).
      IF (TIME.EQ.0.0)  ETOT = A(5)
      ERROR = (A(5) - ETOT)/A(5)
      ETOT = A(5)
      YEARS = TIME/TWOPI
C          OBTAIN ELAPSED CPU TIME & UPDATE TOTAL SINCE LAST OUTPUT/RESTART.
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
      WRITE (6,10)  YEARS,TCOMP,NSTEPN(1),NSTEPN(2),NSTEPN(4),
     &                                      NSTEPN(9),NSTEPN(10),N,ERROR
   10 FORMAT (//,'  T =',F9.1,'  CPU =',F7.1,'  # =',I10,'  MOD =',I7,
     &            '  MERG =',I4,'  FRAG =',I4,'  CRAT =',I4,'  N =',I4,
     &                                                '  DE/E =',1PE9.1)
C
      ZH = 0.0
      ZHM = 0.0
      ZIMAX = 0.0
      VZ = 0.0
      VZM = 0.0
      EZ = 0.0
      EZM = 0.0
      ZMASS = 0.0
      DO 20 I = 1,N
      ZH = ZH + X(3,I)**2
      ZHM = ZHM + BODY(I)*X(3,I)**2
      ZIMAX = MAX (DABS (X(3,I)),ZIMAX)
      VZ = VZ + XDOT(3,I)**2
      VZM = VZM + BODY(I)*XDOT(3,I)**2
      EZ = EZ + 0.5*(X(3,I)**2 + XDOT(3,I)**2)
      EZM = EZM + 0.5*BODY(I)*(X(3,I)**2 + XDOT(3,I)**2)
      ZMASS = ZMASS + BODY(I)
   20 CONTINUE
      ZN = FLOAT (N)
      ZH = DSQRT (ZH/ZN)
      ZHM = DSQRT (ZHM/ZMASS)
      VZ = DSQRT (VZ/ZN)
      VZM = DSQRT (VZM/ZMASS)
      EZ = DSQRT (EZ/ZN)
      EZM = DSQRT (EZM/ZMASS)
      WRITE (6,30)  ZH,ZHM,ZIMAX,VZ,VZM,EZ,EZM
   30 FORMAT (/'  <Z> =',1PE8.1,'  <ZM> =',E8.1,'  ZMAX =',E8.1,
     &                 '  <VZ> =',E8.1,'  <VZM> =',E8.1,'  <EZ> =',E8.1,
     &                 '  <EZM> =',E8.1)
      IF (KZ(7).GT.0)  CALL HISTOG(AZ)
      IF (KZ(6).EQ.0)  GO TO 160
      DUMMY1 = 0.0
      DUMMY2 = 0.0
      WRITE (3,50)  N, YEARS, SCALE, DUMMY1, DUMMY2
   50 FORMAT(/'N T SCALE DUMMY',I10, F15.0, E10.1,2F8.1)
      DO 55  J=1,N
                WRITE (3,60) NAME(J),BODY(J),(X(K,J),K=1,3),
     &        (XDOT(K,J),K=1,3), SEMI(J), ECC(J)
   60 FORMAT(/'M R V', I10, 7E10.1, 2F10.3)  
   55 CONTINUE
  160 ISCALE = 1
      IF (KZ(12).EQ.1)  ISCALE = FLOAT (NSTEPN(4) + N)/(0.999*FLOAT(N))
C          OUTPUT INTERVAL INCREASED BY INTEGER PART OF N0/N IF KZ(12) = 1.
      TNEXT = TNEXT + FLOAT (ISCALE)*DELTAT
      IF (KZ(12).GT.1)  DELTAT = DELTAT*(1.0 + 0.01*KZ(12))
C          OUTPUT TIME INTERVAL INCREASED BY LOGARITHMIC FACTOR IF KZ(12) > 1.
      NSTABL = SQRT (FLOAT (N))
C          SET MAXIMUM RADIUS FOR COLLISION SEARCH PROCEDURE.
      RMAX = 0.0
      IF (N.GE.20)  THEN
          J3 = N/3
      ELSE
          J3 = N
      END IF
      DO 170 J = 1,J3
      IF (RMAX.LT.R(J))  RMAX = R(J)
  170 CONTINUE
C
      IF (KZ(2).GT.1.)  CALL COMMON(1,2)
      IF (TIME.LT.TCRIT)  RETURN
C
      IF (KZ(2).NE.0)  CALL COMMON(1,1)
      WRITE (6,200)
  200 FORMAT (//,10X,'END OF RUN')
C
      STOP
C
      END
