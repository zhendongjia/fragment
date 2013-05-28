      SUBROUTINE INPUT
C
C
C          INITIAL CONDITIONS.
C          -------------------
C
      INCLUDE 'COMMONP.FOR'
      REAL*4  RAN2,LSTAR,VSTAR,MSTAR
C
C
C                  INPUT PARAMETERS
C                  ****************
C
C          ------------------------------------------------------------------
C          KSTART  VALUES 1,2,3 FOR NEW CASE, RESTART, OR MODIFIED RESTART.
C          TCOMP   MAXIMUM COMPUTING TIME IN MINUTES.
C
C          N       TOTAL NUMBER OF PROTO-PLANETS.
C          NBTOT   TOTAL NUMBER OF NEIGHBOUR BINS.
C          NBPERT  PERTURBER BINS ON EITHER SIDE (TOTAL USED = 2*NBPERT + 1).
C          NRPERT  RADIAL PERTURBER BINS ON EITHER SIDE.
C          NRAND   RANDOM NUMBER SEQUENCE SKIP.
C          NRUN    RUN IDENTIFICATION INDEX.
C
C          ETA     ACCURACY FACTOR FOR TIME-STEPS (COMPOSITE CRITERION).
C          DELTAT  OUTPUT INTERVAL IN YEARS.
C          TCRIT   TERMINATION TIME IN YEARS.
C          DFMIN   RELATIVE PERTURBATION FOR SWITCHING OFF STABILIZATION.
C          RIN     RADIUS OF INNER BOUNDARY (UNITS OF AU).
C          ROUT    RADIUS OF OUTER BOUNDARY.
C          SCALE   MASS FACTOR FOR SCALING TO OTHER INITIAL MASSES AND RADII.
C          ZMAX    MAXIMUM HALF-THICKNESS (IN ROCHE RADII).
C
C          KZ      NON-ZERO OPTIONS FOR ALTERNATIVE PATHS (SEE TABLE).
C
C          CI      COEFFICIENT OF RESTITUTION.
C          CII     MODIFIED COEFFICIENT OF RESTITUTION.
C          VC      MAXIMUM REBOUND VELOCITY.
C          ZKM     MASS EXCAVATED COEFFICIENT.
C          SI      IMPACT STRENGTH.
C          CEJ     EJECTA VELOCITY COEFFICIENT.
C
C          ZM1     LOWER MASS LIMIT (POWER LAW WITH OPTION 13).
C          ZM2     UPPER MASS LIMIT (UNITS OF ZMOON).
C
C          RSCALE  INITIAL SCALING FACTOR FOR RADIUS (OPTION 15).
C          ------------------------------------------------------------------
C
C
C                  SECONDARY PARAMETERS
C                  ********************
C
C          ------------------------------------------------------------------
C          ALPHA1  INVERSE ANGLE OF A PERTURBER BIN IN RADIANS (NBTOT/TWOPI).
C          DEGREE  ONE RADIAN IN DEGREES.
C          DFMIN2  SQUARE OF STABILIZATION PERTURBATION DFMIN.
C          DFMAX2  SQUARE OF THE SCALED COLLISION FORCE DFMAX.
C          EMBRYO  MASS OF EMBRYO (4*ZMOON).
C          NBMAX   TOTAL NUMBER OF PERTURBER BINS (2*NBPERT + 1).
C          NEMB    NUMBER OF EMBRYOS.
C          ROCHE   ROCHE RADIUS (AU).
C          ROUT2   SQUARE RADIUS OF OUTER BOUNDARY.
C          TWOPI   TWO PI (8.0*DATAN (1.0D0)).
C          ZMOON   LUNAR MASS IN SOLAR MASSES (SCALING INCLUDED).
C          ------------------------------------------------------------------
C
C
C                  OUTPUT COUNTERS NSTEPN(J)
C                  *************************
C
C          ------------------------------------------------------------------
C          1  INTEGRATION STEPS.
C          2  STABILIZATION TERMINATIONS.
C          3  CLOSE ENCOUNTERS INSIDE 3*R.
C          4  MERGERS.
C          5  ELLIPTIC MERGERS.
C          7  LARGE PERTURBATION WITH NO BODY INSIDE 3*R.
C          8  MERGER OF CLOSE BINARY ON SECOND TIME CHECK TRY.
C          9  FRAGMENTATION. 
C         10  CRATERING.
C          ------------------------------------------------------------------
C
C
C                  OPTIONS KZ(J)
C                  *************
C
C          ------------------------------------------------------------------
C          1  BAUMGARTE-STIEFEL STABILIZATION.
C          2  COMMON SAVE ON UNIT 1 (=2: EVERY 100000 STEPS & OUTPUT ON # 2).
C          3  READ DATA FOR PERTURBING PLANET.
C          4  INITIAL DENSITY (=0: CONST; =1: -1/2; =2: -3/2; =3: -2).
C          5  MERGER OUTPUT (= 2: COLLISION; = 3: CLOSE ENCOUNTERS).
C          6  DATA STORED ON UNIT 3 EVERY OUTPUT.
C          7  DATA ANALYSIS IN ROUTINE HISTOG.
C          8  OUTPUT OF ECCENTRICITY VS MASS & ECCENTRICITY DISTRIBUTION.
C          9  OUTPUT OF INDIVIDUAL BODIES (UP TO I = KZ(9)).
C         10  OUTPUT OF DISPERSIONS AND THE TWO INTEGRALS OF MOTION.
C         11  ALL PERTURBATIONS INCLUDED (= 0: FULL N SUMMATION FOR EMBRYOS).
C         12  EXPAND OUTPUT TIME (=1: BY INT (N0/N); >1: BY (1+0.01*KZ(12)).
C         13  INITIAL MASS FUNCTION N(M) = M**(-3/2) IN RANGE (ZM1,ZM2).
C         14  NON-CIRCULAR ORBIT FOR INITIAL PROTOPLANETS (0.01*KZ(14) USED).
C         15  SCALING OF INITIAL PARTICLE SIZES (READ SCALING FACTOR).
C         16  ITERATION PROCEDURE & OUTPUT FOR Q & L DISTRIBUTION.
C         17  ALWAYS CONSIDER THE IMPACT OF THE PERTURBING PLANET.
C          ------------------------------------------------------------------
C
C
C          INITIALIZE PARAMETERS, COUNTERS & SET USEFUL CONSTANTS.
      TIME = 0.0
      TNEXT = 0.0
      CPUTOT = 0.0
      BINARY = 0.0
      NSTEPS = 0
      NTIMER = 0
      NEMB = 0
      LISTC(1) = 0
      DO 1 K = 1,20
    1 NSTEPN(K) = 0
      DO 2 K = 1,NMAX
    2 LIST(1,K) = 0
C
      ONE3 = 1.0/3.0D0
      ONE6 = 1.0/6.0D0
      ONE9 = 1.0/9.0D0
      ONE12 = 1.0/12.0D0
      TWOPI = 8.0D0*DATAN (1.0D0)
      DEGREE = 360.0D0/TWOPI
C
C          READ & PRINT THE MAIN INPUT PARAMETERS.
      READ (5,*)  N,NBTOT,NBPERT,NRPERT,NRAND,NRUN
      READ (5,*)  ETA,DELTAT,TCRIT,DFMIN,RIN,ROUT,SCALE,ZMAX
      READ (5,*)  (KZ(J),J=1,20)
      READ (5,*)  CI,CII,VC,ZKM,SI,CEJ
C
      WRITE (6,5)
    5 FORMAT (/////,7X,'N  NBTOT  NBPERT  NRPERT  NRAND  NRUN')
      WRITE (6,6)  N,NBTOT,NBPERT,NRPERT,NRAND,NRUN
    6 FORMAT (/,I8,I7,2I8,I7,I6)
      WRITE (6,8)
    8 FORMAT (//,7X,'ETA        DELTAT     TCRIT   DFMIN   RIN   ROUT',
     &                                             '   SCALE     ZMAX')
      WRITE (6,9)  ETA,DELTAT,TCRIT,DFMIN,RIN,ROUT,SCALE,ZMAX
    9 FORMAT (/,F13.4,F11.1,F10.0,F8.3,F6.1,F7.1,F8.1,F9.1)
      WRITE (6,10)  (J,J=1,20)
   10 FORMAT (//,7X,'OPTIONS   ',20I4)
      WRITE (6,11)  (KZ(J),J=1,20)
   11 FORMAT (/,17X,20I4)
      WRITE (6,12)
   12 FORMAT (//,7X,'CI   CII   VC        ZKM       SI        CEJ')
      WRITE (6,13)  CI,CII,VC,ZKM,SI,CEJ
   13 FORMAT (/,F9.1,F6.1,1P4E10.1,/)
C
      LSTAR = 1.5E+13
      VSTAR = 3.0E+06
      MSTAR = 2.0E+33
      BSTAR = 2.0*TWOPI/3.0*(LSTAR/VSTAR)**2*(LSTAR/MSTAR)
      IF (KZ(15).NE.0)  THEN
          READ  (5,*)  RSCALE
          WRITE (6,14)  RSCALE
   14     FORMAT (/,5X,'INITIAL SCALING FACTOR FOR RADII =',F6.2,/)
      END IF
C
      DO 3 I = 1, N
         READ (5,*) BODY(I), (X(K,I), K = 1,3), (XDOT(K,I), K = 1, 3)
 3    CONTINUE 
C
      ALPHA = TWOPI/FLOAT (NBTOT)
      ALPHA1 = 1.0/ALPHA
      NBMAX = 2*NBPERT + 1
      RIN2 = RIN**2
      ROUT2 = ROUT**2
      DFMIN2 = DFMIN**2
C          DFMIN IS A DIMENSIONLESS PARAMETER WHICH SHOULD NOT BE SCALED.
      ZMOON = SCALE/(81.0*330000.0)
      ZMOON = 161.0*ZMOON/DFLOAT (N)
C          MASS OF THE MOON IN UNITS OF THE SOLAR MASS.
      EMBRYO = 4.0*ZMOON
      ROCHE = (ZMOON/3.0)**0.33333
      ZMAX = ZMAX*ROCHE
C          MAXIMUM HALF-THICKNESS IN AU (SCALING BY ROCHE RADIUS).
      RM = 1.15E-05*(161.0*SCALE/DFLOAT (N))**0.3333
C          RADIUS OF THE MOON IN AU.
      DFMAX = ZMOON/(9.0*RM**2)
      DFMAX = DFMAX*ZMOON
C          THE COLLISION FORCE CRITERION IS SCALED BY THE INITIAL MASS.
      DFMAX2 = DFMAX**2
      IF (KZ(4).EQ.0)  DELTA = (ROUT2 - RIN2)/FLOAT (N)
      IF (KZ(4).EQ.1)  DELTA = (DSQRT(ROUT) - DSQRT(RIN))/FLOAT (N)
      IF (KZ(4).EQ.2)  DELTA = (ROUT**1.5 - RIN**1.5)/FLOAT (N)
      IF (KZ(4).EQ.3)  DELTA = DLOG (ROUT/RIN)/FLOAT (N)
C
C          INITIALIZE RANDOM NUMBER GENERATOR AND SKIP FIRST NRAND MEMBERS.
      KKK = -1
      RAN0 = RAN2(KKK)
      DO 15 I = 1,NRAND
      RAN0 = RAN2(KKK)
   15 CONTINUE
      R1 = RIN
C
C          GENERATE INITIAL CONDITIONS.
      DO 20 I = 1,N
      ISTAB(I) = KZ(1)
C          INITIAL STABILIZATION MODE.
      R(I) = RM
C          RADIUS OF THE MOON IN UNITS OF AU.
      R(I) = R(I)*SCALE**0.3333
      IF (KZ(15).NE.0)  R(I) = RSCALE*R(I)
C          SCALING OF INITIAL PARTICLE RADIUS BY INPUT PARAMETER.
      SPIN(I) = 0.0
      RADIUS = R1
      IF (KZ(4).EQ.0)  R1 = DSQRT (R1**2 + DELTA)
      IF (KZ(4).EQ.1)  R1 = (DSQRT (R1) + DELTA)**2
      IF (KZ(4).EQ.2)  R1 = (R1**1.5 + DELTA)**0.66667
      IF (KZ(4).EQ.3)  R1 = RIN*EXP (FLOAT (I - 1)*DELTA)
 16   IF (KZ(14).NE.0)  THEN
          KZ14 = KZ(14)
          XDOT(1,I) = XDOT(1,I)*(1.0 + 0.01*FLOAT (KZ14))
          XDOT(2,I) = XDOT(2,I)*(1.0 + 0.01*FLOAT (KZ14))
      END IF
C          CHECK FOR INITIALLY BOUND PAIRS (LAST THREE BODIES ONLY).
      I1 = I - 3
      IF (I1.LT.1)  I1 = 1
      DO 18 J = I1,I
      IF (J.EQ.I)  GO TO 18
      RIJ2 = (X(1,I)-X(1,J))**2 + (X(2,I)-X(2,J))**2 +(X(3,I)-X(3,J))**2
      IF (RIJ2.GT.0.0025)  GO TO 18
C          QUICK TEST ON SQUARE DISTANCE TO SAVE TIME.
      A1 = (XDOT(1,I) - XDOT(1,J))**2 + (XDOT(2,I) - XDOT(2,J))**2
     &                                + (XDOT(3,I) - XDOT(3,J))**2
      EREL = 0.5D0*A1 - (BODY(I) + BODY(J))/DSQRT (RIJ2)
      IF (EREL.GT.0.0)  GO TO 18
      WRITE (6,17)  I,J,RIJ2,EREL
   17 FORMAT (/,10X,'INITIAL BOUND PAIR REJECTED',2I6,1P2E10.1)
      GO TO 16
   18 CONTINUE
C
C          SET INITIAL PHASE ANGLE BIN & INCLUDE IT IN THE PERTURBER LIST.
      IBIN = 1 + THETA/ALPHA
      ILIST(I) = IBIN
      NNB = LIST(1,IBIN) + 1
      LIST(NNB+1,IBIN) = I
      LIST(1,IBIN) = NNB
      NAME(I) = I
      IBINR = ROUT2/RADIUS**2
      LISTR(I) = IBINR
      MPERT(I) = 1 + SQRT (BODY(I)/EMBRYO)*NBPERT
   20 CONTINUE
      LASTNAME = N
C
      MARS = 0
      IF (KZ(13).EQ.0)  GO TO 28
C
C          GENERATE INITIAL MASS FUNCTION N(M) = M**(-3/2) IN ZM2 TO ZM1.
      READ (5,*)  ZM1,ZM2
      ZMASS = 0.0
      DO 22 I = 1,N
   21 ZM = ZM1/RAN2(KKK)**2
      IF (ZM.GT.ZM2)  GO TO 21
      BODY(I) = ZMOON*ZM
      ZMASS = ZMASS + BODY(I)
   22 CONTINUE
C
      ZM = 1.0
      ZH = 0.0
C          RESCALE TO TOTAL MASS = N*ZMOON AND SET CORRECT SIZE.
      DO 24 I = 1,N
      BODY(I) = BODY(I)*FLOAT (N)*ZMOON/ZMASS
      R(I) = R(I)*(BODY(I)/ZMOON)**0.3333
      IF (BODY(I).GT.ZH)  ZH = BODY(I)
      IF (BODY(I).GT.ZM)  GO TO 24
      ZM = BODY(I)
      IM = I
      MPERT(I) = 1 + SQRT (BODY(I)/EMBRYO)*NBPERT
   24 CONTINUE
C
      DFMAX = ZM**2/(9.0*R(IM)**2)
      DFMAX2 = DFMAX**2
      WRITE (6,25)  ZM/ZMOON,ZH/ZMOON,DFMAX
   25 FORMAT (//,10X,'IMF   MIN =',F7.3,'  MAX =',F7.3,
     &                                          '  DFMAX =',1PE10.1)
C 
      JUPITER = 0
   28 IF (KZ(3).NE.1)  GO TO 30
C          INCLUDE ONE INITIAL BODY WITH SPECIFIED MASS & ORBITAL PARAMETERS.
      N = N + 1
      LASTNAME = LASTNAME + 1
      JUPITER = N
      READ  (5,*)  BODY(N),SEMI(N),ECC(N)
      WRITE (6,29)  BODY(N)/ZMOON,SEMI(N),ECC(N)
   29 FORMAT (//,5X,'PERTURBING PLANET AT PERICENTRE',F12.4,2F10.4)

      R(N) = R(1)*(BODY(N)/ZMOON)**0.3333
      X(1,N) = SEMI(N)*(1.0 - ECC(N))
C          PERTURBING PLANET AT PERICENTRE.
      X(2,N) = 0.0
      X(3,N) = 0.0
      XDOT(1,N) = 0.0
      SUNPL = 1.0 + BODY(N)
      XDOT(2,N) = DSQRT (SUNPL*(1.0D0 + ECC(N)) / 
     &                   (SEMI(N)*(1.0D0 - ECC(N))))
      XDOT(3,N) = 0.0
      ISTAB(N) = KZ(1)
      IBIN = 1
      ILIST(N) = IBIN
      NNB = LIST(1,IBIN) + 1
      LIST(NNB+1,IBIN) = N
      LIST(1,IBIN) = NNB
      NAME(N) = LASTNAME
      MARS = N
      IBINR = ROUT2/(X(1,N)**2 + X(2,N)**2)
      LISTR(N) = IBINR
      SPIN(N) = 0.0
      MPERT(N) = 1 + SQRT (BODY(N)/EMBRYO)*NBPERT
C
C          SET OUTPUT INTERVAL AND TERMINATION TIME IN SCALED TIME UNITS.
   30 DELTAT = DELTAT*TWOPI
      TCRIT = TCRIT*TWOPI
C
C          SET NEW FORCE DIFFERENCES AND TIME-STEPS BY THE EXPLICIT METHOD.
      DO 40 I = 1,N
      CALL UPDATE_ORBIT(I)
      CALL FPOLY(I)
   40 CONTINUE
C
      I = 1 + 0.25*FLOAT (N)
      DTLIST = STEP(I)
      IF (N.EQ.1)  DTLIST = TCRIT
   45 TLIST = DTLIST
      NNB = 0
C
C          INCLUDE ALL BODIES IN TIME-STEP LIST WITH STEP(J) < DTLIST.
      DO 50 J = 1,N
      IF (T0(J) + STEP(J).GE.TLIST)  GO TO 50
      NNB = NNB + 1
      NLIST(NNB+1) = J
   50 CONTINUE
C
      IF (NNB.EQ.0)  DTLIST = 1.2*DTLIST
      IF (NNB.EQ.0)  GO TO 45
      NLIST(1) = NNB
      IF (NNB.LT.100)  GO TO 60
      DTLIST = 0.9*DTLIST
      GO TO 45
C
   60 RETURN
C
      END
