C
C                         F R A G M E N T
C                         ***************
C
C          PROGRAM FOR COLLIDING & FRAGMENTING PLANETESIMALS.
C          --------------------------------------------------
C
C          DEVELOPED BY SVERRE AARSETH & CHRIS BEAUGE, CAMBRIDGE.
C          ......................................................
C
C
      INCLUDE 'COMMONP.FOR'
C
C
C          INITIALIZE THE TIMER.
      CALL CPUTIM(CPU0)
C
      READ  (5,*)  KSTART,TCOMP
      IF (KSTART.EQ.1)  GO TO 2
C
C          READ PREVIOUSLY SAVED COMMON VARIABLES.
      CALL COMMON(0,1)
C          CHECK SAFETY INDICATOR FOR TROUBLESOME BINARY.
      IF (BINARY.GT.0.0)  STOP
C
C          INITIALIZE THE RANDOM NUMBER GENERATOR AFTER RESTART.
      J = -1
      DUM1 = RAN2(J)
C
C          SET TIME LIMIT AND CONTINUE (UNLESS NEW RESTART PARAMETERS).
      CPU = TCOMP
      CPU0 = 0.0
      IF (KSTART.EQ.2)  GO TO 10
C
C          READ MODIFIED PARAMETERS (KSTART = 3 OR 4).
      IF (KSTART.GE.3)  THEN
          READ (5,*)  TSCALE,TCRIT,J,K
          IF (TSCALE.GT.0.0)  DELTAT = TSCALE*DELTAT
          IF (TNEXT.GT.TIME + DELTAT)  TNEXT = TIME
          TCRIT = TWOPI*TCRIT
          IF (J.GT.0)  KZ(J) = K
      END IF
      IF (KSTART.EQ.4)  THEN
          NBP0 = NBPERT
          READ (5,*)  ETA,NBPERT,NRPERT
          NBMAX = 2*NBPERT + 1
C          SEE WHETHER TO UPDATE INDIVIDUAL MASS PERTURBER BINS.
          IF (NBPERT.NE.NBP0)  THEN
              DO 9 I = 1,N
              MPERT(I) = 1 + SQRT (BODY(I)/EMBRYO)*NBPERT
    9         CONTINUE
          END IF
      END IF
C
      WRITE (6,1)  NBPERT,NRPERT,ETA,DELTAT,TCRIT,J,K
    1 FORMAT (//,5X,'NEW PARAMETERS ',2I5,F12.4,2F10.1,2I5)
      GO TO 10
C
    2 CPU = TCOMP
      CALL INPUT
C
    5 CALL OUTPUT
      IF (TIME.GT.0.0)  GO TO 10
C
   10 CALL FLUSH(6)
      CALL INTGRT
C
      GO TO 5
C
      END
