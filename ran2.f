      FUNCTION RAN2(IDUM)
C
C
C          RANDOM NUMBER GENERATOR (PRESS P. 195).
C          ---------------------------------------
C
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1./M)
      COMMON/RAND/  IY,IFF,IR(97) 
CC      DATA IFF /0/
C
      IF (IDUM.LT.0.OR.IFF.EQ.0)  THEN
         IFF = 1
         IDUM = MOD(IC-IDUM,M)
         DO 11 J = 1,97
            IDUM = MOD(IA*IDUM+IC,M)
            IR(J) = IDUM
   11       CONTINUE
         IDUM = MOD(IA*IDUM+IC,M)
         IY = IDUM
      END IF
      J = 1 + (97*IY)/M
      IF (J.GT.97.OR.J.LT.1)  write (6,12) J,IDUM
   12 FORMAT (/,'  TROUBLES IN RAN2   J IDUM ',2I12)
      IY = IR(J)
      RAN2 = IY*RM
      IDUM = MOD(IA*IDUM+IC,M)
      IR(J) = IDUM
      RETURN
      END 
