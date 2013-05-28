      SUBROUTINE MERGE(ICOMP,JCOMP)
C
C
C          MERGING OF COLLIDING BODIES.
C          ----------------------------
C
      INCLUDE 'COMMONP.FOR'
C
      WRITE (6,1) NAME(ICOMP), NAME(JCOMP), NAME(ICOMP)
 1    FORMAT (5X,'MERGE:',I5, I5, " -> ", I5)
C
C          FIRST EVALUATE THE OLD SOLAR INTERACTION TERMS.
      RI = DSQRT (X(1,ICOMP)**2 + X(2,ICOMP)**2 + X(3,ICOMP)**2)
      RJ = DSQRT (X(1,JCOMP)**2 + X(2,JCOMP)**2 + X(3,JCOMP)**2)
      POTS = BODY(ICOMP)/RI + BODY(JCOMP)/RJ
C
C          SET NEW VARIABLES IN ICOMP FOR THE MERGED BODY.
      BCM = BODY(ICOMP) + BODY(JCOMP)
      SPIN(ICOMP) = SPIN(ICOMP) + SPIN(JCOMP) + ((X(1,ICOMP)-X(1,JCOMP))
     &     *(XDOT(2,ICOMP) - XDOT(2,JCOMP)) - (X(2,ICOMP) - X(2,JCOMP))*
     &      (XDOT(1,ICOMP) - XDOT(1,JCOMP)))*BODY(ICOMP)*BODY(JCOMP)/BCM
C
      DO 5 K = 1,3
      X(K,ICOMP) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/BCM
      XDOT(K,ICOMP) = (BODY(ICOMP)*XDOT(K,ICOMP) +
     &                                    BODY(JCOMP)*XDOT(K,JCOMP))/BCM
    5 CONTINUE
C
      R(ICOMP) = R(ICOMP)*(BCM/BODY(ICOMP))**0.3333
      BODY(ICOMP) = BCM
      MPERT(ICOMP) = 1 + SQRT (BODY(ICOMP)/EMBRYO)*NBPERT
C
C          REMEMBER NAME OF HEAVY PERTURBER IF IT COLLIDES.
      IF (KZ(3).GT.0.AND.NAME(JCOMP).EQ.MARS)  MARS = NAME(ICOMP)
      IF (KZ(3).GT.0.AND.JCOMP.EQ.JUPITER) JUPITER = ICOMP
C
C          OBTAIN THE C.M. SOLAR INTERACTION TERM.
      POTCM = BODY(ICOMP)/DSQRT (X(1,ICOMP)**2 + X(2,ICOMP)**2 +
     &                                           X(3,ICOMP)**2)
C          CORRECT THE TOTAL ENERGY FOR CHANGE IN INTERACTION TERMS.
      ETOT = ETOT + POTS - POTCM
C
C          REMOVE JCOMP FROM THE PERTURBER LIST OF JCOMP OR ICOMP.
      JBIN = ILIST(JCOMP)
      IBIN = ILIST(ICOMP)
   10 NNB = LIST(1,JBIN)
      L = 1
   20 IF (L.GT.NNB)  GO TO 24
      L = L + 1
      IF (LIST(L,JBIN).NE.JCOMP)  GO TO 20
      IF (L.GT.NNB) GO TO 23
      DO 22 K = L,NNB
   22 LIST(K,JBIN) = LIST(K+1,JBIN)
   23 LIST(1,JBIN) = NNB - 1
      JBIN = -1
   24 IF (JBIN.EQ.-1.OR.JBIN.EQ.IBIN)  GO TO 25
      JBIN = IBIN
      GO TO 10
C
C          FIRST SEE WHETHER ANY EMBRYO TABLES SHOULD BE REMOVED.
   25 KEMB = 1
   30 IF (NEMB.EQ.0)  GO TO 60
      DO 35 I = 1,N
      IF (NAME(I).EQ.LISTE(KEMB))  GO TO 50
   35 CONTINUE
C
C          REDUCE MEMBERSHIP AND MOVE UP EMBRYO TABLES (UNLESS KEMB IS LAST).
      NEMB = NEMB - 1
      IF (KEMB.GT.NEMB)  GO TO 50
C
      DO 46 L = KEMB,NEMB
      LISTE(L) = LISTE(L+1)
      LISTC(L) = LISTC(L+1)
      NNB = LISTC(L)
      DO 44 M = 1,NNB
      DO 42 K = 1,6
   42 IEMB(K,M,L) = IEMB(K,M,L+1)
   44 CONTINUE
   46 CONTINUE
C
C          CHECK SAME LOCATION AGAIN BECAUSE OF UPDATING.
      GO TO 30
C
C          CONSIDER THE NEXT EMBRYO (IF ANY).
   50 KEMB = KEMB + 1
      IF (KEMB.LE.NEMB)  GO TO 30
C
C          SET CURRENT EMBRYO INDEX (= 0 IF NEW EMBRYO OR TOO SMALL MASS).
   60 KEMB = 0
      NAMEJ = NAME(JCOMP)
      IF (BODY(ICOMP).LT.EMBRYO)  GO TO 70
C
      DO 65 K = 1,NEMB
      IF (LISTC(K).EQ.NAME(ICOMP))  KEMB = K
   65 CONTINUE
C
C          INCREASE EMBRYO INDEX UNLESS IDENTIFIED ABOVE (LIMIT IS MBMAX).
      IF (KEMB.EQ.0.AND.NEMB.LT.MBMAX)  THEN
          NEMB = NEMB + 1
          KEMB = NEMB
          LISTE(KEMB) = NAME(ICOMP)
          LISTC(KEMB) = 0
      END IF
C
C          UPDATE ALL COMMON VARIABLES.
   70 N = N - 1
      IF (JUPITER.GT.JCOMP) JUPITER = JUPITER - 1
      IF (JCOMP.GT.N)  GO TO 90
C
      DO 80 J = JCOMP,N
      J1 = J + 1
      T0(J) = T0(J1)
      T1(J) = T1(J1)
      T2(J) = T2(J1)
      T3(J) = T3(J1)
      BODY(J) = BODY(J1)
      R(J) = R(J1)
      SPIN(J) = SPIN(J1)
      A0(J) = A0(J1)
      STEP(J) = STEP(J1)
      H(J) = H(J1)
      HDOT(J) = HDOT(J1)
      D1HDOT(J) = D1HDOT(J1)
      D2HDOT(J) = D2HDOT(J1)
      D3HDOT(J) = D3HDOT(J1)
      ILIST(J) = ILIST(J1)
      ISTAB(J) = ISTAB(J1)
      LISTR(J) = LISTR(J1)
      NAME(J) = NAME(J1)
      MPERT(J) = MPERT(J1)
C
      DO 75 K = 1,3
      X(K,J) = X(K,J1)
      X0(K,J) = X0(K,J1)
      XDOT(K,J) = XDOT(K,J1)
      X0DOT(K,J) = X0DOT(K,J1)
      F(K,J) = F(K,J1)
      FDOT(K,J) = FDOT(K,J1)
      FI(K,J) = FI(K,J1)
      D1(K,J) = D1(K,J1)
      D2(K,J) = D2(K,J1)
      D3(K,J) = D3(K,J1)
   75 CONTINUE
   80 CONTINUE
C
C          REDUCE HIGHER LOCATIONS IN THE PERTURBER LIST BY ONE.
   90 DO 100 J = 1,NBTOT
      NNB = LIST(1,J) + 1
      IF (NNB.EQ.1)  GO TO 100
      DO 95 L = 2,NNB
      IF (LIST(L,J).GT.JCOMP)  LIST(L,J) = LIST(L,J) - 1
   95 CONTINUE
  100 CONTINUE
C
C           MODIFY THE TIME-STEP LIST DUE TO REMOVAL OF JCOMP.
      NNB = NLIST(1)
      L = 2
  170 IF (NLIST(L).NE.JCOMP)  GO TO 175
      IF (L.GT.NNB)  GO TO 173
C          MOVE UP THE SUBSEQUENT LIST MEMBERS AND REDUCE MEMBERSHIP BY ONE.
      DO 172 K = L,NNB
  172 NLIST(K) = NLIST(K+1)
  173 NLIST(1) = NLIST(1) - 1
C          REDUCE HIGHER PARTICLE LOCATIONS BY ONE.
  175 IF (NLIST(L).GT.JCOMP)  NLIST(L) = NLIST(L) - 1
      L = L + 1
      IF (L.LE.NLIST(1) + 1)  GO TO 170
C
C          SET AN ARBITRARY PARTICLE IN FIRST LOCATION IF NLIST IS EMPTY.
      IF (NLIST(1).EQ.0)  THEN
          NLIST(2) = 1
          NLIST(1) = 1
      END IF
C
C          OBTAIN NEW FORCE DIFFERENCES AND TIME-STEP.
      CALL FPOLY(ICOMP)
C
C          FORM NEW ELEMENTS WITH RESPECT TO THE SUN.
      SEMI = -0.5*(1.0 + BCM)/H(ICOMP)
      RDOT = X(1,ICOMP)*XDOT(1,ICOMP) + X(2,ICOMP)*XDOT(2,ICOMP) +
     &                                  X(3,ICOMP)*XDOT(3,ICOMP)
      ECC = DSQRT ((1.0 - RI/SEMI)**2 + RDOT**2/((1.0 + BCM)*SEMI))
C
      IF (KEMB.EQ.0)  GO TO 176
      IF (LISTC(KEMB).GE.NMAX)  GO TO 176
C
C          UPDATE EMBRYO VARIABLES (EVENT #, TIME, MASS, SEMI, E, EZ & NAME).
      IZ = (X(3,ICOMP)**2 + XDOT(3,ICOMP)**2)/ROCHE**2
C          INCREASE COLLISION COUNTER AND SAVE INTEGER SCALED VARIABLES.
      LISTC(KEMB) = LISTC(KEMB) + 1
      L = LISTC(KEMB)
      IEMB(1,L,KEMB) = TIME/(1000.0*TWOPI)
      IEMB(2,L,KEMB) = BCM/ZMOON
      IEMB(3,L,KEMB) = 10000.0*SEMI
      IEMB(4,L,KEMB) = 10000.0*ECC
      IEMB(5,L,KEMB) = 1000.0*MIN (IZ,30)
      IEMB(6,L,KEMB) = NAMEJ
C
C          CHECK OPTION FOR DIAGNOSTIC OUTPUT.
  176 IF (KZ(5).EQ.0)  GO TO 180
      WRITE (6,178)  ICOMP,BODY(ICOMP)/ZMOON,SEMI,ECC,RI,X(3,ICOMP)
  178 FORMAT (5X,'MERGER :',I5,F9.4,3F8.3,F9.4)
C
  180 NSTEPN(4) = NSTEPN(4) + 1
      RETURN
C
      END
