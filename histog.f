      SUBROUTINE HISTOG(AZ)
C
C
C          DATA ANALYSIS AND OUTPUT.
C          -------------------------
C
      INCLUDE 'COMMONP.FOR'
      REAL *4  WK,DISP,FM,ZN,ZMASS,EBAR
      DIMENSION  A(20),WK(10),DISP(20),FM(NMAX),EBAR(10)
      INTEGER*2  MASS(50),IHIST(NMAX),JHIST(25),NEBAR(10),IAHIST(35),
     &		 ISHIST(10)
C
C
      DO 5 J = 1,20
    5 DISP(J) = 0.0
      KMAX = 0
      KSEMIMAX = 0
      KKMAX = 0
      DO 6 J = 1,25
      JHIST(J) = 0
      IAHIST(J) = 0
      IHIST(J) = 0
    6 CONTINUE
      DO 8 J = 1,10
      EBAR(J) = 0.0
      ISHIST(J) = 0
      NEBAR(J) = 0
    8 CONTINUE
      ZMAX1 = 0.0
      ZMIN = 1.0
      JUP = 0
      IF (KZ(3).GT.0)  THEN
          DO 10 I = 1,N
          IF (BODY(I).GT.1000.0*ZMOON)  JUP = I
   10     CONTINUE
      END IF
      DO 12 I = 1,N
      IF (BODY(I).LT.ZMIN)  ZMIN = BODY(I)
      IF (BODY(I).LE.ZMAX1)  GO TO 12
      IF (I.EQ.JUP)  GO TO 12
      ZMAX1 = BODY(I)
      IMAX = I
   12 CONTINUE
      M1 = ZMAX1/ZMOON
      M0 = ZMIN/ZMOON
      IF (M0.LT.1)  THEN
          M0 = 1
          ZMIN = ZMOON
      END IF
      NB1 = 0
      NB2 = 0
      NDIR = 0
      NRETRO = 0
      ZMASS = 0.0
      EMAX = 0.0
      IGONE = 0
      DO 20 I = 1,N
      IF (I.EQ.JUP)  GO TO 20
      RI = DSQRT (X0(1,I)**2 + X0(2,I)**2 + X0(3,I)**2)
      SUNPL = 1.0 + BODY(I)
      ENERGY = 0.5D0*(X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2) -
     &                                                          SUNPL/RI
C          BINDING ENERGY PER UNIT MASS WITH RESPECT TO THE SUN.
      SEMI = -0.5D0*SUNPL/ENERGY
      RDOT = X0(1,I)*X0DOT(1,I) + X0(2,I)*X0DOT(2,I) +X0(3,I)*X0DOT(3,I)
      ECC = (1.0 - RI/SEMI)**2 + RDOT**2/(SEMI*SUNPL)
      ECC = DSQRT (ECC)
      IF (ECC.GE.1.0.OR.SEMI.GT.4.*ROUT/3.)  THEN
          IGONE = I
          AGONE = SEMI
          EGONE = ECC
      END IF
      IF (ECC.GT.EMAX)  EMAX = ECC
      IF (I.EQ.IMAX)  THEN
           EIMAX = ECC
           AIMAX = SEMI
      END IF
      K = 1 + ECC/0.01
C          INTERVAL CHANGED FROM 0.005 TO 0.01 ON 26 OCTOBER 1982.
      IF (K.GT.25)  K = 25
      IF (K.GT.KMAX)  KMAX = K
      JHIST(K) = JHIST(K) + 1
      IF (SEMI.LE.0.0)  GO TO 13
      KSEMI = INT (1 + SEMI/0.1)
      IF (KSEMI.GT.KSEMIMAX.AND.KSEMI.LE.20)  KSEMIMAX = KSEMI
      IAHIST(KSEMI) = IAHIST(KSEMI) + INT (BODY(I)/ZMOON)
   13 WK(1) = BODY(I)/ZMOON
      IF (WK(1).LE.1.0)  WK(1) = 1.0
      KK = 1 + ALOG10 (WK(1))/ALOG10 (1.999)
      IF (KK.GT.KKMAX)  KKMAX = KK
      IF (KKMAX.GT.7)  KKMAX = 7
      EBAR(KK) = EBAR(KK) + ECC
      NEBAR(KK) = NEBAR(KK) + 1
      IF (KZ(10).EQ.0)  GO TO 14
      ZMASS = ZMASS + BODY(I)
      DISP(1) = DISP(1) + ECC
      DISP(2) = DISP(2) + ECC**2
      DISP(3) = DISP(3) + BODY(I)*ECC
      DISP(4) = DISP(4) + BODY(I)*ECC**2
      DISP(5) = DISP(5) + SEMI
      DISP(6) = DISP(6) + SEMI**2
      DISP(7) = DISP(7) + BODY(I)*SEMI
      DISP(8) = DISP(8) + BODY(I)*SEMI**2
      DISP(9) = DISP(9) + ECC**3
      DISP(10) = DISP(10) + ECC**4
      DISP(14) = DISP(14) + BODY(I)**2
      DISP(18) = DISP(18) + BODY(I)**3
      DISP(19) = DISP(19) + BODY(I)**4
      IF (SEMI*(1.0 - ECC).LT.RIN)  NB1 = NB1 + 1
      IF (SEMI*(1.0 + ECC).GT.ROUT)  NB2 = NB2  +  1
   14 A(1) = BODY(I)
      IF (A(1).EQ.0.0)  A(1) = 1.0
      ROT = SPIN(I)/(0.4*A(1)*R(I)**2)
      DISP(15) = DISP(15) + DABS (ROT)
      DISP(16) = DISP(16) + ROT**2
      DISP(17) = DISP(17) + SPIN(I)
      IF (SPIN(I).GT.0.0)  NDIR = NDIR + 1
      IF (SPIN(I).LT.0.0)  NRETRO = NRETRO + 1
      IF (I.GT.KZ(9))  GO TO 20
      ERROR = (ENERGY - H(I))/H(I)
C          H HOLDS THE VALUE OF ORBITAL ENERGY AT BEGINNING OF STEP.
      THETA = DEGREE*DATAN2 (X(2,I),X(1,I))
C          THE ANGLE IS NOW IN DEGREES.
      IF (THETA.LT.0.0)  THETA = THETA + 360.0
      IMASS = BODY(I)/ZMOON
      TROT = 0.0
      IF (ROT.NE.0.0)  TROT = 365.0/ROT
      WRITE (6,16)  I,THETA,STEP(I),RI,ENERGY,SEMI,ECC,ERROR,HDOT(I),
     &                                                BODY(I),IMASS,TROT
   16 FORMAT (I9,F10.3,F10.4,F10.5,2F15.10,F12.8,1P3E10.1,0P,I6,F8.2)
   20 CONTINUE
      IF (KZ(8).EQ.0)  GO TO 40
      DO 25 J = 1,KKMAX
      IF (NEBAR(J).EQ.0)  GO TO 25
      WK(1) = NEBAR(J)
      EBAR(J) = EBAR(J)/WK(1)
   25 CONTINUE
   40 IF (KZ(10).EQ.0)  GO TO 140
      IF (DISP(16).GT.0.0)  DISP(16) = SQRT(DISP(16)/FLOAT(NDIR+NRETRO))
      IF (DISP(16).GT.0.0)  DISP(16) = 365.0/DISP(16)
      ZN = FLOAT (N)
      IF (KZ(3).GT.0)  ZN = ZN - 1.0
C          CALCULATION OF SAFRONOV'S NUMBER.
      DISP(1) = DISP(1)/ZN
      RMED = 1.15D-05*(ZMAX1/ZMOON)**.3333
      SF = ZMAX1/(RMED*DISP(1)**2)
      IF (SF.GT.99.99)  SF = 99.99
      WRITE (6,50)  NB1,NB2,NDIR,NRETRO,DISP(16),AZ,ETOT,DISP(17),SF
   50 FORMAT (/,'  IN =',I3,'  OUT =',I3,'  DIRECT =',I3,
     &              '  RETRO =',I3,'  <TROT> =',F6.2,'  AZ =',1PE12.4,
     &               '  ETOT =',E13.4,'  SPIN =',E9.1,'  SF =',0PF6.2)
      DISP(2) = DISP(2)/ZN - DISP(1)**2
      IF (DISP(2).LT.0.0)  DISP(2) = 0.0
      DISP(3) = DISP(3)/ZMASS
      DISP(4) = DISP(4)/ZMASS - DISP(3)**2
      IF (DISP(4).LT.0.0)  DISP(4) = 0.0
      DISP(4) = SQRT (DISP(4))
      DISP(5) = DISP(5)/ZN
      DISP(6) = SQRT (DISP(6)/ZN - DISP(5)**2)
      DISP(7) = DISP(7)/ZMASS
      DISP(8) = SQRT (DISP(8)/ZMASS - DISP(7)**2)
      DISP(9) = DISP(9)/ZN - DISP(1)**3 - 3.0*DISP(1)*DISP(2)
      DISP(10) = DISP(10)/ZN - DISP(1)**4 - 4.0*DISP(1)*DISP(9) -
     &                                      6.0*DISP(1)**2*DISP(2)
      IF (DISP(9).LT.0.0)  DISP(9) = 0.0
      IF (DISP(10).LT.0.0)  DISP(10) = 0.0
      DISP(9) = DISP(9)**0.3333
      DISP(10) = DISP(10)**0.25
      DISP(2) = SQRT (DISP(2))
      DISP(13) = ZMASS/ZN
      ZMBAR = DISP(13)/ZMIN
C          MEAN MASS IN UNITS OF THE MINIMUM MASS IS USED FOR Q ITERATION.
      DISP(14) = DISP(14)/ZN - DISP(13)**2
      DISP(18) = DISP(18)/ZN - DISP(13)**3 - 3.0*DISP(13)*DISP(14)
      DISP(19) = DISP(19)/ZN - DISP(13)**4 - 4.0*DISP(13)*DISP(18) -
     &                                          6.0*DISP(13)**2*DISP(14)
      IF (DISP(14).LT.0.0)  DISP(14) = 0.0
      IF (DISP(18).LT.0.0)  DISP(18) = 0.0
      IF (DISP(19).LT.0.0)  DISP(19) = 0.0
      DISP(14) = SQRT (DISP(14))/ZMOON
      DISP(13) = DISP(13)/ZMOON
      DISP(18) = DISP(18)**0.3333/ZMOON
      DISP(19) = DISP(19)**0.25/ZMOON
      IBIN = 0
      IMASS = 0
      IF (KZ(3).EQ.0)  GO TO 70
      DO 60 I = 1,N
      IF (I.EQ.JUP)  THEN
          ZMARS = BODY(I)
          GO TO 60
      END IF
      IF (BODY(I).GT.ZMAX1)  ZMAX1 = BODY(I)
   60 CONTINUE
      ZMASS = ZMASS + ZMARS
      IBINE = 0
   70 ICOUNT = 0
      IF (KZ(13).EQ.0)  GO TO 74
      Z1 = 0.1*ZMOON
      IF (Z1.GT.ZMIN)  Z1 = ZMIN
      Z2 = 2.0*Z1
      IMASS = 0
   71 IBIN = IBIN + 1
      DO 72 I = 1,N
      IF (BODY(I).GT.Z2.OR.BODY(I).LT.Z1)  GO TO 72
      ICOUNT = ICOUNT + 1
      IMASS = IMASS + 1
   72 CONTINUE
      FM(IBIN) = Z1
      IHIST(IBIN) = ICOUNT
      IF (IMASS + KZ(3).GE.N)  GO TO 78
      ICOUNT = 0
      Z1 = Z2
      Z2 = 2.0*Z2
      GO TO 71
   74 DO 75 I = 1,N
      IF (DABS (BODY(I) - ZMAX1).LT.0.0001*ZMAX1)  ICOUNT = ICOUNT + 1
   75 CONTINUE
      IBIN = IBIN + 1
      FM(IBIN) = ZMAX1
      IHIST(IBIN) = ICOUNT
      IMASS = IMASS + ICOUNT
      IF (IMASS + KZ(3).GE.N)  GO TO 78
      ZMAX2 = ZMAX1
      ZMAX1 = 0.0
      DO 77 I = 1,N
      IF (BODY(I).GT.ZMAX1.AND.BODY(I).LT.0.9999*ZMAX2)  ZMAX1 = BODY(I)
   77 CONTINUE
      IF (NSTEPN(9).EQ.0)  IBINE = IBIN
      IF (ZMAX1.LT.ZMOON.AND.IBINE.EQ.0)  IBINE = IBIN
      GO TO 70
C          HISTOGRAM OF MASS DISTRIBUTION FOR SMALL BODIES (M < ZMOON).
   78 COEF = 1.0/LOG(2.0)
      ISMAX = 0
      DO 79 I = 1,N
      IF (BODY(I).GE.ZMOON)  GO TO 79
      ISBIN = 1 - INT (COEF*LOG (BODY(I)/ZMOON))
      IF (ISBIN.GT.ISMAX) ISMAX = ISBIN
      ISHIST(ISBIN) = ISHIST(ISBIN) + 1
   79 CONTINUE
      JGROUP = 1
      KGROUP = 1
      KBIN = 0
      IMASS = 0
      IMASS1 = 0.5001*FLOAT (N + 1)
      IMASS2 = 0.5001*FLOAT (N + 2)
      NNMAX = 1
      SUM = 0.0
   80 KBIN = KBIN + 1
      K = 0
   81 SUM = SUM + FM(KBIN)
      IMASS = IMASS + 1
      IF (IMASS.EQ.IMASS1)  KBIN1 = KBIN
      IF (IMASS.EQ.IMASS2)  KBIN2 = KBIN
   82 IF (SUM.LT.0.09999*ZMASS*FLOAT (JGROUP))  GO TO 83
      MASS(JGROUP) = IMASS
      IF (IMASS.EQ.N)  GO TO 83
      JGROUP = JGROUP + 1
      GO TO 82
   83 IF (IMASS.EQ.NNMAX)  A(KGROUP) = SUM/ZMASS
      IF (IMASS + KZ(3).EQ.N)  GO TO 85
      IF (IMASS.LT.NNMAX)  GO TO 84
      NNMAX = 2*NNMAX
      IF (NNMAX + KZ(3).GT.N)  GO TO 84
      KGROUP = KGROUP + 1
   84 K = K + 1
      IF (K - IHIST(KBIN))  81,80,80
   85 DISP(15) = 0.5*(FM(KBIN1) + FM(KBIN2))/ZMOON
      DO 90 J = 1,IBINE
   90 FM(J) = FM(J)/ZMOON
      WRITE (6,999) IBINE
  999 FORMAT (5X,'IBINE =',I5)
      IF (IBINE.EQ.0)  IBINE = 1
      N0 = IHIST(IBINE)
      IF (KZ(13).NE.0)  N0 = IHIST(1)
      IF (M0.GT.1)  N0 = 0
      IF (N0.LT.100)  WRITE (6,92)  (FM(J),IHIST(J),J=1,IBINE)
   92  FORMAT ('  F(M) ',13(F6.1,I3),/,8X,10(F6.1,I3))
      IF (N0.GE.100)  WRITE (6,93)  (FM(J),IHIST(J),J=1,IBINE)
   93 FORMAT ('  F(M) ',12(F6.1,I4),/,8X,11(F6.1,I4))
      WRITE (6,94)  (ISHIST(J),J=1,ISBIN)
   94 FORMAT ('  F(SM) ',10I5)
      WRITE (6,95)  DISP(1),DISP(2),DISP(9),DISP(10),EMAX,
     &                                             (JHIST(K),K=1,KMAX)
   95 FORMAT ('  <E> =',F6.3,'  DISP(E) =',F6.3,'  S =',F6.3,
     &         '  T =',F6.3,'  EMAX =',F6.3,'  N(E) = ',14I4,/,5X,20I4)
      WRITE (6,998)  (IAHIST(KSEMI),KSEMI=4,KSEMIMAX)
  998 FORMAT ('  N(A;.3,.1) =',25I4)
      WRITE (6,96)  DISP(13),DISP(14),DISP(18),DISP(19),DISP(15),M1,
     &                                           EIMAX, AIMAX,MASS(5),N0
   96 FORMAT ('  <M> =',F6.2,'  DISP(M) =',F6.2,'  S =',F6.2,
     &    '  T =',F6.2,'  MED(M) =',F5.1,'  M(1) =',I4,'  E(1) =',F6.3,
     &                   '  A(1) =',F5.2,'  N(50%) =',I3,'  N(0) =',I3)
      IF (KZ(16).EQ.0)  GO TO 140
C          ITERATION PROCEDURE FOR Q AND L.
      IT = 0
      IT1 = 0
      Q = 1.0 + 1.0/(ZMBAR - 1.0 + 1.0E-03)
      ZLAM = 1.0/FLOAT (M1)
      DO 100 K = 1,30
  100 A(K) = 0.0
      IF (M1.EQ.1)  GO TO 120
      ZM1 = FLOAT (M1)/FLOAT (M0)
      ZLOG = DLOG (ZM1)
  110 Q1 = 1.0 - Q
      Q2 = 2.0 - Q
      A(11) = ZM1**Q1
      A(12) = ZM1**Q2
      ZM = Q1*(A(12) - 1.0)/(Q2*(A(11) - 1.0))
      DM = 1.0/Q2 - 1.0/Q1 + ZLOG*(A(11)/(A(11) - 1.0) -
     &                                             A(12)/(A(12) - 1.0))
      DQ = (ZMBAR - ZM)/(ZM*DM)
      IT = IT + 1
      Q = Q + DQ
      IF (DABS (DQ).GT.1.0D-06.AND.IT.LT.20)  GO TO 110
      Q1 = 1.0 - Q
      Q2 = 2.0 - Q
      Q3 = 3.0 - Q
      Q4 = 4.0 - Q
      Q5 = 5.0 - Q
      A(11) = Q1/(ZM1**Q1 - 1.0)
      A(7) = A(11)*(ZM1**Q2 - 1.0)/Q2
      A(8) = A(11)*(ZM1**Q3 - 1.0)/Q3
      A(9) = A(11)*(ZM1**Q4 - 1.0)/Q4
      A(10) = A(11)*(ZM1**Q5 - 1.0)/Q5
      A(1) = A(8) - A(7)**2
      A(2) = A(9) - A(7)**3 - 3.0*A(7)*A(1)
      A(3) = A(10) - A(7)**4 - 4.0*A(7)*A(2) - 6.0*A(7)**2*A(1)
      A(1) = M0*DSQRT (A(1))
      A(2) = M0*A(2)**0.3333
      A(3) = M0*A(3)**0.25
      A(18) = A(1)/A(7)
      A(19) = A(2)/A(7)
      A(20) = A(3)/A(7)
      ZMBAR = ZMBAR*M0
C          THE LAMBDA ITERATION REQUIRES ACTUAL MEAN MASS IN NATURAL UNITS.
  115 ZL = ZLAM*(M1 - M0)
      A(11) = 1.0/(1.0 - DEXP(-ZL))
      ZM = (1.0 + ZLAM*M0 - ZL*DEXP(-ZL)*A(11))/ZLAM
      DM = - (1.0 - ZL**2*DEXP(-ZL)*A(11)**2)/ZLAM**2
      DLAM = (ZMBAR - ZM)/DM
      IT1 = IT1 + 1
      ZLAM = ZLAM + DLAM
      IF (DABS (DLAM).GT.1.0D-06.AND.IT1.LT.20)  GO TO 115
      ZL = ZLAM*(M1 - M0)
      ZX0 = ZLAM*M0
      A(10) = 1.0 - DEXP(-ZL)
      A(11) = A(10) - ZL*DEXP(-ZL)
      A(12) = 2.0*A(10) - (2.0*ZL + ZL**2)*DEXP(-ZL)
      A(13) = 6.0*A(10) - (6.0*ZL + 3.0*ZL**2 + ZL**3)*DEXP(-ZL)
      A(14) = 24.0*A(10) - (24.0*ZL + 12.0*ZL**2 + 4.0*ZL**3 + ZL**4)*
     &                                                        DEXP(-ZL)
      A(7) = (A(11) + ZX0*A(10))/(ZL*A(10))
      A(8) = (A(12) + 2.0*ZX0*A(11) + ZX0**2*A(10))/(ZL**2*A(10))
      A(9) = (A(13) + 3.0*ZX0*A(12) + 3.0*ZX0**2*A(11) + ZX0**3*A(10))/
     &                                                    (ZL**3*A(10))
      A(10) = (A(14) + 4.0*ZX0*A(13) + 6.0*ZX0**2*A(12) + 
     &                   4.0*ZX0**3*A(11) + ZX0**4*A(10))/(ZL**4*A(10))
      A(4) = A(8) - A(7)**2
      A(5) = A(9) - A(7)**3 - 3.0*A(7)*A(4)
      A(6) = A(10) - A(7)**4 - 4.0*A(7)*A(5) - 6.0*A(7)**2*A(4)
      A(4) = DSQRT (A(4))
      A(5) = A(5)**0.3333
      A(6) = A(6)**0.25
      A(8) = A(5)/A(7)
      A(9) = A(6)/A(7)
      A(7) = A(4)/A(7)
  120 WRITE (6,125)  Q,(A(K),K=1,3),(A(K),K=18,20),IT
  125 FORMAT ('  Q =',F7.2,'  DISP(Q) =',F6.2,'  S =',F6.2,
     &                   '  T =',F6.2,'  Q-RATIOS =',3F6.2,'  IT =',I3)
      WRITE (6,130)  ZLAM,(A(K),K=4,9),IT1,1.0/ZLAM
  130 FORMAT ('  L =',F7.2,'  DISP(L) =',F6.2,'  S =',F6.2,
     &  '  T =',F6.2,'  L-RATIOS =',3F6.2,'  IT =',I3,'  1/L =',F7.3,/)
C
  140 IF (IGONE.LE.1)  GO TO 150
      WRITE (6,145) IGONE,BODY(IGONE)/ZMOON,AGONE,EGONE
  145 FORMAT (3X,'  GONER  :',I5,F9.4,2F8.3)
C          GET RID OF EXPELLED BODY.
      BODY(IGONE) = 0.0
      CALL MERGE(IGONE-1,IGONE)
C
  150 RETURN
C
      END
