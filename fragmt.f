      SUBROUTINE FRAGMT(I)
C
C
C          GENERATION OF FRAGMENTS.
C          ------------------------
C
      INCLUDE 'COMMONP.FOR'
      COMMON/CM/  XCM(3),VCM(3),BCM,ZMEJ,XF(NFMAX),VF(NFMAX),BF(NFMAX),
     & 		  AI,NF,IF(NFMAX),ICOMP,JCOMP
      COMMON/COLL/  VR2,VESC,VREB,EGRAV
      REAL*8  RAN2,WK(2) 
      PARAMETER(LF=3)
      REAL*8  MTOT, BETA, BMAX, BSMAX, VLR(3), A, DELT_V,
     &     MREM, S, V_ESC, PHI, NM(NFMAX), B, XR2,XDR2, ITA, 
     &     MF(LF), VMF(LF), M(LF), BODYI, XDOTI(3), RHO1, M_SUN, RC1, RI
      INTEGER LK, Q, L, NLR, NSLR, NMF(LF), TAR
C
C         DETERMINE THE TARGET
      ITA = BODY(JCOMP)/BODY(ICOMP)
      IF (ITA.GT.1)  THEN
         TAR = JCOMP
      ELSE 
         TAR = ICOMP
      END IF 
      BODYI = BODY(TAR)
      RI = R(TAR)
      DO 1 K = 1, 3 
            XDOTI(K) = XDOT(K,TAR)
 1       CONTINUE
C
      MTOT = BODY(ICOMP) + BODY(JCOMP)
      NLR = 1
      NSLR = 2
      BETA = 2.85
      BMAX = GET_LARGEST_REMNANT()
      BSMAX = MTOT*(3-BETA)*(1-NLR*BMAX/MTOT)/(NSLR*BETA)
      CONST = NSLR*BETA*((3-BETA)*(MTOT-NLR*BMAX)/(NSLR*BETA))**(BETA/3)
C
C          SPECIFY FRAGMENT MASS FUNCTION FOR BODY #I.
      NMF(1) = NLR
      NMF(2) = NSLR
      MF(1) = BMAX
      MF(2) = BSMAX
      DO 2 L = 3, LF
         MF(L) = MF(L-1)/2
         NMF(L) = (CONST/BETA)*MF(L)**(-BETA/3)
 2    CONTINUE
C     
      DO 3 L = 1, LF
         IF (L.EQ.1) THEN 
            DO 4 Q = 1, NMF(L)
               BF(Q) = MF(L)
 4          CONTINUE
         ELSE
            DO 5 Q = NMF(L-1)+1, NMF(L)
               BF(Q) = MF(L)
 5          CONTINUE
         END IF
 3    CONTINUE
C     
C     REDUCE MASS & RADIUS OF BODY #I #J AND UPDATE GRAVITATIONAL ENERGY.
      ZMEJ = MF(1)*NMF(1)
      DO 10 L =2, LF
         ZMEJ = ZMEJ + MF(L)*(NMF(L) - NMF(L-1))
 10   CONTINUE
      BODY(ICOMP) = (MTOT - ZMEJ)/2
      BODY(JCOMP) = (MTOT - ZMEJ)/2
      CALL UPDATE_ORBIT(ICOMP)
      CALL UPDATE_ORBIT(JCOMP)
      R(ICOMP) = RI*(BODY(ICOMP)/BODYI)**0.3333
      R(JCOMP) = RI*(BODY(JCOMP)/BODYI)**0.3333
      EGRAV = EGRAV + 0.6*BODY(ICOMP)**2/R(ICOMP) 
     &     + 0.6*BODY(JCOMP)**2/R(JCOMP)

C          ASSIGN NEW LOCATIONS FOR THE FRAGMENTS AND INCREASE N.
      IFIRSTFRAG = LASTNAME + 1
      DO 20 L = 1,NMF(LF)
         NF = NF + 1
         J = N + 1
         IF (L.EQ.1 .AND. (IS_JUPITER(ICOMP).EQ.1.OR.IS_JUPITER(JCOMP).EQ.1)) THEN
            IS_JUPITER(J) = 1
            NJUPITER = NJUPITER + 1
            JUPITER(NJUPITER) = J
         ELSE
            IS_JUPITER(J) = 0
         END IF
         IF(NF) = J			
         N = N + 1
         LASTNAME = LASTNAME + 1
         NAME(N) = LASTNAME
 20   CONTINUE
      WRITE(6,15) NAME(TAR), NAME(ICOMP+JCOMP-TAR), 
     &     IFIRSTFRAG, LASTNAME
 15   FORMAT(5X,"FRAGMT:", I5, I5, " -> ", I5, " ... ", I5 )
C
C          GENERATE  VELOCITIES OF THE LARGEST REMNANT
C     VCM2 = VCM(1)**2  + VCM(2)**2 + VCM(3)**2
C      VI2 = XDOT(1,ICOMP)**2 + XDOT(2,ICOMP)**2 + XDOT(3.ICOMP)**2
      B = (X(1,ICOMP) - X(1,JCOMP))*(XDOT(1,ICOMP) - XDOT(1,JCOMP))
     &     + (X(2,ICOMP) - X(2,JCOMP))*(XDOT(2,ICOMP) - XDOT(2,JCOMP))
     &     + (X(3,ICOMP) - X(3,JCOMP))*(XDOT(3,ICOMP) - XDOT(3,JCOMP))
C
      XR2 = (X(1,ICOMP) - X(1,JCOMP))**2 + (X(2,ICOMP) - X(2,JCOMP))**2
     &     + (X(3,ICOMP) - X(3,JCOMP))**2
C
      XDR2 = (XDOT(1,ICOMP) - XDOT(1,JCOMP))**2 
     &     +   (XDOT(2,ICOMP) - XDOT(2,JCOMP))**2   
     &     +   (XDOT(3,ICOMP) - XDOT(3,JCOMP))**2
C
      B = B/(SQRT(XR2)*SQRT(XDR2))
      IF(B.LE.0.7) THEN
         DO 30 K=1,3
      		 VLR(K) = VCM(K)
 30           CONTINUE
      ELSE 
         DO 35 K=1,3
      		  VLR(K) = B*(XDOTI(K)-VCM(K))/0.7 + VCM(K)
 35               CONTINUE
      END IF
C
C         GENERATE THE VELOCITY OF SMALLER REMNANT
      M_SUN = 1.989E30
      RHO1 = 1000
      RC1 = (MTOT*M_SUN/((4/3)*3.14*RHO1))**0.333/1.5E11
      A = -0.3*BMAX/MTOT + 0.3
      DELT_V = 1
      MREM = MTOT - BMAX
      S = 10**A/(LOG(10.)*DELT_V*(MREM/MTOT))
      V_ESC = SQRT(2*MTOT/RC1)
      VMF(1) = SQRT(VLR(1)**2 + VLR(2)**2 + VLR(3)**2)
      M(1) = MF(1)*NMF(1)
      DO 40 L = 2, LF
         M(L) = MF(L)*(NMF(L)-NMF(L-1)) + M(L-1)
         VMF(L) = (A-LOG10(10**A-LOG(10.)*S*DELT_V*(M(L)-M(1))/MTOT))/S
 40      CONTINUE
C     
C       GENERATE COORDINATES OF FRAGMENTS RESPECTIVE TO THE CENTER OF MASS AND SET DIRECTION FOR VELOCITY
      RESC = 10.0*RI*(BCM/BODYI)**0.3333
      KKK = 1
C
      DO 41 L = 1, LF
         IF (L.EQ.1) THEN
            DO 42 Q =1, NMF(L)
               PHI = TWOPI*(FLOAT(Q) - 0.5 + (FLOAT(NF) - 4.0)/8.0)/4.0	
               LK = 3*(Q-1)
               XF(LK+1) = RESC*DCOS (PHI)*(1.0 + 0.1*RAN2(KKK))
               XF(LK+2) = RESC*DSIN (PHI)*(1.0 + 0.1*RAN2(KKK))
               XF(LK+3) = 0
               VF(LK+1) = VLR(1)
               VF(LK+2) = VLR(2)
               VF(LK+3) = VLR(3)
 42         CONTINUE
         ELSE
            DO 43 Q = NMF(L-1)+1, NMF(L)
               PHI = TWOPI*(FLOAT(Q) - 0.5 + (FLOAT(NF) - 4.0)/8.0)/
     &                  FLOAT(NMF(LF))	
               LK = 3*(Q-1) 
               XF(LK+1) = RESC*DCOS (PHI)*(1.0 + 0.1*RAN2(KKK))
               XF(LK+2) = RESC*DSIN (PHI)*(1.0 + 0.1*RAN2(KKK))
               XF(LK+3) = 0
               VF(LK+1) = VMF(L)*V_ESC *DCOS(PHI)*(1 + 0.1*RAN2(KKK))
     &                        + VLR(1)
               VF(LK+2) = VMF(L)*V_ESC*DSIN(PHI)*(1 + 0.1*RAN2(KKK))
     &                         + VLR(2)
               VF(LK+3) = VLR(3)
 43         CONTINUE
         END IF         
 41   CONTINUE
C
C        DO 39 L =1, LF
C            WRITE (6, 38) NMF(L), MF(L)/MTOT, VMF(L)
C 38         FORMAT (3X, 'FRAGM. NMV:',I5,F9.4,F9.4)
C 39         CONTINUE
C
C          SPECIFY GLOBAL VARIABLES FOR THE FRAGMENTS (NO Z-DISPERSION).
      DO 50 Q = 1,NMF(LF)
      J = IF(Q)
C
      DO 45 K = 1,3
         LK = 3*(Q - 1) + K
         X(K,J) = XCM(K) + XF(LK)
         XDOT(K,J) = VF(LK)
   45 CONTINUE
C 
      BODY(J) = BF(Q)
      CALL UPDATE_ORBIT(J)
      ISTAB(J) = 0
C          ALLOCATE FRACTIONAL SPIN.
      SPIN(J) = SPIN(I)*BODY(J)/BODYI
      R(J) = RI*(BODY(J)/BODYI)**0.3333
      EGRAV = EGRAV + 0.6*BODY(J)**2/R(J)
      LISTR(J) = ROUT2/(X(1,J)**2 + X(2,J)**2)
      WK(1) = X(1,J)
      WK(2) = X(2,J)
      THETA = ATAN2 (WK(2),WK(1))
      IF (THETA.LT.0.0)  THETA = THETA + TWOPI
      IBIN = 1 + THETA/ALPHA
      ILIST(J) = IBIN
      NNB = LIST(1,IBIN) + 1
      LIST(NNB+1,IBIN) = J
      LIST(1,IBIN) = NNB
      MPERT(J) = 1 + SQRT (BODY(J)/EMBRYO)*NBPERT
   50 CONTINUE
C
C          RETAIN FRACTIONAL PART OF INTRINSIC REMNANT SPIN.
      SPIN(ICOMP) = SPIN(ICOMP)*BODY(ICOMP)/BODYI
      SPIN(JCOMP) = SPIN(JCOMP)*BODY(JCOMP)/BODYI
      NSTEPN(9) = NSTEPN(9) + 1
C
      RETURN
C
      END
