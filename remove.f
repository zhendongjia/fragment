      SUBROUTINE REMOVE(ID)

      INCLUDE 'COMMONP.FOR'

      WRITE (6, '(A,I10)') 'REMOVE:',  NAME(ID)

C     UPDATE ALL COMMON VARIABLES.
      N = N - 1
      IF (IS_JUPITER(ID).EQ.1) THEN
         LOCATION = 0
         DO K=1,NJUPITER
            IF (JUPITER(K).EQ.ID) THEN
               LOCATION = K
            END IF
         END DO
         NJUPITER = NJUPITER - 1
         DO K = LOCATION, NJUPITER
            JUPITER(K) = JUPITER(K+1)
         END DO
      END IF
      DO K=1, NJUPITER
         IF (JUPITER(K).GT.ID) JUPITER(K)=JUPITER(K)-1
      END DO

      DO J = ID,N
         J1 = J + 1
         T0(J) = T0(J1)
         DT01(J) = DT01(J1)
         DT02(J) = DT02(J1)
         DT03(J) = DT03(J1)
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
         SEMI(J) = SEMI(J1)
         ECC(J) = ECC(J1)     
         IS_JUPITER(J) = IS_JUPITER(J1)
         DO K = 1,3
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
         END DO
      END DO

C     REDUCE HIGHER LOCATIONS IN THE PERTURBER LIST BY ONE.
      DO J = 1,NBTOT
         NNB = LIST(1,J)
         LOCATION = 0
         DO L = 2,NNB+1
            IF (LIST(L,J).EQ.ID) LOCATION = L
            IF (LIST(L,J).GT.ID)  LIST(L,J) = LIST(L,J) - 1
         END DO
         IF (LOCATION.NE.0) THEN
            DO L = LOCATION, NNB
               LIST(L,J) = LIST(L+1,J)
            END DO
            LIST(1,J) = NNB - 1
         END IF
      END DO

C     MODIFY THE TIME-STEP LIST DUE TO REMOVAL OF JCOMP.
      NNB = NLIST(1)
      LOCATION  = 0
      DO L = 2, NNB+1
         IF (NLIST(L).EQ.ID) LOCATION = L
         IF (NLIST(L).GT.ID) NLIST(L) = NLIST(L) - 1
      END DO
      IF (LOCATION.NE.0) THEN
         DO L = LOCATION, NNB
            NLIST(L) = NLIST(L+1)
         END DO
         NLIST(1) = NNB - 1
      END IF

C     SET AN ARBITRARY PARTICLE IN FIRST LOCATION IF NLIST IS EMPTY.
      IF (NLIST(1).EQ.0)  THEN
         NLIST(2) = 1
         NLIST(1) = 1
      END IF

      END SUBROUTINE REMOVE
