



      SUBROUTINE WRITE_OBJECT(J, FOUT, PREFIX)
      INCLUDE 'COMMONP.FOR'
      INTEGER J, FOUT
      CHARACTER(LEN=*) PREFIX
      CALL UPDATE_ORBIT(J)
      CALL GET_PRECESSION(J, W(J))
      WRITE (FOUT, '(A,I5,12E12.4)') PREFIX, NAME(J), BODY(J)/M_EARTH,
     &     SEMI(J), ECC(J), (X(K,J),K=1,3), (XDOT(K,J),K=1,3),
     &     T_TIDAL1(J), T_TIDAL2(J), W(J)
      END SUBROUTINE
