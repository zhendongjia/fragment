      SUBROUTINE COMMON(I,J)
C
C
C          COMMON SAVE OR READ WITH TAPE/DISC.
C          -----------------------------------
C
      PARAMETER  (NMAX=250,LMAX=20,MBMAX=20)
      PARAMETER  (NA=375*NMAX,NB=10*NMAX)
      PARAMETER  (NC=(LMAX*NMAX+7*NMAX+MBMAX*(5*NMAX+LMAX+1))/2+52)
      COMMON/NBODY/  A(NA)
      COMMON/ORBIT/  B(NB)
      COMMON/NAMES/  C(NC)
      COMMON/PARAMS/ D(90)
C
C
      REWIND  J
      IF (I.NE.0)  GO TO 1
      READ (J)  A,B,C,D
      RETURN
C
    1 WRITE (J)  A,B,C,D
      END FILE J
C
      RETURN
C
      END
