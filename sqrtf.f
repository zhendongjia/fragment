      SUBROUTINE SQRTF
C
C          LOOK-UP TABLE FOR 1/R**3.
C          -------------------------
C
      REAL*8  R2
      COMMON/TABLE/  RSMIN2,RSMAX2,R3INV(0:10000)
C
C          SET FIRST ELEMENT TO ZERO.
      R3INV(0) = 0.0
C          INITIALIZE THE TABLE ON AN EXPANDING SCALE.
      DO 1 I = 1,10000
      R2 = RSMAX2/DFLOAT (I)
      R3INV(I) = 1.0/(R2*DSQRT (R2))
    1 CONTINUE
C
      RETURN
      END
