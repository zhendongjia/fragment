      SUBROUTINE GET_PRECESSION(I, WI)
      INCLUDE 'COMMONP.FOR'
      REAL*8 R12,HI,RDOT12,SIN_F,COS_F,SIN_W_F,COS_W_F,COS_W,SIN_W
      REAL*8 W_F,TA_F,WI
C
      R12 = SQRT(X(1,I)**2 + X(2,I)**2)
      HI = X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I)
      RDOT12 = (X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I))/R12
C
      SIN_F = SEMI(I)*(1-ECC(I)**2)*RDOT12/(HI*ECC(I))
      COS_F = SEMI(I)*(1-ECC(I)**2)/(R12*ECC(I)) - 1/ECC(I)

      SIN_W_F = X(2,I)/R12
      COS_W_F = X(1,I)/R12 
      
*     W_F = ASIN(SIN_W_F)
*     TA_F = ASIN(SIN_F)

      COS_W = COS_W_F*COS_F + SIN_W_F*SIN_F
      SIN_W = SIN_W_F*COS_F - COS_W_F*SIN_F
      
      W_F = ATAN2(SIN_W_F,COS_W_F)
      W_F = W_F*180/ACOS(-1.)
      TA_F = ATAN2(SIN_F,COS_F)
      TA_F = TA_F*180/ACOS(-1.)
      WI = ATAN2(SIN_W,COS_W)
      WI = WI*180/ACOS(-1.)
*     TA_F = TA_F*180/ACOS(-1.)
*     W_F = W_F*180/ACOS(-1.)
*     W = W_F - TA_F


      END
