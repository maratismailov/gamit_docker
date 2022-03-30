      SUBROUTINE TF2REF( X,XP,SCALE,T,R )
C
C     Perform a 7-parameter transformation on a set of Cartesian
C     coordinates using the linearized model:
C
C             X' = (1+SCALE)*R*X  + T
C
C     where T is a translation vector, SCALE is a scalar, and R
C     is the (linearized) rotatation matrix about axes 1, 2, and 3.
C     The units of X and T are meters; of the angles, arcseconds
C
C     R.King   3 February 1987
C
      implicit none
      integer*4 i,j
      real*8 x,xp,t,r,rot,convds,scale
      DIMENSION T(3),X(3),XP(3),R(3),ROT(3,3)
C
C       Arcseconds to radians
      CONVDS= DATAN(1.D0)*4.D0/180.D0/3600.D0
C
C       Calculate the rotation matrix
C
      ROT(1,1)= 1.D0
      ROT(1,2)=  R(3) * CONVDS
      ROT(1,3)= -R(2) * CONVDS
      ROT(2,1)= -R(3) * CONVDS
      ROT(2,2)=  1.D0
      ROT(2,3)=  R(1) * CONVDS
      ROT(3,1)=  R(2) * CONVDS
      ROT(3,2)= -R(1) * CONVDS
      ROT(3,3)=  1.D0
C
      DO 10 I=1,3
      XP(I)= 0.D0
      DO 10 J=1,3
   10 XP(I)= XP(I) + (1.D0 + SCALE)*ROT(I,J)*X(J)
      DO 20 I=1,3
   20 XP(I)= XP(I) + T(I)
C
      RETURN
      END
