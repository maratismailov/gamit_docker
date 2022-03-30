Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE DISTNC(COORDS,IS1,IS2,DIST)
C
C     DETERMINE DISTANCE FOR A PARTICULAR BASELINE

      implicit none
 
      integer*4 is1,is2,indx
      real*8 dist,x1,y1,z1,x2,y2,z2,dx,dy,dz

      REAL*8 COORDS(*)
C
      INDX=3*(IS1-1)
      X1=COORDS(INDX+1)
      Y1=COORDS(INDX+2)
      Z1=COORDS(INDX+3)
      INDX=3*(IS2-1)
      X2=COORDS(INDX+1)
      Y2=COORDS(INDX+2)
      Z2=COORDS(INDX+3)
C
      DX=X2-X1
      DY=Y2-Y1
      DZ=Z2-Z1
C
      DIST=DSQRT(DX*DX+DY*DY+DZ*DZ)
C
      RETURN
      END
