Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE fillj(g,gg,sinlat,coslat,sinlon,coslon,radius,ind,
     .                 ilive,mode)
c
c     fill Jacobian matrix
c     mode = 1: d(x,y,z)/d(lat,lon,rad)
c     mode = 2: d(u,v,w)/d(x,y,z)
c
      implicit none

      real*8 G(3,6),gg(3,3)
      integer*4 ilive(3) 
      real*8 sinlat,coslat,sinlon,coslon,radius
      integer*4 ind,mode

      if (mode.eq.1) then
          IF(ILIVE(1).EQ.0) GO TO 10
c dx/dlat
          G(1,IND+1)=-RADIUS*SINLAT*COSLON
c dx/dlont
          G(1,IND+2)=-RADIUS*COSLAT*SINLON
c dx/drad
          G(1,IND+3)=COSLAT*COSLON
C
 10       IF(ILIVE(2).EQ.0) GO TO 20
c dy/dlat
          G(2,IND+1)=-RADIUS*SINLAT*SINLON
c dy/dlon
          G(2,IND+2)=RADIUS*COSLAT*COSLON
c dy/drad
          G(2,IND+3)=COSLAT*SINLON
C
 20       IF(ILIVE(3).EQ.0) GO TO 100
c dz/dlat
          G(3,IND+1)=RADIUS*COSLAT
c dz/dlon
C         G(3,IND+2)=0.D0
c dz/drad
          G(3,IND+3)=SINLAT
       endif

      if (mode.eq.2) then
C
C PROPAGATE LOCAL SYSTEM (NORTH,EAST,UP) COVARIANCE MATRIX
C  FORMULATE JACOBIAN MATRIX
c du/dx
         GG(1,1)=-SINLAT*COSLON
c du/dy
         GG(1,2)=-SINLAT*SINLON
c du/dz
         GG(1,3)=COSLAT
c dv/dx
         GG(2,1)=-SINLON
c dv/dy
         GG(2,2)=COSLON
c dv/dz
         GG(2,3)=0.D0
c dw/dx
         GG(3,1)=COSLAT*COSLON
c dw/dy
         GG(3,2)=COSLAT*SINLON
c dw/dz
         GG(3,3)=SINLAT
      endif

 100  continue

      return
      end
