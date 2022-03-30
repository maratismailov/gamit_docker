Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine rotmat(angle,iaxis,r)
C
C Written by Yehuda Bock
C
c this routine computes a rotation matrix
c input: angle  - rotation angle
c        iaxis  - rotation axis
c output : r    - rotation matrix
c
      implicit none

      integer*4 iaxis,j,k
      real*8 angle,r,sinang,cosang

      dimension r(3,3)
c
      j=MOD(iaxis,3) + 1
      k=MOD(j,3) + 1
      r(iaxis,iaxis)=1.d0
      r(iaxis,j)=0.d0
      r(j,iaxis)=0.d0
      r(iaxis,k)=0.d0
      r(k,iaxis)=0.d0
      sinang=dsin(angle)
      cosang=dcos(angle)
      r(j,j)=cosang
      r(k,k)=cosang
      r(j,k)=sinang
      r(k,j)=-sinang
c
      return
      end
