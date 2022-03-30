Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine SNP(prec,rnut,srot,snpmat)
C
C Written by Yehuda Bock
C
C Form SNP transformation matrix from inertial to earth fixed system

      implicit none

      real*8 np,prec,srot,rnut,snpmat
C
      dimension prec(3,3),rnut(3,3),srot(3,3),snpmat(3,3),np(3,3)
C
      call matmpy(rnut,prec,np,3,3,3)
      call matmpy(srot,np,snpmat,3,3,3)
C
      return
      end
