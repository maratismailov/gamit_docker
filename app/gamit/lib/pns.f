Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine PNS(prec,rnut,srot,pnst)
C
C Written by Yehuda Bock
C
C Form PNS transformation matrix from earth-fixed to inertial system

      implicit none

      integer*4 i,j

      real*8 nt,nst,prec,rnut,srot,pt,st,pnst
C
      dimension prec(3,3),rnut(3,3),srot(3,3),
     1          PT(3,3),NT(3,3),ST(3,3),NST(3,3),PNST(3,3)
C
      DO 10 I=1,3
      DO 10 J=1,3
      PT(I,J)=0.D0
      NT(I,J)=0.D0
      ST(I,J)=0.D0
      NST(I,J)=0.D0
      PNST(I,J)=0.D0
  10  CONTINUE
C
      call transp(PREC,PT,3,3)
      call transp(RNUT,NT,3,3)
      call transp(SROT,ST,3,3)
c     write(*,703) ((PT(i,j),j=1,3),i=1,3)
c 703 format(' PT-rotation matrix :',/,3(1x,3D22.14,/))
c     write(*,704) ((NT(i,j),j=1,3),i=1,3)
c 704 format(' NT-rotation matrix :',/,3(1x,3D22.14,/))
c     write(*,705) ((ST(i,j),j=1,3),i=1,3)
c 705 format(' ST-rotation matrix :',/,3(1x,3D22.14,/))
      call matmpy(NT,ST,NST,3,3,3)
      call matmpy(PT,NST,PNST,3,3,3)
C
      return
      end
