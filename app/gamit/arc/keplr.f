Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
       SUBROUTINE KEPLR(X,A,E,OMLIT,OMBIG,FINC,ANOMM)
C
C CONVERTS CARTESIAN ELEMENTS TO KEPLERIAN ELEMENTS
C      
      implicit none  

      include '../includes/dimpar.h'   
      include '../includes/arc.h'


       real*8 x,a,e,omlit,ombig,finc,anomm,h,ev,vnu,eta,zv,ecan,r,dot
       integer*4 i

       DIMENSION X(6),H(3),EV(3),VNU(3),ETA(3)
       DIMENSION ZV(3)     


c     no longer necessary to define gm here (FMU), as it is done init and
c     past in through the common const. FMU replaced everywhere with gm(1)
C	   IT IS NOT TOO IMPORTANT WHAT THE VALUE OF FMU (THE GRAVITATIONAL
C	   MASS OF THE EARTH) IS HERE SINCE IT IS ONLY TO USED TO COMPUTE THE
C	   KEPLERIAN ELEMENTS WHICH ARE THEN ONLY VISUALLY EXAMINED.  IN
C	   FACT FMU MAY NOT BE THE SAME AS GM(1) GIVEN IN SUBROUTINE INIT
C	DATA FMU/398600.64D0/
C      DATA FMU/398600.8D0/
      DATA ZV/0.D0,0.D0,1.D0/

      CALL CROSS(X(1),X(4),H(1))
      CALL CROSS(X(4),H,EV)
      R=DSQRT(DOT(X(1),X(1)))
      DO 710 I=1,3
  710 EV(I)=EV(I)/gm(1)-X(I)/R
C
      E=DSQRT(DOT(EV,EV))
      A=1.D0/(2.D0/R-DOT(X(4),X(4))/gm(1))
C
      CALL CROSS(ZV,H,VNU)
      CALL CROSS(H,VNU,ETA)
C
      OMLIT=DATAN2(DOT(EV,ETA)/DSQRT(DOT(H,H)),DOT(EV,VNU))
      OMBIG=DATAN2(VNU(2),VNU(1))
      IF(OMLIT .LT. 0.D0) OMLIT=OMLIT+TWOPI
      IF(OMBIG .LT. 0.D0) OMBIG=OMBIG+TWOPI
C
      FINC=DATAN2(DSQRT(H(1)**2+H(2)**2),H(3))
      IF(FINC .LT. 0.D0) FINC=FINC+TWOPI
C
      ECAN=DATAN2(DOT(X(1),X(4))/DSQRT(gm(1)*A) ,
     $     R*DOT(X(4),X(4))/gm(1)-1.D0)
C
      ANOMM=ECAN-E*DSIN(ECAN)
      IF(ANOMM.LT.0.D0) ANOMM=ANOMM+TWOPI
C
      RETURN
      END
