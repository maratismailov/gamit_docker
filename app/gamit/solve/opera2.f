Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine OPERA2( phi,iobs,iones,gearf,elev )
    
c     Generate the weight matrix for the DD ionospheric constraints
c     Calls CPMAT2 to fill the covariance matrix cphik(maxwm1)   (in solve.h)



      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
c
      character*4 icall
C
      real*8 DUM(MAXOBS),PHI(MAXSIT,MAXSAT)
     .     , elev(maxsit,maxsat),SCALE(MAXPRM),gearf,small
     .     , gears,gsq,gf,rcond

      integer*4 IWORK(MAXPRM),iobs,iones,iones2,iobs2,ier,ijk,i

      equivalence  (isigma,iwork)
      small = 1.0d-11   
     

c      print *,'OPERA2 irowdt '
c      write(*,'(10i7)') (irowdt(i),i=1,50)
C
C
CD      WRITE(6,100) (D(I),I=1,4*IOBS)
CD 100 FORMAT(2(1X,12F5.1))
CD      WRITE(6,101) (IPNTD(I),I=1,4*IOBS)
CD 101 FORMAT(2(1X,12I5))
CD      WRITE(6,*) (IROWD(I),I=1,IOBS+2)
C
      IONES2=IONES*(IONES+1)/2
C
      icall = 'iono'
C FORM COVARIANCE MATRIX FOR IONOSPHERIC PARAMETERS
      CALL CPMAT2( icall,phi,elev,cphik )        
c       print *,'OPERA2 akappa bkappa ',akappa,bkappa
c       print *,'  CphiK(1-2) ',cphik(1),cphik(2)    
c       print *,'  D(1-3) ',d(1),d(2),d(3) 
CD     CALL MATPRT(CPHIK,WORK,IONES,1)
      IOBS2=IOBS*(IOBS+1)/2
      CALL ATPA(D,IPNTD,IROWD,IOBS,IONES,CPHIK,DPDTK,WORK)       
c      print *,'  DD dpdtk(1-2) ',dpdtk(1),dpdtk(2)

CD     CALL MATPRT(DPDT,WORK,IOBS,1)
CD     CALL MATPRT(DPDTK,WORK,IOBS,1)
C
C ADD COVARIANCE MATRICES
      GEARS=(1.D0+gearf*gearf)/(gearf*gearf)
      DO 50 I=1,IOBS2
CD     WRITE(6,*) DPDT(I),DPDTK(I)
CD     WRITE(16,*) DPDT(I),DPDTK(I)
      DPDTK(I)=DPDT(I)+GEARS*DPDTK(I)
   50 CONTINUE
C INVERT COVARIANCE MATRIX
      CALL NRMSC2(DPDTK,DUM,SCALE,IOBS,0,0)
      CALL INVER2(DPDTK,DUM,1,IOBS,rcond,IER)
      if( ier.ne.0 )call report_stat('WARNING','SOLVE','opera2',' '
     .    , 'Bad inverse of DPDTK matrix',0)
      CALL NRMSC2(DPDTK,DUM,SCALE,IOBS,1,0)
CD     CALL MATPRT(DPDTK,WORK,IOBS,2)
C
      CALL ATPA(DT,IPNTDT,IROWDT,IONES,IOBS,DPDTK,CPHIK,WORK) 
c      print *,'  DD Nk(1-2) ',cphik(1),cphik(2)

C
C GLOBAL MULTIPLICATIVE FACTOR
      GSQ=gearf*gearf
      GF=1.D0+GSQ
C SCALE WEIGHT MATRICES BY GF FACTOR (SEE SCHAFFRIN&BOCK BG PAPER)
      DO 123 IJK=1,IONES2
      CPHIK(IJK)=CPHIK(IJK)/GF
      if (dabs(cphik(ijk)).lt.small) cphik(ijk)=0.0d0
  123 CONTINUE            
c      print *,' After 1/gf  DD Nk(1-2) ',cphik(1),cphik(2)
c      stop

C
CD     CALL MATPRT(CPHIK,WORK,IONES,2)
C
C FOR ZERO IONOSPHERIC CONSTRAINT CPHI=CPHIK
CD     CALL MATPRT(CPHI,WORK,IONES,2)
C
       RETURN
       END
