Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine OPERA( phi,iobs,iones,gearf,elev )
         
c     Generate the weight matrix for the DD observations
c     Calls CPMAT2 to fill the covariance matrix cphi(maxwm1)  (in solve.h)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      character*4 icall

      real*8  DUM(MAXOBS),PHI(MAXSIT,MAXSAT),elev(maxsit,maxsat)
     .     ,  SCALE(MAXPRM),gearf,small,gsq,gf,rcond

      integer IWORK(MAXPRM),iobs,iones,iones2
     .      , iobs2,ier,ijk
             
      logical debug/.false./   
c     for debug only: 
      integer jstat,jsat,ifac 

      equivalence  (isigma,iwork)
      small = 1.0d-11             

      if( debug) then 
        print *,'OPERA elev(2,4) ',elev(2,4) 
        print *,'OPERA nsite nsat phi ',nsite,nsat
        write(*,'(1x,6f12.3)') 
     .    ((phi(jstat,jsat),jsat=1,nsat),jstat=1,nsite)
      endif 

      IONES2=IONES*(IONES+1)/2
C                                        
      if(debug) print *,'OPERA aphi bphi ',aphi,bphi
      icall = 'data'
      if(debug) print *,'OPERA calling CPMAT2 elev(2,4) ',elev(2,4)
      call cpmat2(icall,phi, elev, cphi ) 
      IOBS2=IOBS*(IOBS+1)/2   
      if(debug) print *,'OPERA cphi(1-2) ',cphi(1),cphi(2)    
      if(debug)  print *,'  D(1-3) ',d(1),d(2),d(3)  
      CALL ATPA(D,IPNTD,IROWD,IOBS,IONES,CPHI,DPDT,WORK)
      if(debug) print *,'OPERA iobs2 dpdt',iobs2,(dpdt(ijk),ijk=1,iobs2)
      call copy1d(1,iobs2,0,dpdt,dpdtk)
      if(debug) print *,' dpdtk iobs ',iobs,(dpdtk(ijk),ijk=1,iobs)
      CALL NRMSC2(DPDTK,DUM,SCALE,IOBS,0,0)         
      if(debug) print *,' iobs rcond scaled dpdtk '
     .     ,iobs,rcond,(dpdtk(ijk),ijk=1,iobs)
      CALL INVER2(DPDTK,DUM,1,IOBS,rcond,IER)
      if( ier.ne.0 )call report_stat('WARNING','SOLVE','opera',' '
     .    , 'Bad inverse of DPDTK matrix',0)
      CALL NRMSC2(DPDTK,DUM,SCALE,IOBS,1,0) 
C     TRANSPOSE OPERATOR
      CALL TRANS(D,DT,IPNTD,IPNTDT,IROWD,IROWDT,IOBS,IONES,IWORK)
 
      if( debug) then  
        ifac = 4                                
        print *,'OPERA ifac iobs ',ifac,iobs
        write(*,'(2(1X,12F5.1))') (dt(ijk),ijk=1,ifac*iobs)
        write(*,'(2(1X,12I5))') (irowdt(ijk),ijk=1,iones+2)
      endif 

      CALL ATPA(DT,IPNTDT,IROWDT,IONES,IOBS,DPDTK,CPHI,WORK)
      if(debug)  print *,'OPERA  DD N(1-2) ',cphi(1),cphi(2)
C
C GLOBAL MULTIPLICATIVE FACTOR
      GSQ=gearf*gearf
      GF=1.D0+GSQ
C SCALE WEIGHT MATRICES BY GF FACTOR (SEE SCHAFFRIN&BOCK BG PAPER)
      DO 123 IJK=1,IONES2
      CPHI(IJK)=CPHI(IJK)/GF  
      if (dabs(cphi(ijk)).lt.small) cphi(ijk)=0.0d0
  123 CONTINUE
c      print *,' After 1/gf  DD N(1-2) ',cphi(1),cphi(2)

      RETURN
      END
