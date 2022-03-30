Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine init( isat )
c
c     initialize parameters necessary for gps orbit integrator
c     Rick Abbot - november 1984
c     modified Mmarch 1987 by R.King and Y.Bock (from R.Abbot update)
c       to allow option of WGS72, WGS84, or MERIT (GEM-l2) constants.
c     modified May 1994 by R. King and Y. Bock to allow IERS92 constants
c     modified April 2017 by R. King to remove WGS72 and MERIT constants
c
      implicit none

      include '../includes/dimpar.h'  
      include '../includes/global.h'
      include '../includes/arc.h'


      integer*4 nhmn,nhmx,isat,nhcout,nhc,k,i

      character*5 upperc
      character*256 message

c  Function
      real*8 taiutc

      logical debug/.false./   

c  Get the offset between the integration time type and TDT for interpolating the Sun, Moon
c     (formerly in lib/ghdred, now obsolete)

      if( time_type.eq.'UTC ' ) then
        tdtoff= 32.184d0 + taiutc(jde)
      elseif( time_type.eq.'GPST') then
        tdtoff= 32.184d0 + 19.0d0
      endif
     
   
c  Get the stop, stop times for the integration

c    Compute the actual start start and stop times for the integration by 
c    expanding the observation times by 11 tabular intervals and rounding
c    to the nearest tabular interval. Originally 7 tabular intervals
c    was enough, but now we need to be able to interpolate the tfile
c    up to 1 hour before the start of a session to handle satellites
c    that begin the session in eclipse. McClusky 950405   
c    Occasional peculiar eclipse sequence means that we need 13 tabular 
c    intervals. Matt King/Bob King 010615/010717  
c    Another case found where 13 is not enough; increase to 14. R King 030721
      if( isat.eq.1 ) then 
         if(debug) print *,'INIT jdb tb jde te jdf tf '
     .                       ,   jdb,tb,jde,te,jdf,tf 
         trun0= (jdb-jde)*86400.d0 + (tb-te)
         trunf= (jdf-jde)*86400.d0 + (tf-te)
         if(debug) print *,'INIT trun0 trunf delt ',trun0,trunf,delt
* MOD TAH 210202: Changed limit to make model's expectation (plus one extra
*        tabular point; 5*sdelt in model).  Was 14 delt here)
*        Original code say 14*900s= which is 3.5 hrs and far excceds the 
*        additional 1hr needed mentiones in teh comments above.  Why is
*        not clear.  Code tested with GPS 2020 day 5 with G026 beta 
*        angle -0.81deg.  The 3600.d0 seconds should be a multiple of 
*        delt (Tabular interval).
C        trun0= trun0 - dmod(trun0,delt) - 14.d0*delt
C        trunf= trunf - dmod(trunf,delt) + 14.d0*delt  
         trun0= trun0 - dmod(trun0,delt) - 3600.d0 - 6*delt
         trunf= trunf - dmod(trunf,delt) + 3600.d0 + 6*delt  
         if(debug) print *,' INIT trun0 trunf ',trun0,trunf
         jdb= jde
         tb = te
         call timinc( jdb, tb, trun0 )
         jdf= jde
         tf = te
         call timinc( jdf, tf, trunf )
         if(debug) print *,'INIT jdb tb jdf tf trun0 trunf '
     .                    ,jdb,tb,jdf,tf,trun0,trunf 
      endif


c  Direction of integration

      if (trunf.gt.0.d0) nsign=+1
      if (trunf.lt.0.d0) nsign=-1

       
c  Set the integration step parameters

c       integration step-size
      stdint=168.75
      if (diint.eq.0.0) diint = stdint
      if (delt.eq.0.0) delt = 1350.d0     
c        = positive, step size is nh seconds
c       <= zero, step size is 2**nh seconds
      nhc = 12
      hc = 2.d0**(-nhc)*86400.d0
      nhc = 8
      hc = diint/nhc
c       minimum step size for integrator = hmn
      nhmn = 32
      hmn = diint/nhmn
c       maximum step size for integrator = hmx
      nhmx = delt/diint
      hmx = delt/nhmx
      dinthmx = hmx
c        = positive,tabular interval is inthmx seconds
c       <= zero, tabular interval is 2**inthmx seconds
      dnh = hmx
      npstep=0
      nrec=0
      iptr3=3
      l1=0 
c        integration steps to take before tabular output = nstout
c      nstout=2**(nhmx-nhcout)
       nstout = delt/hmx
                 
      if(debug) print *,'INIT nsign hmx hmn hc nstout '
     .       ,nsign,hmx,hmn,hc,nstout

c Get the number of first order differential equations to integrate 
      neq=6  
      if ( apar.eq."N" .or. apar.eq."n" ) goto 502
      if (nics.ne.15 .and.nics.ne.19 ) then
       write(message,'(a,i3,a)') 'nics = ',nics,'; must be 15 or 19 '
        call report_stat('FATAL','ARC','init',' ', message,0)
      endif
c     neq = 96 assumes 6 for coords and 6*(15 partials)  ECOM1(formerly BERNE), UCLR1, UCLR2 models
c     neq = 120 adds the addional 4 partials for the ECOM2/ECOMC model
c     to generalize, must add parameter to the control file to designate
c     the partials to be integrated                                      
      neq = 60
      if (modrad.eq.1.or.modrad.eq.7.or.modrad.eq.8) neq = 96
      if( modrad.eq.2) neq = 120 
  502 kount=neq/6-1      


c Stability and accuracy criterion (units=km)
      epsi=1.0d-05
c       number of differential equations to include in the stability and accuracy tests 
      meq=6
             

c  Order for predictor-corrector

      ncoret=11
      npredt=11    

cd      print *,npredt,ncoret,epsi,neq,kount,dnh,hmx,nhmx,nhmn,nhc,trun0
cd     .      ,trunf,diint,delt,nsign,tb,tf,tdtoff,jde
                  

c  Initialize the  radiation pressure coefficients

      do i=1,13
       radcon(i) = 0.d0
      enddo
      write(iarh,2990) srpmod
 2990 format(/,1x,a5,': Solar radiation pressure model coefficients ')
c     Set 
      radcon(1)=satics(7)
      radcon(2)=satics(8)
      radcon(3)=satics(9)
c     all current models (ECOM1, ECOM2, UCLR1, UCLR2 use at least constant and once-per-rev direct, Y, B
      radcon(4)=satics(10)
      radcon(5)=satics(11)
      radcon(6)=satics(12)
      radcon(7)=satics(13)
      radcon(8)=satics(14)
      radcon(9)=satics(15)
      if(modrad.ne.2) then
        write (iarh,3000) sname(1:10),(icsnam(i),i=7,15),radcon
 3000   format(1x,a7,3x,9(a4,8x),/,7x,9(d11.4,1x),/)
      else
c       ECOM2 has four extra adjustable parameters
        radcon(10)=satics(16)
        radcon(11)=satics(17)
        radcon(12)=satics(18)
        radcon(13)=satics(19)         
        write (iarh,3001) sname(1:10),(icsnam(i),i=7,19),radcon
 3001   format(1x,a7,3x,13(a4,8x),/,7x,13(d11.4,1x),/)
      endif


c Initialize the equations of motion and partials

      do 1000 k=7,maxyt2
 1000 satics(k)=0.d0
c        dx/dx0
      satics(7)=1.d0
c        dy/dy0
      satics(14)=1.d0
c        dz/dz0
      satics(21)=1.d0
c        dxdot/dxdot0
      satics(28)=1.d0
c        dy/dydot0
      satics(35)=1.d0
c        dz/dzdot0
      satics(42)=1.d0
c    
cd      print *,'INIT satics: ',(satics(i),i=1,6)

      return
      end
