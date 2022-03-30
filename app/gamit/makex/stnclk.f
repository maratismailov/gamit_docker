       Subroutine stnclk ( debug,gnss,iwkntag,sowtag,iprn,prgl1
     .                   , coords,xsat,svclock
     .                   , rclock,iserr,warnings ) 
c
c     Compute the station receiver clock offsets by using
c     the broadcast ephemeris and apriori station coordinates
c
c     Written by Peter Morgan for the apollo, February 1987
c     based on a routine of same name by R. King
c
c     Modified to pass arrays better June 1987, K. Feigl

c     Modified for GNSS and separate navigation-file evaluation
c     into getnav.  R. King December 2015,
c
c     The meanings of the transfered parameters can be found in
c     these print statments.           
        
      implicit none       

      include '../includes/makex.h'


c     Input

c       debug   logical
c       gnss    char*1   GNSS G R C E J I (only G and R coded thus far)
c       iwkntag integer  GPS week number
c       sowtag  real*8   GPS seconds of week
c       iprn    integer  PRN number
c       prgl1   real*8   Observed L1 pseudorange (meters)
c       coords  real*8   Site coordinates (meters)
c       xsat    real*8   SV coordinates (meters0
c       svclock real*8   SV clock offset

      integer*4 iwkntag,iprn
      real*8 sowtag,prgl1,coords(3),xsat(3),svclock
      character*1 gnss
      logical debug     

c     Output
c       rclock  real*8    Computed receiver clock correction 
c       iserr   integer   Error encountered 
c       warnings logical Bad clocks detected 

      integer*4 iserr
       real*8 rclock 
       logical warnings

c     All times now GPST (Glonass conversion from UTC made in subroutine reade)

c     return iserr = 1 on error

c   Local     
      character*256 message
      integer*4 iter,iyr,mon,iday,imin,i
      real*8 tau,ltvel,gm,erate,t,xn,a
     .     , xsatdot(3),tprgl1,dt,testpr,range,sec,utcoff

      data ltvel   /2.99792458d+08/

c     start OK
      iserr = 0

c   Compute the theoretical delay
c
      tau=dsqrt( (xsat(1)-coords(1))**2 + (xsat(2)-coords(2))**2 
     .          +(xsat(3)-coords(3))**2 ) /ltvel 
cd      print*,' sat vector magnitude = '
cd     .    ,dsqrt(xsat(1)**2 + xsat(2)**2 + xsat(3)**2)
cd       print*,' site vector magnitude = ',dsqrt(coords(1)**2 + 
cd    .     coords(2)**2 + coords(3)**2)

c   Difference the theoretical delay and observed pseudorange

      tprgl1=prgl1/ltvel
c      print*,' in stnclk prgl1, ltvel, tprgl1 ',prgl1, ltvel, tprgl1,tau
c     compute the clock correction
      rclock= tprgl1 - tau + svclock

c     if debug needed print the following parameters
c
      if(debug) then
         print 900,iwkntag,sowtag,xsat(1),xsat(2),xsat(3)
     .           , coords
     .           ,tprgl1, tau, svclock, rclock
 900     format(/,1x,'at time tag wkntag= ',i4, ' and second count of ',
     1   f12.3,/,1x,'SV coords =   ',1p3e22.15, ' m,'/,1x,
     5   'site coords  = ',1p3e22.15,' m',/,1x,
     6   'tprgl1 observed range in seconds of time  = ',1pe22.15,/,1x,
     7   'tau, theoretical range in seconds of time = ',1pe22.15,/,1x,
     8   'diff of svclock from true time is           ',1pe22.15,/,1x,
     9   'diff of receiver clock from true time is  = ',1pe22.15)
      endif

c     to see that both observed and predicted distance to sat are
c     reasonable values (RWK 160817: allow synchronous SVs)
      range = tau*ltvel 
      if (range.lt. 1.1d7 .or. range.gt. 5.0d7) then
        warnings = .true.
        write(message,*) 'STNCLK: warning: calc range (m) = ',range
        if(uinfor.gt.0)  write(uinfor,*) message
        if(debug) write(*,*) message
        iserr = 1
      endif                  
      testpr = prgl1
      if( gnss.eq.'R') then    
c       get GPST-UTC
        testpr = prgl1 - utcoff*ltvel
      endif 
      if (testpr.lt.1.0d7.or.testpr.gt.4.0d7) then
        warnings = .true.
        write(message,*) 'STNCLK: warning: obs range (m) = ',testpr
        if(uinfor.gt.0 ) write(uinfor,*) message
        if(debug) write(*,*) message
* MOD TAH 930513: Problem with this warning when Trimble data with out
*        clock jumps is processed.  In these cases the ranges can become
*        very large but are still valid.  Fix at the moment is to not
*        set the error flag .i.e., only warning given.
C        iserr = 1  
* MOD TAH/RWK 160831: Still trap bad pseudoranges that are very small  
*         AMC2 2016 093
         if( dabs(testpr).lt.1000d0 ) iserr = 1  
c RWK 940107: print warnings only in the info file, with one message to
c        screen in main program
      endif   
* MOD RWK 010727: We're not successfully trapping all the bad values from TI 4100s
*        in lib/clkera (called by FIXDRV), so add back at least a crude check
*        which is too large to cause a problem with Trimbles.
* MOD RWK 151209: Skip this for Glonass, which will have a leap-sec-sized offset
      if( dabs(rclock).gt.1.d0 .and.gnss.ne.'R' ) then 
         if(debug) print *,'STNCLK: dabs(clock) > 1. ',rclock 
         iserr = 1 
      endif

      return

      end
