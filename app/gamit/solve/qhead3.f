c
      Subroutine QHEAD3( jparm,symbl,sigtmp,kcoord,coslat,mode)

c    Print the adjustments, converting units when necessary


c     mode = 1: Parameter fixed, not sigmas printed
c     mode = 2: Parameter adjusted, print sigmas
c

C       Table of units

c   Parameter      Internal       Pre/Postfit        Adjust/Sigma
c  ---------------------------------------------------------------

c Site latitude       radians       deg min sec          meters
c     longitude       radians       deg min sec          meters
c      radius         kilometers    kilometers           meters

c Atm zenith delay    meters        meters               meters  
c     gradient        meters        meters               meters (N-S or E-w diff at 10 deg)

c Orbital elements    kilometers    kilometers           meters  

c SV antenna offset   meters        meters               meters


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      character*1 symbl
      character*16 buf16a,buf16b

      integer kcoord,mode,jparm

      real*8 sigtmp,erad,preprt,adjprt,sigprt,postprt
     .     , coslat,m2cyc,at10deg,fract

      data erad/6378137.d0/
c
      m2cyc = 1/(299792458.d0/1.57542d9)
c Changed as the G77 compiler does not support intrinsic trig functions working with degrees 
c      at10deg = 1/(sind(10.d0)*tand(10.d0)+0.003)
      at10deg = 1/(sin(10.d0/180.d0*pi)*tan(10.d0/180.d0*pi)+0.003)
      
      if( jparm.eq.210 ) then
c        print *,'QHEAD3 sigtmp ',sigtmp 
      endif

c     Convert the units and compute the fractional adjustment for printout

c      print *,'QHEAD3 jparm idms: ',jparm,idms(jparm)
c      print *,'  mode lpart ',mode,lpart
      preprt = preval(jparm)
      adjprt = adjust(jparm)
      sigprt = sigtmp
      postprt = postvl(jparm)
c     Convert the adjustments to meters
      if( rlabel(jparm)(6:10).eq.'RADIU' .or. 
     .     rlabel(jparm)(1:5).eq.'ORBIT' )
     .   then      
c         print *,'QHEAD3 jparm sigprt ',jparm,sigprt
         adjprt = adjprt * 1.d3
         sigprt = sigprt * 1.d3
      endif
      if( rlabel(jparm)(6:8).eq.'ATM') then
         preprt = preprt
         adjprt = adjprt
         sigprt = sigprt
         postprt = postprt 
c        for piecewise-linear adjustments, zero out the prefit to
c        avoid confusion (adjust is wrt average not wrt a priori)
         if( rlabel(jparm)(18:20).ne.'   ') then
             preprt = 0.d0
             postprt = adjprt
         endif 
      endif
      if(rlabel(jparm)(6:8).eq.'N/S'.or.
     .   rlabel(jparm)(6:8).eq.'E/W') then
         preprt = preprt*at10deg/m2cyc
         adjprt = adjprt*at10deg/m2cyc
         sigprt = sigprt*at10deg/m2cyc
         postprt = postprt*at10deg/m2cyc
      endif
c      print *,'QHEAD3 jparm idms erad: ',jparm,idms(jparm),erad
      if( idms(jparm).ne.0 ) then
c        lat, long adjustments and sigmas were in radians
         adjprt = adjprt * erad
         sigprt = sigprt * erad  
         if( rlabel(jparm)(11:13).eq.'LON') then
           adjprt = adjprt * coslat
           sigprt = sigprt * coslat
         endif
      endif
      if( mode.eq.2 ) fract = adjust(jparm) / sigtmp


      if (mode.eq.1) then

c        Parameters not in degrees, minutes, seconds
         if(idms(jparm).eq.0) then

c           revert to D-format if F-format field exceeded
            if( dabs(adjprt).le.1.d4 ) then
              if( logprt ) write(6,20)
     .        jparm,symbl,rlabel(jparm),preprt,adjprt
            else
              if( logprt ) write(6,21)
     .        jparm,symbl,rlabel(jparm),preprt,adjprt
            endif

            if (iqflag.eq.1) then
            if( dabs(adjprt).le.1.d4 ) then
              write(10,20)
     .        jparm,symbl,rlabel(jparm),preprt,adjprt
            else
              write(10,21)
     .        jparm,symbl,rlabel(jparm),preprt,adjprt
            endif
            endif

            if (ioflag.eq.1) write(15,21)
     .        jparm,symbl,rlabel(jparm),preprt,adjprt

         else
c           latitude or longitude: write ddd:mm:ss.sssss string
c            print *,' kcoord preprt ',kcoord,preprt
            call wdms(kcoord,preprt,buf16a)
            if( dabs(adjprt).le.1.d4 ) then
              if( logprt ) 
     .           write(6,30)  jparm,symbl,rlabel(jparm),buf16a,adjprt
            else
              if( logprt ) 
     .           write(6,31)  jparm,symbl,rlabel(jparm),buf16a,adjprt
            endif

            if (iqflag.eq.1) then
            if ( dabs(adjprt).le.1.d4 ) then
               write (10,30) jparm,symbl,rlabel(jparm),buf16a,adjprt
            else
               write (10,31) jparm,symbl,rlabel(jparm),buf16a,adjprt
            endif
            endif

            if (ioflag.eq.1)
     .         write (15,31) jparm,symbl,rlabel(jparm),buf16a,adjprt

         endif
      endif


      if (mode.eq.2) then

c       non-bias parameters:  check idms and print fract
        if( jparm.le.lpart ) then

          if(idms(jparm).eq.0) then

            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then
              if( logprt ) write(6,40) jparm,symbl,rlabel(jparm),
     .                   preprt,adjprt,sigprt,fract,postprt
              else
              if( logprt ) write(6,41) jparm,symbl,rlabel(jparm),
     .                   preprt,adjprt,sigprt,fract,postprt
              endif

            if (iqflag.eq.1) then
            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then
              write (10,40) jparm,symbl,rlabel(jparm),preprt
     .                        , adjprt,sigprt,fract,postprt
            else
              write (10,41) jparm,symbl,rlabel(jparm),preprt
     .                        , adjprt,sigprt,fract,postprt
            endif
            endif

            if (ioflag.eq.1) then
              write (15,41) jparm,symbl,rlabel(jparm),preprt
     .                    , adjprt,sigprt,fract,postprt
              endif


         else

            call wdms(kcoord,preprt,buf16a)
            call wdms(kcoord,postprt,buf16b)

            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then
              if( logprt ) write(6,50) jparm,symbl,rlabel(jparm)
     .                  , buf16a,adjprt,sigprt,fract,buf16b
            else
              if( logprt ) write(6,51) jparm,symbl,rlabel(jparm)
     .                  , buf16a,adjprt,sigprt,fract,buf16b
            endif

            if (iqflag.eq.1) then
            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then
              write(10,50) jparm,symbl,rlabel(jparm)
     .                       , buf16a,adjprt,sigprt,fract,buf16b
            else
              write(10,51) jparm,symbl,rlabel(jparm)
     .                     , buf16a,adjprt,sigprt,fract,buf16b
            endif
            endif

            if (ioflag.eq.1) then
            write(15,51) jparm,symbl,rlabel(jparm)
     .                 , buf16a,adjprt,sigprt,fract,buf16b
            endif

         endif

       else
c      bias parameter (jparm.gt.lpart):  don't print fract

            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then     
 
c        if( jparm.eq.210 ) then
c          print *,'QHEAD3 2 sigprt ',sigprt 
c        endif

        
              if( logprt ) write(6,60) jparm,symbl,rlabel(jparm),
     .                   preprt,adjprt,sigprt,postprt
              else
              if( logprt ) write(6,61) jparm,symbl,rlabel(jparm),
     .                   preprt,adjprt,sigprt,postprt
              endif

            if (iqflag.eq.1) then
            if( dabs(adjprt).le.1.d4.and.dabs(sigprt).le.1.d4 ) then
              write (10,60) jparm,symbl,rlabel(jparm),preprt
     .                        , adjprt,sigprt,postprt
            else
              write (10,61) jparm,symbl,rlabel(jparm),preprt
     .                        , adjprt,sigprt,postprt
            endif
            endif

            if (ioflag.eq.1) then
              write (15,61) jparm,symbl,rlabel(jparm),preprt
     .                    , adjprt,sigprt,postprt
              endif
       endif

      endif
c
c20   format (1x,i4,a1,a20,f17.10,1x,d12.5)
 20   format (1x,i4,a1,a20,f17.10,f11.4)
 21   format (1x,i4,a1,a20,f17.10,d11.4)
c30   format (1x,i4,a1,a20,1x,a16,1x,d12.5)
 30   format (1x,i4,a1,a20,1x,a16,f11.4)
 31   format (1x,i4,a1,a20,1x,a16,d11.4)
c40   format (1x,i4,a1,a20,f17.10,1x,d12.5,1x,d10.3,1x,a8,1x,f16.8)
 40   format (1x,i4,a1,a20,f17.10,2f11.4,f6.1,1x,f16.8)
 41   format (1x,i4,a1,a20,f17.10,2d11.4,f6.1,1x,f16.8)
c50   format (1x,i4,a1,a20,1x,a16,1x,d12.5,1x,d10.3,1x,a8,2x,a16)
 50   format (1x,i4,a1,a20,1x,a16,2f11.4,f6.1,1x,a16)
 51   format (1x,i4,a1,a20,1x,a16,2d11.4,f6.1,1x,a16)
 60   format (1x,i4,a1,a20,f17.10,2f11.4,7x,f16.8)
 61   format (1x,i4,a1,a20,f17.10,2d11.4,7x,f16.8)

c
      return
      end


