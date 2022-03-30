      subroutine avgclk ( nsites, iprn, ngoodf, lgoodf, outliers, c
     .                  , c1, c2, debug )

c     Average satellite clock rate estimates, discarding outliers
c     Currently hard-wired for a third-order polynomial
c     See Fiegl et al, 1991

c     R. King 25 June 92  from K. Feigl avclck (MODEL)

      implicit none

      include '../includes/dimpar.h'

c     Input--
c       nsites     : number of sets of estimated clock values to be considered
c       c(3,maxsit): third order polynomial coefficients estimated from nsites sites
c                      c0, c1, c2 in Feigl et al Eq. 6.
c                      Units are cycles, Hz, Hz/s
c     Input and output:
c       ngoodf      : number of sets of clock values actually used
c       lgoodf(maxsit): (T/F) whether for good values for site i (this epoch and PRN)

c     Output--
c       c1, c2     : averaged values of first- and second-order coefficients
c                      units are Hz, Hz/s
c                      (c0 is never averaged nor output)
c       outlier(maxsit,maxsat) : accumulated count of outliers

c     Must change in AVGCLK, J_FROM_C, RATE_EST, and ALLANV
      integer mm
      parameter (mm=3)

      logical lgoodf(maxsit),debug

      integer*4 nsites,iprn,ngoodf,outliers(maxsit,maxsat),i,j

      real*8 c(mm,maxsit),err,tol,rsum1,rsum2,rms1,rms2,c1,c2,xn,xn1

c     tolerance for rate deviation in Hz (0.15 Hz = 1 part in 10**10 ==>
c       very generous for any atomic oscillator
      data tol/0.15d0/


      if( debug) then
        print *,' AVGCLK:  nsites, ngoodf, c(1-3,j): ',nsites,ngoodf
        do j=1,nsites
         print *,' ',(c(i,j),i=1,3)
         enddo
       endif

c       Perform an initial mean and rms for outlier detection

      if ( ngoodf.ge.2 ) then
         xn  = ngoodf
         xn1 = ngoodf - 1
      else
         xn  = 1.d0
         xn1 = 1.d0
      endif
      rsum1 = 0.d0
      rsum2 = 0.d0
      rms1 = 0.d0
      rms2 = 0.d0
      do i=1,nsites
         if( lgoodf(i) ) then
           rsum1 = rsum1 + c(2,i)
           rsum2 = rsum2 + c(3,i)
         endif
      enddo
      c1 = rsum1/xn
      c2 = rsum2/xn
      do i=1,nsites
         if( lgoodf(i) ) then
           rms1 = rms1 + (c(2,i)-c1)**2
           rms2 = rms2 + (c(3,i)-c2)**2
           endif
      enddo
      rms1 = dsqrt(rms1/xn1)
      rms2 =  dsqrt(rms2/xn1)


c        Now compare the individual estimates with the mean and rms
c          Use the combination of the rate and acceleration terms (c1 and c2),
c          assuming a 1 second time-tag offset; i.e., c1 + c2*(1.)

      if( debug ) print *,' AVGCLK initial rms1, rms2: ',rms1,rms2
      ngoodf = 0
      do i=1,nsites
         if( lgoodf(i) ) then
            err = ( (c(2,i)+c(3,i)) - (c1 + c2) )
            if (debug ) print *,'  i, err : ',i,err
            if ( dabs(err).lt.tol ) then
               ngoodf = ngoodf + 1
               lgoodf(i) = .true.
            else
               outliers(i,iprn) = outliers(i,iprn) + 1
               lgoodf(i) = .false.
            endif
          endif
      enddo


c         Compute a new average and rms after removing outliers


      if ( ngoodf.gt.0 ) then

         if ( ngoodf.ge.2 ) then
            xn  = ngoodf
            xn1 = ngoodf - 1
         else
            xn  = 1.d0
            xn1 = 1.d0
         endif
         rsum1 = 0.d0
         rsum2 = 0.d0
         rms1 = 0.d0
         rms2 = 0.d0
         do i=1,nsites
            if( lgoodf(i) ) then
              rsum1 = rsum1 + c(2,i)
              rsum2 = rsum2 + c(3,i)
              endif
         enddo
         if ( debug ) print *,'  lgoodf: ',lgoodf
         c1 = rsum1/xn
         c2 = rsum2/xn
         do i=1,nsites
            if( lgoodf(i) ) then
              rms1 = rms1 + (c(2,i)-c1)**2
              rms2 = rms2 + (c(3,i)-c2)**2
            endif
         enddo
         rms1 = dsqrt(rms1/xn1)
         rms2 = dsqrt(rms2/xn1)
         if(debug) print *,'  Final ngoodf,rms1,rms2: ',ngoodf,rms1,rms2

      else

         if( debug) print *,'  Less than 1 good value in AVGCLK'

      endif

      return
      end
