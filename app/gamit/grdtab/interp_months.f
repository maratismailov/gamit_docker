      Subroutine INTERP_MONTHS( ngval,tval,sdoy,atides )

c   Interpolate to the requested date the S1 and S2 tidal coefficients
c   from an array (dimesion 144) containing values in the order specified 
c   by the LU/MIT grid file for 12 months of the year.  
c
c   R. King  31 January 2013
                         
c   The order of the cos/sin coefficents for U N E from the grid (tval vector) is
c        cS1u sS1u cS2u sS2u   cSs1n sS1n cS2n sS2n   cS1e sS1e cS2e sS2
c           x 12 months
                     
c   The order expected for the atides(2,6) array in grdtab and model is
c         cS1u sS1u  cS1n sS1n  cS1e sS1e  
c         cS2u sS2u  cS2n sS2n  cS2e sS2e   
c     

      implicit none
               
      integer*4  maxtid,maxepc
      parameter (maxtid=12,maxepc=12)

c       ngval: number of coefficients for each lat/lon grid point (record of file)
c       tval(ngval): values from grid (12 coeffs x 12 months) 
c              coefficients sin/cos S1 U, sin/cos S1 N, sin/cos S1 E
c                           sin/cos S2 U, sin/cos S2 N, sin/cos S2 E
c       ntid: number of tidal coefficients for each epoch (usually ngval/12) 
c       t:  time in fraction of year from 1 January
c       s1(3): U N E S1 tide in mm
c       s2(3): U N E S2 tide in mm

   
      integer*4 ngval,ntid,nepc,im1,im2,i,j,ii
      integer*4 j1,j2 
                                
c     for now assume ngval is 144, w/ 12 tide values at 12 monthly epochs
      integer*4 sdoy
      real*4 tval(ngval),coeff(maxtid,maxepc),ttval(maxtid),dm
     .     , s1(3),s2(3),tyr,tm,fract,slope,twopi 
     .     , atides(2,6)
      logical debug/.false./
                                                         
      data twopi/ 6.283185/ 
   
c Move the full line of coefficients into a 2-D array 
      if( ngval.ne.144) then
        write(*,'(i5,a)')  ngval,
     .    'ngval ne 144, need to reset dimensions and logic'
        stop
       endif
      ntid = 12
      nepc = 12   
      ii = 0  
      do j=1,nepc
        do i=1,ntid
          ii = ii + 1
cd          if(ii.eq.1)
cd     .      write(*,'(a,f7.4)') 'Interpolated Jan 1st coeff ',tval(1)
          coeff(i,j) = tval(ii)
        enddo
      enddo
                       
cd      print *,'INTERP_MONTHS tval '
cd      do i=1,12                                      
cd        j1= 12*(i-1)+1
cd        j2= 12*(i-1)+12
cd        write(*,'(12f7.4))') (tval(j),j=j1,j2)  
cd      enddo       
c Linearly interpolate the monthly averages to the date needed
c    --averages are assumed to best represent mid-month
c    --assume monthly, hence dm = 1.0
      dm = 1.0  
      tyr = float(sdoy)/365.
      tm = tyr*12.   
      if( tm.le.0.5) then 
c       if before 15 January, set the beginning month to December
        im1 = 12
        im2 =  1
      elseif( tm.ge.11.5) then
c       if after 15 December, set the end month to January
        im1 = 12
        im2 =  1
      else      
        im1 = ifix(aint(tm+0.5))
        im2 = im1 + 1
      endif
      fract = amod(tm+0.5,1.0)
      do i=1,ntid              
        if(i.eq.1) print *,'coeff 21 ',coeff(i,im2),coeff(i,im1)
        slope = (coeff(i,im2)-coeff(i,im1))/dm  
        ttval(i) = coeff(i,im1) + slope*fract
        if( debug ) 
     .     print *,'i im1 im2 dm fract slope coeff(1,im2) coeff(1,im1) '
     .   ,'ttval(1) ',i,im1,im2,dm,fract,slope,coeff(1,im2),coeff(1,im1)
     .   , ttval(1)
cd        stop
      enddo 
      if( debug ) then
        print *,'sdoy tyr tm im1 im2 slope fract '
     .   , sdoy,tyr,tm,im1,im2,slope,fract
        write(*,'(a,12f6.3)') 'tval ',(ttval(i),i=1,12)
      endif

c Fill the atides array 
       
c     The order of the cos/sin coefficents for U N E from the grid (tval vector) is
c        cS1u sS1u cS2u sS2u   cSs1n sS1n cS2n sS2n   cS1e sS1e cS2e sS2
c           x 12 months
                     
c     The order expected for the atides(2,6) array in grdtab and model is
c         cS1u sS1u  cS1n sS1n  cS1e sS1e  
c         cS2u sS2u  cS2n sS2n  cS2e sS2e   

      atides(1,1) = ttval(1)
      atides(1,2) = ttval(2)
      atides(1,3) = ttval(5)
      atides(1,4) = ttval(6)
      atides(1,5) = ttval(9)
      atides(1,6) = ttval(10)
      atides(2,1) = ttval(3)
      atides(2,2) = ttval(4)
      atides(2,3) = ttval(7)
      atides(2,4) = ttval(8)
      atides(2,5) = ttval(11)
      atides(2,6) = ttval(12)

      return
      end


      


