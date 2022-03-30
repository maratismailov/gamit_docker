      Subroutine interp_azel( zen,az,nzen,naz,dzen,dazi,zenmax,table
     .                      , pcvout )

c     Interpolate from a 2-d table of pcv values for elevation and azimuth

      include '../includes/dimpar.h'

      integer*4 nzen,naz,ix1,ix2,iy1,iy2,rcpar,len

      real*8 zen,az,dzen,dazi,table(maxel,maxaz),zen11,az11
     .     , val11,val12,val21,val22,pcvout,zenmax 
                             
      character*80 prog_name
      character*256 message
    
c Get calling program name and X-file name for report_stat
      len = rcpar(0,prog_name)
  

c  Find the corners of the box to be interpolated

      ix1 = int(zen/dzen) + 1
      ix2 = ix1 + 1  
c** rwk 190827:  Need this to work for zen2 < 90
c        if( dabs(zen-90.d0).lt..001d0 ) then
      if( dabs(zen-zenmax).lt..001d0 ) then 
        ix2 = nzen
        ix1 = nzen -1 
      elseif (dabs(zen).lt.001d0 ) then
        ix1 = 1
        ix2 = 2
      endif
      if( ix1.lt.1.or.ix2.gt.nzen ) then  
         write(message,'(a,f6.2,f5.1,3i3)') 
     .   'Error in bilin interp, zen dzen ix1 ix2 nzen '
     .      ,zen,dzen,ix1,ix2,nzen
         call report_stat('FATAL',prog_name,'lib/interp_azel',' '
     .       ,message,0)  
      endif   
      zen11 = dzen*dfloat(ix1-1)   
      iy1 = int(az/dazi) + 1
      iy2 = iy1 + 1      
      if( dabs(az-360.d0).lt..001d0 ) then
        iy2 = naz
        iy1 = naz -1
      elseif (dabs(az).lt.001d0 ) then
        iy1 = 1
        iy2 = 2
      endif
      if( iy1.lt.1.or.iy2.gt.naz ) then  
         write(message,'(a,f6.2,f5.1,3i3)') 
     .     'Error in bilin interp, az dazi iy1 iy2 naz'
     .     , az,dazi,iy1,iy2,naz
         call report_stat('FATAL',prog_name,'lib/interp_azel'
     .                   ,' ',message,0)     
      endif
      az11 = dazi*dfloat(iy1-1)
      val11 = table(ix1,iy1)
      val12 = table(ix1,iy2)
      val21 = table(ix2,iy1)
      val22 = table(ix2,iy2) 

c Do the interpolation

      call bilin8( zen,az,zen11,az11,dzen,dazi
     .          , val11,val21,val22,val12,pcvout )  
c      print *,'zen az zen11 az11 dzen dazi crds val pcvout ',
c     . zen,az,zen11,az11,dzen,dazi,ix1,iy1,ix2,iy2,val11,val12
c     . ,val21,val22,pcvout
c      if( zen.lt.81.d0 .and.az.gt.31.d0 ) stop
      return
      end

      subroutine bilin8(x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val)

c   computes a bilinear interpolation of the value of a single
c   point within a rectangle. I found this formula via google
c   and it seems to work!
c
c   INPUT:
c         x,y  :  the coordinates of the point to be interpolated
c        x0,y0 :  the coords of the lower left corner
c     xstep,ystep : the grid step sizes
c         u1   : value at lower left corner
c         u2   : value at lower right corner
c         u3   : value at upper right corner
c         u4   : value at upper left corner
c
c    OUTPUT:
c        val   : the interpolated value at the requested point
c
c   P. Tregoning
c   7 January 2004

      implicit none

c  argument variables
      real*8 x,y,x0,y0,xstep,ystep,u1,u2,u3,u4,val

c  local variables
      real*8 dy,dx,dy1,dx1

          dy  = abs(y-y0)/ystep
          dy1 = (ystep - dy*ystep)/ystep
          dx  = (x - x0)/xstep
          dx1 = (xstep - dx*xstep)/xstep
c then
          val = dx1*dy1*u1+dx*dy1*u2+dx1*dy*u4+dx*dy*u3

      return
      end








