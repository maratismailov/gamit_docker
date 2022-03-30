      Subroutine RADDMS( rad,dirflg,ideg,min,sec )

c     Convert radians to deg/min/sec
c     R. King 10 May 2003

c     Note: the direction flag is set '+' or '-' since we don't
c           know if this is lat or lon.  The calling progam 
c           should change this to 'N', 'S', 'E', or 'W'.

      implicit none

      integer*4 ideg,min

      real*8 deg,dmin,sec,pi,rad
      
      character*1 dirflg      

         
      pi = 4.d0*datan(1.d0)

      if( rad.lt.0.d0 ) then
        dirflg = '-'
      else
        dirflg = '+'  
      endif  
      deg = dabs(rad)*180.d0/pi
      ideg = int(deg)
      dmin = (deg-dble(ideg))*60.d0
      min = int(dmin)
      sec =  (dmin-dble(min))*60.d0
      return
      end


