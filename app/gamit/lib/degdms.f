      Subroutine DEGDMS( deg,dirflg,ideg,min,sec )

c     Convert decimal degrees to deg/min/sec
c     New version by R. King 10 May 2003

c     Note: the direction flag is set '+' or '-' since we don't
c           know if this is lat or lon.  The calling progam 
c           should change this to 'N', 'S', 'E', or 'W'.

      implicit none

      integer*4 ideg,min

      real*8 deg,dmin,sec
      
      character*1 dirflg
                 
      if( deg.lt.0.d0 ) then
        dirflg = '-'
        deg = -deg
      else
        dirflg = '+'
      endif 
      ideg = int(deg)    
      dmin = (deg-dble(ideg))*60.d0     
      min = int(dmin)   
      sec =  (dmin-dble(min))*60.d0  
      return
      end


