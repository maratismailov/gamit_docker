      Subroutine DMSRAD( dirflg,ideg,min,sec,rad )

c  Converts degrees, minutes, seconds to radians.
c  R. King 10 May 2003

      implicit none
                    
      character*1 dirflg 
      character*80 prog_name

      integer*4 ideg,min,rcpar,len    

      real*8 sec,deg,rad,pi
          
      pi = 4.d0*datan(1.d0)

c     get calling module name for report_stat
      len = rcpar(0,prog_name)
                     
      deg = dble(ideg) + (dble(min)+sec/60.d0)/60.d0   
      if( deg.lt.0.d0 ) 
     .   call report_stat('FATAL',prog_name,'lib/dmsrad',' '
     .     ,'Negative ideg in DMS conversion--fix code',0)
      rad = deg*pi/180.d0

      if( dirflg.eq.'W'.or.dirflg.eq.'S' ) rad = -rad      
      
      return
      end  
