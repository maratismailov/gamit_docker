      Subroutine DMSDEG( dirflg,ideg,min,sec,deg )

c  Converts degrees, minutes, seconds to decimal degrees.
c  New version R. King 10 May 2003

      implicit none
                    
      character*1 dirflg 
      character*80 prog_name

      integer*4 ideg,min,rcpar,len

      real*8 sec,deg
       
c     get calling module name for report_stat
      len = rcpar(0,prog_name)
      
      if( ideg.lt.0.d0 ) 
     .   call report_stat('FATAL',prog_name,'lib/dmsdeg',' '
     .     ,'Negative deg in DMS conversion--fix code',0)
 
      deg = dble(ideg) + (dble(min)+sec/60.d0)/60.d0

      if( dirflg.eq.'W'.or.dirflg.eq.'S' ) deg = -deg      
      
      return
      end  
