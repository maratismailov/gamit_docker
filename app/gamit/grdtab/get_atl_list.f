      Subroutine get_atl_list ( site,slat,slon,lstn )

c     Read a station list file to get coefficients for atmospheric tidal loading 
c     a particular station.  

c     Currently dummy.  

      implicit none
   
      include '../includes/dimpar.h'
      include '../includes/grdtab.h'
      include '../includes/model.h'

      
      character*4 site
      real*4 slat,slon
      logical lstn

c     temporary dummy statement to avoid compiler warnings:
      if( lstn ) then
        print *,'DUMMY ',site,slat,slon,lstn
      endif

      call report_stat('FATAL','GRDTAB','get_atl_list',' '
     . ,'Reading an atl list file not yet coded',0)

      return
      end



