      Subroutine get_atml_list ( syear,sdoy,sdur,sitcod,slat,slon,lstn)
                         
c     Extract atmospheric (non-tidal) loading values from a station list file

c     Currently dummy

      implicit none
                                      
      include '../includes/dimpar.h'
      include '../includes/grdtab.h'
      include '../includes/model.h'


      integer*4 syear,sdoy
      character*4 sitcod
      real*4 sdur,slat,slon
      logical lstn 
                         

      call report_stat('FATAL','GRDTAB','get_atml_list',' '
     . ,'Reading an atml list file not yet coded',0)

      return
      end


