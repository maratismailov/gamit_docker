      Subroutine rd_met_grid

c     Read the header of grid file for met values.
c     Currently dummy.


      implicit none  

      include '../includes/grdtab.h'

      call report_stat('FATAL','GRDTAB','rd_met_grid',' '
     . ,'Reading a met grid file not yet coded',0)

      return
      end


