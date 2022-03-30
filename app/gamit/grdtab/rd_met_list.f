      Subroutine rd_met_list

c     Read the header of the station list file for met values.
c     Currently dummy.


      implicit none     

      include '../includes/grdtab.h'

      call report_stat('FATAL','GRDTAB','rd_met_list',' '
     . ,'Reading a met list file not yet coded',0)

      return
      end


