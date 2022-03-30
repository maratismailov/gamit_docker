      Subroutine get_met_grid( syear,sdoy,sdur,slat,slon )

c     Interpolate a binary grid file to get met values for zenith delay calculations

c     Currently dummy until we write the files; pattern after get_map_grid.f

c     R. King 10 August 2006

      implicit none  
    
      include '../includes/dimpar.h'
      include '../includes/grdtab.h'
      include '../includes/model.h'
            
      integer*4 syear,sdoy
      real*4 sdur,slat,slon
                  

      call report_stat('FATAL','GRDTAB','get_met_grid',' '
     . ,'Reading a met grid file not yet coded',0)
          
c     Temporary dummy statement to avoid compiler warnings
      if( syear.gt.0 ) then
        print *,'DUMMY ',syear,sdoy,sdur,slat,slon   
      endif
      
  

      return
      end

    



