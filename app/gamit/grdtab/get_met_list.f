      Subroutine get_met_list( syear,sdoy,sdur,sitcod,slat,slon,lstn)
  

c     Read met values from a station (list) file 

c     Currently dummy until we write the files; pattern after get_map_list.f

c     R. King 10 August 2006


      implicit none
         
      include '../includes/dimpar.h'
      include '../includes/grdtab.h' 
      include '../includes/model.h'

      integer*4 syear,sdoy
      character*4 sitcod
      real*4 sdur,slat,slon    
      logical lstn

c     Temporary dummy statement to avoid compiler warnings
      if( lstn ) then
        print *,'DUMMY ',syear,sdoy,sdur,sitcod,slat,slon,lstn   
      endif
      
      call report_stat('FATAL','GRDTAB','get_met_list',' '
     . ,'Reading a met list file not yet coded',0)

      return
      end

