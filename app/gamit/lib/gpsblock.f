c-----Temporary function to get block number for radiation pressure and yaw models

      integer*4 function gpsblock(antbody)
      
      implicit none
                
      integer*4 rcpar,len

      character*20 antbody 
      character*80 prog_name
      character*128 message
                   
c     get calling program name for report_stat
      len = rcpar(0,prog_name)

      gpsblock = 0       
      if( antbody(1:9).eq.'BLOCK I  ' ) then
         gpsblock = 1
      elseif( antbody(1:9).eq.'BLOCK II ' ) then
         gpsblock = 2
      elseif( antbody(1:9).eq.'BLOCK IIA' ) then
         gpsblock = 3
      elseif( antbody(1:11).eq.'BLOCK IIR-A' ) then
         gpsblock = 4
      elseif( antbody(1:11).eq.'BLOCK IIR-B' ) then
         gpsblock = 5
      elseif( antbody(1:11).eq.'BLOCK IIR-M' ) then
         gpsblock = 6
      elseif( antbody(1:9).eq.'BLOCK IIF' ) then
         gpsblock = 7
      elseif( antbody(1:10).eq.'BLOCK IIIA' ) then
         gpsblock = 7  ! DUMMY to get runs through
      else
        write(message,'(a,a20,a)') 'SV body type ('
     .   ,antbody,') not known' 
        call report_stat('FATAL',prog_name,'lib/gpsblock',' '
     .                   ,message,0)
      endif
      return
      end

