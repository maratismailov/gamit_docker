      Subroutine rd_map_grid

c     Read the header of an ANU-style (ATML) grid file for mapping functions and ZHD

c     Written by R. King  9 August 2006
c     using P. Tregoning code old utils/grdtab.f.  
                      

      implicit none  

      include '../includes/grdtab.h'

         
c Input unit numbers from grdtab.h                                        
c   luprnt   print file  
c   lumapg   Mapping function grid file

c Output to grdtab.h
c   mapgver    r*4      version number of grid file
c   mapgswap  logical  true if the grid file needs to be byte-swapped
c   mapgmod   c*8      model name
c   mapglat   i*4      number of latitude values on grid
c   mapglon   i*4      number of longitude values on grid 
c   nmapg     i*4      number of components available from file (3 for VMF1  grid including ZHD)
c   mapg_int  R*4      sampling interval (days) of the binary grid file
                             
c Function

      integer*4 nydays

c Local variables

      integer*2  nversion,imapg,isamp
      integer*4 grecl,ioerr
      character*256 message                
      real*4 start_time,stop_time

c Open the grid file

      grecl =  12
      open(unit=lumapg,file='map.grid',form='unformatted'
     .    ,access='direct',status='old',recl=grecl,iostat=ioerr)
      if(ioerr.ne.0) then
        call report_stat('FATAL','GRDTAB','rd_map_grid'
     .  ,'map.grid','Error opening mapping function grid file:',ioerr)
      else
        call report_stat('STATUS','GRDTAB','rd_map_grid'
     .   ,'map.grid','Opened mapping function grid file:',ioerr)
      endif

                    
c Read the first header record to get the version number
      mapgmod = ' '
      read(lumapg,rec=1,iostat=ioerr) nversion,mapgmod(1:4),mapgmod(7:8)
c      print*,'map.grid line 1: ',nversion,mapgmod
c check if it is a byte-swapped file
      if(nversion.lt.10.or.nversion.gt.100)then
        call swap_bytes(2,nversion,1)
        if(nversion.ge.10.and.nversion.lt.100)then
            call report_stat('STATUS','GRDTAB','rd_map_grid'
     .         ,'map.grid', 'Byte-swapping input grid file',0)
          mapgswap = .true.
        else
          call report_stat('FATAL','GRDTAB','vmftoasc'
     .      ,'map.grid','Invalid values read from binary header',0)
        endif
      else
        mapgswap = .false.
      endif

c Set the 5th character of the mapping function model to "G" to indicate
c that it comes from the gridded values (hence needs height correction in MODEL)
      mapgmod(5:5) = "G"

c Read the first header record to get array limits

      read(lumapg,rec=2,iostat=ioerr) isamp,imapg,mapglon,mapglat
      if( mapgswap )then
        call swap_bytes(2,isamp,1)
        call swap_bytes(2,imapg,1)
        call swap_bytes(4,mapglon,1)
        call swap_bytes(4,mapglat,1) 
      end if       
c      print*,'map.grid Line 2: ',isamp, imapg, mapglon, mapglat
c convert the I*2 header values (sampling interval in hours, number of variables) into
c the relevant types and units for grdtab.h variables
      mapg_int = isamp*1.0/24.0
      nmapg = imapg

      if( nmapg.ne.3 ) then
        write(message,'(a,i6,a)') 'Number values on map.grid ('
     .     ,nmapg,') not equal 3: problem with swapping or format?'
        call report_stat('STATUS','GRDTAB','rd_map_grid'
     .      ,'map.grid',message,0)
      endif      
      
c Third header line contains start and stop times (decimal years) of the file
      read(lumapg,rec=3)start_time,stop_time
      if(mapgswap)then
        call swap_bytes(4,start_time,1)
        call swap_bytes(4,stop_time,1)
      endif
cd     print*,'map.grid Line 3: ',start_time,stop_time
      mapg_start(1) = aint(start_time)
      mapg_start(2) = amod(start_time,1.0)*
     .     float(nydays(int(mapg_start(1))))+1.0
      call round6h(mapg_start(2)) 
      mapg_end(1) = aint(stop_time)          
      mapg_end(2) = amod(stop_time,1.0)*
     .     float(nydays(int(mapg_end(1))))+1.0
      call round6h(mapg_end(2)) 
cd      print *,'converted for grdtab.h: ',mapg_start,mapg_end

c  Write the header information to the print file ('grdtab.out')

      mapgver = nversion
      mapgver = mapgver/10. 
      write(luprnt,'(/,a,f5.1,a,3i6)') 'Opened map.grid version '
     .     ,mapgver,'   # val lat lon ',nmapg,mapglat,mapglon
      write(luprnt,'(a,i2,a)') 'Sampling Interval: ',isamp,' hours.'
      write(luprnt,'(a,1x,a8)') 
     .   'Mapping function and source: ',mapgmod

      return
      end




