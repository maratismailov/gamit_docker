      Subroutine rd_atml_grid

c     Read the header of an ANU-style grid file for atmospheric loading

c     Written by R. King  9 August 2006
c     using P. Tregoning code old utils/grdtab.f.  
                      

      implicit none   

      include '../includes/grdtab.h'

         
c Input unit numbers from grdtab.h                                        
c   luprnt   print file  
c   luatmlg  ATML grid file

c Output to grdtab.h
c   atmlver    r*4      version number of grid file
c   atmlgswap  logical  true if the grid file needs to be byte-swapped
c   atmlgmod   c*8      model name
c   atmlglat   i*4      number of latitude values on grid
c   atmlglon   i*4      number of longitude values on grid 
c   natmlg     i*4      number of components available from file (always 3)
c   atmlg_int  R*4      epoch interval of atml grid file

c Local variables
      integer*2 nversion
      integer*4 nrecl,ioerr
      character*2 spare
      character*5 press_src
      character*3 ref_frame
      character*256 message                
      integer*2 i2val,i2lon,i2lat
      integer*2 yr_start,yr_end
      integer*4 yr_start4,yr_end4
      real*4 doy_start,doy_end


c Open the grid file

      nrecl = 12
      open(unit=luatmlg,file='atml.grid',form='unformatted'
     .    ,access='direct',status='old',recl=nrecl,iostat=ioerr)
      if(ioerr.ne.0) then
        call report_stat('FATAL','GRDTAB','rd_atml_grid'
     .  ,'atml.grid','Error opening atm loading grid file:',ioerr)
      else
        call report_stat('STATUS','GRDTAB','rd_atml_grid'
     .  ,'atml.grid','Opened atm loading grid file:',ioerr)
      endif

c Read the first header record

      read(luatmlg,rec=1,iostat=ioerr) 
     .     nversion,spare,press_src,ref_frame
      if( nversion.gt.100.or.nversion.lt.0) then
c        If the grid size does not seem be correct, swap the bytes and see if OK.
        call swap_bytes(2, nversion, 1)
c         Now see if OK size.
        if( nversion.gt.100.or.nversion.lt.0) then
c          Still a problem so kill.
          call report_stat('FATAL','GRDTAB','rd_amtl_grid','atml_grid',
     .              'Invalid version number of atm loading grid file',0)
        else
*         File seems to be swapped.  Set the swapped status and continue.
*        (We now need to swap every record read)
         call report_stat('STATUS','GRDTAB','rd_atml_grid'
     .      ,'atml.grid','Byte-swapping atm loading  grid file',0)
                atmlgswap = .true.      
        endif 
      endif   
      atmlgmod(1:5) = press_src
      atmlgmod(6:8) = ref_frame 

c  Trap old-format (v1.0) files (repair by running update_atmlg)  
      if( nversion.lt.20 ) 
     .   call report_stat('FATAL','GRDTAB','rd_atml_grid','atml.grid'
     .  ,'Old-version of ATML grid file, fix with update_atmlg',0)

c  Read the next record to get array limits
c  PT060824: line 2 is now 
c     interval  ngval  nglon  nglat  (R*4  I*2  I*2  I*2 ) 

      read(luatmlg,rec=2,iostat=ioerr)atmlg_int,i2val,i2lon,i2lat
      if( atmlgswap )then
        call swap_bytes(4,atmlg_int,1)
        call swap_bytes(2,i2val,1)
        call swap_bytes(2,i2lon,1)
        call swap_bytes(2,i2lat,1) 
      end if            
      natmlg   = i2val
      atmlglon = i2lon
      atmlglat = i2lat

      if( natmlg.ne.3 ) then
        write(message,'(a,i6,a)') 'Number values on atm loading grid ('
     .     ,natmlg,') not equal 3: problem with swapping or format?'
        call report_stat('STATUS','GRDTAB','rd_atml_grid'
     .      ,'atml.grid',message,0)
      endif      

c PT060825: read the start/stop indicators on the 3rd header line
       read(luatmlg,rec=3)yr_start,doy_start,yr_end,doy_end
       if(atmlgswap)then
            call swap_bytes(2,yr_start,1)
            call swap_bytes(4,doy_start,1)
            call swap_bytes(2,yr_end,1)
            call swap_bytes(4,doy_end,1)
       endif
            

c  Write the header information to the print file ('grdtab.out')

      atmlgver = nversion
      atmlgver = atmlgver/10. 
      write(luprnt,'(/,a,f5.1,a,3i6)') 'Opened atml.grid version '
     .     ,atmlgver,'   # val lat lon ',natmlg,atmlglat,atmlglon
      write(luprnt,'(a,f8.4)') 'Interval (decimal days)',atmlg_int
      write(luprnt,'(a,1x,a5,1x,a3)') 
     .   'Pressure source  frame ',press_src,ref_frame      
      yr_start4 = yr_start
      call fix_y2k(yr_start4)
      yr_end4 = yr_end
      call fix_y2k(yr_end4)
      write(luprnt,'(a,2(i5,f9.5))')
     .   'Start,stop times',yr_start4,doy_start,yr_end4,doy_end
      return
      end

