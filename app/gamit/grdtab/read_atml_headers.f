      Subroutine read_atml_headers(lun,nhead,swapped,press_src,ref_frame
     .                            , interval,ngval,nglon,nglat
     .                            , yr_start,doy_start,yr_end,doy_end )

c     Read the headers of an ATML file
c     R King August 2010 from code in atmtoasc/atmfilt         
  
      implicit none
                                  
      integer*2 version,nhead,i2val,i2lon,i2lat
     .        , yr_start,yr_end

      integer*4 lun,ngval,nglon,nglat,ioerr

      real*4 interval,doy_start,doy_end,fversion

      character*2 spare
      character*3 ref_frame
      character*5 press_src

      logical swapped


      write(*,'(/,a,i3)') 'Reading headers for unit ',lun

c  Record 1:   
c
c  Version 1.0 or 2.0       version  spare pressure-source  reference-frame
c                             I*2     C*2     C*5              C*3
c  
c  Version 2.1 (same except:          nhead
c                                     I*2

c  First read just the version to check for byte swapping
              
      read(lun,rec=1,iostat=ioerr) version 
c      print*,'Initial read ',version
      if(version.gt.100.or.version.lt.0) then
*         If the grid size does not seem be correct, swap the bytes
*         and see if OK.   
        write(*,'(a,i4,a)') '**version =',version,' swapping**'
        call swap_bytes(2,version,1)
*         Now see if OK size.
        if(version.gt.100.or.version.lt.0) then 
           write(*,'(a,i4,a)') 'Invalid version (',version
     .            ,') of ATML grid file '
*             Still a problem so kill.   
           stop
        else
*         File seems to be swapped.  Set the swapped status and continue.
*         (We now need to swap every record read)  
          write(*,'(a)') 'Byte-swapping input grid file '
          swapped = .true.    
        end if
      endif 

c  Now read all of the first record

      if( version.eq.10 ) then  
        read(lun,rec=1,iostat=ioerr) version,spare,press_src,ref_frame  
        if( swapped ) call swap_bytes(2,version,1)
        nhead= 2    
        write(*,'(a,1x,i4,1x,a5,1x,a3)') 
     .    'Record 1 Version 1.0  version pressure-source frame: '
     .            ,version,press_src,ref_frame 
        write(*,'(a)') '--assume nhead = 2'
      elseif( version.eq.20 ) then   
        read(lun,rec=1,iostat=ioerr) version,spare,press_src,ref_frame    
       if( swapped ) call swap_bytes(2,version,1)
        nhead = 3      
        write(*,'(a,1x,i4,1x,a5,1x,a3)') 
     .    'Record 1 Version 2.0  version pressure-source frame: '
     .            ,version,press_src,ref_frame 
        write(*,'(a)') '--assume nhead = 3'
      elseif( version.eq.21 ) then 
        read(lun,rec=1,iostat=ioerr) version,nhead,press_src,ref_frame 
        if( swapped ) then          
           call swap_bytes(2,version,1)
           call swap_bytes(2,nhead,1)
        endif 
        write(*,'(a,2i3,1x,a5,1x,a3)') 
     . 'Record 1 Version 2.1  version nhead pressure-source frame: '
     .            ,version,nhead,press_src,ref_frame  
      else
         write(*,'(a,i2)') 'Unknown version ',version
         stop
      endif
                              

c  Record 2:  array limits  

c     Version 1.0 :  ngval   nglon   nglat
c                     i4       i4     i4

c     Version 2.0 :  interval  ngval   nglon   nglat
c                      r4        i2      i2      i2

      if( version.eq.10 ) then 
         read(lun,rec=2,iostat=ioerr) ngval,nglon,nglat
         if(swapped)then
            call swap_bytes(4,ngval,1)
            call swap_bytes(4,nglon,1)
            call swap_bytes(4,nglat,1)
        endif
        write(*,'(a,3i5,/,a)') 
     .       'Record 2  Version 1.0 ngval nglon nglat '
     .            ,ngval,nglon,nglat
     .           ,'--assume interval = 0.25 d'
        interval = 0.25  

      elseif( version.ge.20 ) then  
        read(lun,rec=2,iostat=ioerr) interval,i2val,i2lon,i2lat  
        if( ioerr.ne.0 ) then  
           write(*,'(a)') 'Error reading record 2 of grid file '
           stop    
        endif
        if(swapped)then
            call swap_bytes(4,interval,1)
            call swap_bytes(2,i2val,1)
            call swap_bytes(2,i2lon,1)
            call swap_bytes(2,i2lat,1)   
        endif    
        ngval = i2val
        nglon = i2lon
        nglat = i2lat
        fversion = float(version)/10.d0                 
        write(*,'(a,f3.1,a,1x,f4.2,1x,3i5)')
     .   'Record 2  Version ',fversion,' interval ngval nglon nglat'
     .                     ,interval,ngval,nglon,nglat        
      endif   
            

c  Record 3 (version 2.0ff only): start/stop times 

c     Version 2.0 : yr_start  doy_start  yr_end  doy_end
c                      i2       r4        i2      r4

      if( version.ge.20 ) then
        read(lun,rec=3,iostat=ioerr) yr_start,doy_start,yr_end,doy_end
        if( ioerr.ne.0 ) then
           write(*,'(a,i4)') 
     .             'Error reading record 3 of grid file, iostat=',ioerr
           stop   
        endif
        if( swapped ) then
cd           print *,'swapped '
           call swap_bytes(2,yr_start,1)
           call swap_bytes(4,doy_start,1)
           call swap_bytes(2,yr_end,1)
           call swap_bytes(4,doy_end,1) 
        endif                                
        write(*,'(a,1x,i2,1x,f6.2,2x,i2,1x,f6.2)') 
     .        'Record 3: yr/doy start yr/doy end'
     .       , yr_start,doy_start,yr_end,doy_end
      else  
        write(*,'(a,1x,i2,1x,f6.2,2x,i2,1x,f6.2)') 
     .        'Version 1.0, no Record 3; assume start/stop 1.0 365.75'
        doy_start = 1.0
        doy_end = 365.
      endif   

      return
      end



