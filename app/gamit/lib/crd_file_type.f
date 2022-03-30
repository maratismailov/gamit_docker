      Subroutine crd_file_type( fname, kfflg )

c     Read the coordinate file and determine if it's a spherical (fixed-format)
c     or cartesian (globk free format) l-file
c     R. King 030214

      implicit none

      
      character*(*) fname 
      character*1 alat,alon
      character*10 arad
      character*80 prog_name

      integer*4 kfflg,ioerr,len,rcpar,iul  

      real*8 rad

      logical file_open,unit_open

c      type of file
c      kfflg = 0 : lfile
c            = 1 : apr file 
       kfflg = -1  
       iul = 0
        
c     get calling program for report_stat
      len = rcpar(0,prog_name)

c     see if the file is open 

      inquire(file=fname,number=iul,opened=file_open,iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .  ,'lib/crd_file_type',fname
     .  ,'Error inquiring status of coordinate file',ioerr)  
c      print *,'DEBUG inquire lfile open ',file_open,iul
      if( .not.file_open ) then  
c       find an available unit number between 60 and 99
        iul = 60 
        unit_open = .true.
        do while (unit_open .and. iul.lt.100 ) 
          inquire(unit=iul,opened=unit_open)
          if( unit_open ) iul = iul + 1
        enddo 
        if(unit_open) call report_stat('FATAL',prog_name
     .       ,'lib/crd_file_type',fname
     .       ,'Cannot find unit to open for coordinate file',0)  
        open(file=fname,unit=iul,iostat=ioerr)  
        if( ioerr.ne.0) 
     .      call report_stat('FATAL',prog_name,'lib/crd_file_type'
     .               ,fname,'Error opening coordinate file',ioerr)    
c         print *,'DEBUG lfile temp open iul ',iul
      else
        rewind(unit=iul,iostat=ioerr) 
        if( ioerr.ne.0) 
     .      call report_stat('FATAL',prog_name,'lib/crd_file_type'
     .               ,fname,'Error rewinding coordinate file',ioerr)    
      endif

c     now read the file looking for a valid spherical l-file line
                                           
      do while ( kfflg.lt.0 )
c       look for N or S or blank in column 18, W or E or blank in col 34 and a valid radius
        read(iul,'(17x,a1,15x,a1,16x,a10)',iostat=ioerr) alat,alon,arad   
c        print *,'DEBUG read alat along arad ',alat,alon,arad  
c        write(*,'(1x,a1,1x,a1,1x,a10)') alat,alon,arad
        if( ioerr.eq.-1 ) then
c          eof reached without finding a spherical entry--must be apr file
           kfflg = 1   
        elseif ( ioerr.ne.0 ) then
c         non-eof error reading file 
          call report_stat('FATAL',prog_name,'lib/crd_file_type'
     .               ,fname,'Error reading coordinate file',ioerr)  
        else
          if( alat.eq.'N'.or.alat.eq.'S'.or.alat.eq.' ') then
            if(alon.eq.'E'.or.alon.eq.'W'.or.alon.eq.' ') then
              read(arad,'(f10.0)',iostat=ioerr) rad
               if(ioerr.eq.0.and.
     .           (rad.gt.6.3d6.and.rad.lt.6.4d6) ) then
                 kfflg = 0 
c                print *,'DEBUG line items ',n,alat,alon,rad,kfflg   
               endif
            endif
          endif
        endif  
      enddo

c     leave the file in the right disposition
      if( file_open ) then  
        rewind(unit=iul)
      else
        close(unit=iul)
      endif     
c** debug 
c      inquire(file=fname,number=iul,opened=file_open,iostat=ioerr)
c      print *,'DEBUG lfile now open, iul, kfflg ',file_open,iul,kfflg

      return
      end


 
