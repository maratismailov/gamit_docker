      Subroutine rdsvclk( iusvclk,iprnt,jd,sod,nsat,gnss,itsat,svclk )

c     Read the SV clock corrections from a RINEX-style clock file
c     R. W. King  March 2000    
                      
c     Input
c        iusvclk       logical unit number of RINEX-style clock file
c        iprnt         logical unit number of ttongs.out (print) file
c        jd            PEP JD of requested epoch
c        sod           seconds-of-day of requested epoch
c        iprn          PRN # of requested SV
c        gnss          GNSS system to get clock for (TAH 201211).
c     Output
c        svclk         SV clock offsets

      implicit none
             
      include '../includes/dimpar.h'

      logical started,time_found,end_block,header
      character*(*) gnss   ! GNSS system (GREC) from t-file.

      integer*4 iusvclk,iprnt,jd,jdc,nsat,itsat(maxsat),ioerr,len,rcpar
     .        , iyear,month,iday,ihr,min,nval,iprn,julday,i,iloop

      real*8 sod,sodc,svclk(maxsat),timdif,test_time,tol,sec

      character*1 ignss !  GNSS system read from clock file.

      character*80 prog_name,buff80
      character*256 message
     
c     tolerance for matching clock epochs set to 1 sec (= 1 ps delay = 0.3 mm if clocks .001 ppb)
      data started/.false./,header/.true./,tol/1.d0/
      save started
           

c  Get calling program name  report_stat
                  
      len = rcpar(0,prog_name)
      if( .not.started) then
c     read headers down to the first data record
        write(iprnt,'(//,a)') '----------------------------------------'
        write(iprnt,'(/,a,/)') '** CLOCK FILE HEADER INFORMATION **'
        read(iusvclk,'(a)',iostat=ioerr) buff80 
        if (ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .      ,'orbits/rdsvclk',' '
     .      ,'Error reading first record of clock file',ioerr)
        if( buff80(21:30).ne.'CLOCK DATA') call report_stat(
     .      'FATAL',prog_name,'orbits/rdsvclk',' '
     .      ,'Clock file not recognized as RINEX format',ioerr)
        started = .true. 
        write(iprnt,'(a)') buff80 
        do while( header )
          read(iusvclk,'(a)',iostat=ioerr) buff80 
          if( ioerr.eq.-1 ) then
             call report_stat('FATAL',prog_name,'orbits/rdsvclk'
     .          ,' ','No END OF HEADER in clock file',ioerr)
          elseif( ioerr.ne.0 ) then  
            call report_stat('FATAL',prog_name,'orbits/rdsvclk'
     .             ,' ','No END OF HEADER in clock file',ioerr)
          endif 
          write(iprnt,'(a)') buff80 
          if( buff80(61:73).eq.'END OF HEADER' ) header=.false.
        enddo
      
c     end of header read
      endif
                    
c     initialize the SV clock values -- all nines indicates value is missing
* MOD TAH 030214: Initialization done in sdtrit.f so that same value
*      is used when no clock file is given.
c      do i=1,nsat
c       svclk(i) = 999999.999999d-6
c     enddo

c     find and read the next block (epoch) of SV clock values 
      time_found = .false.         
      end_block = .false.      
      iloop = 1
      do while (.not.end_block )
        read(iusvclk,'(a)',iostat=ioerr,end=99) buff80 
        if( buff80(1:2).eq.'AS' ) then
*         MOD TAH 201211: Added ignss to call to make sure it matches
*         gnss system being read from t-file. (rx->3x,a1)
          read(buff80,'(3x,a1,i2,2x,i4,4i3,f10.6,i3)',iostat=ioerr) 
     .         ignss, iprn,iyear,month,iday,ihr,min,sec,nval 
          if( ignss.eq.' ' ) ignss = 'G'   ! Default to GPS.   
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .            ,'orbits/rdsvclk',' ','Error clock file record',ioerr)
* MOD TAH 030318: Allow upto two values in the clock since the sigma
*         may also be given (as it is in the mit clock files).
          if( nval.gt.2 ) call report_stat('FATAL',prog_name
     .       ,'orbits/rdsvclk',' ','Number of clock values not 1',ioerr)
          jdc =julday( month,iday,iyear )                                             
          sodc = ihr*3600.d0 + min*60.d0 + sec
          test_time = timdif(jdc,sodc,jd,sod)           
          if( sod.gt.86340.d0 ) then
             print *,'jd sod jdc sodc test_time '
     .           ,jd,sod,jdc,sodc,test_time
          endif
* MOD TAH 201211: Also test if correct GNSS system.
          if( dabs(test_time).lt.tol .and. ignss.eq.gnss ) then   
c           times match, save the value    
            time_found = .true.
            do i=1,nsat 
              if(iprn.eq.itsat(i)) then
                read(buff80,'(40x,2e19.12)',iostat=ioerr) svclk(i)
                if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .          ,'orbits/rdsvclk',' ','Error reading clock value',ioerr)
              endif
            enddo                     
          elseif (test_time.gt.0.d0 ) then
c           clock file time is late, stop reading and return
            end_block = .true.
c           backspace to avoid missing a value on the next read
            backspace(iusvclk)
          else
c           time neither matches nor late, must be early, keep reading
          endif
        else
          if( time_found ) end_block = .true.
        endif   
        iloop = iloop + 1
        if( iloop.gt.50000 ) then  
           write(message,'(a,2(i8,f7.0))') 'Possible infinite loop at'
     .          ,jd,sod,jdc,sodc
           call report_stat('FATAL',prog_name,'orbits/rdsvclk',' '
     .          ,message,0)
        endif 
      enddo

   99 return
      end   

