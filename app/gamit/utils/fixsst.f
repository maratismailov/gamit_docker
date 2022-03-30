      Program FIXSST
                                                                           
c     This is a quick-and-dirty utility to correct a RINEX files for Trimble SST 
c     receivers when it was incorrectly translated as an SSE.  The program removes
c     P2 from the observable list, sets the L2 wavelength factor to 1, sets 
c     the # obs = 3. and changes the firmware version blank (read as 0.0 by MAKEX).   
c     It can be called for multiple files using sh_fixsst in /com.  

c     **Warning:  The receiver type is echoed to the screen, but there is not 
c       trap for non-SST receivers.

c     R. King 5 May 2001 
c       modified 17 June 2002 to blank out the firmware version
      
      implicit none
                   
      logical eof

      integer*4 urinex,urinex2,ioerr,linecnt,ifirmware    
                  
      character*12 frinex
      character*16 frinex2
      character*80 line             
      character*20 buff20
          
      write (*,*) 'RINEX file name?'
      read  (*,'(a)') frinex

      urinex = 10
      open (unit = urinex,
     .      file = frinex,
     .      iostat = ioerr,
     .      status = 'old')  
                         
      frinex2 = frinex // '.new'
      urinex2 = 11         
      open (unit = urinex2,
     .      file = frinex2,
     .      iostat = ioerr,
     .      status = 'new')  

      write(*,*) ' '
      write(*,*) '*************************'
      write(*,*) ' ' 
      write(*,*) 'Opening: ',frinex
 
      linecnt = 0
      eof = .false.

      do while (.not.eof)
 
        linecnt = linecnt + 1

        read(urinex,'(a)',end=99,iostat=ioerr ) line  
        
        if(ioerr.eq.-1 ) then  
          write(*,*) 'EOF at line ',linecnt
          write (*,*) 'LINE:',line
          eof = .true.
          goto 99
        elseif( ioerr.ne.0 ) then
          write(*,*) 'Error reading input RINEX file ',ioerr
        endif

        buff20= line(61:80)

        if ( buff20(1:19) .eq. 'REC # / TYPE / VERS' )  then
           read(line(41:41),'(i1)',iostat=ioerr) ifirmware 
           if (ioerr.eq.0.and.ifirmware.gt.4 ) then
              line(41:60) = '                   '  
              write(*,*) 'Replacing firmware version on line: ',line 
           endif
        endif                        
        if(  buff20 .eq. 'WAVELENGTH FACT L1/2' ) then
            line(12:12) = '2' 
            write(*,*) line
        endif

        if( buff20(1:19) .eq. '# / TYPES OF OBSERV' ) then
             line(6:6) = '3'
             line(29:30) = '  '
             write(*,*) line
        endif

        write( urinex2,'(a)') line

      enddo

   99 write(*,*) 'Wrote ',linecnt,' lines to ',frinex2
      write(*,*) ' '

      stop
      end
      

