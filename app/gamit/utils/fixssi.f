      Program FIXSSI
                                                                           
c     This is a quick-and-dirty utility to correct a RINEX files for Trimble SSI 
c     receivers when the header incorrectly listed P1 as an observable.
c
c     R. King 23 October 2007 
      
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


        if( buff20(1:19) .eq. '# / TYPES OF OBSERV' ) then      
          if( line(6:36).eq.'5    L1    L2    C1    P2    P1')  
     .        line(6:36) =  '4    L1    L2    C1    P2      '
          if( line(6:48).eq.     
     .          '7    C1    L1    L2    P2    P1    S1    S2')  
     .        line(6:48) = '6    C1    L1    L2    P2    S1    S2      '
          write(*,*) line
        endif

        write( urinex2,'(a)') line

      enddo

   99 write(*,*) 'Wrote ',linecnt,' lines to ',frinex2
      write(*,*) ' '

      stop
      end
      

