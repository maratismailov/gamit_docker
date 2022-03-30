      Program FIXX
                                                                           
c     This is a quick-and-dirty utility to correct the X-files for TREX 1-18  (December 1986
c     - March 19908), which have binary characters in the comments section and hence are
c     not handled correctly by grep in sh_gamit.  TREX 11 (March 1988) also has missing
c     values of LAMBDA.  The LAMBDA code here is skipped unless lflag = .true.
c     This program can be called for multiple files using sh_fixx in /com.

c     R. King 10 July 2001
      
      implicit none
                   
      logical eof,eoh,lambda

      integer*4 uxfile,uxfile2,ioerr,linecnt    
                  
      character*10 fxfile
      character*14 fxfile2
      character*116 line       
  
      data lambda/.false./      
          
      write (*,*) '***EXECUTING NEW FIXX**'
      write (*,*) 'X-file name?'
      read  (*,'(a)') fxfile

      uxfile = 10
      open (unit = uxfile,
     .      file = fxfile,
     .      iostat = ioerr,
     .      status = 'old')      
      if( ioerr.ne.0 ) then
        write(*,*) 'Error opening ',fxfile 
        stop
      else
        write(*,*) ' ' 
        write(*,*) 'Opened: ',fxfile
      endif
      fxfile2 = fxfile // '.new'
      uxfile2 = 11         
      open (unit = uxfile2,
     .      file = fxfile2,
     .      iostat = ioerr,
     .      status = 'new')  

      if( ioerr.ne.0 ) then
        write(*,*) 'Error opening ',fxfile2 
        stop
      else
        write(*,*) ' '
        write(*,*) 'Opened: ',fxfile2
      endif 

      linecnt = 0
      eof = .false.
      eoh = .false.
                
c     these files were written by ctox, which has a prescribed first line:
c      X-File written from C-File :  CALGOA.076  
c    
c     but when these were created also put binary junk at the end that
c     causes a problem with grep (in sh_preproc): Linux will handle them only
c     with grep -a, but Solaris and HP-UX don't recognize the -a.  Our solution
c     for now is to remove the junk here.  If the C-file passed through multiple
c     edits, there may be more than one of these lines in the comment section
c     of the header, terminated by END in the first three columns
        
      do while (.not.eoh) 

        read(uxfile,'(a)',end=99,iostat=ioerr ) line   
        if( ioerr.ne.0 ) then
           write(*,*) ' Error reading comment line of x-file ',ioerr
        elseif( line(2:27).eq.'X-File written from C-File' ) then
             line(42:80) = '                                      '    
        elseif( line(2:9).eq.'Receiver' ) then 
             line(2:40) = 'Receiver        : TI4100              ' 
        elseif( line(17:30).eq.'antenna height' ) then 
             line(45:80) = '                                   '
        elseif ( line(1:3) .eq. 'END' ) then  
          write( uxfile2,'(a)') ' '
          if( lambda ) then   
            write( uxfile2,'(a)') 
     .         'Rewritten using fixx to correct missing LAMBDA  '
          else
             write( uxfile2,'(a)')
     .         'Rewritten using fixx to remove binary characters'  
          endif 
          write( uxfile2,'(a)') ' '   
          linecnt = linecnt + 3
          eoh = .true.  
        endif
        write( uxfile2,'(a)') line
      enddo 

c     now do the rest of the header and all the data lines
      do while (.not.eof)
        linecnt = linecnt + 1
        read(uxfile,'(a)',end=99,iostat=ioerr ) line  
        if(ioerr.eq.-1 ) then  
          write(*,*) 'EOF at line ',linecnt
          write (*,*) 'LINE:',line
          eof = .true.
          goto 99
        elseif( ioerr.ne.0 ) then
          write(*,*) 'Error reading input x-file ',ioerr
        endif
        if( lambda ) then
          if( line(5:14).eq.'SATELLITES' ) then
            line(28:52) = '4 DATA TYPES:  1  2  3  4'  
          endif
          if( line(2:8) .eq. 'CHANNEL' ) then
            line(35:52) = 'LAMBDA -1 -1  1  1'
          endif
          if( line(2:7) .eq. 'YR DAY' ) then
            line(39:60) = 'DATA INTERVAL  SESSION'   
            write( uxfile2,'(a)') line 
            linecnt = linecnt + 1
            read(uxfile,'(a)',end=99,iostat=ioerr ) line  
c           set the data interval equal to the x-file sampling interval and the session =1
            line(42:44) = line(25:27)       
            line(57:57) = '1'
          endif 
        endif  
        write( uxfile2,'(a)') line
      enddo

   99 write(*,*) 'Wrote ',linecnt,' lines to ',fxfile2
      write(*,*) ' '

      stop
      end
      

