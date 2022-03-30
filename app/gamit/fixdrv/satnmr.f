      Subroutine satnmr( lgfil, ngsat,igsat )

c     Determine the number and PRN #'s of SVs on the G-file
c     R. King 920911 (modification of S. Shimada SATNMR)

      include '../includes/dimpar.h'

      character*80 line   
      character*256 message

      integer*4 lgfil,ngsat,igsat(maxsat),nparam,maxcom
     .        , ioerr,i,j
                           
      logical eof

c     maximum number of comments on the G-file
      parameter (maxcom=100)
c
c        Read the header lines

      read( lgfil,'(a80)')  line
      read( lgfil,'(i2)') nparam
      do i=1,maxcom
           read( lgfil, '(a80)' ) line
           if( line(1:3) .eq. 'END' ) goto 10
      enddo
      call report_stat('FATAL','FIXDRV','satnmr',' '
     .       ,'Something wrong with G-file header',0)

c        Loop over all the SVs

   10 ngsat = 0         
      eof = .false.
      do while (.not.eof)

         read( lgfil,'(a80)',iostat=ioerr ) line
         if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV'
     .     ,'satnmr',' ','Error reading g-file ',ioerr)
         if( line(1:3).eq.'END' ) then
           eof =.true.
         else
           if (line(1:3).ne.'PRN') then
              call report_stat('FATAL','FIXDRV','satnmr',' '
     .               ,'SV names in G-file must be PRN numbers',0)
           else                  
             ngsat = ngsat + 1
             if( ngsat.gt.maxsat ) then   
                write(message,'(a,i3,a,i3,a)') 
     .            'Number of SVs on G-file (',ngsat
     .             ,') exceeds MAXSAT (',maxsat,')'
                 call report_stat('FATAL','FIXDRV','satnmr',' '
     .            ,message,0)
             endif
             read( line(5:6),'(i2)') igsat(ngsat)
           endif
           do j=1,nparam
             read( lgfil,'(a80)') line
           enddo
         endif
      enddo

      end




