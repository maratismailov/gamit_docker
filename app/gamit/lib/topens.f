      subroutine topens (tfname,stat,lunit,ioerr)
c
c     Open a T-file named tfname on logical unit lunit

      character*(*) tfname,stat
      integer lunit
      integer*4 ioerr,len,rcpar,ite(3)

      character*80 prog_name,head
      character*256 message

      logical old

c     get calling program for report_stat
      len = rcpar(0,prog_name)
 

c*    call lowers(stat) causes problem on DEC *
      if( stat(1:3).eq.'old' .or. stat(1:3).eq.'OLD' ) then
         old = .true.
      else
         old = .false.
      endif
      
      open (unit     = lunit,
     .      file     = tfname,
c    .      recl     = 8000,
     .      iostat   = ioerr,
     .      form     = 'unformatted',
     .      access   = 'sequential',
     .      status   = stat)
                   
      if (ioerr .ne. 0) 
     .   call report_stat('FATAL',prog_name,'lib/topens',tfname
     .       ,'Error opening T-file',ioerr)
            
      if( old ) then
c       read part of the first record to check binary compatibility  
c       gFortran writes extra bytes so the best test is the first
c       part of the date (month)
        read(lunit,iostat=ioerr) head,ite  
        if (ioerr .ne. 0)  then
          call report_stat('FATAL',prog_name,'lib/topens',tfname,
     .   'Error reading 1st record to check binary compatibility',ioerr)
        else       
c         check year
          if( ite(1).lt.0.or. ite(1).gt.12 ) then
            write(message,'(a,a)')
     .       'Binary T-file not compatible with machine'
     .      ,' architecture.  Was T-file made on a different machine?'
            call report_stat('FATAL',prog_name,'lib/topens',tfname
     .                    ,message,ioerr)
          endif
          rewind(lunit)
        endif
      endif

      return
      end

