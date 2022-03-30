      subroutine mopens (mfname,stat,lunit,ioerr)
c
c purpose:
c     open a M-file named mfname on logical unit lunit
c
      implicit none

      character*(*) mfname,stat
      character*16 prog_name
      integer lunit
      integer*4 ioerr,len,rcpar

c     get calling module name for report_stat
      len =  rcpar(0,prog_name)

      open (unit     = lunit,
     .      file     = mfname,
     .      iostat   = ioerr,
     .      err      = 100,
     .      form     = 'unformatted',
     .      access   = 'sequential',
     .      status   = stat)

  100 continue
      if (ioerr .ne. 0) then
        call report_stat('WARNING',prog_name,'lib/mopens',' '
     .          , 'Error opening M-file',ioerr)
        call ferror (ioerr,6)
      endif

      return
      end

