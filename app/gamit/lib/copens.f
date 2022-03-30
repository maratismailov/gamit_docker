      subroutine copens (cfname,stat,lunit,ioerr)

c     open a C-file named cfname on logical unit lunit
c     return ioerr error code

c     K. Feigl June 89
          
      include '../includes/dimpar.h'

      character*(*) cfname,stat
      character*80  prog_name
      integer lunit
      integer*4 ioerr,len,rcpar,nblen

c     get calling program for report_stat
      len = rcpar(0,prog_name)
              
      open (unit     = lunit,
     .      file     = cfname(1:nblen(cfname)),
     .      iostat   = ioerr,
     .      err      = 100,
     .      form     = 'unformatted',
     .      access   = 'sequential',
     .      status   = stat)

  100 continue
      if (ioerr .ne. 0) then
         call report_stat('WARNING',prog_name,'/lib/copens',cfname
     .          , 'Error opening C-file: ',0)
         call ferror (ioerr,6)
      endif

      return
      end


