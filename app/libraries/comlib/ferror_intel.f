      subroutine ferror(ios,lu)

      implicit none

c     print the translation for error code

c     Kurt Feigl for sun

      integer ios,lu,len,rcpar
      character*4 buf4
      character*80 prog_name
      character*127 string
                                    
c  Peng Fang 940428
c  orig      call gerror(string)
c  replaced gerror with a unix generic call perror working on Apo/Sun/HP. rwk 960226
      call for_perror (string)
c If you are using the intel fortran compiler prior to version 8 use the perror call below      
c     call perror (string)
  
c     get calling program name for report_stat
      len = rcpar(0,prog_name)
      call report_stat('WARNING',prog_name,'lib/ferror',' ',string,ios)

c  rwk 960229: 'lu' is the output file, usually the screen or log file.
c              With report_error, no longer need this variable in ferror.
c              To avoid a compiler warning, though, do something innocuous with it:
c
c     old code:
c              write (lu,*) string
c     new code:
         write(buf4,'(i4)') lu

      return
      end

