      subroutine wwrite(s1)
c     This routine replaces any Fortran Write to standard output.

      character*(*) s1
c     call cwrite(s1,len(s1))

* MOD TAH 920302: Changed routine to simply write the string to the
*     screen since curses is not being used anymore
      write(*,'(a)') s1

      return
      end

