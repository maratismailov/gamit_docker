      character*1 function upper1(letter)
c
c  Purpose:
c     Return an upper case letter.
c
c     Kurt Feigl April 3, 1989
c     Mark Murray June 19, 1989  Generalized
c     Kurt Feigl Dec. 27, 1991   De-generalized!
c
c  Input:
c     letter   -  character*1
c                 (unchanged on exit).
c  Output:
c     upper1   -  character*(*)
c
c  N.B. This version of the routine is only correct for ASCII code.
c       Installers must modify the routine for other character-codes.
c       For EBCDIC systems the constant IOFF must be changed to -64.
c
c  Functions and Subroutines:
c     FORTRAN  ichar,min,len,lle,lge,char
c

      integer ioff
      character*1 letter
      parameter (ioff = -32)

      if (lge(letter,'a') .and. lle(letter,'z')) then
         upper1=char(ichar(letter)+ioff)
      else
         upper1 = letter
      endif

      return
      end
