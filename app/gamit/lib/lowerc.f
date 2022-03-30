      character*(*) function lowerc(string)
c
c  Purpose:
c     Return a lower case string the same length as the minimum
c     of the lengths of string and lowerc declared in calling
c     routine.
c
c     Kurt Feigl April 3, 1989
c     Mark Murray June 19, 1989  Generalized
c
c  Input:
c     string   -  character*(*)
c                 String to be converted to lower case
c                 (unchanged on exit).
c  Output:
c     lowerc   -  character*(*)
c                 Converted string (length minimum of lengths of
c                 string and lowerc declared in calling routine)
c
c  N.B. This version of the routine is only correct for ASCII code.
c       Installers must modify the routine for other character-codes.
c       For EBCDIC systems the constant IOFF must be changed to +64.
c
c  Functions and Subroutines:
c     FORTRAN  ichar,min,len,lle,lge,char
c

      integer i,ilen,ioff
      character*(*) string
      parameter (ioff = 32)

      ilen = min(len(string),len(lowerc))
      do 10 i=1,ilen
         if (lge(string(i:i),'A') .and. lle(string(i:i),'Z')) then
            lowerc(i:i)=char(ichar(string(i:i))+ioff)
         else
            lowerc(i:i) = string(i:i)
         endif
  10  continue

c     null the rest of lowerc
      do 20 i=ilen+1,len(lowerc)
         lowerc(i:i) = char(0)
  20  continue

      return
      end

