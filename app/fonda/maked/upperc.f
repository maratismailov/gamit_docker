      character*(*) function upperc(string)
c
c  Purpose:
c     Return a upper case string the same length as the minimum
c     of the lengths of string and upperc declared in calling
c     routine.
c
c     Kurt Feigl April 3, 1989
c     Mark Murray June 19, 1989  Generalized       
c
c  Input:
c     string   -  character*(*)
c                 String to be converted to upper case
c                 (unchanged on exit).
c  Output:
c     upperc   -  character*(*)
c                 Converted string (length minimum of lengths of 
c                 string and upperc declared in calling routine)
c
c  N.B. This version of the routine is only correct for ASCII code.
c       Installers must modify the routine for other character-codes.
c       For EBCDIC systems the constant IOFF must be changed to -64.
c
c  Functions and Subroutines:
c     FORTRAN  ichar,min,len,lle,lge,char
c

      integer i,j,ilen,ioff
      character*(*) string 
      parameter (ioff = -32) 
      
      ilen = min(len(string),len(upperc))
      do 10 i=1,ilen
         if (lge(string(i:i),'a') .and. lle(string(i:i),'z')) then
            upperc(i:i)=char(ichar(string(i:i))+ioff)
         else
            upperc(i:i) = string(i:i)
         endif
  10  continue

c     null the rest of upperc
      do 20 i=ilen+1,len(upperc)
         upperc(i:i) = char(0)
  20  continue
                
      return
      end

