ctitle LASTEX

      integer*4 function lastex( string, sub )

      implicit none 

*     Routine to return the position of the last occurenence
*     of 'sub' in 'string'.  This is like a reverse order of
*     index.

      character*(*) string, sub

* LOCAL VARIABLES

*     i   - Loop counter
*     lenstring - Length of string
*     lensub    - Length of sub
*     trimlen   - returns length of used portion of string

      integer*4 i, lenstring, lensub, trimlen

****  Get the length of string
      lenstring = trimlen(string)
      lensub    = trimlen(sub)

      i = lenstring - lensub + 1
      lastex = 0
      do while ( i.gt.1 ) 
         if( string(i:i+lensub-1).eq.sub(1:lensub) ) then
             lastex = i
             i = 0
         end if
         i = i - 1
      end do

***** Thats all
      return
      end
