CTITLE SUB_CHAR

      subroutine sub_char(string, old, new )

      implicit none

*     routine to replace all occurences of character OLD with
*     character NEW

*  string  - string to be manipilated
*  old     - old character to be replaced
*  new     - new character to be used

      character*(*) string
      character*1 old, new

* LOCAL VARIABLES

*  i    - Pointer to finding character

      integer*4 i

****  Loop over string
      i = 1
      do while ( i.gt.0 )
         i = index ( string, old )
         if( i.gt.0 ) then
             string(i:i) = new
         end if
      end do

****  Thats all
      return
      end


