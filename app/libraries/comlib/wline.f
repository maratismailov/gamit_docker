ctitle WLINE
 
      subroutine wline(i)

      implicit none
 
*     This routine will write out [#] and then move the cursor
*     to be beginning of the line.
 
*   i       - Value to be written
 
      integer*4 i
 
*   cr      - Contains carrier return
 
      character*1 cr
 
      cr = char(13)
      write(*,100) i, cr
 100  format('[',i4.4,']',a,$)
 
****  Thats all
      return
      end
 
