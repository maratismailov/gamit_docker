
CTITLE CLEARBUFFER

      subroutine clearbuffer( buff, dim )

      implicit none

*     routine to clear and I*4 dimensioned to dim .

      integer*4 dim, buff(dim)

*     Local variables
      integer*4 i

      do i = 1 , dim
         buff(i) = 0
      end do

****  Thats all
      return
      end
