CTITLE WMOVE

      subroutine wmove(from,to,count)

      implicit none

*     Routine to move count i*4 words from 'from' to 'to'

      integer*4 from(1),to(1)
      integer*4 count
 
c
      integer*4 i
c
      do i=1,count
        to(i)=from(i)
      enddo
c
      return
      end
