      subroutine movewords(from,to,count)

      implicit none
      integer*2 from(1),to(1)
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
