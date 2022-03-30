c
      subroutine clnbuf(istart,iend,array)
c
c     initialize character array
c
c
      integer i,istart,iend
      character*(*) array
c
      do 20 i = istart,iend
         array(i:i) = ' '
 20   continue
c
      return
      end
