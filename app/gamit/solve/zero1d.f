c
      subroutine zero1d(istart,iend,array)

c initialize 1-d array

      implicit none

      integer istart,iend,i

      real*8 array(iend)

      do  i = istart,iend
         array(i) = 0.0d0
      enddo

      return
      end
