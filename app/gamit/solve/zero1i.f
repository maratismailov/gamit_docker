c
      subroutine zero1i(istart,iend,irray)

c initialize 1-d integer array

      implicit none

      integer iend,irray,istart,i

      dimension irray(iend)

      do i = istart,iend
         irray(i) = 0
      enddo

      return
      end
