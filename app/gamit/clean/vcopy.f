      subroutine vcopy (iq1,a1,iq2,a2,n)

c     copy n-vectors from 1 to 2

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'

      integer n,i
      real*8 a1(maxepc)
      real*8 a2(maxepc)
      integer iq1(maxepc)
      integer iq2(maxepc)

c     These loops are separated to avoid bouncing between
c     different sections of memory

      do 10 i=1,n
         iq2(i) = iq1(i)
 10   continue

      do 20 i=1,n
         a2(i) = a1(i)
 20   continue

      return
      end

