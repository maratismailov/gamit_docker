c
      subroutine zero2d(irow,icolmn,ndat1,ndat2,array)

c     initialize 2-d array

      implicit none

      integer irow,icolmn,ndat1,ndat2,i,j

      real*8 array(irow,icolmn)

      do  i = 1,ndat2
         do  j = 1,ndat1
            array(j,i) = 0.0d0
         enddo
      enddo

      return
      end
