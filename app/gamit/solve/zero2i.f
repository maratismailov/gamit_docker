c
      subroutine zero2i(irow,icolmn,ndat1,ndat2,irray)
c
c     initialize 2-d integer array

      implicit none

      integer irow,icolmn,irray,ndat1,ndat2,i,j
      dimension irray(irow,icolmn)


      do  i = 1,ndat2
         do  j = 1,ndat1
            irray(j,i) = 0
         enddo
      enddo

      return
      end
