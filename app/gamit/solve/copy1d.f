c
      subroutine copy1d(istart,iend,ishift,array1,array2)
c
c copy 1-d array from array1 to array2

      implicit none

      integer istart,iend,ishift,i

      real*8 array1(iend),array2(iend+ishift)
c
      do 20 i = istart,iend
         array2(i+ishift) = array1(i)
 20   continue
c
      return
      end
