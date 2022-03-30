c
      subroutine add1(istart,iend,array1,array2)
c
c copy 1-d array from array1 to array2

      implicit none

      integer*4 istart,iend,i
      real*8 array1(iend),array2(iend)
c
      do 20 i = istart,iend
         array2(i) = array1(i) + array2(i)
 20   continue
c
      return
      end
