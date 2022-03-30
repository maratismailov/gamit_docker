c
      subroutine zero1d(istart,iend,array)
c
c initialize 1-d array
c
c
      real*8 array(iend)
      integer istart,iend,i
c
      do 20 i = istart,iend
         array(i) = 0.0d0
 20   continue
c
      return
      end
