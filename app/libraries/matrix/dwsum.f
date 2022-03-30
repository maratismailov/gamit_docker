

CTITLE 'dwsum'
 
      subroutine dwsum(r, v1, inc1, iter)

      implicit none
c
c     routine to sum the elements in V1 and return the result in R
c
*   inc1            - increment on V1
*   i               - loop counter
*   i1              - position in V1
*   iter            - number of values to be searched

      integer*4 inc1, i, i1, iter
*
c
c   v1(inc1,1)                ! vectors to be dotted 
c   r   -  result of dot product

      real*8 v1(*), r    
*
c
****  Loop swapping values
      i1 = 1
      r = 0.d0
      do i = 1, iter
         r = r + v1(i1)
         i1 = i1 + inc1
      end do
c
***** Thats all
      return
      END
 
