

CTITLE 'dwint'
 
      subroutine dwint(scalar, v1, inc1,  iter)

      implicit none
c
*     Routine to initialize v1 with the value scale
c
*   inc1,inc2,inc3  - increment on V1, v2, and v3
*   i               - loop counter
*   i1,i2,i3        - position in V1, v2  and v3
*   iter            - number of values to be searched

      integer*4 inc1, i, i1, iter
*
c
c   v1(inc1,1),v2(inc2,1), v3(inc3,1)     ! vectors to be dotted 

      real*8 scalar, v1(*)
*
****  Loop adding   values
      i1 = 1
      do i = 1, iter
          v1(i1) = scalar
          i1 = i1 + inc1
      end do
c
***** Thats all
      return
      END
 
