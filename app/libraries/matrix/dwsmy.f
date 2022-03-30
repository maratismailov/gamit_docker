CTITLE 'DWSMY'
 
      subroutine dwsmy(scalar, v1, inc1, v2, inc2, iter)

      implicit none
c
c
*     Scalar multiply routine  v2 = scalar*v1
c
*          inc1, inc2      - Increments for vectors
*   iter            - number of elements
*   i1, i2          - Indices for vector one and two
*   i               - Loop counter
      integer*4 inc1, inc2, iter, i1, i2, i
*
c
C     real*8
c    .    v1(inc1,1), v2(inc2,1)    ! Input vectors
*   scalar          - multipling scalar
      real*8 v1(1), v2(1), scalar
*
c
c
****  Do the multipy
      i1 = 1
      i2 = 1
      do i = 1, iter
         v2(i2) = scalar*v1(i1)
c        v2(1,i) = scalar*v2(1,i)
         i1 = i1 + inc1
         i2 = i2 + inc2
      end do
c
***** Thats all
      return
      END
 
