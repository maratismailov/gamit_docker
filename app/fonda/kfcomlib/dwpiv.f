CTITLE 'DWPIV'
 
      subroutine dwpiv(scalar, v1,inc1, v2,inc2, v3,inc3, iter)
c
c
*     Does v3 = scalar*v1 + v2  -- Pivot operationa
c
*         i                     - Loop counter
*   inc1,inc2,inc3        - increments fpor each vector
*   i1,i2,i3              - pointers for each vector
*   iter                  - number of iterations
 
      integer*4 i, inc1,inc2,inc3, i1,i2,i3, iter
*
C     real*8
c    .    v1(inc1,1), v2(inc2,1), v3(inc3,1)   ! Vectors to be operated on
*   scalar                - Scaler multiplier
      real*8 v1(*), v2(*), v3(*), scalar
*
c
c
***** Set up pointers and start
      i1 = 1
      i2 = 1
      i3 = 1
      do i = 1 ,iter
         v3(i3) = scalar*v1(i1) + v2(i2)
c        v3(1,i) = scalar*v1(1,i) + v2(1,i)
         i1 = i1 + inc1
         i2 = i2 + inc2
         i3 = i3 + inc3
      end do
 
****  Thats all
      return
      end
