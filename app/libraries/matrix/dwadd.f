

CTITLE 'dwadd'
 
      subroutine dwadd(v1, inc1, v2,inc2, v3, inc3, iter)

      implicit none
c
*     Routine to add v1 to v2 and return result in v3
c
*   inc1,inc2,inc3  - increment on V1, v2, and v3
*   i               - loop counter
*   i1,i2,i3        - position in V1, v2  and v3
*   iter            - number of values to be searched

      integer*4 inc1,inc2, inc3, i, i1,i2,i3, iter
*
c
c   v1(inc1,1),v2(inc2,1), v3(inc3,1)     ! vectors to be dotted 

      real*8 v1(*), v2(*), v3(*)
*
****  Loop adding   values
      i1 = 1
      i2 = 1
      i3 = 1
      do i = 1, iter
         v3(i3) =  v1(i1) + v2(i2)
         i1 = i1 + inc1
         i2 = i2 + inc2
         i3 = i3 + inc3
      end do
c
***** Thats all
      return
      END
 
