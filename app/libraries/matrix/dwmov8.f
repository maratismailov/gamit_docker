
CTITLE 'dwmov8'
 
      subroutine dwmov8(v1, inc1, v2,inc2, iter)

      implicit none
c
c MOD TAH 190520: Added new version that allows iter to be I*8
c     to handle matrices bigger than 32767x32767.
c
*     Routine to move elements from v1 to v2.
c
*   inc1,inc2       - increment on V1, and v2
*   i               - loop counter
*   i1,i2           - position in V1 and v2
*   iter            - number of values to be searched

      integer*4 inc1,inc2
*
      integer*8 iter

      integer*8 i1,i2, i
c
      real*8 v1(*), v2(*)
*
c
c
****  Loop moving values
      i1 = 1
      i2 = 1
      do i = 1, iter
         v2(i2) = v1(i1)
         i1 = i1 + inc1
         i2 = i2 + inc2
      end do
c
***** Thats all
      return
      END
 
