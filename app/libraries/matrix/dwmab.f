CTITLE 'DWMAB'
 
      subroutine dwmab(scalar, v1, inc1, iter)

      implicit none
c
c
*     Returns index of largest absolute value in v1
c
*           inc1       - increment on V1
*   i          - loop counter
*   i1         - position in V1
*   iter       - number of values to be searched
*   scalar     - index of maximum value
      integer*4 inc1, i, i1, iter, scalar
*
c
C     real*8
c    .    v1(inc1,1)      ! vector to be checked
*   valmax     - maxvalue found
      real*8 v1(*), valmax
*
c
c
***** Find maximum value
* MOD TAH 980609: Initialized valmax to -1 to ensure that a
*     value is found.
      valmax = -1
      i1 = 1
      do i = 1, iter
         if( abs(v1(i1)).gt.valmax ) then
             valmax = abs(v1(i1))
             scalar = i
         end if
* MOD TAH 130420: Moved increment out of test.  Bug noted by
*        Qingping WANG, Fuzhou,Fujian,CHINA.
         i1 = i1 + inc1 
c        if( abs(v1(1,i)).gt.valmax ) then
c            valmax = abs(v1(1,i))
c            scalar = i
c        end if
      end do
c
***** Thats all
      return
      END
 
