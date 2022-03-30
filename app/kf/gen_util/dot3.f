CTITLE DOT3
 
      real*8 function dot3( a,b )

      implicit none 
 
*     routine for dotting two 3, element vectors together
 
* Variables
*     a  - vector 1 (three elements)
*     b  - vector 2 (three elements)
 
      real*8 a(3), b(3)
 
* Local Variables
*     i  - loop counter

      integer*4 i
 
****  Do the dot product
      dot3 = 0.d0
      do i = 1,3
         dot3 = dot3 + a(i)*b(i)
      end do
 
***** That's all
      return
      end
 
