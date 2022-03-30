CTITLE ESORT
 
      subroutine esort( num, ilist)

      implicit none 

*     This routine uses an exchande sort algormithm to sort
*     the list ilist into ascending order.  There are num values
*     in ilist.
 
*   num     - Number of values to be sorted
*   ilist(num)  - List to be sorted in to ascending order.
 
      integer*4 num, ilist(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - Smallest integer in current pass.
*   iswap   - Value used to swap integers
 
 
      integer*4 i,j, smallest_one, iswap
 
****  Start loop using exchange sort
 
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( ilist(j).lt. ilist(smallest_one) ) then
                  smallest_one = j
              end if
          end do 
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              iswap = ilist(smallest_one)
              ilist(smallest_one) = ilist(i)
              ilist(i) = iswap
          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end
 
