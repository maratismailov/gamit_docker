CTITLE 'MOVECHARS'
 
      subroutine MoveChars(from,to,num_strings)

      implicit none
 
*     J.L. Davis 870413
 
*     Move character strings from arrays
 
*         i                         - Loop counter
*   ,   num_strings                 - The number of strings to transfer
 
      integer*4 i, num_strings
 
*             from(1)               - The array to transfer FROM
*   ,   to(1)                       - The array to transfer to
 
      character*(*) from(1), to(1)
 
***** Transfer
      do i = 1, num_strings
 
          to(i) = from(i)
 
      end do
 
      end
 
