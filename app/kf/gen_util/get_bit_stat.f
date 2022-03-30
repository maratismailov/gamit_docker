CTITLE GET_BIT_STATE
 
      subroutine get_bit_state( next_cont, bit_state )

      implicit none 
 
 
*     Routine to check the first character in next_cont and set
*     bit_state accordingly.  If the first character is - then
*     bit_state is set to zero, if + or any other character then
*     bit_state is set to one.  The + or - sign is stripped from
*     next_cont if either value is present
 
*   bit_state   - Bit state to be set depending on + or -
*               - at start of string
*   i           - Loop counter
*   string_len  - length of the next_cont string
 
      integer*4 bit_state, i, string_len
 
 
      character*(*) next_cont
 
***** Initialize bit_state, and then see if first character
*     is minus.  Change bit state if it is
 
      bit_state = 1
      if( next_cont(1:1).eq.'-' ) then
          bit_state = 0
      end if
 
*     Strip + or - from string
      if( next_cont(1:1).eq.'-' .or. next_cont(1:1).eq.'+' ) then
 
          string_len = len(next_cont)
          do i = 2, string_len
              next_cont(i-1:i-1) = next_cont(i:i)
          end do
      end if
 
***** Thats all
      return
      end
 
