CTITLE GET_NAME_FROM_BUFFER
 
      subroutine get_name_from_buffer( buffer, pos, record, start,
     .                                 next_name)
 

      implicit none
 
*     Routine to get the next name from a buffer.  The name we are
*     looking for may not be in this buffer, and if this is the
*     case then a blank name is returned.  The buffer contains names
*     between start+(record-1)*8
 
*   act_pos     - actual position of item in the current buffer
*               - which contains only sixteen of the names
*   pos         - The name number we are looking for.
*   record      - Current record number in the list of names being
*               - read from disk.
*   start       - the first position minus one of the name type
*               - we looking for.
*   true_pos    - Postion of item we are looking for in a buffer
*               - which contained all of the names
 
 
      integer*4 act_pos, pos, record, start, true_pos
 
*   buffer(1)   - Buffer containing upto sixteen names we may be
*               - looking for.
*   next_name   - Name corresponding to pos extracted from the
*               - buffer.  If the name is not in this buffer, then
*               - a blank name is returned
 
      character*(*) buffer(1), next_name
 
***** Clear the name
      next_name = ' '
 
*     See if the name we are looking for is in this buffer
 
      true_pos = start + pos
*                                         ! appears after current buffer
      if( true_pos.gt.record*64 ) RETURN
      if( true_pos.le.(record-1)*64) RETURN
 
*     Get the actual position in buffer
      act_pos = mod( true_pos-1, 64 ) + 1
      next_name = buffer(act_pos)
 
****  Thats all
      return
      end
 
