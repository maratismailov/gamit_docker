CTITLE 'CHECK_RESET'
 
 
      logical function check_reset ( line, jndx )
 

      implicit none
 
*     Routine to check to see if a RESET is first 5 characters of line
 
*         indx      - Pointer to current position in line
*         jndx      - Copy of the pointer which is returned
*                   - unchanged.
 
      integer*4 indx, jndx
 
*              line - line to be decoded.
 
      character*(*) line
 
*           first_word  - First word of line (after blanks)
 
      character*5 first_word
 
****  Get the first word of the line
      indx = jndx
      call GetWord( line, first_word, indx)
      call casefold( first_word )
 
***** Check to see if RESET
      if( first_word(1:5).eq.'RESET' ) then
          check_reset = .true.
      else
          check_reset = .false.
      end if
 
***** That all
      return
      end
 
 
