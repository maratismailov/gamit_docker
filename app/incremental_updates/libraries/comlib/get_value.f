CTITLE 'GET_VALUE'
 
      subroutine get_value( line, indx, type, field, last_field, value)

      implicit none
 
 
*     Routine to get next 'type' value from line with pointer in
*     position indx.  If last field is < field, then the next
*     set of values is read, otherwize the 'line' is decoded from
*     the beginning.
*     WARNING: When the routine is first called: Last_field should
*              greater than field to ensure that the pointers are
*              reset to the correct position in the string, i.e.,
*              unless 'indx' actually matches the correct column in
*              the line, then last should be set large to ensure the
*              allignemt takes place.
 
*         value         - Dummy declaration of value
 
      integer*2 value
 
*           indx        - Current pointer in line
*       field           - Column number to be read from line
*       last_field      - Last column read
*       i               - Loop counter
*       ierr            - IOSTAT error from readline
*       skip            - Number of columns to skip
 
      integer*4 indx, field, last_field, i, ierr,  skip

* MOD TAH 200623: Moved dum from I*4 to I*2 (gfortran 10.01 error during
*     compilation)
*       dum             - Dummy argument for read_line
      integer*2 dum
 
*             line      - line to be decoded
*       type            - type of value to be read
 
      character*(*) line, type
 
*           cdum        - Dummy for read_line
 
      character*1 cdum
 
****  See if we are just processing through the string
      if( last_field.lt. field ) then
          skip = field - last_field - 1
      else
          skip = field - 1
          indx = 1
      end if
 
****  Now skip over the columns
      do i = 1, skip
          call read_line(line,indx,'CH', ierr, dum, cdum)
      end do
 
****  Now get the value we need
      call read_line(line,indx, type, ierr, value, cdum)
 
***** Thats all
      return
      end
 
