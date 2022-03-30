CTITLE EPOSN
 
      subroutine eposn(dcb, ierr, num, mode )
 
      implicit none

*     Routine to emulate the positioning suboutines from the
*     HP1000.  This routine simply updates the record pointer
*     in the dcb buffer
 
* Variables
*   dcb(16)     - DCB buffer, element 3 contains the record
*               - number. Element 5 is checked to make sure this
*               - is a direct access file.
*   ierr        - Error return.  -11 if not direct access.
*             -                -12 not valid mode
*   num         - Number of records to move if mode is 0, and
*               - absolute record number if mode is 1.
*   mode        - Mode, 1 for absolute, 0 for relative.
 
 
      integer*4 dcb(16), ierr, num, mode
 
****  See if file is direct access.
      if( dcb(5).gt.2 ) then
          ierr  = -11
          return
      end if
 
      ierr = -12
 
      if( mode.eq.0 ) then
          dcb(3) = dcb(3) + num
          ierr = 0
      end if
      if( mode.eq.1 ) then
*         Set dcb to 1 less than recpord past, since it will be 
*         updated in READF or WRITF
          dcb(3) = num - 1
        ierr = 0
      end if
 
****  Thats all
      return
      end
 
