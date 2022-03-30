CTITLE IFBRK
 
 
      logical function ifbrk(err)
 
*     Routine to emulate the HP1000 If break function.  If any key
*     is pressed, then this routine will return true.  Otherwize
*     if is flase.
 
*   err         - Error return. NOT USED.
 
      integer*4 err
 
*     LOCAL VARIABLES
 
*   inkey       - Function to return the next character in buffer
*   key         - The key in the buffer
 
      integer*4 inkey, key
               
c     this dummy statement avoids a gfortran compile warning
      err = 0

      if( inkey(key).gt.0 ) then
          ifbrk = .true.
      else
          ifbrk = .false.
      end if
 
 
****  Thats all
      return
      end
