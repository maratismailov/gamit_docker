CTITLE DECIMALTOINT
 
 
      integer*4 function decimaltoint( string, ierr )

      implicit none
 
*     Routine to emulate the same function on HP1000.
*     This routine will get the integer value stored in string;
*     here we do this by simply reading string and rteurning
*     IOSTAT error
 
 
      character*(*) string
 
*         ierr      - IOSTAT error on read
 
      integer*4 ierr
 
*         ivalue    - Value read from string, copied to
*                   - decimaltoint if there is no error, otherwize
*                   - return set to zero.
 
 
      integer*4 ivalue
 
****  Simply read the value
      read(string,*, IOSTAT=ierr ) ivalue
 
      if( ierr.eq.0 ) then
          Decimaltoint = ivalue
      else
          Decimaltoint = 0
      end if
 
****  Thats all
      return
      end
