CTITLE CHECK_CSIZE
 
      subroutine check_csize( dim, max, type, jerr )

      implicit none 
 
*     This routine will check that the dimension of arrays read
*     from a cfile and make sure that we have exceeded dimensions.
*     If we have and error will be returned.
 
* INCLUDE FILES
*     None
* PASSED VARIABLES
 
*   dim     - dimension read from cfile
*   max     - Maximum allowed in the users program
*   jerr        - Error return, which is set if there
*           - is an error.
 
      integer*4 dim, max, jerr
 
*   type        - Type of quanitity being checked.
 
      character*(*) type
 
* LOCAL VARIABLES
*     None
 
***** Check the size the array
      if( dim.gt.max ) then
 
*         Tell user we have a problem
          write(*,100) type, dim, max
 100      format(' ** ERROR READING CFILE ** Array ',a,' sized',
     .        ' with ',i4,' elements.  Max allowed ',i4)
          if( jerr.eq.0 ) jerr = -1001
      end if
 
****  Thats all
      return
      end
 
