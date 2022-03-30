CTITLE CHECK_CFLAG
 
      subroutine check_cflag( flag, expected, jerr )

      implicit none 
 
*     This routine will check that the flag read from a cfile
*     record matches that expected for the read.
 
* INCLUDE FILES
*     None
* PASSED VARIABLES
 
*   flag        - Flag read from cfile
*   expected    - Expected flag
*   jerr        - Error return, which is set if there
*           - is an error.
 
      integer*4 flag, expected, jerr
 
* LOCAL VARIABLES
*     None
 
***** Check the flag
      if( flag.ne.expected ) then
 
*         Tell user we have a problem
          write(*,100) expected, flag
 100      format(' ** ERROR READING CFILE ** Record type ',i2,
     .        ' expected and ',i5,' found for type')
          jerr = -1000
      else
 
          jerr = 0
 
      end if
 
****  Thats all
      return
      end
 
