CTITLE SET_DD_RET
 
      subroutine set_dd_ret(type, user_val, max_ret, max_max_ret)

      implicit none
 
*     Routine to set the values for maximum return values.  The
*     values are checked and the maximum used if the user value
*     is too large
 
* PASSED VARIABLES
 
*   user_val    - User value selected
*   max_ret     - The maximum value set to be used in the
*               - run
*   max_max_ret - Maximum value allowed.
 
      integer*4 user_val, max_ret, max_max_ret
 
*   type        - Type of maxixum being set
 
      character*(*) type
 
***** See if user value is too large, if not save values,
*     otherwise print an error message and use max allowed.
      if( user_val.gt.max_max_ret ) then
          write(*,100) type, user_val, max_max_ret
100       format(' ** WARNING ** User value for ',a,' (',I5,')',
     .            ' too large.  Using ',i5)
          max_ret = max_max_ret
      else
          if( user_val.lt.5 ) then
              write(*,110) type, user_val
110           format(' **WARNING** User value for ',a,' (',I5,')',
     .               ' too small.  At least 5 points required')
              max_ret = 5
          else 
              max_ret = user_val
          end if
      end if
 
****  Thats all
      return
      end
 
