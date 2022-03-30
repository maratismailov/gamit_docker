CTITLE CHECK_APR
 
      subroutine check_apr( old, new, check, type )

      implicit none 
 
 
*     Routine to check of the apriori values match what is expected
*     CHECK is one if the values are to be checked.  (Only values
*     which were fixed as apriori in the SOLVK solution should
*     be checked)
 
*   check   - If 1 then we should check
*   type    - Type of parameter.  Passed so that orbits are not
*             checked since these are expected to change.
 
      integer*4 check, type
 
*   old     - Old value of apriori (if zero then we assume we
*           - we have no value yet)
*   new     - New value for apriori
 
      real*8 old, new
 
 
***** If we should then check values

*     If the apriori has already been saved then copy it back
*     apr.  This way we will keep the first apriori value.  This
*     will stop the messgade below ever being printed.
      if( old.ne.0.d0 ) new = old
 
      if( check.eq.1 .and. type.ne.51 ) then
 
*         Now see if we have a value
*                                     ! we have an old value
          if( old.ne.0.d0 ) then
 
*             Check new against old
*                                     ! Values do not match, Report
C             if( old.ne.new ) then
              if( abs((old-new)/old).gt.1.d-6 ) then
                  write(*,100) type, old, new
  100             format(' **** WARNING **** Apriori values',
     .                   ' for parameter type ',i4,' do not match ',/,
     .              19x, ' Old value is ',d20.12,' New is ',d20.12)
              end if
          end if
      end if
 
***** Thats all
      return
      end
 
