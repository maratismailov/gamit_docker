CTITLE SAVE_TOOFEW
 
      subroutine save_toofew(num,  epochs, dd_edit_num, dd_edit_ep,
     .                    max_dd_edit)

      implicit none
 
*     Routine to see if we have some double or single differences
*     but not enough for patching (i.e, more than zero but less than
*     5).  If this is the case then the data are edited.  Gerenerally
*     this occurrs when we have nearby bias flags on different satelltes.
 
* PASSED VARIABLES
 
*   num     - Number of double or single differneces
*   epochs(num) - Epochs of the single or oduble differences
*   max_dd_edit - Maximum number of edits allowed.
*   dd_edit_num - number of edits in the list of edits
*   dd_edit_ep(max_dd_edit) - List of edited epochs
 
      integer*4 num, epochs(num), max_dd_edit, dd_edit_num,
     .    dd_edit_ep(max_dd_edit)
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
*     See if we need to edit.  If we do then added epochs to liist
      if( num.gt.0 .and. num.le.4 ) then
 
*         Need to edit, so add to list.
          if( num+dd_edit_num.gt.max_dd_edit) then
              write(*,*) ' **ERROR** Too many edits, ignoring these'
          else
              write(*,*) ' Editing ', num, ' points, epoch ', epochs(1)
              do i = 1, num
                  dd_edit_ep(i+dd_edit_num) = epochs(i)
              end do
              dd_edit_num = dd_edit_num + num
          end if
      end if
 
****  Thats all
      return
      end
 
 
 
 
