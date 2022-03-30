CTITLE SET_PRE_EDIT
 
      subroutine set_pre_edit( data_flag, cf_iprn, rcv, ep )

      implicit none
 
*     Routine to pre-edit data based on user specified site/sv
*     epoch ranges.
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   cf_iprn         - PRN number (from cfile)
 
      integer*2 cf_iprn
 
*   data_flag       - Data status flag, will have bit 25 set
*                   - if data is pre-edited.
*   rcv             - Cffile number (rcv number)
*   ep              - current epoch
 
      integer*4 data_flag, rcv, ep
 
* LOCAL VARIABLES
 
*   i           - Loop counter
 
      integer*4 i
 
****  Start, see if there are any pre-edits
      if( num_pre_edit.gt. 0 ) then
 
          do i = 1, num_pre_edit
 
*             See if site matches
              if( pre_edit(1,i).eq.rcv .or. pre_edit(1,i).eq.0 ) then
 
*                 See if PRN matches
                  if( pre_edit(2,i).eq. cf_iprn .or.
     .                pre_edit(2,i).eq. 0 ) then
 
*                     See if with in epoch range
                      if( ep.ge. pre_edit(3,i) .and.
     .                    ep.le. pre_edit(4,i)      ) then
 
*                         This point should be editted
                          call sbit(data_flag,25,1)
                      end if
                  end if
              end if
          end do
      end if
 
****  Thats all
      return
      end
