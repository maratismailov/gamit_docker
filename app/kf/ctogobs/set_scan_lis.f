CTITLE
 
      subroutine set_scan_list(uep )

      implicit none
 
*     This routine will set up the scan epochs for check_continuity
*     to resolve where a skip in the one-ways has occurred.
*     This routine will use the double difference arrays to
*     to decide where the slip has occurred.
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   uep(2)     - Epoch at which there is an unflagged slip and its
*                companion for comparison.
 
      integer*4 uep(2)
 
* LOCAL VARIABLES - None

*     Copy the scan epoch overs
 
      num_scan = 2

*     Put the epoch with the break second so that it will be flagged
      scan_ep(2) = uep(1)
      scan_ep(1) = uep(2)
 
***** Thats all
      return
      end
 
 
 
