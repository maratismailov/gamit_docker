CTITLE JIWS
 
      subroutine jiws(device, option, element, num, ilist, rlist)
 
*     Emulation of inquire about work station.
*
* PASSED VARIABLES
 
*   device  - Workstation number (not used)
*   option  - Not used
*   element - Not used
*   num     - Not used
*   ilist(5)    - 5 integer values return (not used)
 
 
      integer*4 device, option, element, num, ilist(5)
 
*   rlist(2)    - Returns 1.0 and aspect ratio fot yaxis.
 
 
 
      real*4 rlist(2)
 
* LOCAL VARIABLES
 
*     None
 
      rlist(1) = 1.0
      rlist(2) = 1.0
 
      return
      end
 
