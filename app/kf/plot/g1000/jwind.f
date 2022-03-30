CTITLE JWIND
 
      subroutine jwind ( wl, wr, wb, wt )
 
*     Routine to set the window.  Uses the SET command with the current
*     world coordinates scales. (if world coordinates are not yet set
*     then it uses dummy values)
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   wl, wr, wb, wt  - worldl coordinates of left, right,
*                   - bottom and top of user area
 
 
      real*4 wl, wr, wb, wt
 
* LOCAL VARIABLES
 
 
*     Save the world boundaries in the common
      gwl = wl
      gwr = wr
      gwb = wb
      gwt = wt
 
*     Set the values
      call set(gvl,gvr,gvb,gvt,  gwl, gwr, gwb, gwt, 1 )
 
****  Thats all
      return
      end
 
