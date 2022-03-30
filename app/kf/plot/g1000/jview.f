CTITLE JVIEW
 
      subroutine jview ( vl, vr, vb, vt )
 
*     Routine to set the view.  Uses the SET command with the current
*     world coordinates scales. (if world coordinates are not yet set
*     then it uses dummy values)
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   vl, vr, vb, vt  - Virtual coordinates of left, right,
*                   - bottom and top of virtual area
 
 
      real*4 vl, vr, vb, vt
 
* LOCAL VARIABLES
 
 
*     Save the view port in the common
      gvl = vl
      gvr = vr
      gvb = vb
      gvt = vt
 
*     Set the values
      call set(gvl,gvr,gvb,gvt,  gwl, gwr, gwb, gwt,1 )
 
****  Thats all
      return
      end
 
