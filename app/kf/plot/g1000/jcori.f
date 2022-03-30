CTITLE JCORI
 
      subroutine jcori(ox,oy,oz, px,py,pz )
 
*     Subroutine to save the orientation of characters wanted by the
*     user.  As far as I can tell we are current restricted to 90 deg
*     changes in orientation.  We could get fancy and invoke the surface
*     routines but this seems to unnecessary.
 
*     The emulation here saves the orientation in degrees for later use
*     in the writing routines.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   ox,oy,oz        - orientation of baseline.  The only values
*               - used here. We do not allow inversion of
*               - chracters
*   px,py,pz        - Orientation of perpendicular.  Not used.
 
 
      real*4 ox,oy,oz, px,py,pz
 
* LOCAL VARIABLES
 
*   ang         - Orientation angle in rads.
 
 
      real*4 ang
 
*     Compute the angle we want as an atan2 function
 
      ang = atan2(oy, ox)
 
      jor_deg = nint(ang*57.29577951)
 
***** Thats all
      return
      end
 
