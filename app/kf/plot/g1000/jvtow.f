CTITLE JVTOW
 
      subroutine jvtow( vx, vy, wx, wy, wz )
 
*     Converts virtual coordinates to world coordinates.  Here we
*     set wz to 0.
 
* PASSED VARIABLES
 
*   vx, vy      - Virtual coorinates (we take as fractional)
*   wx, wy, wz  - World coordinates
 
 
      real*4 vx, vy, wx, wy, wz
 
* LOCAL VARIABLES
 
*   cfux, cfuy  - Functions to convert from fractional to world
 
 
      real*4 cfux, cfuy
 
****  Use the SPPS conversion routines
 
      wx = cfux(vx)
      wy = cfuy(vy)
      wz = 0.0
 
***** Thats all
      return
      end
 
