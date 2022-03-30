CTITLE J2MRK
 
      subroutine j2mrk ( xp, yp, point_type )
 
*     Emulates the point of point type at position xp,yp (user
*     coordinates.  If point type is <0 then a polymarker
*     is output
 
* PASSED VARIABLES
 
*   xp, yp  - Position of point to be plotted
 
 
      real*4 xp, yp
 
*   point_type  - Type of point (positive value)
 
 
      integer*4 point_type
 
* LOCAL VARIABLES
 
*   mod_point -- There seems to be only 5 poly marks.  Therefore take
*                value mod 5
 
      integer*4 mod_point
 
*     Emulation uses the POINTS subroutines from SPPS.
      mod_point = point_type
      if( point_type.lt.0 ) mod_point = mod(point_type,5)
      call points(xp,yp, 1, mod_point , 0 )
*     Here               ^     ^         ^
*                        |     |         | Do not draw line
*                        |     | if negative of GKS Polymarkers
*                        | number of points
 
*     Thats all
      return
      end
 
