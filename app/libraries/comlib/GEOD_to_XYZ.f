 
      subroutine GEOD_to_XYZ(geod_pos, site_pos )

      implicit none
 
 
*     Routine to convert the geodetic co-latitudes, longitude and
*     elipsoidal height to the cartesian XYZ coordinates of
*     the site.  The algorithm used is adopted from:
*     Heiskanen W A and H. Moritz, Physical Geodesy, W.H.Freeman
*     and company, San Francisco, 364, 1967. Chapter 5.3.
*
*     The routine also returns the tranformation between XYZ and
*     local NEU system.
*
*     The following formula is used for the conversion:
*     NOTE: The colatiude is passed.

*     X = (N+h) cos(phi) cos(lambda)
*     Y = (N+h) cos(phi) sin(lambda)
*     Z = [(1-e**2)N+h]  sin(phi)
*
*     where  e**2 = 2f-f**2; f is flattening
*            N**2 = a**/[1-(e*sin(phi))**2]
*     N is East-West Radius of curvature.
*
      include '../includes/const_param.h'
 
*   geod_pos(3) - the geodetic coordinates.  The order of the
*               - values is:
*               - (1) - Geodetic co-latitudes (rads)
*               - (2) - Geodetic longitude (positive east) (rad)
*               - (3) - ellipsoidal height (m)
*   rad_curve    - radius of curvature at the site (N above)
*   site_pos(3) - XYZ coordinates of the site (m)
*   eccsq       - Eccentricity
 
      real*8  geod_pos(3),  rad_curve, site_pos(3), eccsq

***** Start, get the latitudes and height by iteration.
 
      eccsq   = 2.d0*earth_flat - earth_flat**2
 
      rad_curve = earth_rad /
     .            sqrt(1.d0 - eccsq*cos(geod_pos(1))**2 )

      site_pos(1) = (rad_curve+geod_pos(3))*sin(geod_pos(1)) *
     .               cos(geod_pos(2))
      site_pos(2) = (rad_curve+geod_pos(3))*sin(geod_pos(1)) *
     .               sin(geod_pos(2))
      site_pos(3) = ((1.d0-eccsq)*rad_curve+geod_pos(3))*
     .               cos(geod_pos(1))

 
***** Thats all
      return
      end
 
