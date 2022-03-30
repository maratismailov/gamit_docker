 
      subroutine GEOD_to_XYZ(latd,lond,ht,site_pos )
 
 
*     Routine to convert the geodetic co-latitudes, longitude and
*     elipsoidal height to the cartesian XYZ coordinates of
*     the site.   The algorithm used is adopted from:
*     Heiskanen W A and H. Moritz, Physical Geodesy, W.H.Freeman
*     and company, San Francisco, 364, 1967. Chapter 5.3.
c     Subtoutine adapted from tah GEOD_TO_XYZ in kf/gen_util. rwk 080115
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
                             
*       latd : geodetic latitude degrees
*       lond : longitude degrees
*       ht   : ellipsoidal height in meters


*   geod_pos(3) - the geodetic coordinates.  The order of the
*               - values is:
*               - (1) - Geodetic co-latitudes (rads)
*               - (2) - Geodetic longitude (positive east) (rad)
*               - (3) - ellipsoidal height (m)
*   rad_curve    - radius of curvature at the site (N above)
*   site_pos(3) - XYZ coordinates of the site (m)
*   eccsq       - Eccentricity
 
      real*8  geod_pos(3),  rad_curve, site_pos(3), eccsq  
     .     ,  pi,earth_flat,earth_rad,latd,lond,ht

c  routine expects co-latitude, longitude in radians
                            
      pi = 4.d0*datan(1.d0)
      geod_pos(1) = (90.d0 - latd)*pi/180.d0
      geod_pos(2) = lond*pi/180.d0
      geod_pos(3) = ht                  

***** Start, get the latitudes and height by iteration.
         
      earth_flat =  1.d0/298.257223563d0
      eccsq   = 2.d0*earth_flat - earth_flat**2
      earth_rad =  6378137.D0
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
 
