      Subroutine geoc_to_geod( clatd,rad,glatd,ht,semi,finv )
                                                                        
c     Subroutine to tranform geocentric lat,radius to geodetic lat, height.
c     The algorithm used is adopted from Heiskanen W A and H. Moritz,
c     Physical Geodesy, W.H.Freeman & Co., San Francisco, 364, 1967. Chapter 5.3.
c     Adapted from tah XYZ_TO_GEOD in kf/gen_util.  rwk 18 Feb 2000
                                                                        
      implicit none
                                                                        
c     Input
c       clatd  geocentric latitude in degrees
c       rad    geocentric radius in meters
c       semi   semi-major axis of ellipsoid in meters
c       finv   inverse of the flattening (dimensionless)
c     Output
c       glatd  geodetic latitude in degrees
c       ht     height above the ellipsoid in meters
                                                                        
      real*8 clatd,rad,semi,finv,glatd,ht
                                                                        
                                                                        
*   f           - flattening = 1 /finv
*   eccsq       - eccentricity squared computed from flattening
*   equ_rad     - radial distance of site from rotation axis.
*   zcoord      - cartesian z coordinate
*   lat_i       - approximate latitudes used in iteration
*   lat_p       - approximate latitudes from previous iteration
*   long        - longitude of the site.
*   h_i         - approximate height used in iteration
*   h_p         - approximate height from previous iteration
*   rad_lat     - radius to be used in computing the geodetic
*               - latitudes
*   rad_curve    - radius of curvature at the site (N above)
*   tolerance   - toleranace for convergence on geodetic
*               - coordinates.
                                                                        
      real*8 f, eccsq, equ_rad, zcoord, lat_i, lat_p, h_i, h_p, rad_lat
     .     , rad_curve, tolerance, clat, pi
                                                                        
*   niter       - Number of iterations.  EXITS if more than 50
*                 iterations are required.  (Removes problem if
*                 coordinates are NaN in IEEE floating point)
                                                                        
      integer*4 niter

*   calling program

      character*80 prog_name
                                                                        
*   converged   - Indicate values have converged.
                                                                        
      logical converged
                                                                        
      data
*                                     ! Converge to 0.1 mm
     .    tolerance  / 0.0001d0 /
          
	  
*     Get the name of the calling program

      call rcpar(0,prog_name)
                                                                    
***** Start, get the latitudes and height by iteration.
                                                                        
      pi = datan(1.d0) * 4.d0
      clat = clatd * pi / 180.d0
      zcoord = rad*dsin(clat)
      equ_rad = rad * dcos(clat)
      f = 1.d0/finv
      eccsq   = 2.d0*f - f**2
                                                                        
*                                            ! Set previous iteration values
      lat_p = atan2( zcoord, equ_rad)
      h_p   = 0.d0
      converged = .false.
      niter = 0
                                                                        
      do while ( .not. converged )
                                                                        
*         Get radius of curvature using previous latitudes estimate
          rad_curve = semi /
     .               sqrt(1.d0 - eccsq*sin(lat_p)**2 )
          rad_lat  = equ_rad *
     .               ( 1.d0 - eccsq*rad_curve/(rad_curve+h_p) )
                                                                        
          lat_i = atan2( zcoord, rad_lat)
                                                                        
*                                          ! Use cos lat formula
          if( abs(lat_i).lt. pi/4 ) then
               h_i   = equ_rad/cos(lat_i) - rad_curve
*                                           ! Use sin lat formula
           else
               h_i   = zcoord/sin(lat_i) - (1.d0-eccsq)*rad_curve
           end if
                                                                        
c      print *,'rad_curve rad_lat lat_i h_i ',rad_curve,rad_lat,lat_i,h_i
                                                                        
*         Check for convergence
          if( abs(h_i-h_p)              .lt. tolerance .and.
*                                                             ! Converged
     .        abs(lat_i-lat_p)*rad_curve.lt. tolerance ) then
                                                                        
              converged = .true.
          end if
                                                                        
*         Check for too many iterations
          niter = niter + 1
          if( niter.gt.50 ) then
c           print *,'lat_i lat_p h_i h_p ',lat_i,lat_p,h_i,h_p
            call report_stat('FATAL',prog_name,'utils/geoc_to_geod'
     .         ,' ','Geocentric to geodetci failure to converge',0)
          end if
                                                                        
*         Save the latest values
          h_p   = h_i
          lat_p = lat_i
                                                                        
*                     ! iterating for latitudes and height
      end do
                                                                        
                                                                        
***** Save the final values
                                                                        
      glatd = lat_i * 180.d0 / pi
      ht = h_i
                                                                        
      return
      end

