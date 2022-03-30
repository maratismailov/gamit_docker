CTITLE GEOD_to_GEOD

      subroutine geod_to_geod(incrd, outcrd, intype, outtype,
     .                        indatum, outdatum, zone, hemi)

      implicit none 

*     Routine to make general coordinate system transformations 
*     from either XYZ, GEOD or UTM cordinates to another system 
*     and changed the datum at the same time.  The datum choices 
*     defined in kf/gen_util/datum_def.f.

*     Basic scheme is to covert input to XYZ, apply the datum shift and
*     then convert to outut system


      include '../includes/utmut.h'
      include '../includes/const_param.h'

* PASSED VARIABLES
      real*8 incrd(3)  ! Input coordinates, either XYZ, GEOD 
                       ! (Geodetic co-lat, long and height), or UTM
                       ! (Northing, Easting and Up)
      real*8 outcrd(3) ! The output system (same choices as above)

      integer*4 zone   ! Zone for UTM coordinates (passed in with UTM
                       ! coordinates, returned when converted to UTM)

      character*(*) intype, outtype ! Choices for the input and output
                       ! coordinate types (XYY, GEOD, UTM)
      character*(*) indatum, outdatum ! Choices for datum for in and
                       ! out systems (see datum_def.f.  Current choices
                       ! are WGS84, NAD27 and OMAN.
      character*(*) hemi  ! Hemisphere (N or S).  Passed in with UTM or
                       ! returned when converted to UTM.

* LOCAL VARIABLES
      real*8 inxyz(3), outxyz(3)  ! Input and output cartessian coords
      real*8 geod(3)   ! Geodetic co-lat, long and height
      integer*4 lenp, rcpar  ! Length of program name
      integer*4 j     ! Loop counter
      character*128 prog  ! Name of program calling routine


****  Based on indatum, convert the inout coordinates to XYZ
      if( intype(1:3).eq.'XYZ' ) then
          call datum_def(indatum)
          do j = 1,3
             inxyz(j) = incrd(j)
          enddo
      else if ( intype(1:4).eq.'GEOD' ) then
          call g2x(incrd, inxyz, indatum)
      else if ( intype(1:3).eq.'UTM' ) then
          call UTM_to_GEOD(incrd, geod, zone, hemi, indatum)
          call g2x(geod,inxyz, indatum)
      else
          lenp = rcpar(0,prog)
          call report_stat('WARNING',prog,'geod_to_geod',intype,
     .          'Unknown coordinate intype',0)
          return
      endif

****  OK: Now apply any transformation to bring system to ITRF origin
      do j = 1,3 
         outxyz(j) = inxyz(j) - dXYZ_ut(j)
      enddo

****  Now form the output results.  Get the new datum and apply
*     translation. 
      call datum_def(outdatum)

***** Apply any translation needed for the new datum
      do j = 1,3 
         outxyz(j) = outxyz(j) + dXYZ_ut(j)
      enddo

****  Now convert based on outsystem
      if( outtype(1:3).eq.'XYZ' ) then
         do j = 1,3
            outcrd(j) = outxyz(j)
         end do
      else if ( outtype(1:4).eq.'GEOD' ) then
          call x2g(outxyz,outcrd, outdatum)
      else if ( outtype(1:3).eq.'UTM' ) then
          call x2g(outxyz,outcrd, outdatum)
          zone = 0
          call GEOD_to_UTM(outcrd,outcrd, zone, hemi, outdatum)
      else
          lenp = rcpar(0,prog)
          call report_stat('WARNING',prog,'geod_to_geod',outtype,
     .          'Unknown coordoutate out type',0)
          return
      endif

****  Thats all
      return
      end

CTITLE G2X

      subroutine G2X(geod_pos,xyz_pos, indatum)
 
*     Routine based on GEOD_to_XYZ but with datum information passed.
 
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
      include '../includes/utmut.h'
      include '../includes/const_param.h'
 
*   geod_pos(3) - the geodetic coordinates.  The order of the
*               - values is:
*               - (1) - Geodetic co-latitudes (rads)
*               - (2) - Geodetic longitude (positive east) (rad)
*               - (3) - ellipsoidal height (m)
*   rad_curve    - radius of curvature at the site (N above)
*   xyz_pos(3) - XYZ coordinates of the site (m)
*   eccsq       - Eccentricity
 
      real*8  geod_pos(3),  rad_curve, xyz_pos(3), eccsq

      character*(*) indatum  ! Datum definition.

***** Start, get the latitudes and height by iteration.
      call datum_def(indatum)
 
      eccsq   = 2.d0*F_ut - F_ut**2
 
      rad_curve = ER_ut /
     .            sqrt(1.d0 - eccsq*cos(geod_pos(1))**2 )

      xyz_pos(1) = (rad_curve+geod_pos(3))*sin(geod_pos(1)) *
     .               cos(geod_pos(2))
      xyz_pos(2) = (rad_curve+geod_pos(3))*sin(geod_pos(1)) *
     .               sin(geod_pos(2))
      xyz_pos(3) = ((1.d0-eccsq)*rad_curve+geod_pos(3))*
     .               cos(geod_pos(1))

 
***** Thats all
      return
      end

CTTTLE X2G

      subroutine x2g(xyz_pos,geod_pos, datum)

*     Routine based on XYX_to_GEOD but with datum passed.  This
*     routine does not return rotation matrix.
 
*     Routine to compute the geodetic latitudes, longitude and
*     elipsoidal height from the cartesian XYZ coordinates of
*     the site.  The algorithm used is adopted from:
*     Heiskanen W A and H. Moritz, Physical Geodesy, W.H.Freeman
*     and company, San Francisco, 364, 1967. Chapter 5.3.
*
*     The routine also returns the tranformation between XYZ and
*     local NEU system.
*
*     RESTRICTION: User must supply the XYZ coordinates.  If these
*     are not available the can be computed from geod_pos using:
*
*     X = (N+h) cos(phi) cos(lambda)
*     Y = (N+h) cos(phi) sin(lambda)
*     Z = [(1-e**2)N+h]  sin(phi)
*
*     where  e**2 = 2f-f**2; f is flattening
*            N**2 = a**/[1-(e*sin(phi))**2]
*     N is East-West Radius of curvature.
*
      include '../includes/utmut.h'
      include '../includes/const_param.h'
 
 
*   eccsq       - eccentricity squared computed from flattening
*   equ_rad     - radial distance of site from rotation axis.
*   geod_pos(3) - the geodetic coordinates.  The order of the
*               - values is:
*               - (1) - Geodetic co-latitudes (rads)
*               - (2) - Geodetic longitude (positive east) (rad)
*               - (3) - ellipsoidal height (m)
*   lat_i       - approximate latitudes used in iteration
*   lat_p       - approximate latitudes from previous iteration
*   long        - longitude of the site.
*   h_i         - approximate height used in iteration
*   h_p         - approximate height from previous iteration
*   rad_lat     - radius to be used in computing the geodetic
*               - latitudes
*   rad_curve    - radius of curvature at the site (N above)
*   xyz_pos(3) - XYZ coordinates of the site (m)
 
*   tolerance   - toleranace for convergence on geodetic
*               - coordinates.
 
      real*8 eccsq, equ_rad,
     .    geod_pos(3), lat_i, lat_p, long, h_i, h_p, rad_lat,
     .    rad_curve, xyz_pos(3), tolerance

*   niter       - Number of iterations.  EXITS if more than 50 
*                 iterations are required.  (Removes problem if
*                 coordinates are NaN in IEEE floating point)

      integer*4 niter
 
*   converged   - Indicate values have converged.
 
      logical converged

      character*(*) datum
 
      data
*                                     ! Converge to 0.1 mm
     .    tolerance  / 0.0001d0 /
 
***** Start, get the latitudes and height by iteration.

      call datum_def( datum )
 
      equ_rad = sqrt( xyz_pos(1)**2 + xyz_pos(2)**2 )
      eccsq   = 2.d0*F_ut - F_ut**2
 
*                                            ! Set previous iteration values
      lat_p = atan2( xyz_pos(3), equ_rad)
      h_p   = 0.d0
 
      converged = .false.
      niter = 0
 
      do while ( .not. converged )
 
*         Get radius of curvature using previous latitudes estimate
          rad_curve = ER_ut /
     .               sqrt(1.d0 - eccsq*sin(lat_p)**2 )
          rad_lat  =  equ_rad *
     .               ( 1.d0 - eccsq*rad_curve/(rad_curve+h_p) )
 
          lat_i = atan2( xyz_pos(3), rad_lat)
 
*                                          ! Use cos lat formula
          if( abs(lat_i).lt. pi/4 ) then
               h_i   = equ_rad/cos(lat_i) - rad_curve
*                                           ! Use sin lat formula
           else
               h_i   = xyz_pos(3)/sin(lat_i) - (1.d0-eccsq)*rad_curve
           end if
 
 
*         Check for convergence
          if( abs(h_i-h_p)              .lt. tolerance .and.
*                                                             ! Converged
     .        abs(lat_i-lat_p)*rad_curve.lt. tolerance ) then
 
              converged = .true.
          end if

*         Check for two many iterations
          niter = niter + 1
          if( niter.gt.50 ) then
              write(*,'('' XYZ_to_GEOD ERROR: Failure to converge'')')
              converged = .true.
          end if
 
*         Save the latest values
          h_p   = h_i
          lat_p = lat_i
 
*                     ! iterating for latitudes and height
      end do
 
 
***** Save the final values
      long = atan2( xyz_pos(2),xyz_pos(1) )
 
*                                     ! colatitudes
      geod_pos(1) = pi/2.d0 - lat_i
*                                     ! Add 2*pi
      if( long.lt.0 ) then
          geod_pos(2) = 2*pi + long
      else
          geod_pos(2) = long
      end if
 
      geod_pos(3) = h_i
  
***** Thats all
      return
      end
 

 


      
