CTITLE LOC_TO_GEOD
 
      subroutine loc_to_geod( loc_coord, pos_geod )

      implicit none
 
*     Routine to convert local coordinates (co-latitude, longitude,
*     and ellipsiodal height) to geodectic position coordinates.
*     These coodinates are defined as the distance from the equator,
*     distance from greenwich meridian along the small circle at the
*     latitude of the site, and the height above the ellipsiod.
*
* MOD TAH 871016 Changed the computatoion of the small circle radius at
*     the latitude of the site to use the latitude truncated to the
*     nearest 10 asec.  This way we always use the same radius so that
*     the change in the E component represent changes in the longitude.
*     For high latitude site, a large amount of latitude change can get
*     into the E change.
 
      include '../includes/const_param.h'
 
*   colat           - Colatitude (used so that loc_coord and pos_geod
*                   - can be the same arrays)
*   colat_trun      - Coltitude truncated to the nearest 10 asecs
*                   - Used to compute small circle radius
*   loc_coord(3)    - Colatitude (rad), longitude (rad) and
*                   - height (m)
*   pos_geod(3)     - Distance from equator, distance from Greenwich
*                   - meridian, and height above ellipsiod (all m)
 
 
      real*8 colat, colat_trun, loc_coord(3), pos_geod(3)
 
****  Distance from equator
 
      colat = loc_coord(1)
*                                          ! 2.d4 = 1/10" in radians
      colat_trun = anint(colat*2.d4)/ 2.d4
 
      pos_geod(1) = (pi/2 - loc_coord(1))*earth_rad
 
****  Distance from Greenwich meridian along small circle
 
      pos_geod(2) = loc_coord(2) * earth_rad * sin(colat_trun)
 
***** Height above ellipssiod
 
      pos_geod(3) = loc_coord(3)
 
***** Thats all
      return
      end

CTITLE GEOD_TO_LOC
 
      subroutine geod_to_loc(pos_geod, loc_coord, xyz )

      implicit none

*     Routine to invert geodetic positions to local coordinates 
*     (co-latitude, longitude,and ellipsiodal height).  This is the
*     inverse of loc_to_geod

      include '../includes/const_param.h'
 
*   colat           - Colatitude (used so that loc_coord and pos_geod
*                   - can be the same arrays)
*   colat_trun      - Coltitude truncated to the nearest 10 asecs
*                   - Used to compute small circle radius
*   loc_coord(3)    - Colatitude (rad), longitude (rad) and
*                   - height (m)
*   pos_geod(3)     - Distance from equator, distance from Greenwich
*                   - meridian, and height above ellipsiod (all m)
 
 
      real*8 colat, colat_trun, loc_coord(3), pos_geod(3), xyz(3)

      integer*4 zone    ! UTM zone (not used)
      character*1 hemi  ! Hemisphere (not used)


****  Get the colat and loc_coord(1)
      loc_coord(1) = pi/2 - pos_geod(1)/earth_rad 

      colat = loc_coord(1)
*                                          ! 2.d4 = 1/10" in radians
      colat_trun = anint(colat*2.d4)/ 2.d4

*     Now get longitude

      loc_coord(2) = pos_geod(2)/(earth_rad * sin(colat_trun))
 
***** Height above ellipssiod
 
      loc_coord(3) =  pos_geod(3)

****  Now get XYZ
      call geod_to_geod(loc_coord, xyz, 'GEOD', 'XYZ',
     .                  'WGS84','WGS84', zone, hemi)

****  Thats all
      return
      end

CTITLE CONVERT_GTOL

      subroutine convert_gtol( ref_llu, neu, loc_coord, xyz )

 
      implicit none

*     Routine to convert neu values back to geodetic lat,long, height
*     and Cartessian coordinates when the NEU (E specifically) may
*     have been created with a co-latitude different from the ref_neu
*     value (after the 0.1" arc second quanitification).

      include '../includes/const_param.h'

* PASSED VARIABLES
* INPUT
      real*8 ref_llu(3)   ! Reference Lat/Long (rads) and height (m)
      real*8 neu(3)       ! GLOBK NEU values to be converted
* OUTPUT
      real*8 loc_coord(3) ! Geodetiic lat/long height (rad/m)
      real*8 xyz(3)       ! Cartessain coordinates (m)

* LOCAL VARIABLES

      real*8 loc_llu(3)   ! Lat/Long/Up computed from NEU assuming
                          ! no switch in the reference latitude
      real*8 std_llu(3)   ! Lat/Long/Up with latitude replaced by the
                          ! reference value
      real*8 std_neu(3)   ! NEU values computed using reference latitude
                          ! (normally will match NEU value).
      real*8 dEast        ! Difference in East.  Normally should be zero
                          ! but a change reference latitude will make it
                          ! about 800 meters from mid-latitude site

****  First, convert the NEU values to loc_llu values 
      call geod_to_loc(neu,loc_llu,xyz)
      std_llu = loc_llu
*     Replace the latitude with the reference value to test the conversion
      std_llu(1) = ref_llu(1)
*     Now convert the NEU using reference latitude
      call loc_to_geod(std_llu, std_neu )
*     Get the East difference.  Factor corrects for the wrong radius being used
*     to get longitude in first geod_to_loc call.
      dEast = (neu(2)-std_neu(2))/(1-(neu(2)-std_neu(2))/neu(2))
      if( abs(dEast).gt.1.d-4 ) then
*         Then adjust the East back to value correct for correct latitude
          neu(2) = neu(2) + dEast
*         Now compute the correct XYZ values
          call geod_to_loc(neu,loc_llu,xyz)
!         print *,'dEast ',dEast,neu(2),std_neu(2),' XYZ', xyz
      endif

****  Thats all
      return
      end

