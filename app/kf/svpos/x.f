CTITLE ATM_DEL

      real*8 function atm_del( mjd, geod, elev )

      implicit none

      include '../includes/const_param.h' 

*     Function to return a GPT/GMF model for the atmospheric delay
*

* PASSED  
      real*8 mjd       ! MJD at which delay is to be computed
     .,      geod(3)   ! Geodetic coordinates (colat, long, ht; rads/m) 
     .,      elev      ! Elevation angle (rads)

* LOCAL
      real*8 press, temp_C  ! Pressure and temperature
     .,      undu      ! Geoid height (not used).
     .,      dry_zen, wet_zen   ! Dry and Wet zenith delay
     .,      zend      ! Zenith angle (rad)
     .,      dry_map, wet_map   ! Dry and wet mapping function
     .,      lat       ! Latitude (rad)

***** Call GPT to get pressute and temperatue

      lat = pi/2-geod(1)

      call gpt(mjd,lat,geod(2),geod(3), press, temp_C, undu) 

*     Compute dry zenith delay
      call dry_saas_zen( temp_C,lat ,geod(3), press,   dry_zen)

*     Get mapping function
      call gmf(mjd,lat,geod(2),geod(3), pi/2-elev, dry_map, wet_map)

*     Compute the dry delay (set wet to 0.1 meters) 
      wet_zen = 0.1       
      atm_del = dry_zen*dry_map + wet_zen*wet_map

***** That all
      return 
      end

CTITLE SET_GEOD

      subroutine set_geod

      implicit none

*     Routine to compute geodetic positions based on apriori XYZ coords

      include 'svsp3.h'

*     Use XYZ_to_GEOD to get Geodetic co-lat, longitude (rads) and height from
*     XYZ coordinates


* LOCAL
      real*8 rot_mat(3,3)    ! Rotation to topocentric frame

      call XYZ_to_GEOD(rot_mat, site_xyz, site_geod)
      call XYZ_to_GEOD(rot_mat, diff_xyz, diff_geod)

****  Thats all
      return
      end


