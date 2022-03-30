CTITLE DRY_HGT
 
      subroutine dry_hgt_site (sn, dry_zen )
 
 
      implicit none 

*     Rotuine to compute zenith delay based on pressure computed from
*     height alone.  Hydrostatic equilrium is used to get the pressure
*     This rouitine takes the site number.
*
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   sn      - Indicates site number (1-num_sites)
 
      integer*4 sn
 
*   fphih       - Gravity function
 
      real*4 fphih
 
*   axis_press  - Pressure computed at axis of telescope
*   axis_temp   - Temperature at axis (computed assumming
*               - 293.15K at MSL)
*   dry_zen(2)  - Zenith delay and rate.
 
      real*8 axis_press, axis_temp, dry_zen(2)
 
***** Compute temperature and pressure at axis
 
      axis_temp = 293.15d0 - 6.5d-3*(ellip_hgt(sn) +
     .                               barometer_hgt(sn) )
 
      axis_press = 1013.25d0*(axis_temp/293.15d0)**5.26d0
 
      call phi_function(latitudes(sn),ellip_hgt(sn),fphih)
 
*                                                         ! ps
      dry_zen(1) = (0.0022768d0*axis_press)/vel_light/
     .              fphih*1.d12
*                                                         ! fs/s
      dry_zen(2) =  0.00d0
 
***** Thats all
      return
      end
 
