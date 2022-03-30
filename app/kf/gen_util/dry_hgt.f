CTITLE DRY_HGT
 
      subroutine dry_hgt (in, dry_zen )
 
      implicit none 

 
*     Rotuine to compute zenith delay based on pressure computed from
*     height alone.  Hydrostatic equilrium is used to get the pressure
*
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Indicates site in baseline (in=1 or 2)
 
      integer*4 in
 
*   fphih       - Gravity function
 
C     real*4 fphih
 
*   axis_press  - Pressure computed at axis of telescope
*   axis_temp   - Temperature at axis (computed assumming
*               - 293.15K at MSL)
*   dry_zen(2)  - Zenith delay and rate.
 
C     real*8 axis_press, axis_temp, dry_zen(2)
      real*8 dry_zen(2)
 
***** Compute temperature and pressure at axis

      call dry_hgt_site( site(in), dry_zen )
 
C     axis_temp = 293.15d0 - 6.5d-3*(ellip_hgt(site(in)) +
C    .                               barometer_hgt(site(in)) )
 
C     axis_press = 1013.25d0*(axis_temp/293.15d0)**5.26d0
 
C     call phi_function(latitudes(site(in)),ellip_hgt(site(in)),fphih)
 
C     dry_zen(1) = (0.0022768d0*axis_press)/vel_light/
*                                                         ! ps
C    .              fphih*1.d12
*                                                         ! fs/s
C     dry_zen(2) =  0.00d0
 
***** Thats all
      return
      end
 
