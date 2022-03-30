CTITLE DRY_SAAS
 
      subroutine dry_saas( in, dry_zen )
 
 
      implicit none 

*     Routine to compute dry zenith delay using Saastamoinen formula
*     at site in of a baseine (The in index (=1 or 2) is used only to
*     access the pressure values.  The sign of the correction is handled
*     in the add_medium routine.
*
*                                         10:36 PM SUN., 15 Feb., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Indicates site in baseline (either 1 or 2)
 
      integer*4 in
 
*   fphih       - Latitudes and height function for gravity
 
      real*4 fphih
 
*   axis_press  - Pressure at interestion of axes
*   axis_temp   - Temperature at intersection of axes
*   dry_zen(2)  - Delay and rate for the 'dry' zenith delay.
 
      real*8 axis_press, axis_temp, dry_zen(2)
 
*   kbit        - Function to see if bit is turned on
 
      logical kbit
 
***** Before doing any calculations, see if pressure is available
*     If it isnot, then turn off this meduim contribution and default
*     to height calculation.
 
*                                              ! Pressure available
      if( kbit(avail_met(site(in)),1) ) then
 
*         Now see if pressure good fro this particular observation
*                                                 ! Pressure good
          if( .not.kbit(medium_flag(in),1) ) then
 
*****         Compute the gravity function
 
              call phi_function(latitudes(site(in)),ellip_hgt(site(in)),
     .            fphih)
 
              axis_temp  = temp_C(in,1) - barometer_hgt(site(in))*6.5d-3
              axis_press = pressure(in,1)*
     .            ( (273.15+axis_temp)/(273.15+temp_C(in,1)) )**5.26d0
 
              dry_zen(1) = (0.0022768d0*axis_press)/vel_light/
*                                                                 ! ps
     .                      fphih*1.d12
              dry_zen(2) = (0.0022768d0*pressure(in,2))/vel_light/
*                                                                 ! fs/s
     .                      fphih*1.d15
 
*                                          ! Pressure bad for this obs, flag
          else
*                                          ! data_flag
              call sbit(data_flag,5,1)
          end if
 
*                         ! Turn off this calcualtion, and go straight
      else
*                         ! height dependent pressure
          call sbit(cont_medium(site(in)),1, 0)
*                                                ! Turn on constant
          call sbit(cont_medium(site(in)),6, 1)
 
      end if
 
****  Thats all
      return
      end
 
