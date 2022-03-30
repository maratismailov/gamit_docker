CTITLE WET_SAAS
 
      subroutine wet_saas( in, wet_zen )

      implicit none 
 
 
*     Routine to compute wet zenith delay using Saastamoinen formula
*     at site in of a baseine (The in index (=1 or 2) is used only to
*     access the relative humdity and temperature values values.
*     The sign of the correction is handled in the add_medium routine.
*     NOTE: In computing the rate, we currently neglect the temperature
*           derivative
*
*                                         11:15 PM SUN., 15 Feb., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Indicates site in baseline (either 1 or 2)
 
      integer*4 in
 
*   axis_temp   - Temperature at axis of telescope
*   fphih       - Latitudes and height function for gravity
*   wet_pp      - Saturated water vapor partial pressure (mbars)
 
      real*4 axis_temp, fphih, wet_pp
 
*   dwet_dT     - Partial derivative of the wet saturated press.
*               - with respect to temperature
*   wet_zen(2)  - Delay and rate for the 'wet' zenith delay.
 
      real*8 dwet_dT, wet_zen(2)
 
*   kbit        - Function to see if bit is turned on
 
      logical kbit
 
***** Before doing any calculations, see if temp and rh are available
*     If they are not, then turn off this meduim contribution and default
*     to constant value calculation.
 
*                                                    ! Temp_c  available
      if( kbit(avail_met(site(in)),2) .and.
*                                                    ! RH available
     .    kbit(avail_met(site(in)),3)       ) then
 
*         Now see if the values for this observation are good
          if( .not.kbit(medium_flag(in),2) .and.
*                                                      ! Data Ok, compute
     .        .not.kbit(medium_flag(in),3) ) then
 
*****         Compute the gravity function
 
              call phi_function(latitudes(site(in)),ellip_hgt(site(in)),
     .            fphih)
 
              axis_temp  = temp_C(in,1) - barometer_hgt(site(in))*6.5d-3
 
*****         Calculate the saturation pressure of water vapor (mbars)
              call wet_pressure(axis_temp,1.0,wet_pp)
 
*                                    ! Note: 0.11862 = 6.11*5.3/273
              dwet_dT    = (-0.11862*((axis_temp+273.15)/273)**(-6.3)*
     .                 exp(25.2*axis_temp/(axis_temp+273.15)) +
     .                 wet_pp*25.2/(axis_temp+273.15) *
     .                 (1 - axis_temp/(axis_temp+273.15)) ) *
     .                                             temp_C(in,2)
 
              wet_zen(1) = 0.0022768d0*
     .                     (1255/(axis_temp+273.15d0)+0.05d0)*
*                                                                 ! ps
     .                 rel_hum(in,1)*wet_pp/fphih/vel_light*1.d12
 
 
              wet_zen(2) = 0.0022768d0*
     .                    (1255/(axis_temp+273.15d0)+0.05d0)*
     .                 (rel_hum(in,2)*wet_pp + rel_hum(in,1)*dwet_dT)/
*                                                                  ! fs/s
     .                                       fphih/vel_light*1.d15
 
*                                         ! Data bad, flag data
          else
 
              call sbit( data_flag,5,1)
          end if
 
*                         ! Turn off this calcualtion, and go straight
      else
*                         ! user_constant value
          call sbit(cont_medium(site(in)), 7, 0)
*                                                 ! Turn on constant
          call sbit(cont_medium(site(in)),13, 1)
 
      end if
 
****  Thats all
      return
      end
 
