CTITLE DRY_MARI
 
      subroutine dry_mari( in, dry_map, dry_zen, wet_zen)
 
 
      implicit none 

*     Routine to compute the dry marini mapping function.
*                                  1:06 PM  MON., 16  FEB., 1987
*
* MOD TAH 880930 Changed marini formula to the opticial version
*     to see effect.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Site index either 1 or 2
 
      integer*4 in
 
*   fphih       - Gravity term in formulas
 
      real*4 fphih
 
*   A, B        - A and B terms in Marini formula
*   dry_map(2)  - Delay and rate mapping function. (ps/s and (fs/s)/ps)
*   dry_zen(2)  - Dry zenith delay term (part of A in Marini formula)
*   gamma       - Part of the expression (used for convenience)
*   wet_zen(2)  - Wet zenith delay term (part of A in Marini formula)
*   K           - K constant in marini murray formula
*   temp_K      - Temperature in Kelvin
 
      real*8 A, B, dry_map(2), dry_zen(2), gamma, wet_zen(2), K,
     .    temp_K
 
***** Start, get the gravity function
 
      call phi_function(latitudes(site(in)),ellip_hgt(site(in)),fphih)
 
*     Compute some terms we will need
 
*                                                        ! convert to meters
      A = (dry_zen(1) + wet_zen(1))*1.d-12 * vel_light
C     B = 2.644d-3*exp(-0.14372d-3*ellip_hgt(site(in)) )/fphih
C MOD TAH 880929 change definition of B function
      Temp_k = Temp_c(in,1) + 273.15d0
 
      K = 1.163d0 - 0.00968d0*cos(2*latitudes(site(in)))
     .            - 0.00104d0*Temp_K + 0.00001435d0*pressure(in,1)
 
*     The two constants below are the Marini values * 77.6^2/80.343^2
      B = 1.011d-8*pressure(in,1)*Temp_K*K +
     .    4.416d-8*pressure(in,1)**2*2/(Temp_K*(3-1/K))
 
C     gamma = B/(A+B)/ ( sin(elev(in,1)) + 0.015d0 )
C MOD TAH 880930 Changed last constant to Marini value
      gamma = B/(A+B)/ ( sin(elev(in,1)) + 0.010d0 )
 
*     Do the mapping functions
 
      dry_map(1) =  (A+B)/A/( sin( elev(in,1) ) + gamma)
      dry_map(2) = -(A+B)/A/( sin( elev(in,1) ) + gamma)**2 *
     .    ( cos(elev(in,1)) -
     .      B/(A+B)/( sin(elev(in,1))+0.015d0 )**2*cos(elev(in,1)) )*
*                                 ! Make units (fs/s)/ps
     .      elev(in,2)*1000.d0
 
***** Thats all
      return
      end
 
