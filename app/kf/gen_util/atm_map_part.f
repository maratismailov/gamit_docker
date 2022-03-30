
CTITLE ATM_MAP_PARTIAL
 
      subroutine atm_map_partial( type, in, part)

      implicit none 
 
*     This routine computes the atmospheric delay mapping
*     function partial derivative assuming that the error
*     is due to a temperature error. 

      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'

*   in      - Site number in baseline (i.e., 1 or 2)
 
      integer*4 in
 
*   part    - Partial derivative for mapping function (either delay
*             or rate.
*   sign    - Sign of partial (plus for in=2, minus for in=1)
 
      real*4 part, sign
 
* LOCAL VARIABLES 
*   A,B,C,D     - the A,B,C, and D coeffiecents in mit-2.2
*   beta        - Term in mit-2.2
*   cose        - Cos of elevation angle
*   gamma       - Term in mit-2.2
*   sine        - Sine of elevation angle
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith 
*   map_val(2)  - Mapping function +-1 C about nominal
*   dry_zen(2)  - Zenith delay.and its rate of change (latter part 
*                 neglected)
 
      real*8 A,B,C, beta, cose, gamma, sine, topcon, map_val(2), 
     .       dry_zen(2)

      integer*4 i

*    Tc_seas    - Seasonal value of temperature (C)
*    P_seas     - Seasonal Pressure (mbar)
*    rh_seas    - Seasonal relative humidity
*    cosphi     - Cosine of latitude 
*    hs_km      - Height of site in kms.
*    axis_temp  - Temperature at axis of telescope
*    Tb_seas    - Seasonal bias in temperature

      real*4 Tc_seas, P_seas, rh_seas, cosphi, hs_km, axis_temp, tb_seas

 
*   kbit        - Ckeck to see if bit set
 
      logical kbit

*   type    - type of data (either 'delay' or 'rate')

      character*(*) type
 
***** Do some start calculations
 
      if( in.eq.1 ) then
          sign = -1
      else
          sign = +1
      end if
 
***** Start, see if temperature is available or if we are using
*     seasonal values for the met data..
 
      if( .not. kbit(avail_met  (site(in)), 2) .or.
     .          kbit(cont_medium(site(in)),28) ) then

*         Set modifier to use the seasonal temperature values
*         (Save the estimated temperature so that it will used)
          call met_seasonal( Tc_seas, P_seas, rh_seas, tb_seas,
     .         epoch, latitudes(site(in)), ellip_hgt(site(in)))
          temp_C(in,1)   = Tc_seas
          pressure(in,1) = P_seas
*         Force solution to say use seasonal
          call sbit(cont_medium(site(in)),28,1)
      end if

 
***** Get the axis temperature 
*     Test to see if we are using the max temperature value.
      if( kbit(cont_medium(site(in)),27) .and.
     .    .not.  kbit(cont_medium(site(in)),28) ) then
          axis_temp = atm_max_temp(site(in)) - 
     .             barometer_hgt(site(in))*6.5d-3
      else
          axis_temp  = temp_C(in,1) -
     .             barometer_hgt(site(in))*6.5d-3
      end if

*     adjust temperature for surface bias
      axis_temp = axis_temp + tbias_con(site(in))

*     Now compute the coefficients in the mapping function.
      cosphi = cos(latitudes(site(in)))
      hs_km  = ellip_hgt(site(in))/1000.d0

      sine  = sin( elev(in,1) )
      cose  = cos( elev(in,1) )

****  Now compute the values of the mapping function at two 
*     temperatures

      do i = 1,2
         if( i.eq.1 ) axis_temp = axis_temp + 1
         if( i.eq.2 ) axis_temp = axis_temp - 2

*        Compute the coefficienst of the mapping function
         A  = 1.23200d-3 + 0.01391d-3*cosphi - 0.02089d-3*hs_km
     .                   + 0.002154d-3*(axis_temp - 10.d0)
         B  = 3.16116d-3 - 0.16004d-3*cosphi - 0.03306d-3*hs_km
     .                   + 0.002064d-3*(axis_temp - 10.d0)
         C  =71.24372d-3 - 4.29342d-3*cosphi - 0.14908d-3*hs_km
     .                   - 0.002098d-3*(axis_temp - 10.d0)
 
         beta  = B/( sin(elev(in,1)) + C )
         gamma = A/( sin(elev(in,1)) + beta)

         topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))

         if( type(1:2).eq.'de' ) then

             map_val(i) = topcon / ( sine + gamma )
         else
             map_val(i) = -topcon / ( sine + gamma )**2 *
     .                   ( cose - A/( sine + beta)**2 * cose *
     .                   ( 1.d0 - B/( sine + C )**2 ) ) *
     .                   elev(in,2) * 1000.d0
         end if
      end do

*     Now compute the partial.  For the rate partial neglect rate
*     of change of the met quantities.
      call dry_saas( in, dry_zen )
      part = sign*dry_zen(1)*(map_val(1)-map_val(2))/2.d0

***** Thats all
      return
      end
 
