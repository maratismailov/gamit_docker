
CTITLE dry_mtt
 
      subroutine dry_mtt( in, dry_map )
 
      implicit none 

*     Routine to compute the new dry_mit mapping function which
*     depends only on temperature (and position)
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site in baseline (1 or 2)
 
      integer*4 in
 
*   axis_temp   - Temperature at intersection of axes
 
      real*4 axis_temp
 
*   A,B,C,D     - the A,B,C, and D coeffiecents in mit-2.2
*   beta        - Term in mit-2.2
*   cose        - Cos of elevation angle
*   dry_map(2)  - Delay and rate mapping function
*   gamma       - Term in mit-2.2
*   sine        - Sine of elevation angle
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith 

      real*8 A,B,C, beta, cose, dry_map(2), gamma, sine, topcon

*    Tc_seas    - Seasonal value of temperature (C)
*    P_seas     - Seasonal Pressure (mbar)
*    rh_seas    - Seasonal relative humidity
*    tb_seas    - Seasonal bias in temperater
*    cosphi     - Cosine of latitude 
*    hs_km      - Height of site in kms.

      real*4 Tc_seas, P_seas, rh_seas, cosphi, hs_km, tb_seas
 
*   kbit        - Ckeck to see if bit set
 
      logical kbit
 
***** Start, see if temperature is available or if we are using
*     seasonal values for the met data..
 
      if( .not. kbit(avail_met(site(in)),2) .or.
     .          kbit(cont_medium(site(in)),28) ) then

*         Set modifier to use the seasonal temperature values
*         (Save the estimated temperature so that it will used)
          call met_seasonal( Tc_seas, P_seas, rh_seas, tb_seas, 
     .         epoch, latitudes(site(in)), ellip_hgt(site(in)))
          temp_C(in,1) = Tc_seas
*         Force solution to say use seasonal
          call sbit(cont_medium(site(in)),28,1)
      end if

 
*     See if they are good for this observation or if we are using seasonal
*     model
      if( .not.kbit(medium_flag(in),2) .or.
     .         kbit(cont_medium(site(in)),28)   ) then
 
*****     Get the axis temperature 
*         Test to see if we are using the max temperature value.
          if( kbit(cont_medium(site(in)),27) .and.
     .        .not.  kbit(cont_medium(site(in)),28) ) then
              axis_temp = atm_max_temp(site(in)) - 
     .                 barometer_hgt(site(in))*6.5d-3
          else
              axis_temp  = temp_C(in,1) -
     .                 barometer_hgt(site(in))*6.5d-3
          end if

*         adjust temperature for surface bias. Provided we are not using the
*         seasonal model
          if( .not.kbit(cont_medium(site(in)),28) ) then
               axis_temp = axis_temp - tbias_con(site(in))
          end if

*         Now compute the coefficients in the mapping function.
          cosphi = cos(latitudes(site(in)))
          hs_km  = ellip_hgt(site(in))/1000.d0
 
*         Compute the coefficienst of the mapping function
          A  = 1.23200d-3 + 0.01391d-3*cosphi - 0.02089d-3*hs_km
     .                    + 0.002154d-3*(axis_temp - 10.d0)
          B  = 3.16116d-3 - 0.16004d-3*cosphi - 0.03306d-3*hs_km
     .                    + 0.002064d-3*(axis_temp - 10.d0)
          C  =71.24372d-3 - 4.29342d-3*cosphi - 0.14908d-3*hs_km
     .                    - 0.002098d-3*(axis_temp - 10.d0)
 
          beta  = B/( sin(elev(in,1)) + C )
          gamma = A/( sin(elev(in,1)) + beta)
          sine  = sin( elev(in,1) )
          cose  = cos( elev(in,1) )

          topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))

          dry_map(1) = topcon / ( sine + gamma )

          dry_map(2) = -topcon / ( sine + gamma )**2 *
     .                ( cose - A/( sine + beta)**2 * cose *
     .                ( 1.d0 - B/( sine + C )**2 ) ) *
*                                              !Units (fs/s)/s
     .                elev(in,2) * 1000.d0
 
*                         ! Weather data bag, flag data
      else

          call sbit( data_flag,5,1 )
          call dry_chao( in, dry_map)
 
      end if
 
***** Thats all
      return
      end
 
