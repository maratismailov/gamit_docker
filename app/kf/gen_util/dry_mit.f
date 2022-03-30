CTITLE DRY_MIT
 
      subroutine dry_mit( in, dry_map )
 
 
      implicit none 

*     Routine to compute the mit-2.2 Dry mapping function.  If the pressure
*     temperature and humidity are not available, the mapping is set back
*     to CHAO.
* MOD TAH 870730:  Added calculation of mapping function with maximum
*     temperature, rather than observed value.  Max value may be a better
*     representation of the whole atmosphere.
* MOD TAH 881013: Added Temperature dependent lapse rate and Heioght of
*     of tropopause based on Radiosonde data.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site in baseline (1 or 2)
 
      integer*4 in
 
*   axis_temp   - Temperature at intersection of axes
*   wet_pp      - Partial pressure of water vapor (mb)
 
      real*4 axis_temp, wet_pp
 
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
*    Tb_seas    - Seasonal bias in temperature

      real*4 Tc_seas, P_seas, rh_seas, tb_seas

 
*   kbit        - Ckeck to see if bit set
 
      logical kbit
 
***** Start, see if we pressure, temperature and rel.hum., or
*     if we are using a seasonal model
 
      if( .not.(kbit(avail_met(site(in)),1) .and.
     .          kbit(avail_met(site(in)),2) .and.
     .          kbit(avail_met(site(in)),3) ) .or.
     .    kbit( cont_medium(site(in)),28)          ) then

*         Set modifier to use the seasonal temperature values
*         (Save the estimated temperature so that it will used)
          call met_seasonal( Tc_seas, P_seas, rh_seas, tb_seas,
     .         epoch, latitudes(site(in)), ellip_hgt(site(in)))
          temp_C(in,1)   = Tc_seas
          pressure(in,1) = P_seas
          rel_hum(in,1)  = rh_seas 
*         Force solution to say use seasonal
          call sbit(cont_medium(site(in)),28,1)
      end if

*     See if they are good for this observation or if we are using seasonal model
      if( (.not.kbit(medium_flag(in),1) .and.
     .    .not.kbit(medium_flag(in),2) .and.
     .    .not.kbit(medium_flag(in),3)) .or.
     .         kbit(cont_medium(site(in)),28)  )  Then
 
*         Test to see if we are using the max temperature value.
          if( kbit(cont_medium(site(in)),27) .and.
     .        .not.  kbit(cont_medium(site(in)),28) ) then
              axis_temp = atm_max_temp(site(in)) -
     .                 barometer_hgt(site(in))*6.5d-3
          else
              axis_temp  = temp_C(in,1) -
     .                 barometer_hgt(site(in))*6.5d-3
          end if

*         adjust temperature for surface bias provided we are not using seasonal
*         model.
          if( .not. kbit(cont_medium(site(in)),28) ) then
              axis_temp = axis_temp - tbias_con(site(in))
          end if
 
*****     Calculate the partial pressure of water vapor (mbars). Use
*         actual temperature (Neglect effects of axis height)
 
          call wet_pressure(temp_C(in,1) ,rel_hum(in,1),wet_pp)
 
          A = 0.00125003d0* ( 1.d0
     .            + 6.29580d-5  *(pressure(in,1)-1000.d0)
     .            + 1.67337d-5  *(wet_pp)
     .            + 3.36152d-3  *(axis_temp  - 10.d0)
     .            + 1.99915d-2  *(lapse_con(site(in)) + 5.6d0)
     .            - 4.7599 d-3  *(ht_con(site(in))  - 10.d0)  )

          B = 0.00312108d0* ( 1.d0
     .            + 3.84500d-5  *(pressure(in,1)-1000.d0)
     .            + 6.62430d-5  *(wet_pp)
     .            + 4.62404d-3  *(axis_temp  - 10.d0)
     .            + 4.39239d-2  *(lapse_con(site(in)) + 5.6d0)
     .            - 1.58279d-2  *(ht_con(site(in))  - 10.d0)  )
 
 
          C = 0.06945748d0* ( 1.d0
     .            + 3.956  d-5  *(pressure(in,1)-1000.d0)
     .            + 0.0    d-5  *(wet_pp)
     .            + 2.94868d-3  *(axis_temp  - 10.d0)
     .            + 3.02481d-2  *(lapse_con(site(in)) + 5.6d0)
     .            - 1.29281d-2  *(ht_con(site(in))  - 10.d0)  )

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
 
      end if
 
***** Thats all
      return
      end
 
