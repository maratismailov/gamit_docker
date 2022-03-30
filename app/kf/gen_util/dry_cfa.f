CTITLE DRY_CFA
 
      subroutine dry_cfa( in, dry_map )

      implicit none 
 
 
*     Routine to compute the CfA-2.2 Dry mapping function.  If the pressure
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
 
*   A,B,C,D     - the A,B,C, and D coeffiecents in Cfa-2.2
*   beta        - Term in CfA-2.2
*   cose        - Cos of elevation angle
*   dry_map(2)  - Delay and rate mapping function
*   gamma       - Term in CfA-2.2
*   lapse_rate  - Lapse rate in K/km
*   sine        - Sine of elevation angle
*   tropo_hgt   - Height of tropopause (km)
 
      real*8 A,B,C, beta, cose, dry_map(2), gamma, lapse_rate, sine,
     .    tropo_hgt
 
*   kbit        - Ckeck to see if bit set
 
      logical kbit
 
***** Start, see if we pressure, temperature and rel.hum.
 
      if( kbit(avail_met(site(in)),1) .and.
     .    kbit(avail_met(site(in)),2) .and.
*                                                  ! all mets available
     .    kbit(avail_met(site(in)),3)       ) Then
 
*         See if they are good for this observation
          if( .not.kbit(medium_flag(in),1) .and.
     .        .not.kbit(medium_flag(in),2) .and.
*                                                    ! Good for this one
     .        .not.kbit(medium_flag(in),3)  )  Then
 
****          Check which lapse rate to use
*                                                         ! DM_USR
              if( kbit( cont_medium(site(in)),17) ) then
*                                       ! Use defaults for the moment
                  lapse_rate = -6.5
                  tropo_hgt  = 11.231
*                                       ! Use new values
              else
*                                       ! GSFC values
                  lapse_rate = -5.6
*                                       !  ""    ""
                  tropo_hgt  = 10.000
              end if
 
*****         See if we should update tropo_hgt.  MOVED CODE TO
*             AFTER Atm_temp has been computed.
c             if( kbit( cont_medium(site(in)),18 ) ) then
c                 tropo_hgt = ( 220.d0 -
c    .                         (atm_max_temp(site(in))+273.15d0 -
c    .                          ellip_hgt(site(in))*
c    .                           lapse_rate/1000.d0) )/
c    .                          lapse_rate
c
c             end if
 
*****         Check to see if we should use actual temperature or
*             maximum temperature for mapping function
*                                                         ! Use obs. temp
              if( kbit( cont_medium(site(in)),15 ) ) then
                  axis_temp  = temp_C(in,1) -
     .                         barometer_hgt(site(in))*6.5d-3
*                                                         ! Use max. temp
              else
                  axis_temp  = atm_max_temp(site(in)) -
     .                         barometer_hgt(site(in))*6.5d-3
              end if
 
*MOD TAH 881013:  Now see if we should compute lapse rate and tropopuase
*             height. NOTE: We are including a bias of -4 C in the
*             surface temperature
              if( kbit( cont_medium(site(in)),18 ) ) then
                  lapse_rate = -5.439d0 - 0.0508d0*axis_temp
                  tropo_hgt  = 10.342d0 + 0.133d0 *axis_temp
*                 Include 4 C bias in axis temperature
                  axis_temp = axis_temp + 4.d0
              end if
 
*****         Calculate the partial pressure of water vapor (mbars). Use
*             actual temperature (Neglect effects of axis height)
 
              call wet_pressure(temp_C(in,1) ,rel_hum(in,1),wet_pp)
 
              A = 0.001185d0 * ( 1.d0
     .                + 0.6071d-4*(pressure(in,1)-1000.d0)
     .                - 0.1471d-3*(wet_pp)
     .                + 0.3072d-2*(axis_temp  - 20.d0)
     .                + 0.1965d-1*(lapse_rate + 6.5d0)
     .                - 0.5645d-2*(tropo_hgt  - 11.231)  )
 
              B = 0.001144d0 * ( 1.d0
     .                + 0.1164d-4*(pressure(in,1)-1000.d0)
     .                + 0.2795d-3*(wet_pp)
     .                + 0.3109d-2*(axis_temp  - 20.d0)
     .                + 0.3038d-1*(lapse_rate + 6.5d0)
     .                - 0.1217d-1*(tropo_hgt  - 11.231)  )
 
 
              C = -0.0090d0
 
              beta  = B/( sin(elev(in,1)) + C )
              gamma = A/( tan(elev(in,1)) + beta)
              sine  = sin( elev(in,1) )
              cose  = cos( elev(in,1) )
 
              dry_map(1) = 1.d0 / ( sine + gamma )
 
              dry_map(2) = -1.d0/ ( sine + gamma )**2 *
     .                    ( cose - A/(tan(elev(in,1))+beta)**2 *
     .                    ( 1.d0/cose**2 - B *cose/(sine+C)**2 ) ) *
*                                              !Units (fs/s)/s
     .                    elev(in,2) * 1000.d0
 
 
*                         ! Weather data bag, flag data
          else
 
              call sbit( data_flag,5,1 )
 
          end if
 
*                 ! Not enough data, go to CHAO
      ELSE
 
          if ( kbit( cont_medium(site(in)),15 ) ) then
*                                                      ! Turn off CfA-2.2
              call sbit( cont_medium(site(in)),15,0 )
          else
*                                                      ! Turn off CfA-2.2 Max
              call sbit( cont_medium(site(in)),16,0 )
          end if
 
*                                                  ! Turn on  CHAO
          call sbit( cont_medium(site(in)),19,1 )
 
      END IF
 
***** Thats all
      return
      end
 
