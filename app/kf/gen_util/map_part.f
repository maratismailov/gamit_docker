
CTITLE map_part
 
      subroutine map_part( type, in, part )

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
*   i       - Loop counter
 
      integer*4 in, i

*   type    - Type of partial de or rate

      character*(*) type
 
*   axis_temp   - Temperature at intersection of axes
*   wet_pp      - Partial pressure of water vapor (mb)
 
      real*4 axis_temp, wet_pp
 
*   A,B,C,D     - the A,B,C, and D coeffiecents in mit-2.2
*   beta        - Term in mit-2.2
*   cose        - Cos of elevation angle
*   part(2)  - Delay and rate mapping function
*   gamma       - Term in mit-2.2
*   lapse_rate  - Lapse rate in K/km
*   sine        - Sine of elevation angle
*   tropo_hgt   - Height of tropopause (km)
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith 
*   map_val(2)  - Mapping function +-1 C about nominal
*   dry_zen(2)  - Zenith delay.

      real*4 part(2)
 
      real*8 A,B,C, beta, cose, gamma,  sine,
     .     topcon, map_val(2), dry_zen(2)
 
*   kbit        - Ckeck to see if bit set
 
c      logical kbit
 
***** Start, see if we pressure, temperature and rel.hum.
 
***** Get the axis temperature 
      axis_temp  = temp_C(in,1) -
C     axis_temp  = atm_max_temp(site(in)) -
     .             barometer_hgt(site(in))*6.5d-3

*     adjust temperature for surface bias
      axis_temp = axis_temp + tbias_con(site(in))

***** Calculate the partial pressure of water vapor (mbars). Use
*     actual temperature (Neglect effects of axis height)
 
      call wet_pressure(temp_C(in,1) ,rel_hum(in,1),wet_pp)

*     Now loop over the two temperatures

      do i = 1, 2
          if( i.eq.1 ) axis_temp = axis_temp + 1
          if( i.eq.2 ) axis_temp = axis_temp - 2
 
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

          if( type(1:2).eq.'de' ) then

              map_val(i) = topcon / ( sine + gamma )
          else
              map_val(i) = -topcon / ( sine + gamma )**2 *
     .                    ( cose - A/( sine + beta)**2 * cose *
     .                    ( 1.d0 - B/( sine + C )**2 ) ) *
     .                    elev(in,2) * 1000.d0
*                                              !Units (fs/s)/s
          end if
      end do

*     Now compute the partial
      call dry_saas( in, dry_zen )
      part(1) = dry_zen(1)*(map_val(1)-map_val(2))/2.d0
 
 
***** Thats all
      return
      end
 
