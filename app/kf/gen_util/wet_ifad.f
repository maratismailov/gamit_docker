CTITLE WET_IFAD
 
      subroutine wet_ifad(in,wet_map)

      implicit none 
 
*     J.L. Davis                   2:31 PM  FRI., 10  JULY, 1987
*
*     The Ifadis global wet mapping function.  Reference:
*
*     Ifadis, Ioannis, The atmospheric delay of radio waves:
*     Modeling the elevation dependence on a global scale,
*     Chalmers (Sweden) University Technical Report No. 38L,
*     School of  Electrical and Computer Engineering, p. 98,
*     1986.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*       kbit                    - Check if bit set
 
      logical kbit
 
*       in                      - Index to site in baseline (1 or 2)
 
      integer*4 in
 
*       axis_temp               - Temperature at intersection of axes
*   ,   wet_pp                  - Partial pressure of H2O v
 
      real*4 axis_temp, wet_pp
 
*       a1, a2, a3              - Terms in the formula
*   ,   cose                    - Cosine of elevation angle
*   ,   sine                    - Sine of elevation angle
*   ,   wet_map(2)              - Mapping function and air mass
*                               -   rate-of-change
 
      real*8 a1, a2, a3, cose, sine, wet_map(2)
 
***** Start: See if we have pressure, temperature, and rel. hum.
      if (kbit(avail_met(site(in)),1) .and.
     .    kbit(avail_met(site(in)),2) .and.
     .    kbit(avail_met(site(in)),3)      ) then
 
*****     See if they are OK for this obs
          if (.not. kbit(medium_flag(in),1) .and.
     .        .not. kbit(medium_flag(in),2) .and.
     .        .not. kbit(medium_flag(in),3)      ) then
 
*****         Get the temperature at the level of the barometer
              axis_temp = temp_C(in,1) - barometer_hgt(site(in))*6.5D-03
 
*****         Get the partial pressure in mbars
              call wet_pressure(axis_temp,rel_hum(in,1),wet_pp)
 
*****         Determine the a-terms
 
              a1 = 0.5236D-03 + 0.2471D-06 * (pressure(1,in) - 1000.0)
     .                        - 0.1724D-06 * (axis_temp - 15.0)
     .                        + 0.1328D-04 * sqrt(wet_pp)
 
              a2 = 0.1705D-02 + 0.7384D-06 * (pressure(1,in) - 1000.0)
     .                        + 0.3767D-06 * (axis_temp - 15.0)
     .                        + 0.2147D-04 * sqrt(wet_pp)
 
              a3 = 0.5917D-01
 
*****         Calculate some trig functions
              sine = sin(elev(in,1))
              cose = cos(elev(in,1))
 
*****         Calculate the wet mapping function
              wet_map(1) = 1 / (sine + a1 / (sine + a2 / (sine + a3)))
 
*****         Calculate the rate-of-change of air mass
              wet_map(2) = - wet_map(1) ** 2 * elev(in,2) * cose
     .         * (1 - a1 * (1 - a2 / (sine + a3) ** 2)
     .                   / (sine + a2 / (sine + a3)) ** 2)
 
*****         Multiply by 1000 to give units of (fs/s)/ps
              wet_map(2) = 1000 * wet_map(2)
 
          else
 
*****         Flag as having bad WX data
              call sbit(data_flag,5,1)
 
          end if
 
      else
 
*****     Tell program to use Chao
*                                                ! Turn of Ifadis
          call sbit(cont_medium(site(in)),23,0)
*                                                ! Turn on Chao
          call sbit(cont_medium(site(in)),19,1)
 
      end if
 
      END
 
