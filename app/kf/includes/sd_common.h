 
 
*     Include file for the diurnal and semidiurnal corrections for
*     UT1, pole position and tides
 
*   max_sd_ut1      - Maximum number of UT1 entries
*   max_sd_xy       - Maximum number of xy pole entries
*   max_sd_tides    - Maximum number of tidal entries
*   max_sd_nut      - Maximum number of nutation series coefficient
*                     corrections which can be entered.
 
      integer*4 max_sd_ut1, max_sd_xy, max_sd_tides, max_sd_nut
 
      parameter ( max_sd_ut1   = 20 )
      parameter ( max_sd_xy    = 60 )
      parameter ( max_sd_tides =100 )
      parameter ( max_sd_nut   = 40 )
 
 
*     Declarations
 
*   sd_ut1_arg(6,max_sd_ut1)    - Arguments (Browns + gst + pi)
*                               - for ut1
*   sd_xy_arg(6,max_sd_xy)      - Arguments for xy pole (as above)
*   sd_tides_arg(6,max_sd_tides)    - As above for tides
*   sd_nut_arg(6,max_sd_nut)    - Arguments for the nutation series
*                                (the sixth value gst+pi is ignored)
*   sd_ut1_num                  - Number of values in UT1
*   sd_xy_num                   - Numbet of values for xy pole
*   sd_tides_num                - Number of values for tides
*   sd_nut_num                  - Number of values for nutations 
*   sd_tides_site(max_sd_tides) - Site number for each of the tide
*                               - corrections
 
      integer*4 sd_ut1_arg(6,max_sd_ut1), sd_xy_arg(6,max_sd_xy),
     .    sd_tides_arg(6,max_sd_tides), sd_nut_arg(6,max_sd_nut), 
     .    sd_ut1_num, sd_xy_num,  sd_tides_num,sd_nut_num, 
     .    sd_tides_site(max_sd_tides)
 
*   sd_ut1_val(2,max_sd_ut1)    - Values for cosine and sine for
*                               - UT1 (mas)
*   sd_xy_val(2,max_sd_xy)      - Values for cosine and sine for
*                               - x and y (cosine and sine again)
*                               - (mas)
*   sd_tides_val(6,max_sd_tides)    - Values for cosine and sine
*                               - in potential for tides.
*                               - Radial, East and
*                               - North. spherical harmonic depends
*                               - on gst argument.
*   sd_nut_val(4,max_nut_tides) - Values of the corrections to the
*                                 nutation coefficients.  These go
*                                 as long, oblq in phase and then
*                                 out of phase terms

      real*8 sd_ut1_val(2,max_sd_ut1), sd_xy_val(2,max_sd_xy),
     .    sd_tides_val(6,max_sd_tides), sd_nut_val(4,max_sd_nut)

*    sd_ut1_mid(4) - Values of UT1 corrections at mid epoch for the
*                     diurnal, semidiurnal cos and sin terms (mas)
*    sd_xy_mid(6)  - Values of xy pole position corrections at mid
*                    epoch in the same order that the partials are
*                    computed.
*    sd_tides_mid(12,max_sites) - Values of extended Earth tides for
*                    for the diurnal and semidiurnal bands for each
*                    site.
*    sd_nut_mid(2) - Nutation series corrections (as angles) at mid
*                    epoch.

      real*8 sd_ut1_mid(4), sd_xy_mid(6), sd_tides_mid(12,max_sites),
     .       sd_nut_mid(2)

 
* COMMON DECLARATIONS
 
      common / sd_comm / sd_ut1_val, sd_xy_val, sd_tides_val,
     .    sd_nut_val, 
     .    sd_ut1_mid, sd_xy_mid, sd_tides_mid, sd_nut_mid,
     .    sd_ut1_arg, sd_xy_arg, sd_tides_arg, sd_nut_arg,
     .    sd_ut1_num, sd_xy_num, sd_tides_num, sd_nut_num,
     .    sd_tides_site
 
*-------------------------------------------------------------------
