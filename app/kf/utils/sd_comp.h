
* SD_COMP.H   -  Include file for use with the SD_COMP.F
*                suboutine.  This include declares the variables
*                to be used in evaluating the diurnal and semidiurnal
*                variations in pole position and UT1. 
 
* DEFINE THE PARAMETERS OF THE SYSTEM
 
*   num_sdc_ut1      - number of UT1 entries in data statment.
*   num_sdc_xy       - number of xy pole entries in data statement.
 
      integer*4 num_sdc_ut1, num_sdc_xy
 
      parameter ( num_sdc_ut1   = 10 )
      parameter ( num_sdc_xy    = 12 )

* PHYSICAL CONSTANTS NEEDED FOR SD_COMP

*   pi          - Define here to full precision
*   rad_to_deg  - Conversion from radians to degs.
*   DJ2000      - Julian date of J2000
*   sec360      - number of seconds in 360 degreees.

      real*8 pi, rad_to_deg, DJ2000, sec360
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )

*     Computed quanities
      parameter ( rad_to_deg    = 180.d0   /pi         )

 
* DECLARE THE ARRAYS TO CONTAIN THE TIDAL COEFFICIENTS FOR UT1
* AND POLE POSITION
 
*   sdc_ut1_arg(6,num_sdc_ut1)  - Arguments (Browns + (gst+pi))
*                                 for ut1.  Same ordering as IAU
*                                 nutation series (i.e., l,l',F,D,
*                                 Om)
*   sdc_xy_arg(6,num_sdc_xy)    - Arguments for xy pole (as above)
 
      integer*4 sdc_ut1_arg(6,num_sdc_ut1), sdc_xy_arg(6,num_sdc_xy)
 
*   sdc_ut1_val(2,num_sdc_ut1)  - Values for cosine and sine for
*                               - UT1 (mas)
*   sdc_xy_val(2,num_sdc_xy)    - Values for cosine and sine for
*                               - x and y (cosine and sine again)
*                               - (mas)

      real*8 sdc_ut1_val(2,num_sdc_ut1), sdc_xy_val(2,num_sdc_xy)

* COMMON DECLARATIONS
 
      common / sdc_comm / sdc_ut1_val, sdc_xy_val, 
     .                    sdc_ut1_arg, sdc_xy_arg
 
*-------------------------------------------------------------------

