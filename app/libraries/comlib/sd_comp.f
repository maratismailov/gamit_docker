CTITLE SD_COMP

      subroutine sd_comp( jd, dx, dy, dut1 )

      implicit none

*     Routine to compute the diurnal and semidiurnal contributions
*     to the pole position and UT1 using the VLBI derived values
*     for the tidal varitions.

* USAGE:
*     call sd_comp ( jd, dx, dy, dut1 )
*     where <jd>    is a full julian date with fractional part
*                   of the day added (REAL*8 INPUT)
*     and <dx>, <dy>, <dut1> are the tidally coherent contributions
*                   to the x and y pole positions and to UT1.
*                   dx and dy are returned in milliarcseconds,
*                   dut1 in millitimeseconds (REAL*8 OUTPUT)

* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.


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
      external sd_compbd 
*-------------------------------------------------------------------


* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)

* OUTPUT Values
* dx     - Contribution to X pole position (should be added) (mas)
* dy     - Contribution to Y pole position (should be added) (mas)
* dut1   - Contribution to UT1-AT (should be added) (mts)

      real*8 jd, dx, dy, dut1

* LOCAL VARIABLES

*   i      - Loop counters

      integer*4 i,len,rcpar

*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
*   arg         - Angular argument for the correction (rads)
*   fund_arg(6) - Values of the 5 Brown's arguments (l,l',F,D,
*               - Omega) and gst+pi (rads)

      real*8 epoch, arg, fund_arg(6)

      character*80 prog_name
      character*256 message

c     get calling program name for report_stat
      len = rcpar(0,prog_name)

***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(message,100) jd
 100      format('MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          call report_stat('WARNING',prog_name,'lib/sd_comp',' '
     .                    , message,0)
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if

***** Get the fundamental arguments at this epoch

      call tide_angles( epoch, fund_arg)

*     Clear the contributions
      dx   = 0.d0
      dy   = 0.d0
      dut1 = 0.d0

*     Now loop over the UT1 contributions
      do i = 1, num_sdc_ut1

*         Get the argument
          call sdc_arg( sdc_ut1_arg(1,i), fund_arg, arg, 6)

*         Increment the change to UT1
          dut1 = dut1 + sdc_ut1_val(1,i)*cos(arg)
     .                + sdc_ut1_val(2,i)*sin(arg)
      end do

*     Convert dut1 from mas to mts
      dut1 = dut1 / 15.d0

****  Now do polar motion
      do i = 1, num_sdc_xy

*         Get the argument
          call sdc_arg( sdc_xy_arg(1,i), fund_arg, arg, 6)

*         Increment the change to the X anf Y positions

          dx = dx  - sdc_xy_val(1,i)* cos(arg)
     .             + sdc_xy_val(2,i)* sin(arg)
          dy = dy  + sdc_xy_val(1,i)* sin(arg)
     .             + sdc_xy_val(2,i)* cos(arg)

      end do

***** That is all.
      return
      end

