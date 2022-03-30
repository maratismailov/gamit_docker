CTITLE TIDE_ANGLES

      subroutine tide_angles( epoch, fund_arg )

      implicit none

*     Routine to compute the value of the fundamental argument
*     for Brown's arguments.  The sixth entry is returned as GST
*     plus pi.  The additional pi is used for compatability with
*     Doodson's Tide argument.

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

*-------------------------------------------------------------------

* PASSED VARIABLES

* INPUT
* epoch  - Julian date for arguments (jd + fraction of day)

* OUTPUT
* fund_arg(6) -  Brown's arguments plus GST+pi (rads)

      real*8 epoch, fund_arg(6)

* LOCAL VARIABLES
*      cent             - Julian centuries to DJ2000.
*      el,eld           - Mean longitude of moon minus mean
*                       - longitude of moon's perigee (arcsec)
*      elc(5)           - Coefficients for computing el
*      elp,elpd         - Mean longitude of the sun minus mean
*                       - longitude of sun perigee (arcsec)
*      elpc(5)          - Coeffiecents for computing elp
*      f,fd             - Moon's mean longitude minus omega (sec)
*      fc(5)            - Coefficients for computing f
*      d,dd             - Mean elongation of the moon from the
*                       - sun (arcsec)
*      dc(5)            - coefficients for computing d
*      om,omd           - longitude of the ascending node of the
*                       - moon's mean orbit on the elliptic
*                       - measured from the mean equinox of date
*      omc(5)           - Coefficients for computing om.
*      gst              - Greenwich mean sidereal time (rad)

      real*8 cent, el,eld, elc(5), elp, elpd, elpc(5),
     .    f,fd, fc(5), d,dd, dc(5), om,omd, omc(5), gst

*      fract        - fraction of a day from 0:00 hrs UT.
*      Jd_0hr       - Julian date at zero hours UT
*      t_0hr        - Days since DJ2000 at 0:00 hrs UT
*      gstd         - GMST at 0:00 hrs UT1 of day being evaluated
*      diurnv       - Ratio of solar days to sidreal days on
*                     day of evalution.

      real*8 fract, t_0hr, gstd, diurnv, jd_0hr

****  DATA statements for the fundamental arguments.

      data elc    /     0.064d0,    31.310d0,    715922.633d0,
     .             485866.733d0,    1325.0d0 /
      data elpc   /    -0.012d0,    -0.577d0,   1292581.224d0,
     .            1287099.804d0,      99.0d0 /
      data fc     /     0.011d0,   -13.257d0,    295263.137d0,
     .             335778.877d0,    1342.0d0/
      data dc     /     0.019d0,    -6.891d0,    1105601.328d0,
     .            1072261.307d0,    1236.0d0/
      data omc    /     0.008d0,     7.455d0,    -482890.539d0,
     .             450160.280d0,      -5.0d0/

****  Get the number of centuries to current time

      cent = (epoch-dj2000) / 36525.d0

****  Compute angular arguments
      el = elc(1) * cent**3 + elc(2) * cent**2 + elc(3) * cent
     .          + elc(4) + mod( elc(5) * cent, 1.d0 ) * sec360
      el = mod( el, sec360 )
      eld = 3.d0 * elc(1) * cent**2 + 2.d0 * elc(2) * cent + elc(3)
     .      + elc(5) * sec360
c
      elp = elpc(1) * cent**3 + elpc(2) * cent**2 + elpc(3) * cent
     .     + elpc(4) + mod( elpc(5) * cent, 1.d0 ) * sec360
      elp = mod( elp, sec360 )
      elpd = 3.d0 * elpc(1) * cent**2 + 2.d0 * elpc(2) * cent + elpc(3)
     .       + elpc(5) * sec360
c
      f = fc(1) * cent**3 + fc(2) * cent**2 + fc(3) * cent
     .     + fc(4) + mod( fc(5) * cent, 1.d0 ) * sec360
      f = mod( f, sec360 )
      fd = 3.d0 * fc(1) * cent**2 + 2.d0 * fc(2) * cent + fc(3)
     .     + fc(5) * sec360
c
      d = dc(1) * cent**3 + dc(2) * cent**2 + dc(3) * cent
     .     + dc(4) + mod( dc(5) * cent, 1.d0 ) * sec360
      d = mod( d, sec360 )
      dd = 3.d0 * dc(1) * cent**2 + 2.d0 * dc(2) * cent + dc(3)
     .     + dc(5) * sec360
c
      om = omc(1) * cent**3 + omc(2) * cent**2 + omc(3) * cent
     .     + omc(4) + mod( omc(5) * cent, 1.d0 ) * sec360
      om = mod( om, sec360 )
      omd = 3.d0 * omc(1) * cent**2 + 2.d0 * omc(2) * cent + omc(3)
     .      + omc(5) * sec360
c

***** Now compute GMST.  (CALC 7.1 Algorithm)
*     Remove the fractional part of the julian date
*     Get jd at 0:00 UT
      jd_0hr = aint(epoch-0.5d0) + 0.5d0
*                         ! Days since J2000.0
      t_0hr = jd_0hr - dj2000
*                         ! 0:00 hrs at start of day
      cent = t_0hr / 36525.d0

*                         ! Fraction of a day
      fract = epoch - jd_0hr

      diurnv = ( 1.002737909350795d0 + 5.9006d-11*cent
     .                               - 5.9d-15*cent**2 )
C
C**** COMPUTE GST in cycles
      gstd = ( 24110.54841d0  + 8640184.81266d0*cent
     .                        + 0.093104d0*cent**2
     .                        - 6.2d-6*cent**3 ) /86400.d0

      gstd = mod(gstd,1.d0)
*                                             ! Rads
      gst = (gstd + diurnv*fract) * 2.d0*pi

****  Now save the values.  Convert values from arcseconds to radians

      fund_arg(1) = el / (3600.d0*rad_to_deg)
      fund_arg(2) = elp/ (3600.d0*rad_to_deg)
      fund_arg(3) = f  / (3600.d0*rad_to_deg)
      fund_arg(4) = d  / (3600.d0*rad_to_deg)
      fund_arg(5) = om / (3600.d0*rad_to_deg)
      fund_arg(6) = gst + pi

***** Thats all
      return
      end

CTITLE SDC_ARG

      subroutine sdc_arg( sd_mult, fund_arg, arg, num_arg )

*     This routine computes the argument of the tide given the
*     the fundamental argument values and the multipliers to
*     generate the tide

* PHYSICAL CONSTANTS NEEDED FOR SD_COMP

*   pi          - Define here to full precision

      real*8 pi

      parameter ( pi            = 3.1415926535897932D0 )

*-------------------------------------------------------------------

* PASSED VARIABLES

* INPUT
*   sd_mult(6)  - multipliers for Browns arguments and gst+pi
*   num_arg     - Number of arguments to be summed. (Allows the
*                 gst+pi argument to be skipped.)

      integer*4 sd_mult(6), num_arg

*   fund_arg(6) - Values for the fundamental artguments

* OUPUT
*   arg         - Argumnent at this time (rad)

      real*8 fund_arg(6), arg

* LOCAL VRAIABLES

*   i           - Loop counter

      integer*4 i

****  Initialize argument and combine the fundamental angels with the
*     values passed

      arg = 0.d0
      do i = 1, num_arg
          arg = arg + sd_mult(i)*fund_arg(i)
      end do

      arg = mod(arg, 2.d0*pi)

***** Thats all
      return
      end





