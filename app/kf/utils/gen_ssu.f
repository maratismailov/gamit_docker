
      program gen_ssu

      implicit none 
 
*     Program produces a UT1 tidal plot every six hours

      include '../includes/kalman_param.h'
c     include '../includes/const_param.h'
      include '../includes/obs_header.h'
      include '../includes/sd_common.h'

 
*   date_start(5), date_stop(5)     - Start and stop dates
*   date(5)         - date in question
 
      integer*4 date_start(5), date_stop(5), date(5),
     .          rcpar, len_run, indx, ierr
 
*   jd, jd_start, jd_stop           - Julians dates
*   sectag                          - seconds tag
*   dut                             - Tidal contrubution to ut1 (not used)
*   dxy(2)                          - Contributions to x and y pole
*                                     position. (not used)
*   pmu(3) - Test of new routines

*   Values for testing new standard routine.
*   dX, dY, dUT1  - Values returned from sd_comp.
*   cor_x, cor_y, cor_lod - Values from IERS interp routines (arcsec and sec)
 
      real*8 jd, jd_start, jd_stop, sectag, 
     .       step
      integer*4 ijd, njd  ! loop and number of steps

      real*8 dX, dY, dUT1, eop(3)
      real*8 cor_x, cor_y,  cor_lod

*    sd_file  - Name of sd file to use
      character*256 sd_file

*    runstring - Runstring value

      character*256 runstring
 
c      data date_start / 1989, 8,26, 0 , 0 /
c     .    date_stop / 1989,12, 9 , 0 , 0 /
      data date_start / 2006, 11,4, 0 , 0 /
     .    date_stop / 2007,3, 12 , 0 , 0 /
      data step / 0.041666666666667d0 /

*     Decode the runstrinG
      len_run = rcpar(1, sd_file)
      if( len_run.le.0 ) then
          write(*,100)
 100      format(' GEN_SSU: Generate values of SD ut1 terms',/,
     .           ' Runstring:',/,
     .           ' % gen_ssu <sd file> <range>',/,
     .           ' where <range> is given in ymd start, ymd end',
     .           ' and step enclosed in single quotes',/,
     .           ' <sd_file> can be SD_COMP, RAY, INTERP, IERS',/)
           stop 'GEN_SSU: Improper runstring'
      end if
      write(*,'(''* SD UT1 Generated from '',a)') sd_file(1:len_run)
      len_run = rcpar(2,runstring)
      indx = 1
      if( len_run.gt.0 ) then
          call multiread(runstring, indx, 'I4', ierr, date_start,' ',3)
          call multiread(runstring, indx, 'I4', ierr, date_stop ,' ',3)
          call multiread(runstring, indx, 'R8', ierr, step  ,' ',1)
      end if
 
      sectag = 0
      call ymdhms_to_mjd( date_start, sectag, jd_start)
      call ymdhms_to_mjd( date_stop, sectag, jd_stop)

c      call init_sd(sd_file, .true., .true., ' ', 0, 2440000.d0)

      write(*,120) 
 120  format('*     Date          MJD      dx (mas) ',
     .       '  dy (mas)   dut1 (mts) ')

c     do jd = jd_start, jd_stop, step  
      njd = ((jd_stop-jd_start)/step) 
      do ijd = 0, njd
          jd = jd_start + ijd*step   

*         Compute the instaneous values of the pole and ut1 
c          call eval_sd( jd, pmu )
c          pmu(3) = pmu(3)/15.d0

*         Now compute the amplitudes of the coherent components
c          call get_sd_mid ( jd, .true., 0)

*         Convert the ut1 back to mill time seconds
c          do i = 1,4
c             sd_ut1_mid(i) = sd_ut1_mid(i)/15.d0
c          end do

*****     Test new standard routine:
          if( sd_file.eq.'SD_COMP' ) then 
             call sd_comp( jd , dX, dY, dUT1 )
          elseif( sd_file.eq.'RAY' ) then
             call ray(jd, dx, dy, dut1 )
          elseif( sd_file.eq.'INTERP') then
             call PMUT1_OCEANS (jd,dX, dY, dUT1,cor_lod)             
*            Solid Earth terms          
             call PM_GRAVI (jd,cor_x,cor_y)
             dX = (dX + cor_x)*1000
             dY = (dY + cor_y)*1000
             dUT1 = dUT1*1000
         else
             call ortho_eop (jd, eop)
             dX = eop(1)/1000.d0
             dY = eop(2)/1000.d0
             dUT1 = eop(3)/1000.d0
          endif


          call jd_to_ymdhms( jd, date, sectag)
c         write(*,200) date, pmu, sd_ut1_mid, sd_xy_mid
          write(*,200) date, jd, dX, dY, dUT1
 200  format(i5,4i3,1x,F11.4,2x, 3f10.4,4F10.4, 6F10.3)
 
      end do
      end

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
 
      integer*4 i

*   epoch       - Julian date (jd passed in unless the JD 
*                 appears to be an MJD in which case it is 
*                 converted to JD (2 400 000.5d0 added) 
*   arg         - Angular argument for the correction (rads)
*   fund_arg(6) - Values of the 5 Brown's arguments (l,l',F,D,
*               - Omega) and gst+pi (rads)
 
      real*8 epoch, arg, fund_arg(6)

***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
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
      gstd = ( 24110.54841d0  + 8640184.812866d0*cent
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

      implicit none 
 
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
 
CTITLE SD_COMP_BD

      block data sd_comp_bd

      implicit none 

*     This is a block data routine to define the values of
*     the arguments and coefficients of the diurnal and
*     semidiurnal variations in UT1 and pole position.

*     The coefficients in these data statements are those from
*     the GD_920817 Series solutions and are based on the analysis
*     of 1085 VLBI experiments bewteen Jan 1984 and June 1992.
*     A discussion  of these results appears in:
*     Herring, T.A., Diurnal and Semidiurnal Earth rotation 
*        variations, Advances in Space Research, Pergamon Press,
*        New York, in press, 1992.


* DEFINE THE PARAMETERS OF THE SYSTEM
 
*   num_sdc_ut1      - number of UT1 entries in data statment.
*   num_sdc_xy       - number of xy pole entries in data statement.
 
      integer*4 num_sdc_ut1, num_sdc_xy
 
      parameter ( num_sdc_ut1   = 10 )
      parameter ( num_sdc_xy    = 12 )

 
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

*     Fundumental arguments for UT1 variations.  
*     NOTE: These arguments must be kept in the same order as
*           the values of the coefficients below.
*                         l   l'  F   D  Om  GST+pi
      data sdc_ut1_arg /  0,  0,  0,  0,  0, -1, 
     .                    0,  0,  2, -2,  2, -1, 
     .                    0,  0,  2,  0,  2, -1, 
     .                    1,  0,  2,  0,  2, -1, 
     .                    0,  1,  0,  0,  0, -1, 
     .                    1,  0,  0,  0,  0, -1,
     .                    0,  0,  0,  0,  0, -2, 
     .                    0,  0,  2, -2,  2, -2,
     .                    0,  0,  2,  0,  2, -2,
     .                    1,  0,  2,  0,  2, -2 /
*
*     Values of the coefficients for UT1. NOTE: Values
*     are in milliarcseconds (sd_comp converts to millitimesecs)
*                         Cos         Sin
      data sdc_ut1_val / 0.0976d0,  0.2676d0,
     .                  -0.0584d0, -0.0886d0,
     .                  -0.2597d0, -0.2422d0,
     .                  -0.0468d0, -0.0653d0,
     .                   0.0263d0,  0.0189d0,
     .                   0.0276d0,  0.0346d0, 
     .                   0.0113d0,  0.0548d0,
     .                  -0.0014d0,  0.1291d0,
     .                  -0.1622d0,  0.2140d0,
     .                  -0.0242d0,  0.0426d0 /

*     Fundumental arguments for Polar motion variations.  These
*     are expressed as the circular components  (see sd_comp for
*     sign convention to map back to x and y pole postions.)
*     NOTE: These arguments must be kept in the same order as
*           the values of the coefficients below.
*                         l   l'  F   D  Om  GST+pi
      data sdc_xy_arg  /  0,  0,  0,  0,  0,  1, 
     .                    0,  0, -2,  2, -2,  1, 
     .                    0,  0, -2,  0, -2,  1, 
     .                   -1,  0, -2,  0, -2,  1, 
     .                    0,  0,  0,  0,  0, -2, 
     .                    0,  0,  2, -2,  2, -2, 
     .                    0,  0,  2,  0,  2, -2, 
     .                    1,  0,  2,  0,  2, -2, 
     .                    0,  0,  0,  0,  0,  2, 
     .                    0,  0, -2,  2, -2,  2, 
     .                    0,  0, -2,  0, -2,  2, 
     .                   -1,  0, -2,  0, -2,  2  /

*
*     Values of the coefficients for polar motion circular
*     components. NOTE: Values are in milliarcseconds
*                         Cos         Sin
      data sdc_xy_val  /  0.1330d0, -0.0738d0,
     .                   -0.0487d0,  0.0349d0,
     .                   -0.1778d0,  0.0893d0,
     .                   -0.0328d0,  0.0107d0,
     .                   -0.0263d0,  0.0163d0,
     .                   -0.0665d0,  0.0986d0,
     .                   -0.0103d0,  0.2654d0,
     .                   -0.0095d0,  0.0468d0,
     .                    0.0389d0, -0.0045d0,
     .                    0.0023d0, -0.0120d0,
     .                    0.0014d0, -0.0577d0,
     .                    0.0118d0, -0.0120d0  /

      end
      SUBROUTINE RAY (RJD,CORX,CORY,CORT)

      implicit none 

c   This subroutine implements the Ray model for diurnal/subdirunal
c   tides.  It uses the Simon et al. Fundamental Arguments.  The
c   corrections in x and y are in units of millisec of arc and UT1-UTC
c   in millisec of time  (NOTE CHANGE FROM ORIGINAL, to match units of
c   sd_comp.f--rwk ).  These corrections should be added to "average"
c   EOP values to get estimates of the instantaneous values.

* MOD TAH 990311: Cleaned up code by adding explicit declarations
*     and allowing either MJD (as in original or JD to be passed).  
* MOD RWK 990319: Changed output units from arcsec,sec to mas, ms,
*     to match old routine (sd_comb.f)

C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DOUBLE PRECISION
C    .   L,        LPRIME

* INPUT VARIABLES
* RJD  -- Either MJD or JD for evaluating contributions

      real*8 rjd

* OUTPUT VARIABLES
* corx, cory -- Corrections to X and Y pole (arc secs)
* cort       -- Correction to UT1 (time secs)

      real*8  corx, cory, cort

* LOCAL VARIABLES
* t   -- Julian centuries from J2000
* halfpi -- pi/2
* l, lprime, capf, capd, omega, theta  -- Fundamental arguments

      real*8 t, halfpi
      real*8 l, lprime, capf, capd, omega, theta 

* arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8 -- Argument
*     values for the main diurnal and semidiurnal tides
      real*8 arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8


****  START:
      HALFPI = 1.5707963267948966d0

*     Get julian centuries from J2000
      if( rjd.lt. 2 000 000.d0 ) then
          T = (RJD - 51544.5D0)/36525.0D0
      else
          T = (RJD - 2451545.d0)/36525.d0
      end if

*     Compute fundamental arguments
      L = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      L = DMOD(L,1296000d0)
      LPRIME = -0.00001149d0*T**4 + 0.000136d0*T**3 - 0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      LPRIME = DMOD(LPRIME,1296000d0)
      CAPF = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      CAPF = DMOD(CAPF,1296000d0)
      CAPD = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      CAPD = DMOD(CAPD,1296000d0)
      OMEGA = -0.00005939d0*T**4 + 0.007702d0*T**3 + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      OMEGA = DMOD(OMEGA,1296000d0)
      THETA = (67310.54841d0 + 
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0

*     Compute the tidal arguments
      ARG7 = DMOD((-L - 2.0D0*CAPF - 2.0D0*OMEGA + THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG1 = DMOD((-2.0d0*CAPF - 2.0d0*OMEGA + THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG2 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG3 = DMOD(THETA * 3.14159265D0/648000.0D0,6.28318530718D0)
     .     + HALFPI
      ARG4 = DMOD((-L - 2.0d0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG5 = DMOD((-2.0D0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG6 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + 2.0d0*THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG8 = DMOD((2.0d0*THETA)
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)

*     Compute the corrections.
      CORX = - 0.026D0*DSIN(ARG7) + 0.006D0*DCOS(ARG7)
     .       - 0.133D0*DSIN(ARG1) + 0.049D0*DCOS(ARG1)
     .       - 0.050D0*DSIN(ARG2) + 0.025D0*DCOS(ARG2)
     .       - 0.152D0*DSIN(ARG3) + 0.078D0*DCOS(ARG3)
     .       - 0.057D0*DSIN(ARG4) - 0.013D0*DCOS(ARG4)
     .       - 0.330D0*DSIN(ARG5) - 0.028D0*DCOS(ARG5)
     .       - 0.145D0*DSIN(ARG6) + 0.064D0*DCOS(ARG6)
     .       - 0.036D0*DSIN(ARG8) + 0.017D0*DCOS(ARG8)
      CORY = - 0.006D0*DSIN(ARG7) - 0.026D0*DCOS(ARG7)
     .       - 0.049D0*DSIN(ARG1) - 0.133D0*DCOS(ARG1)
     .       - 0.025D0*DSIN(ARG2) - 0.050D0*DCOS(ARG2)
     .       - 0.078D0*DSIN(ARG3) - 0.152D0*DCOS(ARG3)
     .       + 0.011D0*DSIN(ARG4) + 0.033D0*DCOS(ARG4)
     .       + 0.037D0*DSIN(ARG5) + 0.196D0*DCOS(ARG5)
     .       + 0.059D0*DSIN(ARG6) + 0.087D0*DCOS(ARG6)
     .       + 0.018D0*DSIN(ARG8) + 0.022D0*DCOS(ARG8)
      CORT = + 0.0245D0*DSIN(ARG7) + 0.0503D0*DCOS(ARG7)
     .       + 0.1210D0*DSIN(ARG1) + 0.1605D0*DCOS(ARG1)
     .       + 0.0286D0*DSIN(ARG2) + 0.0516D0*DCOS(ARG2)
     .       + 0.0864D0*DSIN(ARG3) + 0.1771D0*DCOS(ARG3)
     .       - 0.0380D0*DSIN(ARG4) - 0.0154D0*DCOS(ARG4)
     .       - 0.1617D0*DSIN(ARG5) - 0.0720D0*DCOS(ARG5)
     .       - 0.0759D0*DSIN(ARG6) - 0.0004D0*DCOS(ARG6)
     .       - 0.0196D0*DSIN(ARG8) - 0.0038D0*DCOS(ARG8)

*     Convert to output units (originally arcsec and seconds of time)
*     now millarcsec and msec of time -- rwk 990319
*      CORX = CORX * 1.0d-3
*      CORY = CORY * 1.0d-3
*      CORT = CORT * 0.1d-3
       cort = cort * 0.1d0
      RETURN
      END

cxxxxxx
c   subroutines to compute the diurnal and semidiurnal variations
c   in the earth orientation from the version of Richard Ray's ocean
c   tide model that was listed in IERS Technical Note 21, July 1996.
c   This code includes the variations from 71 diurnal and semidiurnal
c   terms instead of the 8 that are listed in the report.
cxxxxxx

      subroutine ortho_eop (time, eop)

      implicit none 
c
c...Purpose: to compute the diurnal and semidiurnal variations
c            in EOP (x,y,UT1) from ocean tides
c
c...Coded by: Richard Eanes, UT/CSR, Feb 1997
c
c...Input: time = Modified Julian Date
c
c...Output: eop = (delta_x, delta_y, delta_UT1)
c                 microarcsec for x and y, microsec for UT1
c
      real*8 eop(3), orthow(12,3), h(12), time
      integer*4 k,j 
c
c...diurnal and semidiurnal orthoweights fit to the 8 constituents
c   listed in IERS Technical Note 21, July 1996 which are from the
c   paper "Diurnal and Semidiurnal Variations in the Earth's Rotation
c   Rate Induced by Ocean Tides"  by Ray,R.D., Steinberg,D.J.,
c   Chao,B.F., and Cartwright, D.E., Science, 264, pp. 830-832.
c
      data orthow /
     . -6.77832,-14.86323,  0.47884, -1.45303,  0.16406,  0.42030,
     .  0.09398, 25.73054, -4.77974,  0.28080,  1.94539, -0.73089,
     . 14.86283, -6.77846,  1.45234,  0.47888, -0.42056,  0.16469,
     . 15.30276, -4.30615,  0.07564,  2.28321, -0.45717, -1.62010,
     . -1.76335,  1.03364, -0.27553,  0.34569, -0.12343, -0.10146,
     . -0.47119,  1.28997, -0.19336,  0.02724,  0.08955,  0.04726/
c
c...compute the partials of the tidal variations to the orthoweights
      call cnmtx (time, h)
c
c...compute eop changes
      do k=1,3
         eop(k) = 0.
         do j=1,12
            eop(k) = eop(k) + h(j)*orthow(j,k)
         enddo
      enddo
c
      return
      end



      subroutine cnmtx (dmjd, h)

      implicit none 
c
c...Purpose: To compute the time dependent part of the second degree
c            diurnal and semidiurnal tidal potential from the dominant
c            spectral lines in the Cartwright-Tayler-Edden harmonic
c            decomposition
c
c...Coded by: Richard Eanes, UT/CSR, Feb 1997
c
c...Input: dmjd = modified julian date
c
c...Output: h = vector of length 12 with partials of the tidal
c           variation with respect to the orthoweights
c
      integer*4 nlines
      parameter (nlines=71)
      character*7 numarg(nlines)
      real*8 h(12)
      integer*4 nj(nlines),mj(nlines), m, j,n 
      integer*4 nmax, i , k
      real*8 hs(nlines),phase(nlines),freq(nlines)
      real*8 anm(2:3,0:3,-1:1), bnm(2:3,0:3,-1:1)
      real*8 p(0:2,2),q(0:2,2),sp(6,2)
      real*8 twopi, dt, dt60, d1960
      real*8 dmjd
      real*8 pinm, alpha, ap,am,bp,bm
c
c...the orthotide weight factors
      data ((sp(i,m),i=1,6),m=1,2) /
     . 0.0298,  0.1408, +0.0805,  0.6002, +0.3025,  0.1517,
     . 0.0200,  0.0905, +0.0638,  0.3476, +0.1645,  0.0923/
c
      data twopi /6.2831853071796d0/
      data dt / 2.d0 /
      data nmax /2/
c
c...tidal potential model for 71 diurnal and semidiurnal lines
c
      data d1960/37076.5d0/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=1,15)
     ./2, 1,  -1.94, 9.0899831, 5.18688050, '117.655',
     . 2, 1,  -1.25, 8.8234208, 5.38346657, '125.745',
     . 2, 1,  -6.64,12.1189598, 5.38439079, '125.755',
     . 2, 1,  -1.51, 1.4425700, 5.41398343, '127.545',
     . 2, 1,  -8.02, 4.7381090, 5.41490765, '127.555',
     . 2, 1,  -9.47, 4.4715466, 5.61149372, '135.645',
     . 2, 1, -50.20, 7.7670857, 5.61241794, '135.655',
     . 2, 1,  -1.80,-2.9093042, 5.64201057, '137.445',
     . 2, 1,  -9.54, 0.3862349, 5.64293479, '137.455',
     . 2, 1,   1.52,-3.1758666, 5.83859664, '145.535',
     . 2, 1, -49.45, 0.1196725, 5.83952086, '145.545',
     . 2, 1,-262.21, 3.4152116, 5.84044508, '145.555',
     . 2, 1,   1.70,12.8946194, 5.84433381, '145.755',
     . 2, 1,   3.43, 5.5137686, 5.87485066, '147.555',
     . 2, 1,   1.94, 6.4441883, 6.03795537, '153.655'/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=16,30)
     ./2, 1,   1.37,-4.2322016, 6.06754801, '155.445',
     . 2, 1,   7.41,-0.9366625, 6.06847223, '155.455',
     . 2, 1,  20.62, 8.5427453, 6.07236095, '155.655',
     . 2, 1,   4.14,11.8382843, 6.07328517, '155.665',
     . 2, 1,   3.94, 1.1618945, 6.10287781, '157.455',
     . 2, 1,  -7.14, 5.9693878, 6.24878055, '162.556',
     . 2, 1,   1.37,-1.2032249, 6.26505830, '163.545',
     . 2, 1,-122.03, 2.0923141, 6.26598252, '163.555',
     . 2, 1,   1.02,-1.7847596, 6.28318449, '164.554',
     . 2, 1,   2.89, 8.0679449, 6.28318613, '164.556',
     . 2, 1,  -7.30, 0.8953321, 6.29946388, '165.545',
     . 2, 1, 368.78, 4.1908712, 6.30038810, '165.555',
     . 2, 1,  50.01, 7.4864102, 6.30131232, '165.565',
     . 2, 1,  -1.08,10.7819493, 6.30223654, '165.575',
     . 2, 1,   2.93, 0.3137975, 6.31759007, '166.554'/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=31,45)
     ./2, 1,   5.25, 6.2894282, 6.33479368, '167.555',
     . 2, 1,   3.95, 7.2198478, 6.49789839, '173.655',
     . 2, 1,  20.62,-0.1610030, 6.52841524, '175.455',
     . 2, 1,   4.09, 3.1345361, 6.52933946, '175.465',
     . 2, 1,   3.42, 2.8679737, 6.72592553, '183.555',
     . 2, 1,   1.69,-4.5128771, 6.75644239, '185.355',
     . 2, 1,  11.29, 4.9665307, 6.76033111, '185.555',
     . 2, 1,   7.23, 8.2620698, 6.76125533, '185.565',
     . 2, 1,   1.51,11.5576089, 6.76217955, '185.575',
     . 2, 1,   2.16, 0.6146566, 6.98835826, '195.455',
     . 2, 1,   1.38, 3.9101957, 6.98928248, '195.465',
     . 2, 2,   1.80,20.6617051,11.45675174, '225.855',
     . 2, 2,   4.67,13.2808543,11.48726860, '227.655',
     . 2, 2,  16.01,16.3098310,11.68477889, '235.755',
     . 2, 2,  19.32, 8.9289802,11.71529575, '237.555'/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=46,60)
     ./2, 2,   1.30, 5.0519065,11.73249771, '238.554',
     . 2, 2,  -1.02,15.8350306,11.89560406, '244.656',
     . 2, 2,  -4.51, 8.6624178,11.91188181, '245.645',
     . 2, 2, 120.99,11.9579569,11.91280603, '245.655',
     . 2, 2,   1.13, 8.0808832,11.93000800, '246.654',
     . 2, 2,  22.98, 4.5771061,11.94332289, '247.455',
     . 2, 2,   1.06, 0.7000324,11.96052486, '248.454',
     . 2, 2,  -1.90,14.9869335,12.11031632, '253.755',
     . 2, 2,  -2.18,11.4831564,12.12363121, '254.556',
     . 2, 2, -23.58, 4.3105437,12.13990896, '255.545',
     . 2, 2, 631.92, 7.6060827,12.14083318, '255.555',
     . 2, 2,   1.92, 3.7290090,12.15803515, '256.554',
     . 2, 2,  -4.66,10.6350594,12.33834347, '263.655',
     . 2, 2, -17.86, 3.2542086,12.36886033, '265.455',
     . 2, 2,   4.47,12.7336164,12.37274905, '265.655'/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=61,71)
     ./2, 2,   1.97,16.0291555,12.37367327, '265.665',
     . 2, 2,  17.20,10.1602590,12.54916865, '272.556',
     . 2, 2, 294.00, 6.2831853,12.56637061, '273.555',
     . 2, 2,  -2.46, 2.4061116,12.58357258, '274.554',
     . 2, 2,  -1.02, 5.0862033,12.59985198, '275.545',
     . 2, 2,  79.96, 8.3817423,12.60077620, '275.555',
     . 2, 2,  23.83,11.6772814,12.60170041, '275.565',
     . 2, 2,   2.59,14.9728205,12.60262463, '275.575',
     . 2, 2,   4.47, 4.0298682,12.82880334, '285.455',
     . 2, 2,   1.95, 7.3254073,12.82972756, '285.465',
     . 2, 2,   1.17, 9.1574019,13.06071921, '295.555'/
c
c...compute the time dependent potential matrix
c
      do k=-1,1
         dt60 = (dmjd - k*dt) - d1960
c         anm(2,1:2,k) = 0.
c         bnm(2,1:2,k) = 0.
         anm(2,1,k) = 0.
         bnm(2,1,k) = 0.
         anm(2,2,k) = 0.
         bnm(2,2,k) = 0.
         do j=1,nlines
            n = nj(j)
            m = mj(j)
            pinm = mod(n+m,2)*twopi/4.
            alpha = phase(j) + freq(j)*dt60 - pinm
            alpha = mod(alpha,twopi)
            anm(n,m,k) = anm(n,m,k) + hs(j)*cos(alpha)
            bnm(n,m,k) = bnm(n,m,k) - hs(j)*sin(alpha)
         enddo
      enddo
c
c...orthogonalize the response terms
c
      do m = 1,2
        ap = anm(2,m,1) + anm(2,m,-1)
        am = anm(2,m,1) - anm(2,m,-1)
        bp = bnm(2,m,1) + bnm(2,m,-1)
        bm = bnm(2,m,1) - bnm(2,m,-1)
        p(0,m) = sp(1,m)*anm(2,m,0)
        p(1,m) = sp(2,m)*anm(2,m,0) - sp(3,m)*ap
        p(2,m) = sp(4,m)*anm(2,m,0) - sp(5,m)*ap + sp(6,m)*bm
        q(0,m) = sp(1,m)*bnm(2,m,0)
        q(1,m) = sp(2,m)*bnm(2,m,0) - sp(3,m)*bp
        q(2,m) = sp(4,m)*bnm(2,m,0) - sp(5,m)*bp - sp(6,m)*am
        anm(2,m,-1) = p(0,m)
        anm(2,m, 0) = p(1,m)
        anm(2,m, 1) = p(2,m)
        bnm(2,m,-1) = q(0,m)
        bnm(2,m, 0) = q(1,m)
        bnm(2,m, 1) = q(2,m)
      enddo
c
c...fill partials vector
      j = 1
      do n=2,nmax
         do m = 1,n
            do k = -1,1
               h(j)   = anm(n,m,k)
               h(j+1) = bnm(n,m,k)
               j = j + 2
             enddo
         enddo
      enddo
c
      return
      end 
C----------------------------------------------------------------
C Routines from ftp://hpiers.obspm.fr/eop-pc/models/interp.f
C
C----------------------------------------------------------------
      SUBROUTINE PMUT1_OCEANS (rjd,cor_x,cor_y,cor_ut1,cor_lod)

      implicit none 
C
C    This subroutine provides, in time domain, the diurnal/subdiurnal
C    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
C    listed in the program above, have been extracted from the procedure   
C    ortho_eop.f coed by Eanes in 1997.
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
C     cor_lod  - tidal correction in length of day (sec. of time)
C
C     coded by Ch. Bizouard (2002), initially coded by McCarthy and 
C     D.Gambis(1997) for the 8 prominent tidal waves.  
      
      
      INTEGER nlines
      PARAMETER(nlines=71)
      DOUBLE PRECISION ARG(6),    ! Array of the tidal arguments   
     .                 DARG(6)    ! Array of their time derivative 
      
      REAL*4 XCOS(nlines),XSIN(nlines),
     .YCOS(nlines),YSIN(nlines),UTCOS(nlines),UTSIN(nlines)
     
      REAL*8 t,ag,dag,rjd,halfpi,secrad,
     .       cor_x,cor_y,cor_ut1,cor_lod
      INTEGER NARG(nlines,6),i,j
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),UTSIN(j),UTCOS(j),j=1,nlines)/
     &1,-1, 0,-2,-2,-2,  -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078,
     &1,-2, 0,-2, 0,-1,   0.06,   0.64,  -0.64,   0.06,  0.195, -0.059,
     &1,-2, 0,-2, 0,-2,   0.30,   3.42,  -3.42,   0.30,  1.034, -0.314,
     &1, 0, 0,-2,-2,-1,   0.08,   0.78,  -0.78,   0.08,  0.224, -0.073,
     &1, 0, 0,-2,-2,-2,   0.46,   4.15,  -4.15,   0.45,  1.187, -0.387,
     &1,-1, 0,-2, 0,-1,   1.19,   4.96,  -4.96,   1.19,  0.966, -0.474,
     &1,-1, 0,-2, 0,-2,   6.24,  26.31, -26.31,   6.23,  5.118, -2.499,
     &1, 1, 0,-2,-2,-1,   0.24,   0.94,  -0.94,   0.24,  0.172, -0.090,
     &1, 1, 0,-2,-2,-2,   1.28,   4.99,  -4.99,   1.28,  0.911, -0.475,
     &1, 0, 0,-2, 0, 0,  -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070,
     &1, 0, 0,-2, 0,-1,   9.22,  25.06, -25.06,   9.22,  3.025, -2.280,
     &1, 0, 0,-2, 0,-2,  48.82, 132.91,-132.90,  48.82, 16.020,-12.069,
     &1,-2, 0, 0, 0, 0,  -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078,
     &1, 0, 0, 0,-2, 0,  -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154,
     &1,-1, 0,-2, 2,-2,  -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074,
     &1, 1, 0,-2, 0,-1,  -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050,
     &1, 1, 0,-2, 0,-2,  -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271,
     &1,-1, 0, 0, 0, 0,  -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751,
     &1,-1, 0, 0, 0,-1,  -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151,
     &1, 1, 0, 0,-2, 0,  -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137,
     &1, 0,-1,-2, 2,-2,   1.54,   3.03,  -3.03,   1.54,  0.315, -0.189,
     &1, 0, 0,-2, 2,-1,  -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035,
     &1, 0, 0,-2, 2,-2,  26.13,  51.25, -51.25,  26.13,  5.512, -3.095,
     &1, 0, 1,-2, 2,-2,  -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025,
     &1, 0,-1, 0, 0, 0,  -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070,
     &1, 0, 0, 0, 0, 1,   1.54,   3.00,  -3.00,   1.54,  0.348, -0.171,
     &1, 0, 0, 0, 0, 0, -77.48,-151.74, 151.74, -77.48,-17.620,  8.548,
     &1, 0, 0, 0, 0,-1, -10.52, -20.56,  20.56, -10.52, -2.392,  1.159,
     &1, 0, 0, 0, 0,-2,   0.23,   0.44,  -0.44,   0.23,  0.052, -0.025,
     &1, 0, 1, 0, 0, 0,  -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065,
     &1, 0, 0, 2,-2, 2,  -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111,
     &1,-1, 0, 0, 2, 0,  -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043,
     &1, 1, 0, 0, 0, 0,  -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187,
     &1, 1, 0, 0, 0,-1,  -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037,
     &1, 0, 0, 0, 2, 0,  -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005,
     &1, 2, 0, 0, 0, 0,  -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005,
     &1, 0, 0, 2, 0, 2,  -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037,
     &1, 0, 0, 2, 0, 1,  -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023,
     &1, 0, 0, 2, 0, 0,  -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005,
     &1, 1, 0, 2, 0, 2,  -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024,
     &1, 1, 0, 2, 0, 1,  -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015,
     &2,-3, 0,-2, 0,-2,  -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011,
     &2,-1, 0,-2,-2,-2,  -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032,
     &2,-2, 0,-2, 0,-2,  -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177,
     &2, 0, 0,-2,-2,-2,  -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222,
     &2, 0, 1,-2,-2,-2,  -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015,
     &2,-1,-1,-2, 0,-2,   0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013,
     &2,-1, 0,-2, 0,-1,   2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058,
     &2,-1, 0,-2, 0,-2, -56.87, -12.93,  11.15,  32.88, -3.795, -1.556,
     &2,-1, 1,-2, 0,-2,  -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015,
     &2, 1, 0,-2,-2,-2, -11.01,  -2.40,   1.89,   6.41, -0.698, -0.298,
     &2, 1, 1,-2,-2,-2,  -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014,
     &2,-2, 0,-2, 2,-2,   0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022,
     &2, 0,-1,-2, 0,-2,   1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025,
     &2, 0, 0,-2, 0,-1,  12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266,
     &2, 0, 0,-2, 0,-2,-330.15, -26.96,  37.58, 195.92,-16.195, -7.140,
     &2, 0, 1,-2, 0,-2,  -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021,
     &2,-1, 0,-2, 2,-2,   2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034,
     &2, 1, 0,-2, 0,-2,   9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117,
     &2,-1, 0, 0, 0, 0,  -2.35,   0.37,   0.47,   1.41, -0.106, -0.029,
     &2,-1, 0, 0, 0,-1,  -1.04,   0.17,   0.21,   0.62, -0.047, -0.013,
     &2, 0,-1,-2, 2,-2,  -8.51,   3.50,   3.29,   5.11, -0.437, -0.019,
     &2, 0, 0,-2, 2,-2,-144.13,  63.56,  59.23,  86.56, -7.547, -0.159,
     &2, 0, 1,-2, 2,-2,   1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000,
     &2, 0, 0, 0, 0, 1,   0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001,
     &2, 0, 0, 0, 0, 0, -38.48,  19.14,  17.72,  23.11, -2.104,  0.041,
     &2, 0, 0, 0, 0,-1, -11.44,   5.75,   5.32,   6.87, -0.627,  0.015,
     &2, 0, 0, 0, 0,-2,  -1.24,   0.63,   0.58,   0.75, -0.068,  0.002,
     &2, 1, 0, 0, 0, 0,  -1.77,   1.79,   1.71,   1.04, -0.146,  0.037,
     &2, 1, 0, 0, 0,-1,  -0.77,   0.78,   0.75,   0.45, -0.064,  0.017,
     &2, 0, 0, 2, 0, 2,  -0.33,   0.62,   0.65,   0.19, -0.049,  0.018/

      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   
      DARG(1) = (876600d0*3600d0 + 8640184.812866d0 
     .         + 2.d0 * 0.093104d0 * T - 3.d0 * 6.2d-6*T**2)*15.d0
      DARG(1) = DARG(1)* secrad / 36525.0D0   ! rad/day


      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      
      DARG(2) = -4.d0*0.00024470d0*T**3 + 3.d0*0.051635d0*T**2 
     .  + 2.d0*31.8792d0*T + 1717915923.2178d0 
      DARG(2) = DARG(2)* secrad / 36525.0D0   ! rad/day

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

      DARG(3) = -4.D0*0.00001149d0*T**3 - 3.d0*0.000136d0*T**2
     .  -  2.D0*0.5532d0*T + 129596581.0481d0
      DARG(3) = DARG(3)* secrad / 36525.0D0   ! rad/day
          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

      DARG(4) = 4.d0*0.00000417d0*T**3 - 3.d0*0.001037d0*T**2 
     .- 2.d0 * 12.7512d0*T + 1739527262.8478d0 
      DARG(4) = DARG(4)* secrad / 36525.0D0   ! rad/day
    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

      DARG(5) = -4.d0*0.00003169d0*T**3 + 3.d0*0.006593d0*T**2
     . - 2.d0 * 6.3706d0*T + 1602961601.2090d0
      DARG(5) = DARG(5)* secrad / 36525.0D0   ! rad/day

      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad

      DARG(6) = -4.d0*0.00005939d0*T**3 + 3.d0 * 0.007702d0*T**2
     .  + 2.d0 * 7.4722d0*T - 6962890.2665d0
      DARG(6) = DARG(6)* secrad / 36525.0D0   ! rad/day

C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0
	cor_ut1= 0.d0
	cor_lod= 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
 	dag = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
 		dag = dag + dble(narg(j,i))*DARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x  = cor_x   + dble(XCOS(j)) *dcos(ag) + 
     .                     dble(XSIN(j)) * dsin(ag)
        cor_y  = cor_y   + dble(YCOS(j)) *dcos(ag) + 
     .                     dble(YSIN(j)) * dsin(ag)
        cor_ut1= cor_ut1 + dble(UTCOS(j))*dcos(ag) + 
     .                     dble(UTSIN(j))* dsin(ag)
        cor_lod= cor_lod -(-dble(UTCOS(j)) * dsin(ag) 
     &                    + dble(UTSIN(j)) * dcos(ag) ) * dag   	 

        enddo
  
       cor_x   = cor_x * 1.0d-6   ! arcsecond (")
       cor_y   = cor_y * 1.0d-6   ! arcsecond (")
       cor_ut1 = cor_ut1 * 1.0d-6 ! second (s)
       cor_lod = cor_lod * 1.0d-6 ! second (s)
 
      RETURN
      END
      	
C----------------------------------------------------------------
      SUBROUTINE PM_GRAVI (rjd,cor_x,cor_y)
C
C    This subroutine provides, in time domain, the diurnal
C    lunisolar effet on polar motion (")
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C
C     coded by Ch. Bizouard (2002)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER nlines
      PARAMETER(nlines=10)
      DOUBLE PRECISION ARG(6)    ! Array of the tidal arguments   
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines)
      INTEGER NARG(nlines,6)
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),j=1,nlines)/    
     & 1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44,
     & 1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31,
     & 1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44,
     & 1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14,
     & 1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36,
     & 1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84,
     & 1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76,
     & 1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27,
     & 1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93,
     & 1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76/
 
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   

      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

  
      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad


C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x+dble(XCOS(j))*dcos(ag)+dble(XSIN(j))*dsin(ag)
        cor_y =cor_y+dble(YCOS(j))*dcos(ag)+dble(YSIN(j))*dsin(ag) 

        enddo
  
      cor_x = cor_x * 1.0d-6   ! arcsecond (")
      cor_y = cor_y * 1.0d-6   ! arcsecond (")
 
      RETURN

      END
