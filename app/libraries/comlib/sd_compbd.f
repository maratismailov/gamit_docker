CTITLE SD_COMP_BD

      block data sd_compbd

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

