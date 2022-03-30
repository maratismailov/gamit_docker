CTITLE NUT_COEFF_PART
 
      subroutine nut_coeff_part( deriv, epoch, fcn_period, coeff_deriv )

      implicit none 
 
*     This rotuine will compute the nutation series coefficients
*     for all terms allowed in the GLOBK program.
*     The odd numbered partials are multiplied by deriv(1), (PSI)
*     and the even numbered ones by deriv(2) (EPS)
*
 
*   i,j         - Loop counters
*   x1(6,15)    - Coefficent mulipliers for Nutation arguments
 
      integer*4 i, x1(6,15)
 
*   arg         - Argument of the nutation term
*   cent        - Number of centuries since J2000.0
*   coeff_deriv(2,30)    - Series coefficients for in and
*               - out of phase
*   convds      - Converts cycles to radians
*   deriv(2)    - The PSI ane EPS partial multipliers
*   DJ2000      - JD of epoch J2000.0
*   epoch       - Epoch of this experiment
*   elc(5), elpc(5), fc(5), dc(5), omc(5) - Coefficients for
*               - the lunar and solar arguments
*   el, elp, f, d, om   - Fundamental arguments
*   FCN_arg     - Argument for FCN (from J2000.0)
*   FCN_period  - Fcn period (days)
*   sec360      - Number of seconds in 360 degs
 
      real*8 arg, cent, coeff_deriv(2,30), convds, deriv(2), DJ2000,
     .    epoch, elc(5), elpc(5), fc(5), dc(5), omc(5),
     .    el, elp, f, d, om, FCN_arg, FCN_period, sec360
 
C
C      CONTSTANTS ARE BASED ON VALUES GIVEN IN NUTW IN CALC V5.0
C      pp182-183, with addition of space for FCN.
C
C
C                 MULTIPLE OF
C                L    L'   F    D  OMEGA FCN
      DATA X1 /  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .           0 ,  0 ,  0 ,  0 ,  0 ,  1 ,
     .           0 ,  0 ,  0 ,  0 ,  1 ,  0 ,
     .           0 ,  0 ,  2 , -2 ,  2 ,  0 ,
     .           0 ,  0 ,  2 ,  0 ,  2 ,  0 ,
     .           0 ,  0 ,  0 ,  0 ,  2 ,  0 ,
     .           0 ,  1 ,  0 ,  0 ,  0 ,  0 ,
     .           1 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     .           0 ,  1 ,  2 , -2 ,  2 ,  0 ,
     .           0 ,  0 ,  2 ,  0 ,  1 ,  0 ,
     .           1 ,  0 ,  2 ,  0 ,  2 ,  0 ,
     .           0 , -1 ,  2 , -2 ,  2 ,  0 ,
     .           1 ,  0 ,  0 , -2 ,  0 ,  0 ,
     .           0 ,  0 ,  2 , -2 ,  1 ,  0 ,
     .          -1 ,  0 ,  2 ,  0 ,  2 ,  0  /
 
C**** ARGUMENTS FOR LUNAR AND SOLAR ANGLES
      DATA ELC / 0.064D0, 31.31D0, 715922.633D0, 485866.733D0, 1325.D0/
      DATA ELPC/-0.012D0,-0.577D0,1292581.224D0,1287099.804D0,   99.D0/
      DATA FC  / 0.011D0,-13.257D0,295263.137D0, 335778.877D0, 1342.D0/
      DATA DC  / 0.019D0,-6.891D0,1105601.328D0,1072261.307D0, 1236.D0/
      DATA OMC / 0.008D0, 7.455D0,-482890.539D0, 450160.280D0,   -5.D0/
C
C**** SOME CONSTANTS WHICH NEED VALUES
      DATA SEC360 / 1296000.D0 /,
     .     DJ2000 / 2451545.D0 /,
     .     CONVDS / 4.84813681109536D-6 /
 
 
*     Get the number of Julian centuries since J2000.0
 
      cent = (epoch - dj2000)/36525.d0
 
C**** COMPUTE FUNDAMENTAL ARGUMENTS
      EL = ELC(1)*CENT**3  + ELC(2)*CENT**2  + ELC(3)*CENT
     .   + ELC(4) + MOD(ELC(5)*CENT,1.D0)*SEC360
      EL = MOD(EL,SEC360)
C
      ELP = ELPC(1)*CENT**3 + ELPC(2)*CENT**2 + ELPC(3)*CENT
     .    + ELPC(4) + MOD(ELPC(5)*CENT,1.D0)*SEC360
      ELP = MOD(ELP,SEC360)
C
      F  =  FC(1)*CENT**3  +  FC(2)*CENT**2  +  FC(3)*CENT
     .   +  FC(4) +  MOD( FC(5)*CENT,1.D0)*SEC360
      F  =  MOD(F ,SEC360)
C
      D  =  DC(1)*CENT**3  +  DC(2)*CENT**2  +  DC(3)*CENT
     .   +  DC(4) +  MOD( DC(5)*CENT,1.D0)*SEC360
      D  =  MOD(D,SEC360)
C
      OM =  OMC(1)*CENT**3  + OMC(2)*CENT**2  + OMC(3)*CENT
     .   +  OMC(4) + MOD(OMC(5)*CENT,1.D0)*SEC360
      OM =  MOD(OM,SEC360)
 
****  Get the argument for the FCN
*                                                 ! Argument in seconds
      FCN_arg = (Cent*36525.d0)/FCN_period*sec360
C
      DO I = 1,15
C
C****     COMPUTE ARGUMENT FOR THIS TERM
          ARG = X1(1,I)*EL + X1(2,I)*ELP + X1(3,I)*F + X1(4,I)*D +
     .          X1(5,I)*OM + X1(6,I)*FCN_arg
C
          ARG = DMOD(ARG,SEC360)*CONVDS
C
*                             ! All but the precession constant
          if( i.gt.1 )  then
*                                                      ! In phase PSI
              coeff_deriv(1,2*i-1) = deriv(1)*sin(arg)
*                                                      ! Out of phase PSI
              coeff_deriv(2,2*i-1) = deriv(1)*cos(arg)
 
*                                                      ! In phase EPS
              coeff_deriv(1,2*i  ) = deriv(2)*cos(arg)
*                                                      ! Out of phase EPS
              coeff_deriv(2,2*i  ) = deriv(2)*sin(arg)
*                             ! Do the precession constant
          else
 
*                                                      ! Convert to years
              coeff_deriv(1,1) = deriv(1)*cent*100.d0
*                                                ! Offset for Out-of-phase
              coeff_deriv(2,1) = deriv(1)
 
*                                                      ! Convert to years
              coeff_deriv(1,2) = deriv(2)*cent*100.d0
*                                                ! Offset for Out-of-phase
              coeff_deriv(2,2) = deriv(2)
 
          end if
*                     ! Loop over terms
      end do
 
****  THATS ALL
      return
      end
 
