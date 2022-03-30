CTITLE 'FUNDAMENTAL_ARG'
 
      subroutine fundamental_arg( epoch, iarg, arg, arg_period )
 
      implicit none
 
*     Routine to compute the value of the fundamental argument
*     for the sequence iarg.  The period is also return in sidereal
*     days.
 
      include '../includes/const_param.h'
 
*         iarg(5)       - Fundamental argumemts
 
      integer*4 iarg(5)
 
*      arg              - Argument (radians)
*      argr             - Rate of change of argument (cycles per
*                       - Julian century)
*      arg_period       - Period of the argument in sidereal days
*      cent             - Julian centuries to DJ2000.
*      epoch            - Julian date at which arqument required
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
 
 
      real*8 arg, argr, arg_period, cent, epoch, el,eld, elc(5),
     .    elp,elpd, elpc(5), f,fd, fc(5), d,dd, dc(5), om,omd, omc(5)
 
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
      arg = mod(iarg(1)*el + iarg(2)*elp + iarg(3)*f + iarg(4)*d +
     .          iarg(5)*om, sec360) / (3600.d0 * rad_to_deg)
 
      argr = (iarg(1)*eld + iarg(2)*elpd + iarg(3)*fd + iarg(4)*dd +
*                                    ! Cycles per Julian century
     .        iarg(5)*omd)/sec360
 
      arg_period = (36525.d0 / argr) * solar_to_sidereal
 
***** Thats all
      return
      end
 
