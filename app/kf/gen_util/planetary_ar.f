CTITLE 'PLANETARY_ARG'
 
      subroutine planetary_arg( epoch, iarg, arg, arg_period)
 
 
      implicit none 

*     Routine to compute of planetary arguments for planetary
*     nutation.  The longitudes of the major planets is computed
*     according to:
*     Bretagnon, P., The'orie du mouvement de l'ensemble des
*     plane`tes. Solution VSOP82*, Astron. Astrophys., 114,
*     278--288, 1982. We adopt only the linear term in the rate,
*     high order terms are of order 1.d-6 per century**2.
*
*     Expressions are taken from Tables 4 and 7. Table 7 units have
*     been changed from rads/thousand years to rads/cent.  Values are
*     relative to the fixed dynamic ecliptic of J2000.
*
*     Period is return in sidereal days.
 
      include '../includes/const_param.h'
 
*         iarg(5)   - Arguments for Venus,Earth,Mars,Jupiter,
*                   - Saturn
 
      integer*4 iarg(5)
 
*      arg          - Argument in radians
*      argr         - Rate of change of argument (cycles/Jul.
*                   - cent)
*      arg_period   - Period of argument in sidereal days
*      cent         - Centuries since J2000.
*      epoch        - Julian date for evaluation
*      vl,vld       - Venus long and rate (rads, rads/cent)
*      vlc(2)       - coefficients for computing vl
*      tl,tld       - Earth long and rate (rads, rads/cent)
*      tlc(2)       - coefficients for computing tl
*      ml,mld       - Mars  long and rate (rads, rads/cent)
*      mlc(2)       - coefficients for computing ml
*      jl,jld       - Jupliter long and rate (rads, rads/cent)
*      jlc(2)       - coefficients for computing jl
*      sl,sld       - Saturn long and rate (rads, rads/cent)
*      slc(2)       - coefficients for computing sl
 
      real*8 arg, argr, arg_period, cent, epoch, vl,vld, vlc(2),
     .    tl,tld, tlc(2), ml,mld, mlc(2), jl,jld, jlc(2), sl,sld,
     .    slc(2)
 
      data vlc / 3.176 146 696 89d0, 1021.3285 546 211 00d0 /
      data tlc / 1.753 470 314 35d0,  628.3075 849 180 00d0 /
      data mlc / 6.203 480 913 41d0,  334.0612 431 492 30d0 /
      data jlc / 0.599 546 497 39d0,   52.9690 965 094 60d0 /
      data slc / 0.874 016 756 50d0,   21.3299 095 438 00d0 /
 
***** Get number of Centuries since J2000
 
      cent = (epoch-DJ2000) / 36525.d0
 
*     Compute arguments and rates
      vl  = vlc(1) + vlc(2)*cent
      vld = vlc(2)
 
      tl  = tlc(1) + tlc(2)*cent
      tld = tlc(2)
 
      ml  = mlc(1) + mlc(2)*cent
      mld = mlc(2)
 
      jl  = jlc(1) + jlc(2)*cent
      jld = jlc(2)
 
      sl  = slc(1) + slc(2)*cent
      sld = slc(2)
 
***** Now compute the argument
      arg  = mod( ( iarg(1)*vl + iarg(2)*tl + iarg(3)*ml + iarg(4)*jl
*                                             ! Radians
     .            + iarg(5)*sl ), 2.d0*pi)
      argr = ( iarg(1)*vld + iarg(2)*tld + iarg(3)*mld + iarg(4)*jld
*                                             ! Cycles per Julian cent
     .       + iarg(5)*sld)/ (2.d0*pi)
 
      arg_period = (36525.d0/ argr ) * solar_to_sidereal
 
***** Thats all
      return
      end
 
 
