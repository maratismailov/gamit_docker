      subroutine imfh1p0( a,b,c, height,elev, hmf ) 

c     Received from A. Niell 12 January 2004. 
c
c    Modified by PT 040113: The mapping function coefficients are 
c    calculated from the met information in program GRDTAB. The 
c    appropriate values (not corrected for station height) are 
c    passed into this routine as arguments in the array mapcof.

c    Modified by RWK 060722:  Met information now read from the u-file,
c    written by GRDTAB from a station list or grid.  Changed call to 
c    input a,b,c as separate variables instead of array mapcof, and
c    to return imfh and dimfhdel as an array, to be consistent  with 
c    other IMF and VMF subroutines

*     000519 aen converted from imfh_dz.m
*     000622 aen correct value of pi
*     000623 aen change variables to real*8 except for passed ones.

*     Routine to compute the hydrostatic mapping function imfh which
*     depends on station position (latitude and height above geoid)
*     and height of 200 mb isobar; (use ellipsoid height for now).
*     and derivatives with respect to elevation, dimfh/dz, and 
*     z200, dimfhdz200.

*   input:
c     a, b, c      extended fraction coefficients from u-file
*     latitude     latitude of site (degrees)
*     height       height of site above geoid (meters)
*     z200         height of 200 mb isobar (meters)
*     elev         satellite elevation angle (degrees)

*   output:      
*     hmf(2)       hydrostatic mapping function value and derivative wrt elev (1/rad)
*     (Removed by PT:  dimfhdz200   derivative of imfh with respect to z200  (1/m);

*   define continued fraction coefficients for isobar dependence
*     a         - a coefficient of imf
*     b         - b coefficient of imf
*     c         - c coefficient of imf
*     z200      - height of 200 mb isobar (above sea level, meters)
*     z0        - mean value of z
*     z1        - amplitude of latitude dependence of mean value of z

*     latz0     - bias in latitude of z0
*     a00       - mean value of a
*     a01       - amplitude of latitude dependence of mean value of a
*     dadz0     - mean value of z dependence of a
*     dadz1     - amplitude of latitude dependence of dadh
*     bm0       - median value of b
*     cm0       - mean value of c
*     cm1       - amplitude of latitude dependence of mean value of c
*     lata0, latd0, latc0 - bias in latitude for each coefficient
      
*   a_ht, b_ht, c_ht - parameters for continued fraction for height corr'n.
*   dhcc_del    - derivative of height correction coefficient with elevation
*   dht_corr_del - derivative of height correction with elevation

      implicit none

      real*8 height,   elev
      real*8 imfh, dimfhdel, hmf(2)
      real*8 hs_km,    pi,       deg2rad
      real*8 a_ht, b_ht, c_ht
      real*8 a, b, c, dmda
      real*8 elr,  num,  se, ce, t1, t2, denom
      real*8 ht_corr_coef, ht_corr, dhcc_del, dht_corr_del
      real*8 temp1, temp2

      parameter(a_ht = 2.53e-5)
      parameter(b_ht = 5.49e-3)
      parameter(c_ht = 1.14e-3)
      parameter(pi   = 3.14159265358979)

c debug
c      print*,'a b c  ',a,b,c
c      print*,'height,elev ',height,elev

* calculate conversion factor:
      deg2rad = pi/180.0

* convert height of site in meter to kilometers
      hs_km = height/1000.

c** rwk 060733
c*      a       =  mapcof(1)               
c*      b       =  mapcof(2)
c*      c       =  mapcof(3)                                

c PT040113: the rest of the code below remains intact as per
c           A. Niell's original code.
c debug
cd     write (*,'("hs_km, alat, dadzlat, z200, zm, a, b, c",
cd    .  f15.6, 2e15.6, 2f15.6, 3e15.6)') 
cd    .  hs_km, alat, dadzlat, z200, zm, a, b, c

*  Now the coefficients exist; calculate for the sea level part
*  the mapping function, imfh, and the derivative wrt elevation
*  dimfh/d_el = dimfhdel.

*  To get delay-rate correction d_tau/dt:
*     d_tau/dt = d_tau-zen/dt * imfh + tau-zen * dimfhdel * d_el/dt
*                + tau-zen * dimfh/dz200 * dz200/dt
*     where  dimfhdel=dhmf/d_el

* calculate continued fraction explicitly
      elr=pi/180*elev
      num   = 1+a/(1+b/(1+c))
      se    = sin(elr)
      ce    = cos(elr)
      t1    = se+c
      t2    = se+b/t1
      denom = se+a/t2
* denom = se+a/(se+b/(se+c))
      imfh     = num/denom
      dimfhdel = -num*ce / denom**2
     .           *( 1.d0 - a/t2**2
     .           *( 1.d0 - b/t1**2 ) )

c debug
cd     write (*,'("elr, num, se, ce, t1, t2, denom, imfh, dimfhdel",
cd    .  9e15.6)') 
cd    .  elr, num, se, ce, t1, t2, denom, imfh, dimfhdel

* Calculate partial derivative w.r.t. z200
* dmdz = dmda*dadz where dadz=dadzlat
      temp1 = 1/(1+a+b/(1+c))
      temp2 = 1/(se*t2+a)
      dmda  = imfh*(temp1-temp2)
c dimfhdz is never used. Comment it out
c      dimfhdz  = dmda * dadzlat

c debug
cd     write (*,'("temp1, temp2, dmda, dimfhdz",
cd    .  4e15.6)') 
cd    .  temp1, temp2, dmda, dimfhdz

*  Apply height correction to mapping function and derivative wrt elevation:

*     1) height correction coefficient is 
*        1/sine(elev) - continued fraction(a_ht,b_ht,c_ht).
*     2) height correction is ht_corr_coef times height in km.
*     3) height correction to derivative wrt elevation is (derivative of
*        height correction coefficient wrt elevation)*height in km.

      t1     = b_ht/( se + c_ht )
      t2     = a_ht/( se + t1)
      num    = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
      ht_corr_coef = 1/se - num/(se+t2)
      ht_corr      = ht_corr_coef * hs_km
      imfh         = imfh + ht_corr
      hmf(1) = imfh

c debug
cd     write (*,'("t1, t2, num, ht_corr_coef, ht_corr, imfh",
cd    .  6e15.6)') 
cd    .  t1, t2, num, ht_corr_coef, ht_corr, imfh

*   951129 The derivative of the height correction wrt elevation is added
*   to dimfhdel after Chris Jacobs pointed out the magnitude of the term.

      dhcc_del   = -ce/se**2
     .             +num*ce / ( se + t2)**2
     .             *( 1.d0 - a_ht / ( se + t1  )**2
     .             *( 1.d0 - b_ht / ( se + c_ht)**2) )
      dht_corr_del = dhcc_del * hs_km
      dimfhdel     = dimfhdel + dht_corr_del
      hmf(2) = dimfhdel

c debug
cd     write (*,'("dhcc_del, dht_corr_del, dimfhdel",
cd    .  3e15.6)') 
cd    .  dhcc_del, dht_corr_del, dimfhdel

      return
      end
