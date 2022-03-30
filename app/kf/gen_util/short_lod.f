 
      subroutine short_period_lod( jd, dut1, dlod)
 
 
      implicit none 

*     Routine to compute the short period UT1 and LOD tidal contributions
*     given the Julian date.  Routine based on UT1ZT from CALC which
*     uses the Yoder corrections.
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
*   arg     - Argument for tidal correction
*   dj2000  - Julian date of J2000
*   dut1    - Tidal contribution to UT1 (seconds)
*   dlod    - Tidal contribition to LOD (seconds)
*   jd      - Julian date of dersired correction
*   el, elp, f, d, om   - Fundamental arguments (seconds)
*   seccon  - Converts seconds to radians
*   t       - Julian centuries since J2000
 
 
      real*8 arg, dj2000, dut1, dlod,
     .       jd, el, elp, f, d, om, seccon, t
 
*   x(8,41) - Terms in short period tidal contributions
 
 
      real*4 x(8,41)
 
      data
     .    dj2000 / 2451545.d0 /
     .,   seccon / 206264.8062470964d0 /
 
*               L    L'   F    D  Omega  DUT coeff  DLOD    Domega
      data x  / 1.,  0.,  2.,  2.,  2.,  -0.02,      0.3,    -0.2,
     .          2.,  0.,  2.,  0.,  1.,  -0.04,      0.4,    -0.3,
     .          2.,  0.,  2.,  0.,  2.,  -0.10,      0.9,    -0.8,
     .          0.,  0.,  2.,  2.,  1.,  -0.05,      0.4,    -0.4,
     .          0.,  0.,  2.,  2.,  2.,  -0.12,      1.1,    -0.9,
     .          1.,  0.,  2.,  0.,  0.,  -0.04,      0.3,    -0.2,
     .          1.,  0.,  2.,  0.,  1.,  -0.41,      2.8,    -2.4,
     .          1.,  0.,  2.,  0.,  2.,  -0.99,      6.8,    -5.8,
     .          3.,  0.,  0.,  0.,  0.,  -0.02,      0.1,    -0.1,
     .         -1.,  0.,  2.,  2.,  1.,  -0.08,      0.5,    -0.5,
     .         -1.,  0.,  2.,  2.,  2.,  -0.20,      1.3,    -1.1,
     .          1.,  0.,  0.,  2.,  0.,  -0.08,      0.5,    -0.4,
     .          2.,  0.,  2., -2.,  2.,   0.02,     -0.1,     0.1,
     .          0.,  1.,  2.,  0.,  2.,   0.03,     -0.1,     0.1,
     .          0.,  0.,  2.,  0.,  0.,  -0.30,      1.4,    -1.2,
     .          0.,  0.,  2.,  0.,  1.,  -3.21,     14.8,   -12.5,
     .          0.,  0.,  2.,  0.,  2.,  -7.76,     35.7,   -30.1,
     .          2.,  0.,  0.,  0., -1.,   0.02,     -0.1,     0.1,
     .          2.,  0.,  0.,  0.,  0.,  -0.34,      1.5,    -1.3,
     .          2.,  0.,  0.,  0.,  1.,   0.02,     -0.1,     0.1,
     .          0., -1.,  2.,  0.,  2.,  -0.02,      0.1,    -0.1,
     .          0.,  0.,  0.,  2., -1.,   0.05,     -0.2,     0.2,
     .          0.,  0.,  0.,  2.,  0.,  -0.73,      3.1,    -2.6,
     .          0.,  0.,  0.,  2.,  1.,  -0.05,      0.2,    -0.2,
     .          0., -1.,  0.,  2.,  0.,  -0.05,      0.2,    -0.2,
     .          1.,  0.,  2., -2.,  1.,   0.05,     -0.1,     0.1,
     .          1.,  0.,  2., -2.,  2.,   0.10,     -0.3,     0.2,
     .          1.,  1.,  0.,  0.,  0.,   0.04,     -0.1,     0.1,
     .         -1.,  0.,  2.,  0.,  0.,   0.05,     -0.1,     0.1,
     .         -1.,  0.,  2.,  0.,  1.,   0.18,     -0.4,     0.3,
     .         -1.,  0.,  2.,  0.,  2.,   0.44,     -1.0,     0.9,
     .          1.,  0.,  0.,  0., -1.,   0.53,     -1.2,     1.0,
     .          1.,  0.,  0.,  0.,  0.,  -8.26,     18.8,   -15.9,
     .          1.,  0.,  0.,  0.,  1.,   0.54,     -1.2,     1.0,
     .          0.,  0.,  0.,  1.,  0.,   0.05,     -0.1,     0.1,
     .          1., -1.,  0.,  0.,  0.,  -0.06,      0.1,    -0.1,
     .         -1.,  0.,  0.,  2., -1.,   0.12,     -0.2,     0.2,
     .         -1.,  0.,  0.,  2.,  0.,  -1.82,      3.6,    -3.0,
     .         -1.,  0.,  0.,  2.,  1.,   0.13,     -0.3,     0.2,
     .          1.,  0., -2.,  2., -1.,   0.02,      0.0,     0.0,
     .         -1., -1.,  0.,  2.,  0.,  -0.09,      0.2,    -0.1  /
 
 
***** Compute the fundamental arguments
 
      t = ( jd - dj2000 ) / 36525.d0
 
      el  = ((0.064d0*t + 31.310d0)*t + 715922.633d0)*t +
     .       485866.733d0 + mod(1325.0d0*t,1.d0) * 1296000.d0
      el  = mod(el,1296000.d0)
 
      elp = ((-0.012d0*t - 0.577d0)*t + 1292581.224d0)*t +
     .      1287099.804d0 + mod(99.d0*t,1.d0) * 1296000.d0
      elp = mod(elp,1296000.d0)
 
      f   = ((0.011d0*t - 13.257d0)*t + 295263.137d0)*t +
     .      335778.877d0 + mod(1342.d0*t,1.d0) * 1296000.d0
      f   = mod(f,1296000.d0)
 
      d   = ((0.019d0*t - 6.891d0)*t + 1105601.328d0)*t +
     .      1072261.307d0 + mod(1236.d0*t,1.d0) * 1296000.d0
      d   = mod(d,1296000.d0)
 
      om  = ((0.008d0*t + 7.455d0)*t - 482890.539d0)*t +
     .      450160.280d0 - mod(5.d0*t,1.d0) * 1296000.d0
      om  = mod(om,1296000.d0)
 
*     Now start the summation (do dut1 and dlod)
 
      dut1 = 0.d0
      dlod = 0.d0
 
      do j = 1,41
          i = 42 - j
 
          arg = x(1,i) * el + x(2,i)*elp + x(3,i)*f +
     .          x(4,i) * d  + x(5,i)*om
          arg = mod(arg,1296000.d0) / seccon
 
*         Evaluate the zonal terms
          dut1 = x(6,i) * sin(arg) + dut1
          dlod = x(7,i) * cos(arg) + dlod
      end do
 
*     Convert results to time seconds
      dut1 = dut1 * 1.d-4
      dlod = dlod * 1.d-5
 
      return
      end
