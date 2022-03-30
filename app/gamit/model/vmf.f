cc

      subroutine vmf (ah,aw,dlat,height,zd,dvmfh,dvmfw)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine generates the values for the
c  fast approach of the Vienna mapping function
c
c  input data
c  ----------
c  ah:     hydrostatic coefficient (zero height)
c  aw:     wet coefficient
c  dlat:   latitude in radians
c  height: ellipsoidal height in meters
c  zd:     zenith distance in radians
c
c  output data
c  -----------
c  dvmfh:  hydrostatic mf
c  dvmfw:  wet mf
c
c  jboehm, 2003 Feb 27
c  PT040114: subroutine used exactly as it was provided by J. Boehm
c            except that implicit declarations were removed (to remove
c            nasty-looking (but benign) compiler warnings
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      implicit double precision (a-h,o-z)
      implicit none

      real*8 pi,zd,zd1,bh,cm0,cm1,ch,elev,sine,cose,beta,gamma,topcon
     .      ,dlat,a_ht,b_ht,c_ht,hs_km,height,ah,aw,bw,cw,ht_corr
     .      ,ht_corr_coef,dvmfh,dvmfw
      PI = 3.14159265359d0

c  coefficients for the determination of the mapping functions

      zd1  = zd  *180.0d0/PI


c  imf_h coefficients

      bh  = 0.002905
      cm0 = 0.0634
      cm1 = 0.0014
      ch  = cm0 + cm1*cos(2*dlat)

      elev = 90.d0 - zd1

      sine   = sin(elev * pi/180.d0)
      cose   = cos(elev * pi/180.d0)
      beta   = bh/( sine + ch )
      gamma  = ah/( sine + beta)
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))

      dvmfh = topcon/(sine+gamma)

c  height correction

      a_ht = 2.53e-5
      b_ht = 5.49e-3
      c_ht = 1.14e-3

      hs_km  = height/1000.d0

      beta         = b_ht/( sine + c_ht )
      gamma        = a_ht/( sine + beta)
      topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
      ht_corr_coef = 1/sine - topcon/(sine + gamma)
      ht_corr      = ht_corr_coef * hs_km
      dvmfh        = dvmfh + ht_corr

c  wet mapping function

c coefficients from Niell 96

      bw = 0.00146
      cw = 0.04391

      sine  = sin( elev * pi/180.d0)
      cose  = cos( elev * pi/180.d0)
      beta  = bw/( sine + cw )
      gamma = aw/( sine + beta)
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))

      dvmfw = topcon/(sine+gamma)

      return
      end
