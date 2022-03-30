      subroutine vmf1 ( ah,aw,dmjd,dlat,ht,zd,vmf1h,vmf1w,grdvmf1 )
      
C     !!! This is the version with height correction !!!  
C     !!! It has to be used with the grid !!!    
C
C     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
C     Reference: Boehm, J., B. Werl, H. Schuh (2006), 
C     Troposphere mapping functions for GPS and very long baseline interferometry 
C     from European Centre for Medium-Range Weather Forecasts operational analysis data,
C     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
c
c     PT060601: Johannes' routine has been modified in two ways:
c               1. Explicit declarations of variables added
c               2. An extra argument added to indicate whether
c                  to apply the height correction or not.
C
C     input data
C     ----------
C     ah:   hydrostatic coefficient a (www.hg.tuwien.ac.at/~ecmwf1)
C     aw:   wet coefficient a         (www.hg.tuwien.ac.at/~ecmwf1)  
C     dmjd: modified julian date
C     dlat: latitude in radians
C     ht:   ellipsoidal height in meter
C     zd:   zenith distance in radians
c  grdvmf1: indicates whether vmf1 coefficients have been interpolated
c              from a grid (and hence need height corrections applied)
C
C     output data
C     -----------
C     vmf1h: hydrostatic mapping function
C     vmf1w: wet mapping function
C
C     Johannes Boehm, 2005 October 2
C

c      implicit double precision (a-h,o-z)
c PT060601:    change to explicit declarations for GAMIT  
      implicit none

      real*8 ah,aw,pi,dmjd,dlat,zd,doy,bh,c0h,phh,c11h,c10h,ch,sine
     .     , beta,gamma,topcon,vmf1h,bw,cw,vmf1w
     .     , a_ht,b_ht,c_ht,ht,hs_km,ht_corr,ht_corr_coef
      logical grdvmf1

      pi = 3.14159265359d0
      
C     reference day is 28 January
C     this is taken from Niell (1996) to be consistent
      doy = dmjd  - 44239.d0 + 1 - 28
      
      bh = 0.0029
      c0h = 0.062
      if (dlat.lt.0.d0) then   ! southern hemisphere
          phh  = pi
          c11h = 0.007
          c10h = 0.002
      else                     ! northern hemisphere
          phh  = 0.d0
          c11h = 0.005
          c10h = 0.001
      end if
      ch = c0h + ((dcos(doy/365.25d0*2.d0*pi + phh)+1.d0)*c11h/2.d0 
     .     + c10h)*(1.d0-dcos(dlat))


      sine   = dsin(pi/2.d0 - zd)
      beta   = bh/( sine + ch  )
      gamma  = ah/( sine + beta)
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))
      vmf1h   = topcon/(sine+gamma)
      
C  height correction [Niell, 1996]     
      if(grdvmf1)then
c        print*,'applying ht correction to vmf1 mf'
        a_ht = 2.53e-5
        b_ht = 5.49e-3
        c_ht = 1.14e-3
        hs_km  = ht/1000.d0
        beta         = b_ht/( sine + c_ht)
        gamma        = a_ht/( sine + beta)
        topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
        ht_corr_coef = 1.d0/sine - topcon/(sine + gamma)
        ht_corr      = ht_corr_coef * hs_km
        vmf1h        = vmf1h + ht_corr
      endif

      bw = 0.00146
      cw = 0.04391
      beta   = bw/( sine + cw )
      gamma  = aw/( sine + beta)
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))
      vmf1w   = topcon/(sine+gamma)
      
      end 

