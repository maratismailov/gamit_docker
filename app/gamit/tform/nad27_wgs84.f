      program test

      real*8 radcon,alat,along,dphi,dlam,dhgt

      radcon = datan(1.d0)/45.d0



c  DMA test case
c      alat = 34.d0 + (47.d0 + (08.833d0/60.d0))/60.d0
c      along=273.d0 + (25.d0 + (07.825d0/60.d0))/60.d0


c  Pinon 7256    N33 36 33.18476   W116 27 28.54386   1244.561 m
c                                  E243 32 31.45614

      alat = 33.d0 + (36.d0 + (33.18476d0/60.d0))/60.d0
      along=243.d0 + (32.d0 + (31.45614d0/60.d0))/60.d0

      call nad27_wgs84(alat,along,dphi,dlam,dhgt)
      print*, dphi,dlam,dhgt
      end

      subroutine nad27_wgs84(phi,lam,dphi,dlam,dhgt)
c
c Purpose:
c     Convert nad27 geodetic coordinates to wgs84 geodetic coordinates
c     using multiple regression equations
c
c Input:
c     phi,lam -- geodetic latitude (-90:90), longitude (0:360 east)
c                  (in decimal degrees)
c Output:
c     dphi,dlam -- corrections to convert to WGS84 latitude and longitude
c                  (in seconds)
c     dhgt     --  correction to convert to WGS84 geodetic height above ellipsoid
c                  (in meters)
c
c Source:
c     DMA TR 8350.2-B Section 19
c     RMS differences from doppler derived coordinates:
c        phi: 1.3 m , lam: 1.3 m, hgt: 1.2 m
c
      real*8 phi,lam,dphi,dlam,dhgt
      real*8 k,phi0,lam0,u,v,p(23),l(23),h(20)

      data k,phi0,lam0/0.05235988d0,  37.d0,   265.d0/

      data p/   0.16984d0, -0.76173d0,  0.09585d0,   1.09919d0
     .      ,  -4.57801d0, -1.13239d0,  0.49831d0,  -0.98399d0
     .      ,   0.12415d0,  0.11450d0, 27.05396d0,   2.03449d0
     .      ,   0.73357d0, -0.37548d0, -0.14197d0, -59.96555d0
     .      ,   0.07439d0, -4.76082d0,  0.03385d0,  49.04320d0
     .      ,  -1.30575d0, -0.07653d0,  0.08646d0/

      data l/  -0.88437d0,  2.05061d0,  0.26361d0,  -0.76804d0
     .      ,   0.13374d0, -1.31974d0, -0.52162d0,  -1.05853d0
     .      ,  -0.49211d0,  2.17204d0, -0.06004d0,   0.30139d0
     .      ,   1.88585d0, -0.81162d0, -0.05183d0,  -0.96723d0
     .      ,  -0.12948d0,  3.41827d0,  0.44507d0,   0.18882d0
     .      ,  -0.01444d0,  0.04794d0, -0.59013d0/

      data h/ -36.52600d0,  3.90000d0, -4.72300d0, -21.55300d0
     .      ,   7.29400d0,  8.88600d0, -8.44000d0,  -2.93000d0
     .      ,  56.93700d0,-58.75600d0, -4.06100d0,   4.44700d0
     .      ,   4.90300d0,-55.87300d0,212.00500d0,   3.08100d0
     .      ,-254.51100d0, -0.75600d0, 30.65400d0,  -0.12200d0/

      u = k*(phi - phi0)
      v = k*(lam - lam0)

      dphi = p(1)           + p(2)*u     + p(3)*v     + p(4)*u*u
     .     + p(5)*u**3      + p(6)*u*u*v + p(7)*v**3  + p(8)*u**3*v
     .     + p(9)*u*v**3    + p(10)*v**4 + p(11)*u**5 + p(12)*u**4*v
     .     + p(13)*u*u*v**3 + p(14)*v**5 + p(15)*v**6 + p(16)*u**7
     .     + p(17)*v**7     + p(18)*u**8 + p(19)*v**8 + p(20)*u**9
     .     + p(21)*u**6*v**3 + p(22)*u**3*v**9 + p(23)*u**4*v**9

      dlam = l(1)            + l(2)*v       + l(3)*u*u    + l(4)*u*v
     .     + l(5)*v*v        + l(6)*u**3    + l(7)*u*u*v  + l(8)*u*v*v
     .     + l(9)*u*u*v*v    + l(10)*u*v**3 + l(11)*v**4  + l(12)*u**4*v
     .     + l(13)*u*v**4    + l(14)*u*v**5 + l(15)*v**6  + l(16)*u*v**6
     .     + l(17)*u**3*v**5 + l(18)*u**9  + l(19)*u**8*v + l(20)*u*v**8
     .     + l(21)*v**9      + l(22)*u*v**9 + l(23)*u**9*v**3

      dhgt = h(1)            + h(2)*u       + h(3)*v      + h(4)*u*u
     .     + h(5)*u*v        + h(6)*v*v     + h(7)*u*u*v  + h(8)*u*v*v
     .     + h(9)*u**4       + h(10)*u**3*v + h(11)*v**4  + h(12)*u**4*v
     .     + h(13)*u*u*v**3  + h(14)*u**6  + h(15)*u**5*v + h(16)*v**6
     .     + h(17)*u**7*v    + h(18)*v**8  + h(19)*u**8*v + h(20)*u*v**9

      return
      end
