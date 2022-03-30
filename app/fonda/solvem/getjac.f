	  subroutine getjac(xin,yin,zin,cjaco,job)
c	
c	get Jacob matrix from one coordinate system to another
c	job:
c        1 = d(lon,lat,radius)/d(x,y,z) spherical to Cartesian
c        2 = d(x,y,z)/d(lon,lat,radius) Cartesian to spherical
c        3 = d(lon,lat,height)/d(x,y,z) ellipsoidal to Cartesian
c        4 = d(x,y,z)/d(lon,lat,height) Cartesian to ellipsoidal 
c        5 = d(e,n,u)/d(x,y,z) local Cartesian to global Cartesian
c        6 = d(x,y,z)/d(e,n,u) global Cartesian to local Cartesian
c
c	xin,yin,zin   : components of original coordinate (m or radian)
c	in spherical coordinate: xin->lon, yin->lat, zin->radius
c	in ellipsoidal coordinate: xin->lon, yin->lat, zin->height
c   The reason for such an order is to keep 2-D at the spherical
c   surface
c
      implicit real*8 (a-h,o-z)
      include 'solvem.fti'

      integer job
      real*8 cjaco(9)

      if (job.eq.1) then
         s1 = dsin(xin)
         c1 = dcos(xin)
         s2 = dsin(yin)
         c2 = dcos(yin)
         cjaco(1) = -s1/zin/c2
         cjaco(2) = c1/zin/c2
         cjaco(3) = 0.0d0
         cjaco(4) = -s2*c1/zin
         cjaco(5) = -s2*s1/zin
         cjaco(6) = c2/zin
         cjaco(7) = c2*c1
         cjaco(8) = c2*s1
         cjaco(9) = s2
      endif

      if (job.eq.2) then
         rh = xin**2+yin**2
         r = dsqrt(rh+zin**2)
         rh = dsqrt(rh)
         cjaco(1) = -yin
         cjaco(2) = -xin*zin/rh
         cjaco(3) = xin/r
         cjaco(4) = xin
         cjaco(5) = -yin*zin/rh
         cjaco(6) = yin/r
         cjaco(7) = 0.0d0
         cjaco(8) = rh
         cjaco(9) = zin/r
      endif

      if (job.eq.3) then
         s1 = dsin(xin)
         c1 = dcos(xin)
         s2 = dsin(yin)
         c2 = dcos(yin)
         f = 1.0d0/finv
         e2 = 2.0d0*f-f*f
         fac = 1.0d0/(1.0d0-e2*s2**2)
         zno = radius*dsqrt(fac)
         zh = zno+zin
         fac1 = 1.0d0-e2+e2*zin/zh
         cjaco(1) = -s1/zh/c2
         cjaco(2) = c1/zh/c2
         cjaco(3) = 0.0d0
         cjaco(4) = -s2*c1/zh
         cjaco(5) = -s2*s1/zh
         cjaco(6) = c2/zh*fac1
         cjaco(7) = c2*c1*fac
         cjaco(8) = c2*s1*fac
         cjaco(9) = s2*(1.0d0-e2)*fac
      endif

      if (job.eq.4) then
      call geoxyz(radius,finv,tx,ty,tz,a1,a2,ah,xin,yin,zin,2,rh)
         f = 1.0d0/finv
         e2 = 2.0d0*f-f*f
         fac = 1.0d0-e2+e2*ah/rh
         sf = dsin(a1)
         cf = dcos(a1)
         sl = dsin(a2)
         cl = dcos(a2)
         dnf = (rh-ah)*e2*sf*cf/(1.0d0-e2*sf*sf)
         cjaco(1) = -yin
         cjaco(2) = -rh*sf*cl+dnf*cf*cl
         cjaco(3) = cf*cl
         cjaco(4) = xin
         cjaco(5) = -rh*sf*sl+dnf*cf*sl
         cjaco(6) = cf*sl
         cjaco(7) = 0.0d0
         cjaco(8) = rh*fac*cf+(1.0d0-e2)*dnf*sf
         cjaco(9) = sf
      endif

      if (job.eq.5) then
c        input is lon,lat and radius instead of Xe,Xn and Xu
         sf = dsin(yin)
         cf = dcos(yin)
         sl = dsin(xin)
         cl = dcos(xin)
         cjaco(1) = -sl
         cjaco(2) = cl
         cjaco(3) = 0.0d0
         cjaco(4) = -sf*cl
         cjaco(5) = -sf*sl
         cjaco(6) = cf
         cjaco(7) = cf*cl
         cjaco(8) = cf*sl
         cjaco(9) = sf
      endif

      if (job.eq.6) then
c        input is lon,lat and radius instead of Xe,Xn and Xu
         sf = dsin(yin)
         cf = dcos(yin)
         sl = dsin(xin)
         cl = dcos(xin)
         cjaco(1) = -sl
         cjaco(4) = cl
         cjaco(7) = 0.0d0
         cjaco(2) = -sf*cl
         cjaco(5) = -sf*sl
         cjaco(8) = cf
         cjaco(3) = cf*cl
         cjaco(6) = cf*sl
         cjaco(9) = sf
      endif

      return
      end
