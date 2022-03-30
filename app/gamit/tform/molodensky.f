      program test

      real*8 alat,along,hgt,dlat,dlon,dhgt
      real*8 a,finv,da,df,dx,dy,dz

c     nad27 to wgs84
      data a,finv/6378206.4d0, 294.9786982d0/
      data da,df/-69.4d0, -0.37264639d-04/
c     Contiguous US mean datum shift
      data dx,dy,dz/  -8.d0, 160.d0, 176.d0/

c  DMA test case
c     alat = 42.d0 + (56.d0 + (51.9d0/60.d0))/60.d0
c     along=288.d0 + (22.d0 + (22.6d0/60.d0))/60.d0
c     hgt = 232.d0

c     from charts:
c     dx = -13.d0
c     dy = 165.d0
c     dz = 185.d0

C  Pinon 7256
      alat = 33.d0 + (36.d0 + (33.18476d0/60.d0))/60.d0
      along=243.d0 + (32.d0 + (31.45614d0/60.d0))/60.d0
      hgt = 1244.561

c     from charts:
      dx =  -9.d0
      dy = 158.5d0
      dz = 174.d0

      call  molodensky(a,finv,da,df,dx,dy,dz,
     .                   alat,along,hgt,dlat,dlon,dhgt)

      print*, dlat,dlon,dhgt

      end


      subroutine molodensky(a,finv,da,df,dx,dy,dz,
     .                      lat,lon,hgt,dlat,dlon,dhgt)
c
c Purpose:
c     Convert from one geodetic curvilinear system (1) to another
c     geodetic curvilinear system (2).
c
c     Uses the Standard Molodensky Datum Transformation formulas.
c
c Input:
c     a,finv    -- semimajor axis (metric) and inverse of flattening
c                  for system (1) ellipsoid
c     lat,lon   -- geodetic latitude and longitude in system (1)
c                   (decimal degrees)
c     hgt       -- geodetic height in system (1) (metric)
c     dx,dy,dz  -- shift between centers of systems: (2)-(1) (metric)
c     da,df     -- difference between semimajor axis and flattening
c                  between systems (2)-(1)
c
c Output:
c     dlat,dlon -- latitude and longitude corrections (arc seconds)
c                    lat(2) = lat(1) + dlat, etc.
c     dhgt      -- height correction (metric)
c
c Note:
c     To convert from NAD27 to WGS84 geodetic coordinates use the
c     following values (metric in meters):
c
c     data a,finv/6378206.4d0, 294.9786982d0/
c     data da,df/-69.4d0, -0.37264639d-04/
c     Contiguous US mean datum shift
c     data dx,dy,dz/  -8.d0, 160.d0, 176.d0/
c
c     For superior results consult graphes for local estimates of dx,dy,dz
c     e.g. LA gives: dx = -10, dy = 159, dz = 175
c
      real*8 a,finv,da,df,dx,dy,dz
      real*8 lat,lon,hgt,dlat,dlon,dhgt
      real*8 radcon,onesec,sinone,sinlat,coslat,sinlon,coslon
      real*8 f,boa,esq,N,M,den,num

      radcon = datan(1.d0)/45.d0
      onesec = 1.d0/3600.d0
      sinone = dsin(onesec*radcon)

      sinlat = dsin(lat*radcon)
      coslat = dcos(lat*radcon)
      sinlon = dsin(lon*radcon)
      coslon = dcos(lon*radcon)

c.....flattening of (1)
      f   = 1.d0/finv

c.....ratio of semiminor axis to semimajor axis of (1):  b = a*(1-f)
      boa = 1.d0 - f

c.....square of the first eccentricity of (1)
      esq = 2.d0*f - f*f

c.....radius of curvature of the prime vertical of (1)
      den = dsqrt(1.d0-esq*sinlat*sinlat)
      N = a/den

c.....radius of curvature of the meridian of (1)
      M = a*(1-esq)/den**3

c.....standard molodensky formulas
      num = -dx*sinlat*coslon - dy*sinlat*sinlon + dz*coslat
     .      + da* N/a*esq *sinlat*coslat
     .      + df* (M/boa + N*boa) *sinlat*coslat

      dlat = num / ((M+hgt)*sinone)

      dlon = (-dx*sinlon + dy*coslon) / ((N+hgt)*coslat*sinone)

      dhgt = dx*coslat*coslon + dy*coslat*sinlon + dz*sinlat
     .       -da*a/N + df*boa*N*sinlat*sinlat

      return
      end
