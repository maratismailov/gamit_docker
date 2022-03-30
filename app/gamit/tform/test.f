      program test

      real*8 x0,y0,z0,dx,dy,dz,n,e,u
      real*8 alat,along,hght,radcon

      radcon = datan(1.d0)/45.d0
      alat = 42.d0 + (56.d0 + (51.9d0/60.d0))/60.d0
      along=288.d0 + (22.d0 + (22.6d0/60.d0))/60.d0
      hght = 232.d0
      alat = alat*radcon
      along= along*radcon
      call geoxyz(alat,along,hght,x0,y0,z0,1)
      print*,x0,y0,z0
C      print*,'input x0,y0,z0'
C      read (*,*) x0,y0,z0
C      print*,'input dx,dy,dz'
C      read (*,*) dx,dy,dz
C      call locgeo(x0,y0,z0,dx,dy,dz,n,e,u)
C      print*,n,e,u
      end
c---------------------------------------------------------
      subroutine locgeo(x0,y0,z0,dx,dy,dz,n,e,u)
c
c Purpose:
c     Convert vector from geocentric Cartesian to local geodetic system
c
c Input:
c     (x0,y0,z0) :   Cartesian site coordinate (meters)
c     (dx,dy,dz) :   Cartesian vector to be converted
c
c Output:
c     ( n, e, u) :   Local geodetic vector

      integer iflag
      real*8 x0,y0,z0,dx,dy,dz,n,e,u
      real*8 alat,along,hght,slat,slon,clat,clon

c     convert cartesian site coordinates (x0,y0,z0) to
c        geodetic curvilinear coordinates (alat,alon,hght)

      iflag = 2
      call geoxyz(alat,along,hght,x0,y0,z0,iflag)

c     convert local vector (dx,dy,dz) to (n,e,u)

      slat = dsin(alat)
      slon = dsin(along)
      clat = dcos(alat)
      clon = dcos(along)
c
      n = -slat*clon*dx -slat*slon*dy +clat*dz
      e =      -slon*dx      +clon*dy
      u =  clat*clon*dx +clat*slon*dy +slat*dz

      return
      end
c---------------------------------------------------------------------
      subroutine geoxyz(alat,along,hght,x,y,z,iflag)
c
c Purpose:
c     Convert geodetic curvilinear coordinates to geocentric Cartesian
c        coordinates and vice versa
c
c Input:
c     iflag = 1  convert geodetic coordinates to cartesian coordinates
c           = 2  convert cartesian coordinates to geodetic coordinates
c
c Input/Output:
c     alat,along : geodetic latitude and longitude (radians)
c     hght       : height above the reference ellipsiod (meters)
c     x,y,z      : geocentric Cartesian coordinates (meters)
c
c Notes:
c     Currently assumes the WGS84 reference ellipsoid;
c     Clarke 1866 ellipsoid with approximate translation parameters
c        for NAD27 are commented.
c     Cartesian to geodetic conversion uses an iterative scheme
c        valid to the millimeter level.
c
      integer iflag
      real*8 alat,along,hght,x,y,z
      real*8 semi,finv,tx,ty,tz
      real*8 twopi,f,e2,curvn,sqr,alat0,cutoff
      real*8 sinlat,coslat,sinlon,coslon

c  NAD27: Clarke 1866 ellipsoid, approximate translations
      data semi,finv,tx,ty,tz/ 6378206.4d0, 294.9786982d0,
c     1                         -8.d0, 160.d0, 176.d0/
     1                         -13.d0, 165.d0, 185.d0/
c    1                         -12.010d0, 162.970d0, 189.740d0/
c
c  NAD83 = WGS84 ellipsoid
C      data semi,finv,tx,ty,tz/6378137.d0,298.257222101d0,3*0.d0/
C      data semi,finv,tx,ty,tz/6378137.d0,298.257223563d0,3*0.d0/

      twopi= 8.d0*datan(1.d0)
      f= 1.d0/finv
      e2= 2.d0*f - f*f
      if( iflag.eq.1) then
         sinlat= dsin(alat)
         coslat= dcos(alat)
         sinlon= dsin(along)
         coslon= dcos(along)
         curvn= semi/(dsqrt(1.d0-e2*sinlat*sinlat))
c
         x= (curvn+hght)*coslat*coslon + tx
         y= (curvn+hght)*coslat*sinlon + ty
         z= (curvn*(1.d0-e2)+hght)*sinlat + tz
      else
         x= x - tx
         y= y - ty
         z= z - tz
         along= datan2(y,x)
         if( along.lt.0d0 ) along=along + twopi
c        starting value for latitude iteration
         sqr= dsqrt(x*x+y*y)
         alat0= datan2(z/sqr,1.d0-e2)
         alat= alat0
   40    sinlat= dsin(alat)
         curvn= semi/(dsqrt(1.d0-e2*sinlat*sinlat))
         alat= datan2((z+e2*curvn*sinlat),sqr)
c        iterate to millimeter level
         if( dabs(alat-alat0).lt.1.d-10) goto 30
         alat0= alat
         goto 40
   30    continue
         cutoff= 80.d0*twopi/360.d0
         if(alat.le.cutoff) then
            hght= (sqr/dcos(alat))-curvn
         else
            hght= z/dsin(alat)-curvn+e2*curvn
         endif
      endif
      return
      end
