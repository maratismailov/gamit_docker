      subroutine dxyz_to_neu (xyzvec, sigxyz, lat, long, height, 
     1                         iflag, neuvec, signeu)
c
c     Written by Y. Bock March 1, 2000
c
c     Purpose:
c        Convert Cartesian xyz baselines to neu baselines, given the
c        spherical or geodetic latitude and longitude, Cartesian
c        xyz vector and covariance

C Input
C     xyzvec : baseline vector between site 1 and site 2 (sense 2-1)
C     sigxyz : covariance matrix of xyzvec (symmetric storage)
C     lat, long : geodetic latitude and longitude of site 1 in radians (iflag=1)
C     lat, long : geocentric latitude and longitude of site 1 in radians (iflag=2)
C Output:
C     neuvec : local vector with origin at site 1
C     signeu : covariance matrix of neuvec    
c
      integer i,j, iflag
      real*8 xyzvec(3), neuvec(3), sigxyz(6),signeu(6)
      real*8 lat,long,height
      real*8 slat,slon,clat,clon
      real*8 x1,y1,z1,delx,dely,delz,dist
      real*8 jac(3,3)
      real*8 semi,finv,geodrad
      real*8 temp1(3),temp2(3)

C Hard wire WGS84 ellipsoid parameters (semi-major axis in km)
      data semi/6378.137d0/finv/298.257223563d0/

C Convert to geodetic coordinates for base station
      if(iflag.eq.1) then
C   Convert geocentric to Cartesian
         call sphxyz(lat, long, height, x1, y1, z1, 1)
C   Convert Cartesian to geodetic
         call geoxyz(semi,finv,lat, long, height, geodrad
     2              ,x1, y1, z1 , 2 )
      endif      

      delx = xyzvec(1)
      dely = xyzvec(2)
      delz = xyzvec(3)
      slat = dsin(lat)
      slon = dsin(long)
      clat = dcos(lat)
      clon = dcos(long)

C Set up Jacobian matrix
c     (n,e,u) = J(x,y,z)
      do j = 1, 3
      do i = 1, 3
         jac(i,j) = 0.d0
      enddo
      enddo
C
      jac(1,1) = -slat*clon
      jac(1,2) = -slat*slon
      jac(1,3) =  clat
      jac(2,1) = -slon
      jac(2,2) =  clon
      jac(2,3) =  0.d0
      jac(3,1) =  clat*clon
      jac(3,2) =  clat*slon
      jac(3,3) =  slat

C Calculate local coordinates
      neuvec(1) = jac(1,1)*delx + jac(1,2)*dely + jac(1,3)*delz
      neuvec(2) = jac(2,1)*delx + jac(2,2)*dely + jac(2,3)*delz
      neuvec(3) = jac(3,1)*delx + jac(3,2)*dely + jac(3,3)*delz

C Check distance calculation
      dist = dsqrt(neuvec(1)*neuvec(1)+neuvec(2)*neuvec(2)+
     .             neuvec(3)*neuvec(3))
C      write(6,*) 'dist', dist

C Calculate covariance matrix for local coordinates
      call gpgt(jac,sigxyz,3,3,signeu,temp1,temp2)

      return
      end









