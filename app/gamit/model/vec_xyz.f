C Adapted from /lib/xyztogeod Aug 2007 by EJP
C Compute rotation matrix to convert the magnetic field vector
C from local spherical N,E,down to X,Y,Z

      subroutine vec_xyz(latdeg, londeg, rot_mat)

C In: 	latdeg, londeg	latitude, longitude, (decdeg)
C	
C Out:
C   rot_mat(3,3)    - Transformation matrix NED and XYZ where
C               - N is local North (colatitudes), E is longitude (East) and 
C		  D is down towards centre of (spherical) earth

CCC   site_pos(3) - XYZ coordinates of the site (m) - don't need at the moment

      implicit none

      real*8 sph_pos(3), long, colat, rot_mat(3,3),rot_mat2(3,3),
     . site_pos(3), earth_rad, pi, latdeg, londeg, lonr,latr
       parameter ( pi            = 3.1415926535897932D0 )
C Spherical earth (using IGRF10 value for radius of 6371.2km)
      parameter ( earth_rad     = 6371200.D0           )

C Debug
C      Print*, 'MODEL\vec_xyz: latdeg',latdeg, 'londeg',londeg
C Convert lat, lon to radians
      latr=latdeg*pi/180.d0
      lonr=londeg*pi/180.d0
C Debug
C      Print*,'latr',latr,'lonr',lonr
C Convert to colatitude
      colat = pi/2.d0 - latr
C Debug
C      Print*, 'colat',colat
C      sph_pos(1)=colat
C Make longitude 0-2pi radians
      if( lonr.lt.0 ) then
 	  long = 2.d0*pi + lonr
      else
          long = lonr
      end if
C Debug
C      Print*, 'long 0-2pi',long
C      sph_pos(2)=long
C      sph_pos(3) = earth_rad

C Call subroutine to obtain XYZ cartesian coordinates from lat long radius
C Get X,Y,Z coordinates
C      Call sph2xyz( colat,long,earth_rad,site_pos )

C     Now do rotation matrix (row,column) for converting from N E down to X Y Z
C Try each rotation matrix separately
C first rotation positive by pi-colat around axis(2) - east

C     North component
C      rot_mat(1,1) = dcos(pi-colat)
C      rot_mat(1,2) = 0.d0
C      rot_mat(1,3) = -dsin(pi-colat)
C     Now do EAST component
C      rot_mat(2,1) = 0.d0
C      rot_mat(2,2) =  1.d0
C      rot_mat(2,3) = 0.d0

C     Now do DOWN component
C      rot_mat(3,1) = dsin(pi-colat)
C      rot_mat(3,2) = 0.d0
C      rot_mat(3,3) = dcos(pi-colat)

C second rotation by negative longitude around down/z/axis3
C Axis(1) component
C      rot_mat2(1,1) = dcos(long)
C      rot_mat2(1,2) = -dsin(long)
C      rot_mat2(1,3) = 0.d0
C Axis2 component
C      rot_mat2(2,1) = dsin(long)
C      rot_mat2(2,2) = dcos(long)
C      rot_mat2(2,3) = 0.d0
CAxis3 component
C      rot_mat2(3,1) = 0.d0
C      rot_mat2(3,2) = 0.d0
C      rot_mat2(3,3) = 1.d0

C Both rotations in one combined matrix
C     North component
      rot_mat(1,1) = -dcos(colat)*dcos(long)
      rot_mat(1,2) = -dsin(long)
      rot_mat(1,3) = -dsin(colat)*dcos(long)

C     Now do EAST component
      rot_mat(2,1) = -dcos(colat)*dsin(long)
      rot_mat(2,2) =  dcos(long)
      rot_mat(2,3) = -dsin(colat)*dsin(long)

C     Now do DOWN component
C 
      rot_mat(3,1) = dsin(colat)
      rot_mat(3,2) = 0.d0
      rot_mat(3,3) = -dcos(colat)

C Debug
C      print*, rot_mat(1,1), rot_mat(1,2),rot_mat(1,3)
C      print*, rot_mat(2,1), rot_mat(2,2),rot_mat(2,3)
C      print*, rot_mat(3,1), rot_mat(3,2),rot_mat(3,3)
      return
      end

