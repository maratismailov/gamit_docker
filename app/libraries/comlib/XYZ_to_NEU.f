ctitle
 
      subroutine XYZ_to_NEU(rot_mat,site_pos, loc_coord)

      implicit none
 
c
c     routine to compute the rotation matrix between global X,Y,Z
c     and the local n,e,up coordinates for a site at site_pos (XYZ
c     system)
c
c Variables
c ---------
c rot_mat -- the rotation matrix between global XYZ and local NEU systems.
c site_pos -- the XYZ coordinates of site
c loc_coord -- the local colatitude, longitude and radius
c
 
      real*8 rot_mat(3,3), site_pos(3), loc_coord(3)
 
*   i       - Loop counter
      integer*4 i
 
c
c
c.... Firstly compute the values of the loc_coord
c.... compute the colatitude
      loc_coord(1) = atan2( sqrt(site_pos(1)**2 +
     .   site_pos(2)**2), site_pos(3) )
c
c.... compute the longitude
      loc_coord(2) = atan2( site_pos(2),site_pos(1) )
c
c.... compute radius (geocentric)
      loc_coord(3) = sqrt(site_pos(1)**2 + site_pos(2)**2 +
     .   site_pos(3)**2 )
c
c.... Now compute the rotation matrix
c
c.... latitude -- north component
c
      rot_mat(1,1) = -cos(loc_coord(1))*cos(loc_coord(2))
      rot_mat(1,2) = -cos(loc_coord(1))*sin(loc_coord(2))
      rot_mat(1,3) =  sin(loc_coord(1))
c
c
c.... longitude -- east component
c
      rot_mat(2,1) = -sin(loc_coord(2))
      rot_mat(2,2) =  cos(loc_coord(2))
      rot_mat(2,3) =  0.d0
c
c.... radius    -- up component
c
      do i = 1,3
        rot_mat(3,i) = site_pos(i)/loc_coord(3)
      end do
c
      return
c
      end
 
c......................................................................
