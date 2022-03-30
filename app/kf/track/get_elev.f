CTITLE GET_ELEV
      subroutine get_elev( pos, svs, elev, az)

      implicit none

      include '../includes/const_param.h'

*     Modified from Herring's get_elev. G. Chen
*     Modified back to original correct form TAH 991010.

      real*8 pos(3), svs(3), elev, rot_mat(3,3), loc_coord(3)

      real*8 dpos(3),  az, udp(3)
      real*8 magdpos
      integer*4 i,j

      do i = 1,3
          dpos(i) = svs(i) - pos(i)
      end do

****  Compute the normal to the ellipsod
      call XYZ_to_GEOD(rot_mat, pos, loc_coord)

****  Transform the dpos vector in local frame
      magdpos = sqrt(dpos(1)**2+dpos(2)**2+dpos(3)**2)
      do i = 1,3
         udp(i) = 0.d0
         do j = 1, 3
            udp(i) = udp(i)+rot_mat(i,j)*dpos(j)/magdpos
         enddo
      end do

****  Now get the azimuth and elevation
      az = atan2(udp(2),udp(1))*180/pi
      if( az.lt.0 ) az = az + 360.d0
      elev = atan2(udp(3),sqrt(udp(2)**2+udp(1)**2))*180/pi


      end                                           

CTITLE GET_NADIR

      subroutine get_nadir( pos, svs, nadir)

      implicit none

*     Compute the nadir angle from satellite to ground
      include '../includes/const_param.h'

      real*8 pos(3), svs(3), nadir
      real*8 dpos(3), cnadir 
      integer*4 i
      do i = 1,3
          dpos(i) = pos(i) - svs(i)
      end do
      cnadir = -(dpos(1)*svs(1) + dpos(2)*svs(2) + dpos(3)*svs(3))/
     .          (sqrt(dpos(1)**2+dpos(2)**2+dpos(3)**2)*
     .           sqrt(svs(1)**2 +svs(2)**2 +svs(3)**2))
      nadir = acos(cnadir)*180.d0/pi

      return
      end


