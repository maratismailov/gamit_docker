      program poles

      implicit none 

*     This program uses the frame_to_frame subroutine to compute
*     the poles of rotations of the plates.  The absolue poles in
*     NUVEL-NNR are first given, and then the relative poles of
*     of rotation.


      include '../includes/const_param.h'
 
 
*   num_plates  - Maximum number of frames supported
 
      integer*4 num_plates
 
      parameter ( num_plates = 35 )
 
* PASSED VARIABLES
 
*   rot_vec(3)  - Rotation vector to use
*   est(3)      - Values used for Lat/Long and magnitide
*   xyz_adj(3), neu_adj(3)  _ dummmy values
*   loc_coord(3)            - coLat/Long and magnitude
*   rot_mat(3,3)            - Dummy
 
      real*8 rot_vec(3), xyz_adj(3), neu_adj(3), loc_coord(3),
     .       rot_mat(3,3), est(3)
 
*   i       - Loop counter
 
      integer*4 i, j, k
 
*   frame_names(num_plates) - Names of frames supported
 
      character*8 frame_names(num_plates)
 
      data frame_names / 'PCFC    ', 'COCO    ', 'NAZC    ','CARB    ',
     .                   'SAFD    ', 'ANTA    ', 'INDI    ','AUST    ',
     .                   'AFRC    ', 'ARAB    ', 'EURA    ','NAFD    ',
     .                   'JUAN    ', 'PHIL    ', 'RIVERA  ','SCOTIA  ',
     .                   'NUV-NNR ', 'AM-02   ', 'ITRF93  ','ITRF00  ',
     .                   'GG_PCFC ', 'ANTA_I00', 'AUST_I00','EURA_I00',
     .                   'NOAM_I00', 'PCFC_I00', 'SOAM_I00', 'ANTA_I05',
     .                   'AUST_I05', 'EURA_I05', 'NOAM_I05', 'PCFC_I05',
     .                   'SOAM_I05', 'ARAB_MCC', 'ITRF05  '     / 

      write(*,100)
 100  format(/,' PLATE ROTATION POLES',/,
     .         ' --------------------')

****  Now write the results (last plate is reference)
      write(*, 200)
 200  format(/,' PLATE       Wx (deg/My)       Wy (deg/My)   ',
     .         '       Wz (deg/My)  ')
      do i = 1, num_plates 

         call frame_to_frame(frame_names(i), frame_names(num_plates),
     .                       rot_vec)
         write(*, 220) frame_names(i), (rot_vec(j)*180.d6/pi, j=1,3)
 220     format(1x,a8,1x,3(F10.6,10x))
      end do

****  Now compute the lat, long, and length of each vector
      write(*, 300)
 300  format(/,' PLATE       Lat. (deg)        Long (deg)    ',
     .         '       Mag (deg/My)')

      do i = 1, num_plates 

         call frame_to_frame(frame_names(i), frame_names(num_plates),
     .                       rot_vec)
         call rotat_crd(xyz_adj, neu_adj, 'XYZ', 'NEU',
     .                  rot_vec , loc_coord, rot_mat)

*        Now convert the results:
         est(1) = 90.d0 - loc_coord(1)*180.d0/pi
         est(2) = loc_coord(2)*180.d0/pi
         est(3) = loc_coord(3)*180.d0/pi*1.d6

         write(*, 320) frame_names(i), (est(j), j=1,3)
 320     format(1x,a8,1x,2(F10.3,10x),F10.6)
      end do


****  Now do the "differences between plates"
      if( num_plates.gt.1 ) write(*, 400)
 400  format(/,' PLATE -  PLATE        Wx (deg/my)       ',
     .         'Wy (deg/my)          Wz (deg/My)       ')

      do i = 1, num_plates
         do j = 1, num_plates 

            if( i.ne.j ) then
                call frame_to_frame(frame_names(i), frame_names(j),
     .                           rot_vec)
                write(*,420) frame_names(i), frame_names(j), 
     .                         (rot_vec(k)*180.d0/pi*1.d6, k=1,3)
 420            format(1x,a8,'-',a8,1x,3(F10.6,10x))
            end if
         end do
      end do

****  Now do the "differences between plates".  Write the results as
*     lat/long/and magnitude.  Here we re-form the differences and
*     then convert to lat, long and magnitude
      if( num_plates.gt.1 ) write(*, 500)
 500  format(/,' PLATE -  PLATE        Lat (deg)         ',
     .         'Long (deg)           Mag (deg/My)')

      do i = 1,num_plates
         do j = 1, num_plates

             if( i.ne.j ) then
                 call frame_to_frame(frame_names(i), frame_names(j),
     .                           rot_vec)
                 call rotat_crd(xyz_adj, neu_adj, 'XYZ', 'NEU',
     .                         rot_vec, loc_coord, rot_mat)

*                Now convert the results:
                 est(1) = 90.d0 - loc_coord(1)*180.d0/pi
                 est(2) = loc_coord(2)*180.d0/pi
                 est(3) = loc_coord(3)*180.d0/pi*1.d6
                 write(*,520) frame_names(i), frame_names(j), 
     .                         (est(k), k=1,3)
 520            format(1x,a8,'-',a8,1x,2(F10.3,10x),F10.6)
            end if
         end do
      end do

****  Thats all
      end 

         
