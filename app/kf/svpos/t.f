CTITLE SVS_skymap
 
      subroutine svs_skymap

      implicit none 
 
*     Routine produce sky map of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].skymap

      include '../includes/const_param.h' 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4 ierr, i, j, date(5), trimlen
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
*   unit_xyz(3), unit_loc(3)  - Units vectors from site to satellite
*             in global and local frame
*   rang, hlen - Range to satellite and horizontal length of unit
*              vector
 
      real*8 t, sectag, unit_xyz(3), unit_loc(3), rang, hlen,
     .        svs_plt(3,max_sat)
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.skymap', name)
 
      open(200, file = name, iostat=ierr )
      call report_error('IOSTAT',ierr,'open',name,1,'svs_skymap')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
      do t = start, stop, step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'e')
*             Now convert the vector position of satellite into
*             topocentric frame so that we can AZ and El.
              rang = 0
              do j = 1,3
                 rang = rang + (svs_xyz(j,i)-site_xyz(j))**2
              end do
              rang = sqrt(rang)
              do j = 1,3
                 unit_xyz(j) = (svs_xyz(j,i)-site_xyz(j))/rang
              end do

*             now transform global XYZ into local NEU (loc_rot is
*             computed when the site coords are given in runstring.)
              do j = 1,3
                  call dvdot(unit_loc(j), loc_rot(j,1),3, unit_xyz,1,3)
              end do

*             Compute the horizontal length of the unit vector
              hlen = sqrt(unit_loc(1)**2+unit_loc(2)**2)
              svs_ang(1,i) = atan2(hlen, unit_loc(3))
              svs_ang(2,i) = atan2(unit_loc(2),unit_loc(1))
              svs_ang(3,i) = rang

*             Convert to plotable quanity.  North Up on the page.  This
*             way we can just plot the x y coordinates (i.e., mapping
*             to polar plot is done manually.)
              svs_plt(1,i) = min(svs_ang(1,i),pi/2)*sin(svs_ang(2,i))
              svs_plt(2,i) = min(svs_ang(1,i),pi/2)*cos(svs_ang(2,i))
              svs_plt(3,i) = svs_ang(3,i)
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)

c         write(200,200) date,
c    .        (svs_ang(2,i), min(svs_ang(1,i),pi/2)*180/pi,
c    .         svs_ang(3,i)/1.d6,i=1,num_sat)
          write(200,200) date, (svs_plt(1,i), svs_plt(2,i),
     .         svs_plt(3,i)/1.d6,i=1,num_sat)
 200      format(i5,4i3,32(3(F8.3,1x),2x))
      end do
 
      close(200)
      return
      end

