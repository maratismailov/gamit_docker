          do j = 1, num_sat
              if( omc_OK(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
*                 Compute the transmission time

                  send_epoch = data_epoch -
     .                         (average+p1c(j))/vel_light/86400.d0
                  call eph_to_xyz( send epoch, j, 'E')

                  p1c(j) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)
                  if( p2o(j).gt.0 ) then
                      pc(j) = (p1o(j)*(77.d0/60.d0)**2-p2o(j))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc(j) = p1o(j)
                  end if

****              compute O-C
                  omc(j) = pc(j) -p1c(j)+svclk(j)*vel_light
              end if

*             Now re-compute the average clock offset
              average = 0.d0
              num_av  = 0
              do i = 1, num_chan
                 if( omc_OK(i) ) then
                     average = average + omc(i)
                     num_av  = num_av + 1
                 end if
              end do
 
              if( num_av.gt.0 ) average = average / num_av
          end do
      END DO

