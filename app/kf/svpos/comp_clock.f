CTITLE COMP_CLOCK
 
      subroutine comp_clock ( clk, ep )
 
*     Routine to compute the clock offset in the receiver using
*     the psuedorange data.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
* PASSED Variables
 
*   clk - Receiver clock offset computed here (sec)
 
      real*8 clk
 
*   ep  - Current epoch
 
      integer*4 ep
 
* LOCAL Variables
 
*   i,j - Loop counters
*   nn  - Tabular point immediately after data_epoch
*       - in the SP3 file.
 
      integer*4 i,j, num_av, iter
 
*   average - Average range residual (m)
*   send_epoch  - Time at which signal was transmitted.
*   pc(max_sat) - Dual frequency range (if L2 available)
*   omc(max_sat) - Observed minus computed range (m)
 
 
      real*8 average, send_epoch, pc(max_sat),
     .       omc(max_sat)
 
***** First evaluate the satellie clock corrections at this time
*     Find the SP3 entries that surrond the current data_epoch

      call comp_svs_clk( data_epoch )
       
****  Now start the light time iteration
      do j = 1, num_sat
         omc(j) = 0.d0
      end do
      
      do j = 1, num_sat
          if( omc_OK(j) ) then
 
              call eph_to_xyz( data_epoch-0.066666/86400.d0, j, 'E')
*             Compute a rough range
              p1c(j) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                       (site_xyz(2)-svs_xyz(2,j))**2+
     .                       (site_xyz(3)-svs_xyz(3,j))**2)
 
*             accumulate the average clock offset
              omc(j) =  p1o(j)-p1c(j)+svs_clk(j)*vel_light
           end if
      end do
 
*     Scan the omc values and make sure OK
      call scan_omc(omc, omc_OK, num_sat, ep )
 
*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_sat
         if( omc_OK(i) ) then
             average = average + omc(i)
             num_av  = num_av + 1
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
 
****  Now do the final range computation
 
      DO ITER = 1,2
 
          do j = 1, num_sat
              if( omc_OK(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
*                 Compute the transmission time
 
                  send_epoch = data_epoch -
     .                         (average+p1c(j))/vel_light/86400.d0
                  call eph_to_xyz( send_epoch, j, 'E')
 
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
                  omc(j) = pc(j) -p1c(j)+svs_clk(j)*vel_light
              end if
 
*             Now re-compute the average clock offset
              average = 0.d0
              num_av  = 0
              do i = 1, num_sat
                 if( omc_OK(i) ) then
                     average = average + omc(i)
                     num_av  = num_av + 1
                 end if
              end do
 
              if( num_av.gt.0 ) average = average / num_av
          end do
      END DO
 
****  Finally save the average as the station clock offset at
*     this epoch
      clk = average/vel_light
 
****  Thats all
      return
      end
 
CTITLE SCAN_OMC
 
      subroutine scan_omc( omc, omc_OK,  num_sat, ep )
 
*     Routine to make sure that all the omc's look good before
*     being used in the filter
 
      include '../includes/const_param.h'
 
* num_sat     - number of satellites
* ep          - GAMIT Epoch number
 
 
      integer*4 num_sat, ep
 
* omc(num_sat)  - Observed - computed range (m)
 
      real*8 omc(num_sat)
 
* omc_OK(num_sat) - indicates that OMC is OK.
 
 
      logical omc_OK(num_sat)

 
      real*8 average, biggest_diff

      integer*4 i,j, num, k
 
*     Scan over all combinations.  If a least one combination
*     looks good, then we keep it, else we toss the value.
*     Start the iteration loop throwing out "bad" data
      k = 99
      do while ( k.ne.0 )
          average = 0
          num = 0
          do j = 1, num_sat
                 if( omc_OK(j) ) then
                  average = average + omc(j)
                  num = num + 1
                 end if
          end do
          if( num.gt.0 ) average = average/num
*
*         Now see who is the biggest residual > 300 meters
          biggest_diff = 0.d0 
          k = 0
          do i = 1, num_sat
              if( abs(omc(i)-average).gt.300.d0 .and. omc_OK(i) ) then
                  if( abs(omc(i)-average).gt. biggest_diff ) then
                      biggest_diff = abs(omc(i)-average)
                      k = i
                  end if
              end if
          end do

*****     See if any edits
          if( k.ne.0 ) then
              omc_OK(k) = .false.
              write(*,120) ep,k, omc(k), omc(k)-average
 120          format('BAD RANGE: Epoch ',i4,' PRN',i2.2,
     .            ' OMC (m) ',f12.2,' OMC-Average (m)', f12.2)
          end if
      end do
 
***** Thats all
      return
      end
 
 
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'modear.h'

*  order - Order of the interpolating polynomial.  Must be even 
*          value

      integer*4 order

      parameter ( order = 14 )
 
*   t       - Time for computation (day number)
 
 
 
      real*8 t
 
*   i       - Satellite number
 
 
 
      integer*4 i
 
*   sys     - System for results.
 
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j, bin, ein
 
* int_xyz(3) - Inertial coordinates of satellite.  Copied back to
*              return variable

      real*8 err,  xa(100), ya(100), za(100), int_xyz(3)
 
 
      logical found
 
****  Find out which two points are closest to the epoch
      j = 0
      found = .false.
      do while ( .not.found .and. j.le.num_sp3)
          j = j + 1
          if( sp3_time(j).gt. t  ) then
              found = .true.
          end if
      end do
 
*     Get the start index for interpolation
      bin = j - order/2
      if( bin.le.0 ) bin = 1
      ein = bin + order - 1
      if( ein.gt.num_sp3 ) then
          ein = num_sp3
          bin = ein - (order-1)
      end if
 
***** Now do interporlation
      do j = 0, order-1
          xa(j+1) = sp3_xyz(1,i,bin+j)
          ya(j+1) = sp3_xyz(2,i,bin+j)
          za(j+1) = sp3_xyz(3,i,bin+j)
      end do
      call polint(sp3_time(bin), xa, order, t,
     .            svs_xyz(1,i), err)
      if( abs(err).gt.0.01 ) write(*,500) i,t, 'X', err
      call polint(sp3_time(bin), ya, order, t,
     .            svs_xyz(2,i), err)
      if( abs(err).gt.0.01 ) write(*,500) i,t, 'Y', err
 
      call polint(sp3_time(bin), za, order, t,
     .            svs_xyz(3,i), err)
      if( abs(err).gt.0.01 ) write(*,500) i,t, 'Z', err
 
  500 format(' Interpolation error PRN ',i2,' Epoch ',f12.4,
     .       ' Component ',a1,' Magnitude ', f8.4,' m')
     
***** Check to see if Inertial coordinates wanted.
      if( sys(1:1).eq.'I' ) then
          call earth_to_inert(t, svs_xyz(1,i), int_xyz, 'E','I')

*         Copy the interial coordinates to the return position
*         variable
          do j = 1, 3
             svs_xyz(j,i) = int_xyz(j)
          end do
      end if   
 
****  Thats all
      return
      end
 
CTITLE POLINT
 
      subroutine polint(xa,ya,n,x,y,dy)
 
*     Routine to do an 11 point interpolation of the SP3
*     file entries.  (From Numerical Recipes)
 
      integer*4 nmax, n, i,m, ns
 
      parameter (nmax=21)
 
      real*8 xa(n), ya(n), c(nmax),d(nmax), dif, dift,
     .    den, w, hp, ho, x,y, dy
 
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*) 'POLINT: Interpolation out of bounds'
             den = 1.d-8
          endif
           den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
 
 
CTITLE MIT_DRY
 
      subroutine mit_dry( elev, part )
 
*     Routine to compute the MIT_DRY mapping function under
*     nominal conditions.
 
      include '../includes/const_param.h'
 
      real*8 elev, part, topcon, sine, A, B, C
 
 
      real*8 beta, gamma
 
      A = 0.00125003d0
      B = 0.00312108d0
      C = 0.06945748d0
 
      sine  = sin( elev*pi/180.d0)
      beta  = B/( sine + C )
      gamma = A/( sine + beta)
      topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))
 
      part = topcon / ( sine + gamma )
 
      end
      
CTITLE comp_svs_clk

      subroutine comp_svs_clk( epoch )
      
*     Rouitne to compute satellite clock coorection at time epoch

      include 'modear.h'
 
* epoch - MJD at which clocks are needed

      real*8 epoch
      
* LOCAL VARIABLES
*   i,j - Loop counters
*   nn  - Tabular point immediately after data_epoch
*       - in the SP3 file.
 
      integer*4 i, nn
 
*   dt  - Time between the data_epoch and SP3 entry
*       - immediately before it (days)
*   dt_sp3  - Spacing between SP3 entries (days)
*   dc_sp3  - Difference in the SP3 clocks (sec)

 
      real*8 dt, dt_sp3, dc_sp3

           
      
***** First evaluate the satellie clock corrections at this time
*     Find the SP3 entries that surrond the current data_epoch
      nn = 2
      do while ( sp3_time(nn).lt. epoch .and.
     .        nn.lt.num_sp3 )
          nn = nn + 1
      end do
 
****  The two entries needed are nn-1 and nn.  Warn user if these
*     are actually outside the range of the data
      if( sp3_time(nn-1).gt.epoch .or.
     .    sp3_time(nn).lt. epoch ) then
          write(*,120) data_epoch, nn, num_sp3
 120      format('WARNING: Data epoch ',f12.4,' before or after ',
     .        'SP3 file range. Using point ',i4,' of ',i4)
      end if
 
****  Now get the interpolated satellite time
      dt = epoch - sp3_time(nn-1)
      dt_sp3 = sp3_time(nn) - sp3_time(nn-1)
      do i = 1, num_sat
          dc_sp3 = sp3_clk(i,nn) - sp3_clk(i,nn-1)
          svs_clk(i) = sp3_clk(i,nn-1) + (dc_sp3/dt_sp3)*dt
      end do
      
      
****  Thats all
      return
      end 

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
