 
      program svdiff

      implicit none
 
*     Program to compute the rise and set times of the GPS
*     satelites based on the ephemeris in the rinex navigiation
*     files
 
*     The runstring of the program is
*     % svdiff [navfile] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
 
      include 'svdiff.h'
 
* Main pogram variables

      integer*4 i, ep
 
*   eof     - Indicates end of file.
*   process  - Indicates that we have common data so process
*   read_data  = Set true if we should read the data file
*   read_diff  = set true if we should read difference file
 
      logical eof, process, read_data, read_diff
 
***** Get the runstring runstring
 
      call get_svsrun
 
*     Read in the ephemeris file
      call read_nav
 
*     Now loop over the data and get an estimate of the site position.
 
      call read_data_head
      call read_diff_head
 
*     Now loop over data.
      call init_est

      eof = .false.
      process = .false.
      read_data = .true.
      read_diff = .true.

      if( out_type(1:1).eq.'N' ) then
          write(*,310)
 310      format(/'*     Date          Seconds        dN (m)  +-  ',
     .            '       dE (m) +-         dU (m)  +-         ',
     .            ' dclk (m)  +-  Nsat Nusd  chi**2 SYS')
      else
          write(*,320)
 320      format(/'*     Date          Seconds        dX (m)  +-  ',
     .            '       dY (m) +-         dZ (m)  +-         ',
     .            ' dclk (m)  +-  Nsat Nusd  chi**2 SYS')
      end if

      ep = 0

      do while ( .not.eof )
          if( read_data ) call read_range(eof)
          if( read_diff .and. .not.eof ) call read_rdiff( eof )

          if( data_epoch.gt.diff_epoch+0.1/86400.d0 ) then
              read_data = .false.
              read_diff = .true.
              process = .false.
          end if
          if( data_epoch.lt.diff_epoch-0.1/86400.d0 ) then
              read_data = .true. 
              read_diff = .false.
              process = .false.
          end if
          if( abs(data_epoch -  diff_epoch).lt.0.2/86400.0d0 ) then
              read_data = .true.
              read_diff = .true.
              process = .true. 
          end if

          if( .not.eof .and. process ) then
              ep = ep + 1
              call increment_est(ep)
          end if
 
      end do
 
*     Now write out the resulys
c     call write_est
      write(*,200) site_name, (site_xyz(i)+sol_vec(i),
     .              sqrt(cov_parm(i,i)), i=1,3),
     .              sqrt(chi/nchi)
 200  format('*',/,'* Final position for ',a,' X ',
     .             F13.2,' +- ',f6.2,' Y ', 
     .             F13.2,' +- ',f6.2,' Z ',
     .             F13.2,' +- ',f6.2,/,
     .             '* Final sqrt(chi**2/n) ',f12.3)

 
****  Thats all
      end

CTITLE init_est

      subroutine init_est

      implicit none

*     Routine to initialize the Kalman filter covariance matrices

      include 'svdiff.h'

      integer*4 i,j

***** Initial the covariance matrix
*     Sites +- 1000 meters,
*     clocks +- 10000 meters

      do i = 1,4
         do j = 1,4
            cov_parm(i,j) = 0.d0
         end do
      end do

      do i = 1,3
         cov_parm(i,i) = 1.d6 
      end do
      cov_parm(4,4) = 1.d10

      do i = 1,4
         sol_vec(i) = 0.d0
      end do

      do i = 1, 3
         site_apr(i) = site_xyz(i)
      end do

***** Thats all
      return
      end

      subroutine write_mat( name, mat, nr, nc, dr, dc)

      implicit none

*  nr, nc, dr, dc - size and dimension of rows and cols

      integer*4  nr, nc, dr, dc, i,j
      real*8 mat(dr, dc)
      character*(*) name

c      write(*,*) name
      do i = 1,nr 
          j = i
c         write(*,*) i, (mat(i,j),j=1,nc)
      end do
      return
      end
         
 
CTITLE INCREMENT_EST
 
      subroutine increment_est(ep) 

      implicit none
 
*     Routine to increment the solution for the new data.
 
      include  '../includes/const_param.h'
      include 'svdiff.h'

* sec10 - minutes and seconds in tenths of seconds units.

      integer*4 i,j, k, l, num_av, ipivot(max_sat), date(5), sec10
      integer*4 ep  ! Epoch counter

      real*8 average , dummy(max_sat), scale(max_sat), sectag,
     .       dp(max_sat), dchi, svclk(max_sat), dchiout

      real*8 pc(max_sat)

      real*8 pos_xyz_fin(3), pos_xyz_adj(3), pos_neu_fin(3),
     .       pos_neu_adj(3), loc_coord(3), rot_matrix(3,3),
     .       covar(3,3), neu_covar(3,3), temp_covar(3,3)

* in_sys, out_sys - Strings with XYZ and NEU in them

      character*4 in_sys, out_sys

      data in_sys  / 'XYZ ' /, out_sys / 'NEU ' /
 
****  compute the epheremis position at the measurement
*     time
 
c     write(*,*) 'For epoch ',data_epoch
      pc = 0
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,110) 'Satellite Position at t-66.66 us', ep, 0.0
 110      format(/,'+',1x,a,' Epoch ',i5,' Site clock ',f17.2,' (m)',/,
     .           '+PRN',5x,'Xe (m)',7x,'Ye (m)',7x,'Ze (m)',
     .           5x,'SV Clock (m)')
      end if
 
      do i = 1, num_sat
 

*         Find entry
          k = 1
          do j = 1,num_ent(i)
             if( toe_jd(i,j).le.data_epoch ) k = j
          enddo
 
          call eph_to_xyz( data_epoch-0.066666/86400.d0, i, 'E')
          svclk(i) = (af0(i,k)+
     .                af1(i,k)*(data_epoch-toe_jd(i,k))*86400.d0)*
     .               vel_light
 
*         See if debug wanted:
          if( ep.ge.debug_start .and. ep.le.debug_end ) then
              write(*,120) prn(i), k, (svs_xyz(j,i),j=1,3), svclk(i)
 120          format('+ P ',i2.2, ' E ', i3, 3F13.2, F13.2)
          end if

      end do

      average = 0.0

      do i = 1, num_data
 
          omc_OK_data(i) = .false.
          do j = 1, num_sat
              if( chan_data(i).eq.prn(j) ) then
 
*                 Compute a rough range
                  p1c(i) = sqrt( (site_apr(1)-svs_xyz(1,j))**2+
     .                           (site_apr(2)-svs_xyz(2,j))**2+
     .                           (site_apr(3)-svs_xyz(3,j))**2)
c                 write(*,200) chan_data(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j)         
 200              format(2i5,10F13.3)
 
*                 accumulate the average clock offset
                  omc_data(i) =  p1o_data(i)-p1c(i)+svclk(j)
                  omc_OK_data(i) = .true.
               end if
           end do
      end do 

*     Scan the omc values and make sure OK
C     call scan_omc(omc_data, omc_OK_data, data_noise, num_data, 
C    .              num_used )

*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_data
         if( omc_OK_data(i) ) then
             average = average + omc_data(i)
             num_av  = num_av + 1
         end if
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,105) i,omc_data(i), omc_OK_data(i), num_av, 
     .              p1o_data(i)
 105         format('Data ',i2,' OMC ',f20.2,' OK ',L1,' # ',i4,
     .              ' Obs ',F20.2)
         end if
       end do
 
      if( num_av.gt.0 ) average = average / num_av
 
      do i = 1, num_data
 
          omc_OK_data(i) = .false.
          do j = 1, num_sat
              if( chan_data(i).eq.prn(j) ) then
                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
 
*                 Compute a rough range
                  p1c(i) = sqrt( (site_apr(1)-svs_xyz(1,j))**2+
     .                           (site_apr(2)-svs_xyz(2,j))**2+
     .                           (site_apr(3)-svs_xyz(3,j))**2)
c                 write(*,200) chan_data(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j)         
 
*                 accumulate the average clock offset
                  omc_data(i) =  p1o_data(i)-p1c(i)+svclk(j)
                  omc_OK_data(i) = .true.
               end if
           end do
      end do 


*     Scan the omc values and make sure OK
C     call scan_omc(omc_data, omc_OK_data, data_noise, num_data, 
C    .              num_used )

*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_data
         if( omc_OK_data(i) ) then
             average = average + omc_data(i)
             num_av  = num_av + 1
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
C     average = 0
 
****  Now do the final range computation
      do i = 1, num_data
          do j = 1, num_sat
              if( chan_data(i).eq.prn(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
                  p1c(i) = sqrt( (site_apr(1)-svs_xyz(1,j))**2+
     .                           (site_apr(2)-svs_xyz(2,j))**2+
     .                           (site_apr(3)-svs_xyz(3,j))**2)
                  if( p2o_data(i).gt.0 ) then
                      pc(i) = (p1o_data(i)*
     .                     (77.d0/60.d0)**2-p2o_data(i))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc(i) = p1o_data(i)
                  end if
                  pc(i) = p1o_data(i)
c                 write(*,200) chan(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j),
c    .                       pc    -p1c(i)+svclk(j)           
C                 print *, ' Data ', average, p1c(i), svs_xyz(1,j),
C    .                               svs_xyz(2,j), svs_xyz(3,j), 
C    .                               svclk(j)
 
****              compute O-C
                  omc_data(i) = pc(i)  -p1c(i)+svclk(j)-average 
 
****              Form the partial derivatives
                  apart_data(1,i) = (site_apr(1)-svs_xyz(1,j))/p1c(i)
                  apart_data(2,i) = (site_apr(2)-svs_xyz(2,j))/p1c(i)
                  apart_data(3,i) = (site_apr(3)-svs_xyz(3,j))/p1c(i)
                  apart_data(4,i) = 1
 
              end if
          end do
      end do
*     See if debug for o-minus-c wanted
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,140) 'Range DATA site differences', ep, average
 140      format('+ ',a,' Epoch ',i4,' Clock ',F13.2,' m',/,
     .           '+CHAN PRN',3x,'Obs (m)',6x,'Range (m)',3x,
     .           'Obs-Calc (m)',3x,'OMC Clock (m)')
       
          do i = 1, num_data
              do j = 1, num_sat
                  if( chan_data(i).eq.prn(j) ) then
                      write(*,150) i,prn(j), pc(i), p1c(i), 
     .                             omc_data(i),
     .                             omc_data(i) - average
 150                  format('+',i3,1x,i4, 4F13.2)
                  end if
              end do
          end do
      end if

***** Now do difference site:     

      do i = 1, num_sat


*         Find entry
          k = 1
          do j = 1,num_ent(i)
             if( toe_jd(i,j).le.diff_epoch ) k = j
          enddo
 
          call eph_to_xyz( diff_epoch-0.066666/86400.d0, i, 'E')
          svclk(i) = (af0(i,k)+
     .                af1(i,k)*(diff_epoch-toe_jd(i,k))*86400.d0)*
     .               vel_light
 
 
      end do
      average = 0
      do i = 1, num_diff
 
          omc_OK_diff(i) = .false.
          do j = 1, num_sat
              if( chan_diff(i).eq.prn(j) ) then
 
                  call eph_to_xyz( diff_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
*                 Compute a rough range
                  p1c(i) = sqrt( (diff_xyz(1)-svs_xyz(1,j))**2+
     .                           (diff_xyz(2)-svs_xyz(2,j))**2+
     .                           (diff_xyz(3)-svs_xyz(3,j))**2)
c                 write(*,200) chan(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j)         
 
*                 accumulate the average clock offset
                  omc_diff(i) =  p1o_diff(i)-p1c(i)+svclk(j)
                  omc_OK_diff(i) = .true.
               end if
           end do
      end do 

*     Scan the omc values and make sure OK
C     call scan_omc(omc_diff, omc_OK_diff, data_noise, num_diff, 
C    .              num_used )

*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_diff
         if( omc_OK_diff(i) ) then
             average = average + omc_diff(i)
             num_av  = num_av + 1
             if( ep.ge.debug_start .and. ep.le.debug_end ) then
                write(*,205) i,omc_diff(i), omc_OK_diff(i), num_av, 
     .                 p1o_diff(i)
 205            format('Diff ',i2,' OMC ',f20.2,' OK ',L1,' # ',i4,
     .                 ' Obs ',F20.2)
            end if
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
C     average = 0
 
****  Now do the final range computation
      do i = 1, num_diff
          do j = 1, num_sat
              if( chan_diff(i).eq.prn(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
                  call eph_to_xyz( diff_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
                  p1c(i) = sqrt( (diff_xyz(1)-svs_xyz(1,j))**2+
     .                           (diff_xyz(2)-svs_xyz(2,j))**2+
     .                           (diff_xyz(3)-svs_xyz(3,j))**2)
                  if( p2o_diff(i).gt.0 ) then
                      pc(i) = (p1o_diff(i)*(
     .                       77.d0/60.d0)**2-p2o_diff(i))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc(i) = p1o_diff(i)
                  end if
                  pc(i) = p1o_diff(i)
c                 write(*,200) chan_diff(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j),
c    .                       pc    -p1c(i)+svclk(j)           
c                 print *, ' Diff ', average, p1c(i), svs_xyz(1,j),
c    .                               svs_xyz(2,j), svs_xyz(3,j), 
c    .                               svclk(j)
 
****              compute O-C
                  omc_diff(i) = pc(i)  -p1c(i)+svclk(j)-average 
              end if
          end do
      end do

*     See if debug for o-minus-c wanted
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,240) 'Range REF sites differences', ep, average
 240      format('+ ',a,' Epoch ',i4,' Clock ',F13.2,' m',/,
     .           '+CHAN PRN',3x,'Obs (m)',6x,'Range (m)',3x,
     .           'Obs-Calc (m)',3x,'OMC Clock (m)')
       
          do i = 1, num_diff
              do j = 1, num_sat
                  if( chan_diff(i).eq.prn(j) ) then
                      write(*,250) i,prn(j), pc(i), p1c(i), 
     .                             omc_diff(i),
     .                             omc_diff(i) - average
 250                  format('+',i3,1x,i4, 4F13.2)
                  end if
              end do
          end do
      end if


***** Now try to match the channels
      k = 0
      do i = 1, num_data
         if( omc_OK_data(i) ) then
             do j = 1, num_diff
               if ( chan_data(i).eq. chan_diff(j) .and.
     .              omc_OK_diff(j)   ) then
                  k = k + 1
                  omc(k) = omc_data(i)- omc_diff(j)
                  omc_OK(k) = .true.
C                 print *,' i,j, omc', i,j, k, omc_data(i),
C    .                      omc_diff(j), omc(k)
                  do l = 1,4
                     apart(l,k) = apart_data(l,i)
                  end do
               end if
             end do
          end if
      end do

      num_chan = k
      call scan_omc(omc, omc_OK, data_noise, num_chan, 
     .              num_used )


****  Now increment the filter
      do i = 1,4
          cov_parm(i,i) = cov_parm(i,i) + wn(i)
      end do
 
***** Start computing the Kalman gain.  Get the acat matrix
      call acat(apart, cov_parm, 4, num_chan, acat_mat, temp_gain,
     .         max_sat)
 
****  Now add the data noise (assume 10 meter data noise)
      num_used = 0
      do i = 1, num_chan
          if( omc_OK(i) ) then
              num_used = num_used + 1
              acat_mat(i,i) = acat_mat(i,i) + data_noise**2 
          else
              acat_mat(i,i) = acat_mat(i,i) + 1.d30
          end if
      end do
 
***** Now invert this results
      call invert_vis(acat_mat, dummy, scale, ipivot, 
     .                num_chan, max_sat, 0)


***** Now form the Kalman Gain
      do i = 1,4
          do j = 1, num_chan
              kgain(i,j) = 0
              do k = 1, num_chan
                  kgain(i,j) = kgain(i,j) +
     .                    temp_gain(i,k)*acat_mat(k,j)
              end do
          end do
      end do

* MOD TAH 030120: Added call to dummy routine.  Write statments
*     here are commented out, but call seems to be needed with 
*     GNU Fortran (GCC 3.2.1) 3.2.1 20021119 (release) or else
*     kgain is corrupt.
      call write_mat('kgain', kgain, 4, num_chan, 4, max_sat)

***** Update the parameter estimates
      do j = 1,num_chan
          dp(j) = omc(j) 
          do i = 1,4
              dp(j) = dp(j) -apart(i,j)*sol_vec(i)
          end do
      end do
      do i = 1,4
          dx(i) = 0.d0
          do j = 1, num_chan 
              dx(i) = dx(i) + kgain(i,j)*dp(j)
          end do
      end do

*     Start incrementing the prefit chi**2
      dchi = 0
      do i = 1, num_chan
         do j = 1, num_chan
             dchi = dchi + dp(i)*acat_mat(i,j)*dp(j)
         end do
      end do

      if( num_used.gt.0 ) dchi = dchi / num_used 
* MOD TAH 990628: Skip the rest if dchi too large
      if( dchi.gt.100 ) then
          write(*,220) dchi
 220      format('+ dChi too large, skipping.  Value ',d16.3)
          RETURN
      end if

      chi = chi + dchi*num_used
      nchi = nchi + num_used
      
      do i = 1,4
          sol_vec(i) = sol_vec(i) + dx(i)
      end do
 
***** Now update covariance matrix
      do i = 1,4
          do j = 1,4
              do k = 1, num_chan
                  cov_parm(i,j) = cov_parm(i,j) -
     .                     kgain(i,k)*temp_gain(j,k)
              end do
          end do
      end do
 
****  Write out results
*     See if we need to output
      dchiout = dchi
      if( dchiout.gt.99999.0 ) dchiout = 999.99
      call jd_to_ymdhms(data_epoch, date, sectag)

      sec10 = nint(sectag*10.0)+date(5)*600
      if( sec10-nint(sec10/out_spacing)*out_spacing.eq.0 ) then

*****     Need to output.  See what coordinates
          do i = 1,3
             pos_xyz_fin(i) = site_apr(i) + sol_vec(i)
             pos_xyz_adj(i) = pos_xyz_fin(i) - site_xyz(i)
             do j = 1,3
                covar(i,j) = cov_parm(i,j)
             end do
          end do

          if( out_type(1:1).eq.'N') then
              call rotate_geod(pos_xyz_adj,pos_neu_adj,in_sys,out_sys,
     .            pos_xyz_fin,loc_coord,rot_matrix)
c
c....         convert the local coordinates to adjustments
              call loc_to_geod( loc_coord, pos_neu_fin )
c
c....         Now compute the sigmas of the local coordinates. Firstly
c             save the covariance elements
              call var_comp(rot_matrix,covar,NEU_covar,temp_covar, 
     .                      3,3,1)

              write(*,300) date, sectag, (pos_neu_adj(i),
     ,               sqrt(NEU_covar(i,i)),i=1,3),
     .               sol_vec(4), sqrt(abs(cov_parm(4,4))),
     .               num_used, num_chan, dchiout
 300          format(i5,4i3,1x,f8.3,1x,3(f10.2,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2, 2I3,1x, f8.2,' NEU')
           else
C             write(*,310) date, sectag, (sol_vec(i), 
              write(*,310) date, sectag, (pos_xyz_fin(i)-site_xyz(i),
     .               sqrt(cov_parm(i,i)), i = 1,3), 
     .               sol_vec(4), sqrt(abs(cov_parm(4,4))),
     .               num_used, num_chan, dchiout
 310          format(i5,4i3,1x,f8.3,1x,3(f10.2,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2, 2I3,1x, f8.2,' XYZ')
           end if
      end if

*     Update the apriori so we start with these coordinates
      do i = 1,3

*         If the site is stochastic then update the site position
*         and do not zero sol_vec.  (By not zero'ing sol_vec we are
*         assuming a constant velocity model).
          if( wn(i).ne.0 ) then
              site_apr(i) = site_apr(i) + sol_vec(i)
*             Now re-zero sol_vec since the the new apriori will
*             be used in the filter (assumes zero velocity)
              sol_vec(i) = 0.d0             
          end if
      end do 

****  Thats all
      return
      end

CTITLE SCAN_OMC

      subroutine scan_omc( omc, omc_OK, data_noise, num_chan, num_used ) 

      implicit none

*     Routine to make sure that all the omc's look good before
*     being used in the filter

      include '../includes/const_param.h'

* num_chan    - number of channels at this epoch
* num_used    - number  of channels which seem to good data
     
      integer*4 num_chan, num_used

* omc(num_chan)  - Observed - computed range (m)
* data_noise     - (m)

      real*8 omc(num_chan), data_noise 

* omc_OK(num_chan) - indicates that OMC is OK.

      logical omc_OK(num_chan)

      integer*4 i,j

*     Scan over all combinations.  If a least one combination
*     looks good, then we keep it, else we toss the value.
      do j = 2, num_chan
             omc(j) = omc(j) - nint((omc(j)-omc(1))/
     .                  (vel_light*1.d-3))*vel_light*1.d-3
      end do
      do i = 1, num_chan
         omc_OK(i) = .false.
         do j = 1, num_chan
            if( i.ne.j ) then    ! Check size of differnce 
                if( abs(omc(i)-omc(j)).lt.100*data_noise ) then
                    omc_OK(i) = .true.
                end if
            end if
         end do
      end do

*     See if millisecond biases
      do i = 1, num_chan
         if( .not.omc_ok(i) ) then

*            look for jump
             do j = 1, num_chan
                if( i.ne.j .and. omc_ok(j) ) then
                    
                    omc(i) = omc(i) - nint((omc(i)-omc(j))/
     .                         (vel_light*1.d-3))*vel_light*1.d-3
                end if
             end do
         end if
      end do

    

*     See if we found bad data
      num_used = 0
      do i = 1, num_chan
         if( omc_ok(i) ) num_used = num_used + 1
      end do

*     If we lost data (Output line)
      if( num_used.ne.num_chan ) then
          write(*,100) (i, omc_OK(i), omc(i), i=1,num_chan)
 100      format('*',100(I3,1x,L1,1x,F10.2))
      end if

***** Thats all
      return
      end

CTITLE ACAT
 
      subroutine acat(a,c,np,nc, acat_mat, temp_gain, max_sat)

      implicit none
 
*     Routine to compute a*c*at matrix
 
*   np      - Number of parameters (assummed first dimension of
*           - a)
*   nc      - number of channels
*   max_sat - Maxiumum number of channels possible
 
      integer*4 np, nc, max_sat
 
*   a(np,nc)    - Partials matix
*   c(np,np)    - covariance matrix
*   temp_gain(np, max_sat)  - contains a*c
 
      real*8 a(np,nc), c(np,np), acat_mat(max_sat,max_sat),
     .    temp_gain(np, max_sat)
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
***** First compute temp_gain
      do i = 1, np
          do j = 1, nc
              temp_gain(i,j) = 0.d0
              do k = 1, np
                  temp_gain(i,j) = temp_gain(i,j) + a(k,j)*c(i,k)
              end do
          end do
      end do
 
***** Now complete acat
      do i = 1, nc
          do j = 1, nc
              acat_mat(i,j) = 0.d0
              do k = 1, np
                  acat_mat(i,j) = acat_mat(i,j) + temp_gain(k,i)*a(k,j)
              end do
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE READ_RANGE
 
      subroutine read_range(eof)

      implicit none
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include 'svdiff.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(max_data_types), ierr, 
     .          jerr,i,j, k,id, trimlen,n1,nel
 
      integer*4 nblk   ! Number of lines to read
     .,         nend, nstr ! End and Start num_data entries
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(max_data_types)
 
*   eof     - Indicates end of file
*   noP2    - Set true if no P2 data available

      logical eof, noP2
 
      character*256 line
      character*1 cr
      character*1 stype(max_sat)  ! Type of satellite.  Only G GPS used

* MOD TAH 151021: Replaced code with more general version that 
*     reads any number of channels.  We need to read first line to 
*     get number of channels
      read(101,'(a)', iostat=ierr) line
      call sub_char(line,cr,' ')
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if

*     See if comment block
      id = 2
      do while ( id.eq.2 ) 
          read(line,'(26x,I3,I3)') id, num_data
          if( id.ge.2 .and. id.le.4 ) then
             do i = 1, num_data
                read(101,'(a)', iostat=ierr) line
             end do
             id = 2
             read(101,'(a)', iostat=ierr) line
          endif
      enddo 

      read(line,100,iostat=ierr) date, sectag, id, num_data,
     .            (stype(i),chan_data(i),i=1,min(12,num_data))
 100  format(5i3,f11.7,i3,i3,12(a1,i2))

      if( date(1).gt.50 ) then
          date(1) = date(1) + 1900
      else
          date(1) = date(1) + 2000
      endif
      call ymdhms_to_jd( date, sectag, data_epoch)
    
*     Now read the remainder of the channels
      if( num_data.gt.max_sat ) then 
          call report_stat('FATAL','SVPOS','Too many channels',
     .       data_file,'Too many satellites ', max_sat)
      endif
      nblk = (num_data-1)/12 + 1
      do i = 2, nblk
         read(101,'(a)', iostat=ierr) line
         call sub_char(line,cr,' ')
         if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
             eof = .true.
             RETURN
         end if
         nstr = (i-1)*12+1
         nend = nstr+min(12,num_data-(i-1)*12)
         read(line,110,iostat=ierr) 
     .            (stype(j),chan_data(j),j=nstr,nend)
 110     format(32x,12(a1,i2))

      end do

      do i = 1,num_data
         if( stype(i).eq.' ' ) stype(i) = 'G'
      enddo

 
*     Now loop over the data records.
      noP2 = .true.
      do i = 1, num_data  
* MOD TAH 981113: Changed to the IOSTAT variable to jerr and
*         set the number of data to 0 if error reading.  This is
*         to handle problems marker names inside the data record.
* MOD TAH 130330: Loop reading all data 
* MOD RWK 150514 to fix logic 
          n1 = int((num_data_types-1)/5)+1 
          do k = 1,n1
              read(101,'(a)', iostat=ierr) line
              call sub_char(line,cr,' ')
              nel = min0(5*k,num_data_types)
              read(line,120, iostat=jerr) (vals(j), flags(j), 
     .                j =(k-1)*5+1,nel)
 120          format( 5(f14.3,1x,i1))
          end do
c         write(*,900) I, stype(i), chan_data(i),
c    .                (data_types(j), vals(j), flags(j),
c    .                    j = 1, num_data_types)
c900      format('DA ',i3,1x,a2,1x,I2.2,1x,20(a2,1x,F14.3,1x,I3,1x))
*         Now assign the phase and range measurements
          do j = 1,num_data_types
              if((data_types(j).eq.'C1' .or.
     .            data_types(j).eq.'P1') .and.vals(j).ne.0  )
     .                   p1o_data(i) = vals(j)
              if( data_types(j).eq.'L1' ) l1o_data(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( data_types(j).eq.'L2' ) l2o_data(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( data_types(j).eq.'P2'.and.vals(j).ne.0  ) then
                  p2o_data(i) = vals(j)
                  noP2 = .false.
              endif
              if( data_types(j)(1:2).eq.'C2' .and.vals(j).ne.0 ) then
                  P2o_data(i) = vals(j)
                  noP2 = .false.
              end if
          end do
      end do

*     If there is not P2 data, set the values the same as P1
      if ( noP2 ) then
         do i = 1 ,num_data
            P2o_data(i) = P1o_data(i)
         end do
      end if

* MOD TAH 990624: Now do some simple checking of the data.
      i = 0
      k = 0
      do while ( i.lt.num_data )
          i = i + 1
          k = k + 1
          if( stype(i).ne.'g' .and. stype(i).ne.'G') then
c             write(*,210) date, k, stype(i)
 210          format('DELETING REF DATA AT ',i4,4(1x,i2),
     .               ' Chan ',i3,' Type ',a1)
              do j = i + 1, num_data
                 chan_data(j-1) = chan_data(j)
                 p1o_data(j-1) = p1o_data(j)
                 p2o_data(j-1) = p2o_data(j)
                 stype(j-1) = stype(j)
              end do
              i = i - 1
              num_data = num_data - 1
          else if( abs(p1o_data(i)-p2o_data(i)).gt.200.d0) then
              write(*,220) date, k, abs(p1o_data(i)-p2o_data(i))
 220          format('DELETING REF  DATA AT ',i4,4(1x,i2),
     .               ' Chan ',i3,' P1-P2= ',f12.1,' m different')
* MOD TAH 050705: Move done the data stack
              do j = i + 1, num_data
                 chan_data(j-1) = chan_data(j)
                 p1o_data(j-1) = p1o_data(j)
                 p2o_data(j-1) = p2o_data(j)
                 stype(j-1) = stype(j)
              end do
              i = i - 1
              num_data = num_data - 1
          end if
      end do

* MOD TAH 981113: Changed error return to set num_data equal 0
*     if error.  This handles the case of marker names in files.
      if( jerr.ne.0 ) num_data = 0
 
****  Thats all
      return
      end
 
CTITLE READ_rdiff  
 
      subroutine read_rdiff(eof)

      implicit none
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include 'svdiff.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(max_data_types), ierr, 
     .          jerr,i,j, k,id, trimlen,n1,nel
 
      integer*4 nblk   ! Number of lines to read
     .,         nend, nstr ! End and Start num_chan entries
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(max_data_types)
 
*   eof     - Indicates end of file
*   noP2    - Set true if no P2 data available

      logical eof, noP2
 
      character*256 line
      character*1 cr
      character*1 stype(max_sat)  ! Type of satellite.  Only G GPS used

* MOD TAH 151021: Replaced code with more general version that 
*     reads any number of channels.  We need to read first line to 
*     get number of channels
      read(102,'(a)', iostat=ierr) line
      call sub_char(line,cr,' ')
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if

*     See if comment block
      id = 2
      do while ( id.eq.2 ) 
          read(line,'(26x,I3,I3)') id, num_diff
          if( id.ge.2 .and. id.le.4 ) then
             do i = 1, num_diff
                read(102,'(a)', iostat=ierr) line
             end do
             id = 2
             read(102,'(a)', iostat=ierr) line
          endif
      enddo 
      read(line,100,iostat=ierr) date, sectag, id, num_diff,
     .            (stype(i),chan_diff(i),i=1,min(12,num_diff))
 100  format(5i3,f11.7,i3,i3,12(a1,i2))

      if( date(1).gt.50 ) then
          date(1) = date(1) + 1900
      else
          date(1) = date(1) + 2000
      endif
      call ymdhms_to_jd( date, sectag, diff_epoch)
    
*     Now read the remainder of the channels
      if( num_diff.gt.max_sat ) then 
          call report_stat('FATAL','SVPOS','Too many channels',
     .       data_file,'Too many satellites ', max_sat)
      endif
      nblk = (num_diff-1)/12 + 1
      do i = 2, nblk
         read(102,'(a)', iostat=ierr) line
         call sub_char(line,cr,' ')
         if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
             eof = .true.
             RETURN
         end if
         nstr = (i-1)*12+1
         nend = nstr+min(12,num_diff-(i-1)*12)
         read(line,110,iostat=ierr) 
     .            (stype(j),chan_diff(j),j=nstr,nend)
 110     format(32x,12(a1,i2))

      end do

      do i = 1,num_diff
         if( stype(i).eq.' ' ) stype(i) = 'G'
      enddo

 
*     Now loop over the data records.
      noP2 = .true.
      do i = 1, num_diff  
* MOD TAH 981113: Changed to the IOSTAT variable to jerr and
*         set the number of data to 0 if error reading.  This is
*         to handle problems marker names inside the data record.
* MOD TAH 130330: Loop reading all data 
* MOD RWK 150514 to fix logic 
          n1 = int((num_diff_types-1)/5)+1 
          do k = 1,n1
              read(102,'(a)', iostat=ierr) line
              call sub_char(line,cr,' ')
              nel = min0(5*k,num_diff_types)
              read(line,120, iostat=jerr) (vals(j), flags(j), 
     .                j =(k-1)*5+1,nel)
 120          format( 5(f14.3,1x,i1))
          end do
c         write(*,900) I, stype(i), chan_diff(i),
c    .                (data_types(j), vals(j), flags(j),
c    .                    j = 1, num_diff_types)
c900      format('DA ',i3,1x,a2,1x,I2.2,1x,20(a2,1x,F14.3,1x,I3,1x))
*         Now assign the phase and range measurements
          do j = 1,num_diff_types
              if((diff_types(j).eq.'C1' .or.
     .            diff_types(j).eq.'P1') .and.vals(j).ne.0  )
     .                   p1o_diff(i) = vals(j)
              if( diff_types(j).eq.'L1' ) l1o_diff(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( diff_types(j).eq.'L2' ) l2o_diff(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( diff_types(j).eq.'P2'.and.vals(j).ne.0  ) then
                  p2o_diff(i) = vals(j)
                  noP2 = .false.
              endif
              if( diff_types(j)(1:2).eq.'C2' .and.vals(j).ne.0 ) then
                  P2o_diff(i) = vals(j)
                  noP2 = .false.
              end if
          end do
      end do

*     If there is not P2 data, set the values the same as P1
      if ( noP2 ) then
         do i = 1 ,num_diff
            P2o_diff(i) = P1o_diff(i)
         end do
      end if

* MOD TAH 990624: Now do some simple checking of the data.
      i = 0
      k = 0
      do while ( i.lt.num_diff )
          i = i + 1
          k = k + 1
          if( stype(i).ne.'g' .and. stype(i).ne.'G') then
c             write(*,210) date, k, stype(i)
 210          format('DELETING REF DATA AT ',i4,4(1x,i2),
     .               ' Chan ',i3,' Type ',a1)
              do j = i + 1, num_diff
                 chan_diff(j-1) = chan_diff(j)
                 p1o_diff(j-1) = p1o_diff(j)
                 p2o_diff(j-1) = p2o_diff(j)
                 stype(j-1) = stype(j)
              end do
              i = i - 1
              num_diff = num_diff - 1
          else if( abs(p1o_diff(i)-p2o_diff(i)).gt.200.d0) then
              write(*,220) date, k, abs(p1o_diff(i)-p2o_diff(i))
 220          format('DELETING REF  DATA AT ',i4,4(1x,i2),
     .               ' Chan ',i3,' P1-P2= ',f12.1,' m different')
* MOD TAH 050705: Move done the data stack
              do j = i + 1, num_diff
                 chan_diff(j-1) = chan_diff(j)
                 p1o_diff(j-1) = p1o_diff(j)
                 p2o_diff(j-1) = p2o_diff(j)
                 stype(j-1) = stype(j)
              end do
              i = i - 1
              num_diff = num_diff - 1
          end if
      end do


* MOD TAH 981113: Changed error return to set num_diff equal 0
*     if error.  This handles the case of marker names in files.
      if( jerr.ne.0 ) num_diff = 0
 
****  Thats all
      return
      end

CTITLE READ_DATA_HEAD
 
      subroutine read_data_head
 
      implicit none

*     Routine to read the header from the data files.
*     Returns will be:
*     Approximate site position
*     station name
 
      include 'svdiff.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i

      integer*4 blks   ! Numnber of data type blocks based on num_data_type
     .,         it     ! Counter reading through blocks
     .,         str, rem   ! Start element (9 per line) and
                       ! remaining number on line.
 
*   eoh     - End of header flags
*   rinex   - True if rinex file
 
      logical eoh, rinex
 
*   line    - Line read from file
 
 
      character*256 line
      character*1 cr
 
****  Open the data file
      cr = char(13) 
      open(101,file=data_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',data_file, 1,
     .        'read_data_file/svdiff')
 
*     Start looping over the header
      eoh = .false.
      rinex = .false.
      do while ( .not.eoh )
          read(101,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          call report_error('IOSTAT',ierr,'read', data_file,1,
     .                'End of hedader not found')
          if( trimlen(line).lt.10) eoh = .true.
          if( index(line,'END OF HEAD').gt.0 ) eoh = .true.

*         Get Rinex version
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*) data_rxver
              write(*,120) trim(data_file), data_rxver
 120          format('Data File ',a,' Rinex Version ',F5.2)
              if( data_rxver.ge.3.0 ) then
                 call report_stat('FATAL','SVPOS',
     .                'RINEX 3 not suppported yet',data_file,'',0)
              end if
          end if
 
*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.eq.0 ) then
              indx = index(line,'O')
              if ( indx.ne.21 ) indx = 0
          end if 
          if( indx.gt.0 ) rinex = .true.
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.59 ) site_name = line(1:8)
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.59 .and. site_xyz(1).eq.0 ) then
              read(line,*) site_xyz
          end if
*         See if data types
          indx = index(line,'TYPES OF OBS')
          if( indx .gt.59 ) then
              read(line,'(i6,9a6)', iostat=ierr) num_data_types, 
     .               (data_types(i), i=1,min(9,num_data_types))
* MOD TAH 130330: Allow more data types
              if( num_data_types.gt.max_data_types) then
                  call report_stat('FATAL','SVPOS','Read TYPES OF OBS',
     .               data_file,'Too many types of data',
     .               max_data_types)
              end if
              if( num_data_types.gt.9 ) then
* MOD TAH 151020: Loop over the remaining blocks needed
                  blks = (num_data_types-1)/9+1  ! Includes count for values
                                                 ! already read
                  do it = 2,blks
                     rem = min(9,num_data_types-(it-1)*9)
                     str = (it-1)*9
                     read(101,'(a)', iostat=ierr) line
                     read(line,'(6x,9a6)', iostat=ierr)  
     .                  (data_types(str+i), i=1,rem)
                  enddo
              endif


          end if
      end do
      
      call sub_char( site_name, ' ','_' )
      
*     See if we found rinex file
      if( rinex ) then
          write(*,150) data_file(1:trimlen(data_file)),
     .        site_name, site_xyz
 150      format(/'* For RINEX data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) data_file(1:trimlen(data_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'svdiff: Wrong type data file'
      end if

      write(*,210) num_data_types, (data_types(i), i=1,num_data_types)
 210  format('There are ',i3,' data types: ',100(1x,a2))
 
****  Thats all
      return
      end
 
 
CTITLE READ_diff_HEAD
 
      subroutine read_diff_head

      implicit none
 
*     Routine to read the header from the data files.
*     Returns will be:
*     Approximate site position
*     station name
 
      include 'svdiff.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i
 
      integer*4 blks   ! Numnber of data type blocks based on num_data_type
     .,         it     ! Counter reading through blocks
     .,         str, rem   ! Start element (9 per line) and
                       ! remaining number on line.
 
*   eoh     - End of header flags
*   rinex   - True if rinex file
 
      logical eoh, rinex
 
*   line    - Line read from file
 
 
      character*256 line
      character*1 cr
 
****  Open the data file
      cr = char(13) 
      open(102,file=diff_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',diff_file, 1,
     .        'read_diff_file/svdiff')
 
*     Start looping over the header
      eoh = .false.
      rinex = .false.
      do while ( .not.eoh )
          read(102,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          call report_error('IOSTAT',ierr,'read', diff_file,1,
     .                'End of hedader not found')
          if( trimlen(line).lt.10) eoh = .true.
          if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
 
*         Get Rinex version
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*) diff_rxver
              write(*,120) trim(diff_file), diff_rxver
 120          format('Diff File ',a,' Rinex Version ',F5.2)
              if( diff_rxver.ge.3.0 ) then
                 call report_stat('FATAL','SVPOS',
     .                'RINEX 3 not suppported yet',diff_file,'',0)
              end if
          end if

*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.eq.0 ) then
              indx = index(line,'O')
              if ( indx.ne.21 ) indx = 0
          end if 
          if( indx.gt.0 ) rinex = .true.
 
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.59 ) diff_name = line(1:8)
 
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.59 .and. diff_xyz(1).eq.0 ) then
              read(line,*) diff_xyz
          end if
 
*         See if data types
          indx = index(line,'TYPES OF OBS')
          if( indx .gt.59 ) then
              read(line,'(i6,9a6)', iostat=ierr) num_diff_types, 
     .               (diff_types(i), i=1,min(9,num_diff_types))
* MOD TAH 130330: Allow more data types
              if( num_diff_types.gt.max_data_types) then
                  call report_stat('FATAL','SVPOS','Read TYPES OF OBS',
     .               diff_file,'Too many types of data',
     .               max_data_types)
              end if
              if( num_diff_types.gt.9 ) then
* MOD TAH 151020: Loop over the remaining blocks needed
                 blks = (num_diff_types-1)/9+1  ! Includes count for values
                                                ! already read
                 do it = 2,blks
                    rem = min(9,num_diff_types-(it-1)*9)
                    str = (it-1)*9
                    read(102,'(a)', iostat=ierr) line
                    read(line,'(6x,9a6)', iostat=ierr)  
     .                 (diff_types(str+i), i=1,rem)
                 enddo
              endif
         end if
      end do

*     See if we found rinex file
      if( rinex ) then
          write(*,150) diff_file(1:trimlen(diff_file)),
     .        diff_name, diff_xyz
 150      format(/'* For RINEX data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) diff_file(1:trimlen(diff_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'svdiff: Wrong type data file'
      end if
 
      write(*,210) num_diff_types, (diff_types(i), i=1,num_diff_types)
 210  format('There are ',i3,' data types: ',100(1x,a2))
 

****  Thats all
      return
      end
 
 
CTITLE GET_SVSRUN
 
      subroutine get_svsrun

      implicit none
 
*     Routine to get the cunstring.
 
      include '../includes/const_param.h'
      include 'svdiff.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
 
 
      integer*4 len_run, i, rcpar
 
*  runstring   - Elements of runstring
 
      character*256 runstring
 
*     Get the first runstring parameter
      len_run = rcpar(1,nav_file)
      if( len_run.le.0 ) then
          call proper_runstring('svdiff.hlp','svdiff/nav file',1)
      end if
 
*     See if data difle names passed
      len_run = rcpar(2,data_file)
      if( len_run.le.0 ) then
          call proper_runstring('svdiff.hlp','svdiff/data file',1)
      end if
*     See if diff difle names passed
      len_run = rcpar(3,diff_file)
      if( len_run.le.0 ) then
          call proper_runstring('svdiff.hlp','svdiff/diff file',1)
      end if


*     Get process noise for the position
      len_run = rcpar(4, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) wn(1)
          wn(2) = wn(1)
          wn(3) = wn(1)
      else
          do i = 1,3
             wn(i) = 0.d0
          end do
      end if

*     Get clock process noise
      len_run = rcpar(5, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) wn(4)
      else
          wn(4) = 1.d6
      end if

*     Get data noise
      len_run = rcpar(6, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) data_noise 
      else
          data_noise = 100.d0
      end if

*     Now see if reference positions passed 
      do i = 1,3
         len_run = rcpar(6+i,runstring)
         if(  len_run.gt.0 ) then
             read(runstring,*) diff_xyz(i)
         else
             diff_xyz(i) = 0.d0
         end if
      end do

*     Get data spacing (results output when modulo this value
*     past the hourt
      len_run = rcpar(10, runstring)
      if( len_run.gt.0 ) then
           read(runstring,*) out_spacing
      else
           out_spacing = 0.1
      end if

*     Get output coordinate type
      len_run = rcpar(11, out_type )
      if( len_run.eq.0 ) out_type = 'XYZ'   

* MOD SCM Added apriori site coordinates on command line. 6/1/00
*     Now see if new site positions passed 
      do i = 1,3
         len_run = rcpar(11+i,runstring)
         if(  len_run.gt.0 ) then
             read(runstring,*) site_xyz(i)
         else
             site_xyz(i) = 0.d0
         end if
      end do
*     See if debug wanted
      len_run= rcpar(15, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) debug_start
      else
          debug_start = 0   
      end if

*     See if debug or processing start (+ve if debug)
      if( debug_start.lt.0 ) then
          proc_start = -debug_start 
          debug_start = 0
      else
          proc_start = 0 
      endif

*     See if debug wanted
      len_run= rcpar(16, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*) debug_end   
      else
*         See if end of processing or number of epochs
          if( debug_end.lt.0 ) then
              proc_end = -debug_end
              debug_end = 0
          else
              if( proc_start.gt.0 ) then
                 proc_end = proc_start + debug_end
                 debug_end = 0
              endif
          endif
      end if
      if( proc_start.gt.0 ) then
         write(*,310) proc_start, proc_end
 310     format('Processing data between epoch ',i6, 
     .          ' and ',i6) 
      endif

       if( debug_start.gt.0 ) then
         write(*,320) debug_start, debug_end
 320     format('Debug output between epoch ',i6, 
     .          ' and ',i6) 
      endif

 
****  Thats all
      return
      end


CTITLE READ_NAV
 
      subroutine read_nav

      implicit none
 
*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero)
 
      include 'svdiff.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
 
 
      integer*4 ierr, i, j, k, rp, trimlen, ll, date(5)
      real*4 rinex_version
 
*   sectag      - Seconds tag in date
 
       real*8 sectag
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
*   have            - Set true if we already have this satellite
 
       logical still_header, have
 
*   line            - line read from file
 
       character*256 line

*   cr - Carriage return (for handling dos files)
      character*1 cr

      cr = char(13)
 
****  Open the nav_file
      open(100, file=nav_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',nav_file,1,
     .            'svpos/read_nav')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read
 
      do i = 1, max_sat
          prn(i) = 0
          num_ent(i) = 0
      end do
      rinex_version = 1
 
      still_header = .true.
      do while ( still_header )
          read(100,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          ll = trimlen(line)
          if( ll.eq.0 .or. index(line,'END OF HEAD').gt.0 ) then
              still_header = .false.
          end if
          call report_error('IOSTAT',ierr,'read','NAV FILE HEADER',
     .                      1,'read_nav')
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*, iostat=ierr) rinex_version
              write(*,120) rinex_version
 120          format('* Rinex version ',f3.1,' Nav file found')
          end if
      end do
 
      num_sat = 0
 
*     Now start reading the entries
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          call sub_char(line,cr,' ')
          if( trimlen(line).eq.0 ) ierr = -1
*                                     ! See which PRN
          if( ierr.eq.0 ) then
              read(line,*) rp
 
*             See if we already have this prn
              have = .false.
              do i = 1, num_sat
                  if( rp.eq.prn(i) ) then
                      have = .true.
                      j = num_ent(i) + 1
                      k = i
                  end if
              end do
 
*                                         ! This is first
              if( .not.have ) then
                  num_sat = num_sat + 1
                  rp = num_sat
                  num_ent(rp) = 1
                  j = 1
              else
                  rp = k
                  num_ent(rp) = num_ent(rp) + 1
                  j = num_ent(rp)
              endif
 
              read(line,200) prn(rp), date, sectag,
     .                    af0(rp,j),af1(rp,j), af2(rp,j)
 200          format(i2,5i3,f5.1,3d19.8)
 
              call ymdhms_to_jd( date, sectag, toe_jd(rp,j))
 
*             Read the rest of the entires
              read(100,210) aode(rp,j), crs(rp,j), dn(rp,j), m0(rp,j)
              read(100,210) cuc(rp,j), ecc(rp,j), cus(rp,j), art(rp,j)
              read(100,210) toe(rp,j), cic(rp,j), om0(rp,j), cis(rp,j)
              read(100,210) i0(rp,j) , crc(rp,j), w(rp,j)  , omd(rp,j)
              read(100,210) idt(rp,j), cflg12(rp,j), weekno, 
     .                      pflg12(rp,j)
              read(100,210) svacc, svhealth, tgd, aodc(rp,j)
 210          format(3x,4d19.8)
              if( rinex_version.gt.1 ) read(100,'(a)') line
          end if
      end do
 
      write(*,300) num_sat, nav_file(1:trimlen(nav_file))
 300  format('* ', i5,' satellites found in ',a)
      do j = 1,num_sat
          write(*,310) j,prn(j), num_ent(j), 
     .                 toe_jd(j,1),toe_jd(j,num_ent(j))
 310      format('SV ',i2,' PRN ',i2.2,' Num_ent ',i3,' JDs ',2F14.4)
      enddo
 
      if( num_sat.eq.0 ) stop ' SVPOS: No satellites found'
 
      return
      end

CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'svdiff.h'
 
*   t       - Time for computation (day number)
 
 
      real*8 t
 
*   i       - Satellite number
 
 
      integer*4 i
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
*   k           - Index to time 
 
      integer*4 j, k
 
*   gm      - GM
*   eom     - Earth rotation rate (rads/sec)
*   a       - Semimajor axis
*   n0      - Mean motion
*   tk      - Time of epoch from toe (seconds)
*   n       - COrrected mean motion
*   mk      - Mean anomaly
*   ek      - Eccentric anomaly
*   vk      - true anomaly
*   sinvk, cosvk    - Sin and cos of true anomaly
*   pk      - argument of latitude
*   duk, drk, dik   - Coorections to arg of lat, radius and
*           - inclinations
*   uk      - Argumenr of latitude
*   rk      - radius at time tk
*   ik      - Inclination at tk
*   xpk, ypk    - Inplane coordiantes
*   omk     - Longitude of the asecdning node
*   rot_mat(3,3)    - Rotation matrix from XYZ to NEU
 
 
 
 
      real*8 gm, eom, a, n0, tk, n, mk, ek, vk, sinvk, cosvk, pk,
     .    duk, drk, dik, uk, rk, ik, xpk, ypk, omk, rot_mat(3,3)
 

****  Get the time entry we need (use latest entry before current time)
      k = 1
      do j = 1,num_ent(i)
         if ( toe_jd(i,j).le.t ) then
             k = j
         endif
      end do

      gm = 3.986005d14
      eom = 7.2921151467d-5
 
      a = art(i,k)*art(i,k)
      n0 = sqrt(gm/a**3)
 
      tk = (t-toe_jd(i,k))*86400.0d0
      n = n0 + dn(i,k)
      mk = m0(i,k) + n*tk
 
****  Solve Keplers equation
      ek = mk
      do j = 1, 10
          ek = mk + ecc(i,k)*sin(ek)
      end do
 
***** Get the true anomaly
      sinvk = sqrt(1-ecc(i,k)**2)*sin(ek)/(1 - ecc(i,k)*cos(ek))
      cosvk = (cos(ek)-ecc(i,k))/(1-ecc(i,k)*cos(ek))
 
      vk = atan2(sinvk, cosvk)
 
*     Argument of latitude
      pk = vk + w(i,k)
 
*     Correction terms
      duk = cus(i,k)*sin(2*pk) +cuc(i,k)*cos(2*pk)
      drk = crs(i,k)*sin(2*pk) +crc(i,k)*cos(2*pk)
      dik = cis(i,k)*sin(2*pk) +cic(i,k)*cos(2*pk)
 
      uk = pk + duk
      rk = a*(1-ecc(i,k)*cos(ek)) + drk
      ik = i0(i,k) + dik + idt(i,k)*tk
 
*     Get inplane coordinates
      xpk = rk*cos(uk)
      ypk = rk*sin(uk)
 
*     Compute the longitude of the ascending node
      omk = om0(i,k) + omd(i,k)*tk
 
*     If we are in Earth fixed frame account for rotation of Earth
      if( sys(1:1).eq.'E' .or. sys(1:1).eq.'e') then
          omk = omk - eom*(tk+toe(i,k))
      end if
 
*     Get three_d coordinates
      svs_xyz(1,i) = xpk*cos(omk) - ypk*sin(omk)*cos(ik)
      svs_xyz(2,i) = xpk*sin(omk) + ypk*cos(omk)*cos(ik)
      svs_xyz(3,i) = ypk*sin(ik)
 
*     Now convert to latitude and longitude
      call xyz_to_geod( rot_mat, svs_xyz(1,i), svs_loc(1,i))
 
****  Thats all
      return
      end
 
