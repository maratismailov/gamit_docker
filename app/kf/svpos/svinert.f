 
      program svinert

      implicit none 

*
*     Program to Kalman filter GPS and Inertial navigation data
*      
 
*     The runstring of the program is
*     % svinert [navfile] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
 
      include 'svinert.h'
 
* Main pogram variables

      integer*4 i
 
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

****  Read the set up for the iertial system

      command_file = 'svinert.dat'
      call read_inert
 
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
          write(30,315)
 315      format(/'*     Date          Seconds        vN (m/s) +-  ',
     .            '      vE (m/s) +-         vU (m/s)+-         ')
          write(31,317)
 317      format(/'*     Date          Seconds        aN (m/ss) +-  ',
     .            '      aE (m/ss) +-         aU (m/ss)+-         ')
      else
          write(*,320)
 320      format(/'*     Date          Seconds        dX (m)  +-  ',
     .            '       dY (m) +-         dZ (m)  +-         ',
     .            ' dclk (m)  +-  Nsat Nusd  chi**2 SYS')
          write(30,325)
 325      format(/'*     Date          Seconds        vX (m/s) +-  ',
     .            '       vY (m/s) +-       vZ (m/s)  +-       ')
          write(31,327)
 327      format(/'*     Date          Seconds        aX (m/ss) +-  ',
     .            '       aY (m/ss) +-       aZ (m/ss)  +-       ')
      end if


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

          if( .not.eof .and. process ) call increment_est

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

CTITLE READ_INERT

      subroutine read_inert

      implicit none 

*     Routine to read the inertial setup file

      include 'svinert.h'

* LOCAL VARIABLES
      integer*4 ierr, i
      character*20 dum

****   OPen the input command file
      open(150, file = command_file, status='old', iostat=ierr)
      if( ierr.ne.0 ) write(*,120) ierr, command_file(1:20)
 120  format('IOSTAT error ',i5,' opening ',a)
      if( ierr.ne.0 ) stop 'Error opening command file'

****  Now read the entries.
      read(150,*,iostat=ierr) X_apr_sig
      call report_error('IOSTAT',ierr,'read','X   apr_s',0,'read_inert')     
      read(150,*,iostat=ierr) Xd_apr_sig
      call report_error('IOSTAT',ierr,'read','Xd  apr_s',0,'read_inert')     
      read(150,*,iostat=ierr) Xdd_apr_sig
      call report_error('IOSTAT',ierr,'read','Xdd apr_s',0,'read_inert')     
      read(150,*,iostat=ierr) Xddb_apr_sig
      call report_error('IOSTAT',ierr,'read','Xddbapr_s',0,'read_inert')     
      read(150,*,iostat=ierr) T_apr_sig
      call report_error('IOSTAT',ierr,'read','T   apr_s',0,'read_inert')     
      read(150,*,iostat=ierr) Td_apr_sig
      call report_error('IOSTAT',ierr,'read','Td  apr_s',0,'read_inert')     
      read(150,*,iostat=ierr) Tdb_apr_sig
      call report_error('IOSTAT',ierr,'read','Tdb apr_s',0,'read_inert')     

*     Now read the process noise
      read(150,*,iostat=ierr) wn
      call report_error('IOSTAT',ierr,'read','WN Markov',0,'read_inert')     

*     Now read the rate data sigma
      read(150,*,iostat=ierr) Xdd_sig
      call report_error('IOSTAT',ierr,'read','Td data s',0,'read_inert')     
      read(150,*,iostat=ierr) Td_sig
      call report_error('IOSTAT',ierr,'read','Xdd datas',0,'read_inert')     

*     Read the name of the inertial data file
      read(150,'(a)', iostat=ierr) inert_file
      open(110, file = inert_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',inert_file,0,'read_inert')     
      if( ierr.eq. 0 ) then
         read(110,'(a)') dum
         read(110,'(a)') dum
         inert_OK = .true.
      else
         inert_OK = .false.
         do i = 1,3
            Xdd_sig(i) = 10.
         end do
      end if


****  Thats all
      return
      end

CTITLE init_est

      subroutine init_est

      implicit none 

*     Routine to initialize the Kalman filter covariance matrices

      include 'svinert.h'

      integer*4 i,j

***** Initial the covariance matrix
*     Sites +- 1000 meters,
*     clocks +- 10000 meters

      do i = 1,22
         do j = 1,22
            cov_parm(i,j) = 0.d0
         end do
      end do

      do i = 1,3
         cov_parm(i,i) = X_apr_sig(i)**2
         cov_parm(i+ 4,i+ 4) = Xd_apr_sig(i)**2
         cov_parm(i+ 7,i+ 7) = Xdd_apr_sig(i)**2
         cov_parm(i+10,i+10) = Xddb_apr_sig(i)**2
         cov_parm(i+13,i+13) = T_apr_sig(i)**2
         cov_parm(i+16,i+16) = Td_apr_sig(i)**2
         cov_parm(i+19,i+19) = Tdb_apr_sig(i)**2
      end do
*     Set the clock variance large
      cov_parm(4,4) = 1.d6 

      do i = 1,22
         sol_vec(i) = 0.d0
      end do

      do i = 1, 3
         site_apr(i) = site_xyz(i)
         site_vel(i) = 0
         site_acc(i) = 0
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

      write(*,*) name
      do i = 1,nr 
         write(*,*) i, (mat(i,j),j=1,nc)
      end do
      return
      end
         
 
CTITLE INCREMENT_EST
 
      subroutine increment_est

      implicit none 
 
*     Routine to increment the solution for the new data.
 
      include  '../includes/const_param.h'
      include 'svinert.h'

* sec10 - minutes and seconds in tenths of seconds units.

      integer*4 i,j, k, l, num_av, ipivot(max_sat), date(5), sec10
      real*8 pc, average , dummy(max_sat), scale(max_sat), sectag,
     .       dp(max_sat), dchi, svclk(max_sat), dchiout

      real*8 pos_xyz_fin(3), pos_xyz_adj(3), pos_neu_fin(3),
     .       pos_neu_adj(3), loc_coord(3), rot_matrix(3,3),
     .       covar(3,3), neu_covar(3,3), temp_covar(3,3),
     .       norm_eq(4,4), fvec(4)

      real*8 vel_xyz_fin(3), vel_xyz_adj(3), vel_neu_fin(3),
     .       vel_neu_adj(3)
      real*8 acc_xyz_fin(3), acc_xyz_adj(3), acc_neu_fin(3),
     .       acc_neu_adj(3)

* in_sys, out_sys - Strings with XYZ and NEU in them

      character*4 in_sys, out_sys

      data in_sys  / 'XYZ ' /, out_sys / 'NEU ' /
 

****  Start by stepping the filter forward so that we linearize
*     about the best estimate of the state.  (site_xyz are updated
*     in this loop).

      call step_state

***** Clear the normal eqautions
      do i = 1,4
         fvec(i) = 0
         do j = 1,4
            norm_eq(j,i) = 0
         end do
      end do
 
****  compute the epheremis position at the measurement
*     time
 
c     write(*,*) 'For epoch ',data_epoch
 
      do i = 1, num_sat
 
          call eph_to_xyz( data_epoch-0.066666/86400.d0, i, 'E')
          svclk(i) = (af0(i)+
     .                af1(i)*(data_epoch-toe_jd(i))*86400.d0)*
     .               vel_light
 
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
                      pc = (p1o_data(i)*(77.d0/60.d0)**2-p2o_data(i))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc = p1o_data(i)
                  end if
                  pc = p1o_data(i)
c                 write(*,200) chan(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j),
c    .                       pc    -p1c(i)+svclk(j)           
C                 print *, ' Data ', average, p1c(i), svs_xyz(1,j),
C    .                               svs_xyz(2,j), svs_xyz(3,j), 
C    .                               svclk(j)
 
****              compute O-C
                  omc_data(i) = pc    -p1c(i)+svclk(j)-average 
 
****              Form the partial derivatives
                  apart_data(1,i) = (site_apr(1)-svs_xyz(1,j))/p1c(i)
                  apart_data(2,i) = (site_apr(2)-svs_xyz(2,j))/p1c(i)
                  apart_data(3,i) = (site_apr(3)-svs_xyz(3,j))/p1c(i)
                  apart_data(4,i) = 1
 
              end if
          end do
      end do

***** Now do difference site:     

      do i = 1, num_sat
 
          call eph_to_xyz( data_epoch-0.06666d0/86400.d0, i, 'E')
          svclk(i) = (af0(i)+
     .                af1(i)*(diff_epoch-toe_jd(i))*86400.d0)*
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
                      pc = (p1o_diff(i)*(77.d0/60.d0)**2-p2o_diff(i))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc = p1o_diff(i)
                  end if
                  pc = p1o_diff(i)
c                 write(*,200) chan_diff(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+svclk(j),
c    .                       pc    -p1c(i)+svclk(j)           
c                 print *, ' Diff ', average, p1c(i), svs_xyz(1,j),
c    .                               svs_xyz(2,j), svs_xyz(3,j), 
c    .                               svclk(j)
 
****              compute O-C
                  omc_diff(i) = pc    -p1c(i)+svclk(j)-average 
              end if
          end do
      end do

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

*     Check for "bad" data.
      num_chan = k
C     call scan_omc(omc, omc_OK, data_noise, num_chan, 
C    .              num_used )

****  Now get the gps position estimate
      do i = 1, 4
         fvec(i) = 0
         do j = 1,num_chan
            fvec(i) = fvec(i) + apart(i,j)*omc(j)
         end do
      end do
      do i = 1,4
         do j = 1,4
            do k = 1,num_chan
               norm_eq(i,j) = norm_eq(i,j) + apart(i,k)*apart(j,k)
            end do
         end do
      end do
    
      call invert_vis(norm_eq, fvec  , scale, ipivot, 4,4,1)
      write(20,997 ) fvec
 997  format( 'Fvec ', 4f10.4)

****  Now we have the GPS solution.  Add this to the kalman filter
      do i = 1, 4
         do j = 1,4
            if( i.eq.j ) apart(i,j) = 1.d0
            if( i.ne.j ) apart(i,j) = 0.d0
         end do
      end do 

***** Start computing the Kalman gain.  Get the acat matrix
      call acat(apart,  4, 4       , cov_parm, 22, 
     .        acat_mat, temp_gain, max_sat)
 
****  Now add the data noise (assume 10 meter data noise)
      num_used = 1
      do i = 1, 3
          do j = i,i
             acat_mat(i,j) = acat_mat(i,j) + 
     .                       norm_eq(i,j)*data_noise**2  
          end do
      end do
      acat_mat(4,4) = acat_mat(4,4) + 1.d0
 
***** Now invert this results
      call invert_vis(acat_mat, dummy, scale, ipivot, 
     .                4, max_sat, 0)


***** Now form the Kalman Gain
      do i = 1,22
          do j = 1, 4
              kgain(i,j) = 0
              do k = 1, 4        
                  kgain(i,j) = kgain(i,j) +
     .                    temp_gain(i,k)*acat_mat(k,j)
              end do
          end do
      end do
 
C     call write_mat('kgain', kgain, 4, num_chan, 4, max_sat)
      write(20,998) kgain(1,1), kgain(5,1), kgain(8,1)
 998  format('KGAIN X', 3f10.6)

***** Update the parameter estimates
      do j = 1, 4       
          dp(j) = fvec(j)
      end do
      do i = 1,22
          dx(i) = 0.d0
          do j = 1, 4         
              dx(i) = dx(i) + kgain(i,j)*dp(j)
          end do
      end do

      do i = 1,22
          sol_vec(i) = sol_vec(i) + dx(i)
      end do
 
***** Now update covariance matrix
      do i = 1,22
          do j = 1,22
              do k = 1, 4         
                  cov_parm(i,j) = cov_parm(i,j) -
     .                     kgain(i,k)*temp_gain(j,k)
              end do
          end do
      end do
*
*     read the inert_data file

      call get_inert_data

*     Now increment the estimates for the inertial data
      call increment_inert
 
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
 300          format(i5,4i3,1x,f8.3,1x,3(f12.4,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2, 2I3,1x, f8.2,' NEU')
           else
C             write(*,310) date, sectag, (sol_vec(i), 
              write(*,310) date, sectag, (pos_xyz_fin(i)-site_xyz(i),
     .               sqrt(cov_parm(i,i)), i = 1,3), 
     .               sol_vec(4), sqrt(abs(cov_parm(4,4))),
     .               num_used, num_chan, dchiout
 310          format(i5,4i3,1x,f8.3,1x,3(f12.4,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2, 2I3,1x, f8.2,' XYZ')
           end if

*****     Output the velocity estimates. 
          do i = 1,3
             vel_xyz_fin(i) = sol_vec(i+4)
             vel_xyz_adj(i) = vel_xyz_fin(i) 
             do j = 1,3
                covar(i,j) = cov_parm(i+4,j+4)
             end do
          end do

          if( out_type(1:1).eq.'N') then
              call rotate_geod(vel_xyz_adj,vel_neu_adj,in_sys,out_sys,
     .            pos_xyz_fin,loc_coord,rot_matrix)
c
c....         Now compute the sigmas of the local coordinates. Firstly
c             save the covariance elements
              call var_comp(rot_matrix,covar,NEU_covar,temp_covar, 
     .                      3,3,1)

              write(30,320) date, sectag, (vel_neu_adj(i),
     .               sqrt(NEU_covar(i,i)),i=1,3)
 320          format(i5,4i3,1x,f8.3,1x,3(f11.4,1x,f8.2), 1x,'NEU')
           else
              write(30,330) date, sectag, (vel_xyz_fin(i),
     .               sqrt(cov_parm(i+4,i+4)), i = 1,3)
 330          format(i5,4i3,1x,f8.3,1x,3(f11.4,1x,f8.2), 1x,'XYZ')
           end if

*****     Output the acceleration estimates.
          do i = 1,3
             acc_xyz_fin(i) = sol_vec(i+7)
             acc_xyz_adj(i) = acc_xyz_fin(i) 
             do j = 1,3
                covar(i,j) = cov_parm(i+7,j+7)
             end do
          end do

          if( out_type(1:1).eq.'N') then
              call rotate_geod(acc_xyz_adj,acc_neu_adj,in_sys,out_sys,
     .            pos_xyz_fin,loc_coord,rot_matrix)
c
c....         Now compute the sigmas of the local coordinates. Firstly
c             save the covariance elements
              call var_comp(rot_matrix,covar,NEU_covar,temp_covar, 
     .                      3,3,1)

              write(31,340) date, sectag, (acc_neu_adj(i),
     .               sqrt(NEU_covar(i,i)),i=1,3)
 340          format(i5,4i3,1x,f8.3,1x,3(1x,f16.9,1x,f8.5), 1x,'NEU')
           else
              write(31,350) date, sectag, (acc_xyz_fin(i),
     .               sqrt(cov_parm(i+7,i+7)), i = 1,3)
 350          format(i5,4i3,1x,f8.3,1x,3(1x,f16.9,1x,f8.5), 1x,'XYZ')
           end if
      end if

****  Adopt the current position as the expansion
c     do i = 1, 3
c        site_apr(i) = site_apr(i) + sol_vec(i)
c        sol_vec(i) = 0.d0
c     end do

*     Update the apriori so we start with these coordinates
C     do i = 1,3

*         If the site is stochastic then update the site position
*         and do not zero sol_vec.  (By not zero'ing sol_vec we are
*         assuming a constant velocity model).
C         if( wn(i).ne.0 ) then
C             site_apr(i) = site_apr(i) + sol_vec(i)
*             Now re-zero sol_vec since the the new apriori will
*             be used in the filter (assumes zero velocity)
C             sol_vec(i) = 0.d0             
C         end if
C     end do 

****  Thats all
      return
      end

CTITLE step_state

      subroutine step_state

      implicit none 

*     Rouitne to step the state vector forward one step

      include 'svinert.h'

      integer*4 dim
C     real*8 cov_parm(dim,dim), sol_vec(dim), wn(dim)
      

* STM - State transission matrix
* stm_tmp - Temp array at SCov calcualtion

      integer*4 i,j,k
      real*8 stm(22,22), stm_tmp(22,22), dt, svtmp(22), dsol

***** Set dimension
      dim = 22


****  Form up state transission for a signle 1 second time step
      do i= 1,22
         do j = 1,22
            stm(i,j) = 0.d0
         end do
      end do

***** All elements have 1 transission for themselves
      do i = 1, 22
         stm(i,i) = 1.d0
      end do

****  Now form the invidual elments
*     Xpos, vel and acc.  DeltaTime for Ashtech data is 0.5 seconds
      dt = 0.5d0  
      do i = 1,3
         stm(i,4+i) = dt   
         stm(i,7+i) = dt**2/2
         stm(4+i,7+i) = dt
      end do

*     Now do the rotation rate
      do i = 1,3
         stm(13+i,16+i) = dt
      end do

****  Now propagate the state
      do i = 1,22
         svtmp(i) = sol_vec(i)
      end do

      do i = 1,22
         dsol = 0.d0
         do j = 1,22
            dsol = dsol + stm(i,j)*svtmp(j)
         end do
         sol_vec(i) = dsol
      end do

****  For the position add the updated position to the
*     apriori position.
      do i = 1, 3
         site_apr(i) = site_apr(i) + sol_vec(i)
         sol_vec(i)  = 0
      end do
      sol_vec(4) = 0.

****  Now propagate the covariance matrix of the state.
      do i = 1, 22
         do j = 1, 22
            stm_tmp(i,j) = 0.d0
            do k = 1,22
*              Notice transpose of matrix
               stm_tmp(i,j) = stm_tmp(i,j)+cov_parm(i,k)*stm(j,k)
            end do
         end do
      end do

*     Now complete with S*(VST)
      do i = 1, 22
         do j = 1, 22
            dsol = 0.d0         
            do k = 1,22
*              Notice transpose of matrix
               dsol  = dsol + stm(i,k)*stm_tmp(k,j)
            end do
            cov_parm(i,j) = dsol
         end do
      end do

      write(20,998) 'COV X Pre', cov_parm(1,1), cov_parm(5,5),
     .  cov_parm(8,8), cov_parm(1,5), cov_parm(1,8), cov_parm(5,8)
 998  format( a,6d12.4)

****  Now propagate the process noise as well.   Propagate the
*     process noise according to the simple state transistion.
*     This should make the positions and accelerations constistent.
      do i = 1, 22
         do j = 1, 22
            stm_tmp(i,j) = wn(i)*stm(j,i)
         end do
      end do

*     Now complete with S*(VST)
      do i = 1, 22
         do j = 1, 22
            dsol = 0.d0         
            do k = 1,22
*              Notice transpose of matrix
               dsol  = dsol + stm(i,k)*stm_tmp(k,j)
            end do
            cov_parm(i,j) = cov_parm(i,j) + dsol
         end do
      end do
      write(20,998) 'COV X Pos', cov_parm(1,1), cov_parm(5,5),
     .  cov_parm(8,8), cov_parm(1,5), cov_parm(1,8), cov_parm(5,8)

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
                if( abs(omc(i)-omc(j)).lt.100.d4*data_noise ) then
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
 
      subroutine acat(a,np, nc, c, nv, acat_mat, temp_gain, max_sat)

      implicit none 
 
*     Routine to compute a*c*at matrix
 
*   np      - Number of parameters (assummed first dimension of
*           - a)
*   nc      - number of channels
*   nv      - Dimension of covariance matrix
*   max_sat - Maxiumum number of channels possible
 
      integer*4 np, nc, max_sat, nv
 
*   a(np,nc)    - Partials matix
*   c(np,np)    - covariance matrix
*   temp_gain(np, max_sat)  - contains a*c
 
      real*8 a(np,nc), c(nv,nv), acat_mat(max_sat,max_sat),
     .    temp_gain(nv, max_sat)
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
***** First compute temp_gain
      do i = 1, nv
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
      include 'svinert.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(5), ierr, jerr,i,j, id, trimlen
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(5)
 
*   eof     - Indicates end of file.
 
      logical eof
 
      character*256 line
      character*1 cr

****  Read in the next line
      cr = char(13)
      read(101,'(a)', iostat=ierr) line
      call sub_char(line,cr,' ') 
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if
 
*     Decode the line
      read(line,*,iostat=ierr) date, sectag, id, num_data,
     .            (chan_data(i),i=1,num_data)
      call ymdhms_to_jd( date, sectag, data_epoch)
 
*     Now loop over the data records.
      do i = 1, num_data
          read(101,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ') 
          read(line,120, iostat=ierr) (vals(j), flags(j), j =1,5)
 120      format(10(f14.3,1x,i1))
 
*         Now assign the phase and range measurements
          do j = 1,5
              if((data_types(j).eq.'C1' .or.
     .            data_types(j).eq.'P1') .and.vals(j).ne.0  )
     .                   p1o_data(i) = vals(j)
              if( data_types(j).eq.'L1' ) l1o_data(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( data_types(j).eq.'L2' ) l2o_data(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( data_types(j).eq.'P2' ) p2o_data(i) = vals(j)
          end do
          if( num_data_types.gt.5 ) 
     .        read(101,120, iostat=jerr) (vals(j), flags(j), 
     .                   j =1,num_data_types-5)
      end do

      if( ierr.ne.0 ) eof = .true.
 
****  Thats all
      return
      end
 
CTITLE READ_rdiff  
 
      subroutine read_rdiff(eof)
 
      implicit none 

*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include 'svinert.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(5), ierr,jerr, i,j, id, trimlen
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(5)
 
*   eof     - Indicates end of file.
 
      logical eof
 
      character*256 line
      character*1 cr
 
****  Read in the next line
      cr = char(13)
      read(102,'(a)', iostat=ierr) line
      call sub_char(line,cr,' ') 
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if
 
*     Decode the line
      read(line,*,iostat=ierr) date, sectag, id, num_diff,
     .            (chan_diff(i),i=1,num_diff)
      call ymdhms_to_jd( date, sectag, diff_epoch)
 
*     Now loop over the data records.
      do i = 1, num_diff
          read(102,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ') 
          read(line,120, iostat=ierr) (vals(j), flags(j), j =1,5)
 120      format(10(f14.3,1x,i1))
 
*         Now assign the phase and range measurements
          do j = 1,5
              if((diff_types(j).eq.'C1' .or.
     .            diff_types(j).eq.'P1') .and. vals(j).ne.0 ) 
     .                  p1o_diff(i) = vals(j)
              if( diff_types(j).eq.'L1' ) l1o_diff(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( diff_types(j).eq.'L2' ) l2o_diff(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( diff_types(j).eq.'P2' ) p2o_diff(i) = vals(j)
          end do
          if( num_diff_types.gt.5 ) 
     .        read(102,120, iostat=jerr) (vals(j), flags(j), 
     .                   j =1,num_diff_types-5)
      end do

      if( ierr.ne.0 ) eof = .true.
 
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
 
      include 'svinert.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i
 
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
     .        'read_data_file/svinert')
 
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
 
*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.gt.0 ) rinex = .true.
 
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.0 ) site_name = line(1:8)
 
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.0 .and. site_xyz(1).eq.0 ) then
              read(line,*) site_xyz
          end if
 
*         See if data types
          indx = index(line,'TYPES OF OBS')
          if( indx .gt.0 ) then
              read(line,*) num_data_types, 
     .                    (data_types(i), i = 1, num_data_types)
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
          stop 'svinert: Wrong type data file'
      end if
 
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
 
      include 'svinert.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i
 
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
     .        'read_diff_file/svinert')
 
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
 
*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.gt.0 ) rinex = .true.
 
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.0 ) diff_name = line(1:8)
 
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.0 .and. diff_xyz(1).eq.0 ) then
              read(line,*) diff_xyz
          end if
 
*         See if data types
          indx = index(line,'TYPES OF OBS')
          if( indx .gt.0 ) then
              read(line,*) num_diff_types, 
     .                    (diff_types(i), i=1,num_diff_types)
          end if
      end do
 
*     See if we found rinex file
      if( rinex ) then
          write(*,150) diff_file(1:trimlen(data_file)),
     .        diff_name, diff_xyz
 150      format(/'* For RINEX data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) diff_file(1:trimlen(diff_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'svinert: Wrong type data file'
      end if
 
****  Thats all
      return
      end
 
 
CTITLE GET_SVSRUN
 
      subroutine get_svsrun
 
      implicit none 

*     Routine to get the cunstring.
 
      include '../includes/const_param.h'
      include 'svinert.h'
 
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
          call proper_runstring('svinert.hlp','svinert/nav file',1)
      end if
 
*     See if data difle names passed
      len_run = rcpar(2,data_file)
      if( len_run.le.0 ) then
          call proper_runstring('svinert.hlp','svinert/data file',1)
      end if
*     See if diff difle names passed
      len_run = rcpar(3,diff_file)
      if( len_run.le.0 ) then
          call proper_runstring('svinert.hlp','svinert/diff file',1)
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

*     Now see if positions passed 
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
 
****  Thats all
      return
      end


CTITLE READ_NAV
 
      subroutine read_nav

      implicit none 
 
*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero)
 
      include 'svinert.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
 
 
      integer*4 ierr, i, rp, trimlen, ll, date(5)
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
 120          format(' Rinex version ',i2,' Nav file found')
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
                  end if
              end do
 
*                                         ! This is first
              if( .not.have ) then
*                                             ! on this prn so read
                  num_sat = num_sat + 1
                  rp = num_sat
                  read(line,200) prn(rp), date, sectag,
     .                        af0(rp),af1(rp), af2(rp)
 200              format(i2,5i3,f5.1,3d19.8)
 
                  call ymdhms_to_jd( date, sectag, toe_jd(rp))
 
*                 Read the rest of the entires
                  read(100,210) aode(rp), crs(rp), dn(rp), m0(rp)
                  read(100,210) cuc(rp), ecc(rp), cus(rp), art(rp)
                  read(100,210) toe(rp), cic(rp), om0(rp), cis(rp)
                  read(100,210) i0(rp) , crc(rp), w(rp)  , omd(rp)
                  read(100,210) idt(rp), cflg12(rp), weekno, pflg12(rp)
                  read(100,210) svacc, svhealth, tgd, aodc(rp)
 210              format(3x,4d19.8)
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              else
*                 Skip 6 lines in file
                  do i = 1,6
                      read(100,'(a)') line
                  end do
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              end if
          end if
      end do
 
      write(*,300) num_sat, nav_file(1:trimlen(nav_file))
 300  format('* ', i5,' satellites found in ',a)
 
      if( num_sat.eq.0 ) stop ' SVPOS: No satellites found'
 
      return
      end
 
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none 
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'svinert.h'
 
*   t       - Time for computation (day number)
 
 
      real*8 t
 
*   i       - Satellite number
 
 
      integer*4 i
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j
 
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
 
      gm = 3.986005d14
      eom = 7.2921151467d-5
 
      a = art(i)*art(i)
      n0 = sqrt(gm/a**3)
 
      tk = (t-toe_jd(i))*86400.0d0
      n = n0 + dn(i)
      mk = m0(i) + n*tk
 
****  Solve Keplers equation
      ek = mk
      do j = 1, 10
          ek = mk + ecc(i)*sin(ek)
      end do
 
***** Get the true anomaly
      sinvk = sqrt(1-ecc(i)**2)*sin(ek)/(1 - ecc(i)*cos(ek))
      cosvk = (cos(ek)-ecc(i))/(1-ecc(i)*cos(ek))
 
      vk = atan2(sinvk, cosvk)
 
*     Argument of latitude
      pk = vk + w(i)
 
*     Correction terms
      duk = cus(i)*sin(2*pk) +cuc(i)*cos(2*pk)
      drk = crs(i)*sin(2*pk) +crc(i)*cos(2*pk)
      dik = cis(i)*sin(2*pk) +cic(i)*cos(2*pk)
 
      uk = pk + duk
      rk = a*(1-ecc(i)*cos(ek)) + drk
      ik = i0(i) + dik + idt(i)*tk
 
*     Get inplane coordinates
      xpk = rk*cos(uk)
      ypk = rk*sin(uk)
 
*     Compute the longitude of the ascending node
      omk = om0(i) + omd(i)*tk
 
*     If we are in Earth fixed frame account for rotation of Earth
      if( sys(1:1).eq.'E' .or. sys(1:1).eq.'e') then
          omk = omk - eom*(tk+toe(i))
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
 
 
CTITLE INCREMENT_INERT
 
      subroutine increment_inert

      implicit none 

*     Routine to add the inertial information into the filetr

      include 'svinert.h'

      integer*4 i, j, k, date(5), imn, jmn, imx, jmx
      real*8 tg(22), kg(22), var, sectag, max_corr, min_corr, rho

      call jd_to_ymdhms(data_epoch, date, sectag)
      write(20,100) 'Pre', date, sectag, 
     .              (sol_vec(i), sqrt(cov_parm(i,i)), i = 1,22)
 100  format(/,a,i5,4i3,1x,f8.3,1x,8(1x,d15.4), /,
     .      10(6(1x,d15.4),:/))
      max_corr = -10
      min_corr = +10
      do i = 1,22
         if( cov_parm(i,i).gt.0 ) then
         do j = 1,22
            if( i.ne.j .and. cov_parm(j,j).gt.0 ) then
               rho = cov_parm(i,j)/sqrt(cov_parm(i,i)*cov_parm(j,j))
               if( rho.lt.min_corr ) then
                   min_corr = rho
                   imn = i
                   jmn = j
               end if
               if( rho.gt.max_corr ) then
                   max_corr = rho
                   imx = i
                   jmx = j
               end if
            end if
         end do
         end if
      end do
      if( min_corr.lt.10.d0 ) then
      write(20,120) 'Min Corr ',imn, jmn, min_corr, 
     .  cov_parm(jmn,imn)/sqrt(cov_parm(imn,imn)*cov_parm(jmn,jmn))
 120  format(a,2i5,2d18.10)
      end if
      if( max_corr.gt.-10.d0 ) then
      write(20,120) 'Max Corr ',imx, jmx, max_corr, 
     .  cov_parm(jmx,imx)/sqrt(cov_parm(imx,imx)*cov_parm(jmx,jmx))
      end if

      do i = 1, 22
         do j = 1,22
            if( i.ne.j ) then
                if( abs(cov_parm(i,j)).lt.abs(cov_parm(j,i)) ) then
                    cov_parm(j,i) = cov_parm(i,j)
                else
                    cov_parm(i,j) = cov_parm(j,i)
                end if
            end if
         end do
      end do
            

****  Now assume that we have acceleration data
*     Form CAt(ACAt + sigma)-1 = kgain
*     Do for 1 variable (later need to change to 3 due to
*     correlations from orientation).

      write(20,150) 'Xdd_obs',
     .              (Xdd_obs(i)-(sol_vec(7+i)),i=1,3)
 150  format(a,3f12.8)

      do i = 1,3

*        Only update the filter if the acc data is OK,
         if( xdd_sig(i).lt.1.d0 ) then
*           Partial is unit therefore ACAt = C
            do j = 1,22
               tg(j) = cov_parm(j,7+i)
            end do
            var = tg(7+i) + xdd_sig(i)**2
            do j = 1, 22
               kg(j) = tg(j)/var
            end do

*           Now update the covaraince matrix
            do j = 1, 22
               sol_vec(j) = sol_vec(j) + kg(j)*
     .                     (Xdd_obs(i)-(sol_vec(7+i)))
               do k = 1, 22
                  cov_parm(j,k) = cov_parm(j,k) - kg(j)*tg(k)
               end do
            end do
         end if
      end do
      write(20,100) 'Pos', date, sectag, 
     .              (sol_vec(i), sqrt(cov_parm(i,i)), i = 1,22)
      

      return 
      end

CTITLE GET_INERT_DATA

      implicit none 

      subroutine get_inert_data

*     Routine to read inertdata connected to unit 110

      include 'svinert.h'

      integer*4 ierr, date(5), i
      real*8 sec, dum

*     Read the next record

      if( .not. inert_OK ) Return

      read(110,*,iostat=ierr) date, sec, (Xdd_obs(i), dum,i=1,3)
      call ymdhms_to_jd( date, sec, ep_inert)
      if( abs(ep_inert-data_epoch).gt.1.d-5 ) then
          write(*,120) date, sec
 120      format('Inertial Record out of sync at ',i5,4i3,1x,f5.2)
      end if

      return
      end
