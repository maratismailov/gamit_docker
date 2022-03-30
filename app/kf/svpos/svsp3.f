 
      program svsp3

      implicit none
 
*     Program to compute the rise and set times of the GPS
*     satelites based on the ephemeris in the rinex navigiation
*     files
 
*     The runstring of the program is
*     % svsp3 [sp3 orbit file] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
 
      include '../includes/const_param.h'
      include 'svsp3.h'
 
* Main pogram variables

      integer*4 i, ep
 
*   eof     - Indicates end of file.
*   process  - Indicates that we have common data so process
*   read_data  = Set true if we should read the data file
*   read_diff  = set true if we should read difference file
 
      logical eof, process, read_data, read_diff

      real*8 mean, rms   ! Mean and RMS for satellite depependent
                         ! resiudals
 
***** Get the runstring runstring
 
      call get_svsrun
 
*     Read in the ephemeris file
      if( debug_start.ne.0 ) out_sp3freqs = .true.
      call read_sp3 ( nav_file ) 
 
*     Now loop over the data and get an estimate of the site position.
 
      call read_data_head
      if ( .not. Nodiff ) call read_diff_head
 
*     Now loop over data.
      call init_est

*     Setup geodetic quanties for later use.
      call set_geod 

      eof = .false.
      process = .false.
      read_data = .true.
      read_diff = .true.

      ep = 0
      if( out_type(1:1).eq.'N' ) then
         write(*,120) 
 120     format('*  Yr  M  D  H  M     sec         dN       +- ',
     .          '       dE       +-        dU       +-         ',
     .          ' dClk       +-   dZen1   +-    dZen2   +-   Nu Na',
     .          '     dChi   Ep  SYS')
      else
         write(*,130) 
 130     format('*  Yr  M  D  H  M     sec         dX       +- ',
     .          '       dY       +-        dX       +-         ',
     .          ' dClk       +-   dZen1   +-    dZen2   +-   Nu Na',
     .          '     dChi   Ep  SYS')
      end if

      do while ( .not.eof )
          if( read_data ) call read_range(eof)
          if( .not.Nodiff ) then
             if( read_diff .and. .not.eof ) call read_rdiff( eof )

             if( data_epoch.gt.diff_epoch ) then
                 read_data = .false.
                 read_diff = .true.
                 process = .false.
             end if
             if( data_epoch.lt.diff_epoch ) then
                 read_data = .true. 
                 read_diff = .false.
                 process = .false.
             end if
             if( data_epoch.eq.diff_epoch ) then
                 read_data = .true.
                 read_diff = .true.
                 process = .true. 
             end if
          else
*            Point position set process OK
             process = .true.
          end if
          if( .not.eof .and. process ) then
              ep = ep + 1
              call increment_est(ep)
          end if
 
      end do

****  Output the satellite offsets
      write(*,220) 
 220  format('* SVMEAN CHANNEL MEAN OFFSETS BY PRN',/,
     .       '* SVMEAN PRN  # Ep   Mean (m)      +- (m)   RMS (m)',
     .       '  Mean (ns)    L1 Freq (MHz)' )
      do i = 1, num_sat     
         if( svres_num(i).gt.1 ) then
            mean = svres_sum(i)/svres_num(i)
            rms  = sqrt((svres_var(i)-mean**2*svres_num(i))/
     .                  (svres_num(i)-1))
            write(*,240) prn_sp3(i), svres_num(i), mean, 
     .            sqrt(rms**2/svres_num(i)), rms, mean/vel_light*1.d9,
     .            fL1(i)*1.d-6
 240        format('* SVMEAN ',i3,1x,i5,1x,4(F10.2,1x),1x,F10.4)
         endif
      end do
         
 
*     Now write out the resulys
c     call write_est
      write(*,320) site_name, (site_xyz(i)+sol_vec(i),
     .              sqrt(cov_parm(i,i)), i=1,3),
     .              sqrt(chi/nchi)
 320  format('*',/,'* Final position for ',a,' X ',
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

      include 'svsp3.h'

      integer*4 i,j

***** Initial the covariance matrix
*     Sites  +-  1000 meters,
*     clocks +- 10000 meters
*     Atmospheres +- 1.0 meters

      do i = 1,6
         do j = 1,6
            cov_parm(i,j) = 0.d0
         end do
      end do

      do i = 1,3
         cov_parm(i,i) = 1.d4   ! 100m^2 
      end do
      cov_parm(4,4) = 1.d4      ! 100m^2

*     Currently set the atmospheric delay apriori sigmas to zero
*     so that that 
      cov_parm(5,5) = atm_apr_sig(1)**2
      cov_parm(6,6) = atm_apr_sig(1)**2

      do i = 1,6
         sol_vec(i) = 0.d0
      end do

      svres_sum = 0
      svres_var = 0
      svres_num = 0

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
 
      subroutine increment_est(ep)
 
      implicit none

*     Routine to increment the solution for the new data.
 
      include  '../includes/const_param.h'
      include 'svsp3.h'

* sec10 - minutes and seconds in tenths of seconds units.

      integer*4 i,j, k, l, num_av, ipivot(max_sat), date(5), sec10
      integer*4 ep    ! Epoch counter
      integer*4 omc_prn(max_sat)  ! PRN numbers in diffs
      integer*4 num_org  ! Number of data before cleaning 

      real*8 average , dummy(max_sat), scale(max_sat), sectag,
     .       dp(max_sat), dchi,  dchiout

      real*8 pc(max_sat)

      real*8 pos_xyz_fin(3), pos_xyz_adj(3), pos_neu_fin(3),
     .       pos_neu_adj(3), loc_coord(3), rot_matrix(3,3),
     .       covar(3,3), neu_covar(3,3), temp_covar(3,3), elev

*
      real*8 pcf1, pcf2   ! Factors to ionoshere free observables 
                          ! Computed using specific frequencies for each
                          ! satellite

      real*8 clk_data, clk_diff  ! Clock errors removed from site 
                          ! 1 (data) and reference site 2 (diff).

      real*8 sod   ! Seconds of day to compute output spacing

      real*8 svs_xyz_i(3)  !  Satellite coordinates rotated to account 
                      ! for Earth rotation during light propagatio. time

      real*8 pcobs    ! Function to return the linear combination of the
                      ! observed P1/P2 values.  PC is free from ionospheric
                      ! delay.
 
      real*8 atm_del  ! Function to return atmospheric delay based on GPT/GMF

****  compute the epheremis position at the measurement
*     time
 
c     write(*,*) 'For epoch ',data_epoch
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,110) 'Satellite Position at t-66.66 us', ep, 0.0
 110      format(/,'+',1x,a,' Epoch ',i5,' Site clock ',f17.2,' (m)',/,
     .           '+PRN',5x,'Xe (m)',7x,'Ye (m)',7x,'Ze (m)',
     .           5x,'SV Clock (m)')
      end if
 
 
      do i = 1, num_sat
 
          call eph_to_xyz( data_epoch-0.083/86400.d0, i, 'E')
          call sp3_to_clk( data_epoch-0.083/86400.d0, i) 
*         See if debug wanted:
          if( ep.ge.debug_start .and. ep.le.debug_end ) then
              write(*,120) prn_sp3(i), k, (svs_xyz(j,i),j=1,3),
     .                                     svs_dt(i)
 120          format('+ P ',i3.2, ' E ', i3, 3F13.2, F13.2)
          end if

      end do

      average = 0.0
*     Get first clock estimate using just P1 data (more robust to
*     data quality).
      do i = 1, num_data
 
          omc_OK_data(i) = .false.
          do j = 1, num_sat
              if( chan_data(i).eq.prn_sp3(j) ) then
*                 Get an approximate Range first
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)

                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
                  
*                 Now compute the Earth's rotation contributions
                  call orb_drot(svs_xyz(1,j), svs_xyz_i, p1c(i))

*                 Compute a rough range
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz_i(1))**2+
     .                           (site_xyz(2)-svs_xyz_i(2))**2+
     .                           (site_xyz(3)-svs_xyz_i(3))**2)

c                 write(*,200) chan_data(i), j, p1c(i), p1o_data(i),
c    .                       p1o_data(i)-p1c(i)+sp3_clk(j)         
 200              format(2i5,10F13.3)
 
*                 accumulate the average clock offset
                  omc_data(i) =  p1o_data(i)-p1c(i)+svs_dt(j) 
                  omc_OK_data(i) = .true.
                  if( ep.ge.debug_start .and. ep.le.debug_end ) then
                      write(*,210) i, j, chan_data(i), omc_data(i), 
     .                    omc_OK_data(i), p1o_data(i), p1c(i), svs_dt(j)
 210                  format('DETAIL I, J ',2i4,' PRN ',i3,' omc ',
     .                     F12.2,1x,L1,1x,3(F12.2,1x))
                  endif   
               end if
           end do
           if( .not. omc_OK_data(i) ) then
                write(*,215) i, chan_data(i), omc_data(i), 
     .              omc_OK_data(i), p1o_data(i), p1c(i)
 215            format('BAD I ',i4,' PRN ',i3,' omc ',
     .               F12.2,1x,L1,1x,3(F12.2,1x))
           endif

      end do 

*     Scan the omc values and make sure OK
C     call scan_omc(omc_data, omc_OK_data, data_noise, num_data, 
C    .              num_used )
      num_org = num_data   ! Value of num_data changes in clean_data.
      call clean_data( omc_data, omc_OK_data, num_data, 100.d0, 
     .                 'BASE', ep)

*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_data
         if( omc_OK_data(i) ) then
             average = average + omc_data(i)
             num_av  = num_av + 1
         end if
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,105) i, chan_data(i), omc_data(i), omc_OK_data(i),  
     .              num_av, p1o_data(i)
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
      if( ep.ge.debug_start .and. ep.le.debug_end )
     .    write(*,107) 'P1', num_av, average
 
      do i = 1, num_data
 
          omc_OK_data(i) = .false.
          do j = 1, num_sat
              if( chan_data(i).eq.prn_sp3(j) ) then
*                 Get an approximate Range first
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)

                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')

*                 Now compute the Earth's rotation contributions
                  call orb_drot(svs_xyz(1,j), svs_xyz_i, p1c(i))

*                 Compute a rough range
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz_i(1))**2+
     .                           (site_xyz(2)-svs_xyz_i(2))**2+
     .                           (site_xyz(3)-svs_xyz_i(3))**2)

*                 Get observed value
                  pc(i) = pcobs( p1o_data(i), p2o_data(i), j)

*                 accumulate the average clock offset
                  call get_elev(site_xyz, svs_xyz_i, elev)
                  omc_data(i) =  pc(i)-p1c(i)+svs_dt(j )- 
     .                 atm_del(data_epoch, site_geod, elev) 
                  omc_OK_data(i) = .true.
                  if( ep.ge.debug_start .and. ep.le.debug_end ) then  
                     write(*,103) ep, chan_data(i), elev, 
     .                  atm_del(data_epoch, site_geod, elev)  
 103                 format('ATM DEL EP ',i5,' PRN ',I4,' Elev ',
     .                  F6.2,' deg, DELAY ',F8.3,' m')
                  end if
              
               end if
           end do
      end do 

      num_org = num_data
      call clean_data( omc_data, omc_OK_data, num_data, 100.d0,
     .                'BASE',ep)

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
             write(*,105) i,chan_data(i), omc_data(i), omc_OK_data(i),
     .                    num_av,   pc(i)
 105         format('Data ',i2,' PRN ',I3,' OMC ',f20.2,' OK ',L1,
     .              ' # ',i4, ' Obs ',F20.2)
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
      if( ep.ge.debug_start .and. ep.le.debug_end )
     .    write(*,107) 'PC', num_av, average
 107  format('Clock: Type ',a,' Num ',I3,' Values, Average ',F20.2,' m')

 
****  Now do the final range computation
      do i = 1, num_data
          do j = 1, num_sat
              if( chan_data(i).eq.prn_sp3(j) ) then
*                 Get an approximate Range first
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)

 
*                 Compute the range accounting for the propagation
*                 delay.
                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')

*                 Now compute the Earth's rotation contributions
                  call orb_drot(svs_xyz(1,j), svs_xyz_i, p1c(i))
 
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz_i(1))**2+
     .                           (site_xyz(2)-svs_xyz_i(2))**2+
     .                           (site_xyz(3)-svs_xyz_i(3))**2)

*
*                 Get observed value
                  pc(i) = pcobs( p1o_data(i), p2o_data(i), j)
 
****              compute O-C with atm delay 
                  call get_elev(site_xyz, svs_xyz_i, elev)
                  omc_data(i) = pc(i) -p1c(i)+svs_dt(j)- average - 
     .                 atm_del(data_epoch, site_geod, elev) 
 
****              Form the partial derivatives
                  apart_data(1,i) = (site_xyz(1)-svs_xyz_i(1))/p1c(i)
                  apart_data(2,i) = (site_xyz(2)-svs_xyz_i(2))/p1c(i)
                  apart_data(3,i) = (site_xyz(3)-svs_xyz_i(3))/p1c(i)
                  apart_data(4,i) = 1

****              Get elevation angles for the two sites
                  call mit_dry(elev, apart_data(5,i))
                  if( elev.lt.10.d0 ) omc_OK_data(i) = .false.
                  if( .not.Nodiff ) then
                      call get_elev(diff_xyz, svs_xyz_i, elev)
                      call mit_dry(elev, apart_data(6,i))
                      if( elev.lt.10.d0 ) omc_OK_diff(i) = .false.
                  else
                      apart_data(6,i) = 0
                  endif
                  if( ep.ge.debug_start .and. ep.le.debug_end ) then  
                     write(*,103) ep, chan_data(i), elev, 
     .                  atm_del(data_epoch, site_geod, elev)  
                  end if
            
              end if
         end do
         clk_data = average
      end do

*     See if debug for o-minus-c wanted
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,140) 'Range DATA site differences', ep, average
 140      format('+ ',a,' Epoch ',i4,' Clock ',F13.2,' m',/,
     .           '+CHAN PRN',3x,'Obs (m)',6x,'Range (m)',3x,
     .           'SVS Clk (m)',1x,'Obs-Calc (m)')
       
          do i = 1, num_data
              do j = 1, num_sat
                  if( chan_data(i).eq.prn_sp3(j) ) then
                      write(*,150) i,prn_sp3(j), pc(i), p1c(i), 
     .                             svs_dt(j), omc_data(i)
 150                  format('+',i3,1x,i4, 5F13.2)
                  end if
              end do
          end do
      end if


***** Now do difference site:     
* MOD TAH 151021: Only do this if difference data is being used
      clk_diff = 0.d0
      if( .not. Nodiff ) then
         do i = 1,num_diff
             omc_OK_diff(i) = .true.
         end do

         do i = 1, num_sat
             call eph_to_xyz( data_epoch-0.083/86400.d0, i, 'E')
             call sp3_to_clk( data_epoch-0.083/86400.d0, i) 
         end do
         average = 0
         do i = 1, num_diff
 
            omc_OK_diff(i) = .false.
            do j = 1, num_sat
                 if( chan_diff(i).eq.prn_sp3(j) ) then
*                    Get an approximate Range first 
                     p1c(i) = sqrt( (diff_xyz(1)-svs_xyz(1,j))**2+
     .                              (diff_xyz(2)-svs_xyz(2,j))**2+
     .                              (diff_xyz(3)-svs_xyz(3,j))**2)

                     call eph_to_xyz( diff_epoch-
     .                              (average+p1c(i))/vel_light/86400.d0,
C    .                                average/vel_light/86400.d0,
     .                                j, 'E')
*                    Now compute the Earth's rotation contributions
                     call orb_drot(svs_xyz(1,j), svs_xyz_i, p1c(i))
 
*                    Compute a rough range
                     p1c(i) = sqrt( (diff_xyz(1)-svs_xyz_i(1))**2+
     .                              (diff_xyz(2)-svs_xyz_i(2))**2+
     .                              (diff_xyz(3)-svs_xyz_i(3))**2)

*                    Get observed value

                     pc(i) = pcobs( p1o_diff(i), p2o_diff(i), j)
 
*                    accumulate the average clock offset accounting for
*                    atmospheric delay
                     call get_elev(diff_xyz, svs_xyz_i, elev)

                     omc_diff(i) =  pc(i)-p1c(i)+svs_dt(j) -
     .                       atm_del(diff_epoch, diff_geod, elev)
                     omc_OK_diff(i) = .true.
                  end if
              end do
         end do 

         num_org = num_diff
         call clean_data( omc_diff, omc_OK_diff, num_diff, 100.d0,
     .                    'REF ', ep)

*        Scan the omc values and make sure OK
C        call scan_omc(omc_diff, omc_OK_diff, data_noise, num_diff, 
C    .                 num_used )

*        Now get the average clock offset
         average = 0.d0
         num_av  = 0
         do i = 1, num_diff
            if( omc_OK_diff(i) ) then
                average = average + omc_diff(i)
                num_av  = num_av + 1
            end if
            if( ep.ge.debug_start .and. ep.le.debug_end ) then
                write(*,205) i,chan_diff(i), omc_diff(i), 
     .                  omc_OK_diff(i),num_av, p1o_diff(i)
 205            format('Diff ',i2,' PRN ',I3,' OMC ',f20.2,' OK ',L1,
     .                 ' # ',i4, ' Obs ',F20.2)
            end if
         end do
 
         if( num_av.gt.0 ) average = average / num_av
         if( ep.ge.debug_start .and. ep.le.debug_end )
     .       write(*,207) num_av, average
 207     format('Clock: ',I3,' Values, Average ',F20.2,' m')
 
****     Now do the final range computation
         do i = 1, num_diff
             do j = 1, num_sat
                 if( chan_diff(i).eq.prn_sp3(j) ) then
*                    Get an approximate Range first 
                     p1c(i) = sqrt( (diff_xyz(1)-svs_xyz(1,j))**2+
     .                              (diff_xyz(2)-svs_xyz(2,j))**2+
     .                              (diff_xyz(3)-svs_xyz(3,j))**2)
 
*                    Compute the range accounting for the propagation
*                    delay.
                     call eph_to_xyz( diff_epoch-
     .                             (average+p1c(i))/vel_light/86400.d0,
     .                                j, 'E')

*                    Now compute the Earth's rotation contributions
                     call orb_drot(svs_xyz(1,j), svs_xyz_i, p1c(i))
 
                     p1c(i) = sqrt( (diff_xyz(1)-svs_xyz_i(1))**2+
     .                              (diff_xyz(2)-svs_xyz_i(2))**2+
     .                              (diff_xyz(3)-svs_xyz_i(3))**2)

*                    Get the observed value
                     pc(i) = pcobs(p1o_diff(i), p2o_diff(i), j)

c                    write(*,200) chan_diff(i), j, p1c(i), p1o_data(i),
c    .                          p1o_data(i)-p1c(i)+svs_dt(j),
c    .                          pc    -p1c(i)+svs_dt(j)           
 
****                 compute O-C
                     call get_elev(diff_xyz, svs_xyz_i, elev)

                     omc_diff(i) = pc(i) -p1c(i)+svs_dt(j)-average -
     .                 atm_del(diff_epoch, diff_geod, elev)
                 end if
             end do
         end do
         clk_diff = average 
*        See if debug for o-minus-c wanted
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,240) 'Range REF sites differences', ep, average
 240         format('+ ',a,' Epoch ',i4,' Clock ',F13.2,' m',/,
     .              '+CHAN PRN',3x,'Obs (m)',6x,'Range (m)',3x,
     .              'Obs-Calc (m)',3x,'OMC Clock (m)')
          
             do i = 1, num_diff
                 do j = 1, num_sat
                     if( chan_diff(i).eq.prn_sp3(j) ) then
                         write(*,250) i,prn_sp3(j), pc(i), p1c(i), 
     .                                omc_diff(i),
     .                                omc_diff(i) - average
 250                     format('+',i3,1x,i4, 4F13.2)
                     end if
                 end do
             end do
         end if

      endif


***** Now try to match the channels
      if ( .not.Nodiff ) then 
         k = 0
         do i = 1, num_data
            if( omc_OK_data(i) ) then
                do j = 1, num_diff
                  if ( chan_data(i).eq. chan_diff(j) .and.
     .                 omc_OK_diff(j)   ) then
                     k = k + 1
                     omc(k) = omc_data(i)- omc_diff(j)
                     omc_OK(k) = .true.
                     omc_prn(k) = chan_data(i)
C                    print *,' i,j, omc', i,j, k, omc_data(i),
C    .                         omc_diff(j), omc(k)
                     do l = 1,6
                        apart(l,k) = apart_data(l,i)
                     end do
                  end if
                end do
             end if
         end do
      else
         k = 0
         do i = 1, num_data
            if( omc_OK_data(i) ) then
               k = k + 1
               omc(k) = omc_data(i)
               omc_OK(k) = .true.
               omc_prn(k) = chan_data(i)
               do l = 1,6
                  apart(l,k) = apart_data(l,i)
               end do
            end if
         end do
      endif


  
      num_chan = k
      call scan_omc(omc, omc_OK, data_noise, num_chan, 
     .              num_used )

      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,270) num_chan
 270     format('KF OMC num ',i4)
         do i = 1,num_chan
            write(*,280) i,omc_prn(i), omc(i), apart(1:6,i)
 280        format('KF ',i3,1x,i3,1x,F12.2,' APart ',6(F6.2,1x))
         enddo
      end if  

****  Now increment the filter
      do i = 1,6
          cov_parm(i,i) = cov_parm(i,i) + wn(i)
      end do
 
***** Start computing the Kalman gain.  Get the acat matrix
      call acat(apart, cov_parm, 6, num_chan, acat_mat, temp_gain,
     .         max_sat)
 
****  Now add the data noise (assume 10 meter data noise)
      num_used = 0
      do i = 1, num_chan
          if( omc_OK(i) ) then
              num_used = num_used + 1
              acat_mat(i,i) = acat_mat(i,i) + data_noise**2 
          else
              acat_mat(i,i) = acat_mat(i,i) + 1.d6
          end if
      end do
 
***** Now invert this results
      call invert_vis(acat_mat, dummy, scale, ipivot, 
     .                num_chan, max_sat, 0)


***** Now form the Kalman Gain
      do i = 1,6
          do j = 1, num_chan
              kgain(i,j) = 0
              do k = 1, num_chan
                  kgain(i,j) = kgain(i,j) +
     .                    temp_gain(i,k)*acat_mat(k,j)
              end do
          end do
      end do
 
C     call write_mat('kgain', kgain, 6, num_chan, 6, max_sat)

***** Update the parameter estimates
      do j = 1,num_chan
          dp(j) = omc(j) 
          do i = 1,6
              dp(j) = dp(j) -apart(i,j)*sol_vec(i)
          end do
      end do
      do i = 1,6
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

* MOD TAH 180318: Update the residuals for the new parameter
*     estimates
      do j = 1, num_chan
         do k = 1, 6
            dp(j) = dp(j) - apart(k,j)*dx(k)
         end do
      end do

*     Accumulate mean and RMS
      do j = 1, num_chan
         k = chan_data(j)
         do k = 1, num_sat
            if( prn_sp3(k).eq.omc_prn(j) ) then
               svres_sum(k) = svres_sum(k) + dp(j) 
               svres_var(k) = svres_var(k) + dp(j)**2
               svres_num(k) = svres_num(k) + 1
               if( ep.ge.debug_start .and. ep.le.debug_end ) 
     .         write(*,290) ep, omc_prn(j), dp(j), omc(j)
 290           format('RES EP ',i5,' PRN ',I3,' Res ',2F10.2,' m')
               exit
            end if
         end do
      end do

       
      if( num_used.gt.0 ) dchi = dchi / num_used 
c     if( dchi.gt.20 ) RETURN

      chi = chi + dchi*num_used
      nchi = nchi + num_used
      
      do i = 1,6
          sol_vec(i) = sol_vec(i) + dx(i)
      end do
 
***** Now update covariance matrix
      do i = 1,6
          do j = 1,6
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
      call mjd_to_ymdhms(data_epoch, date, sectag)

C     sec10 = nint(sectag*10.0)+date(5)*600
C     if( abs(sec10-nint(sec10/out_spacing*10)*out_spacing*10).lt.
C    .                                         1.d0  ) then
      sod = date(4)*3600+date(5)*60+sectag 
      if( abs(sod-nint(sod/out_spacing)*out_spacing).lt.1.d-2 ) then

*****     Need to output.  See what coordinates
          do i = 1,3
             pos_xyz_fin(i) = site_xyz(i) + sol_vec(i)
             pos_xyz_adj(i) = sol_vec(i) 
             do j = 1,3
                covar(i,j) = cov_parm(i,j)
             end do
          end do

          if( out_type(1:1).eq.'N') then
              call rotate_geod(pos_xyz_adj,pos_neu_adj,'XYZ','NEU',
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
     .               sol_vec(4)+clk_data-clk_diff, sqrt(cov_parm(4,4)),
     .               sol_vec(5), sqrt(cov_parm(5,5)),
     .               sol_vec(6), sqrt(cov_parm(6,6)),
     .               num_used, num_chan, dchiout, ep
 300          format(i5,4i3,1x,f8.3,1x,3(f10.2,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2,4F7.3, 2I3,1x, f8.2,1x,I5,' NEU')
           else
              write(*,310) date, sectag, (sol_vec(i), 
     .               sqrt(cov_parm(i,i)), i = 1,3),
     .               sol_vec(4)+clk_data-clk_diff, sqrt(cov_parm(4,4)),
     .               sol_vec(5), sqrt(cov_parm(5,5)),
     .               sol_vec(6), sqrt(cov_parm(6,6)),
     .               num_used, num_chan, dchiout, ep
 310          format(i5,4i3,1x,f8.3,1x,3(f10.2,1x,f8.2), 1x,
     .                 f12.2,1x,f9.2, 4F7.3, 2I3,1x, f8.2,1x,I5,' XYZ')
           end if
      end if

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
C     if( num_used.ne.num_chan ) then
C         write(*,100) (i, omc_OK(i), omc(i), i=1,num_chan)
C100      format('*',100(I3,1x,L1,1x,F10.2))
C     end if

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
      include 'svsp3.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(max_data_types), ierr, 
     .          jerr,i,j, k,id, trimlen,n1,nel
 
      integer*4 nblk   ! Number of lines to read
     .,         nend, nstr ! End and Start num_data entries

      integer*4 conoff  ! Function to return constellation offset value
      logical  OK       ! Set true if the are keeping data being read
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(max_data_types)
 
*   eof     - Indicates end of file
*   noP2    - Set true if no P2 data available
*   eph     - Set true if satellite in ephemeris

      logical eof, noP2, eph
 
      integer*4 ep

      character*256 line
      character*1 cr
      character*1 stype(max_sat)  ! Type of satellite.  Only G GPS used

      character*2 saved_dtype(2,max_sat)   ! Obserable saved for L1, L2

      data ep / 0 /

* MOD TAH 151021: Replaced code with more general version that 
*     reads any number of channels.  We need to read first line to 
*     get number of channels
      ep = ep + 1
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
      call ymdhms_to_mjd( date, sectag, data_epoch)
    
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
      P1o_data = 0
      P2o_data = 0
      vals = 0

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
              if( data_types(j).eq.'P1' .and.vals(j).ne.0 ) then
                  p1o_data(i) = vals(j)
                  saved_dtype(1,i) = data_types(j)
              endif
              if( data_types(j).eq.'L1' ) l1o_data(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( data_types(j).eq.'L2' ) l2o_data(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( data_types(j).eq.'P2'.and.vals(j).ne.0  ) then
                  p2o_data(i) = vals(j)
                  saved_dtype(2,i) = data_types(j)
                  noP2 = .false.
              endif
*             See if C1 is available if P1 is not
              if( data_types(j).eq.'C1' .and.vals(j).ne.0 .and.
     .            p1o_data(i).eq.0 ) then
                  p1o_data(i) = vals(j)
                  saved_dtype(1,i) = data_types(j)
              endif

              if( data_types(j)(1:2).eq.'C2' .and.vals(j).ne.0 ) then
*                 Special treatment for Beidou (C2 is the higher frequency)
                  if( stype(i).eq.'C' ) then
                      P1o_data(i) = vals(j)
                      saved_dtype(1,i) = data_types(j)
                  elseif ( P2o_data(i).eq.0 ) then
                      P2o_data(i) = vals(j)
                      saved_dtype(2,i) = data_types(j)
                      noP2 = .false.
                  endif
              end if
              if( data_types(j)(1:2).eq.'C5' .and.vals(j).ne.0 .and.
     .            stype(i).ne.'G' ) then
                  P2o_data(i) = vals(j)
                  saved_dtype(2,i) = data_types(j)
                  noP2 = .false.
              end if
              if( stype(i).eq.'C' ) then
                  if( data_types(j).eq.'C6'.and.vals(j).ne.0  ) then
                     p2o_data(i) = vals(j)
                     saved_dtype(2,i) = data_types(j)
                     noP2 = .false.
                  endif
              endif  

          end do
      end do

      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         do i = 1,num_data
            write(*,150) ep, i,  stype(i), chan_data(i),
     .                   P1o_data(i), saved_dtype(1,i), 
     .                   P2o_data(i), saved_dtype(2,i)
 150        Format('RAWDATA EP ',i3,' C ',I2,1x,a1,1x,I3,1x,
     .             2(F14.3,1x,a,1x))
         enddo
      endif


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
*         See if we are processing these data and if we are
*         set the PRN 
          OK = .false.
          if( index(gnsstouse,stype(i)).gt.0 ) then
*             This type is in list so update the PRN
*             number
              chan_data(i) = chan_data(i)+conoff(stype(i))
              OK = .true.
          endif 

*         Make sure we have an empheris for the data
          eph = .false.
          do j = 1, num_sat
             if( chan_data(i).eq. prn_sp3(j) ) then
                eph = .true.
                exit
             endif
          end do
          OK = OK .and. eph

          if( .not. OK ) then
              do j = i + 1, num_data
                 chan_data(j-1) = chan_data(j)
                 p1o_data(j-1) = p1o_data(j)
                 p2o_data(j-1) = p2o_data(j)
                 stype(j-1) = stype(j)
              end do
              i = i - 1
              num_data = num_data - 1
          else if( abs(p1o_data(i)-p2o_data(i)).gt.200.d0 .and.
     .                 p2o_data(i).ne.0.d0 .and. rep_dels ) then
              write(*,220) date, sectag, chan_data(i), 
     .                     abs(p1o_data(i)-p2o_data(i))
 220          format('DELETING BASE DATA AT ',i4,4(1x,i2),1x,F6.2,
     .               ' PRN ',i3,' P1-P2= ',f12.1,' m different')
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
      include 'svsp3.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(max_data_types), ierr, 
     .          jerr,i,j, k,id, trimlen,n1,nel
 
      integer*4 nblk   ! Number of lines to read
     .,         nend, nstr ! End and Start num_chan entries

      integer*4 conoff  ! Function to return constellation offset value
      integer*4 ep      ! Epoch counter

      logical  OK       ! Set true if the are keeping data being read
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(max_data_types)
 
*   eof     - Indicates end of file
*   noP2    - Set true if no P2 data available

      logical eof, noP2
 
      character*256 line
      character*1 cr
      character*1 stype(max_sat)  ! Type of satellite.  Only G GPS used

      character*2 saved_dtype(2,max_sat)   ! Obserable saved for L1, L2

      data ep / 0 / 

* MOD TAH 151021: Replaced code with more general version that 
*     reads any number of channels.  We need to read first line to 
*     get number of channels
      ep = ep + 1
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
      call ymdhms_to_mjd( date, sectag, diff_epoch)
    
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

              if( diff_types(j).eq.'P1' .and.vals(j).ne.0  ) then
                   p1o_diff(i) = vals(j)
                   saved_dtype(1,i) = diff_types(j)
              end if
              if( diff_types(j).eq.'L1' ) l1o_diff(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( diff_types(j).eq.'L2' ) l2o_diff(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( diff_types(j).eq.'P2'.and.vals(j).ne.0  ) then
                  p2o_diff(i) = vals(j)
                  saved_dtype(2,i) = diff_types(j)
                  noP2 = .false.
              endif

*             See if C1 is available if P1 is not
              if( diff_types(j).eq.'C1' .and.vals(j).ne.0 .and.
     .            p1o_diff(i).eq.0 ) then
                  p1o_diff(i) = vals(j)
                  saved_dtype(1,i) = diff_types(j)
              endif

              if( diff_types(j)(1:2).eq.'C2' .and.vals(j).ne.0 ) then
*                 Special treatment for Beidou (C2 is the higher frequency)
                  if( stype(i).eq.'C' ) then
                      P1o_diff(i) = vals(j)
                      saved_dtype(1,i) = diff_types(j)
                  else
                      P2o_diff(i) = vals(j)
                      saved_dtype(2,i) = diff_types(j)
                      noP2 = .false.
                  endif
              end if
              if( diff_types(j)(1:2).eq.'C5' .and.vals(j).ne.0 .and.
     .            stype(i).ne.'G' ) then
                  P2o_diff(i) = vals(j)
                  saved_dtype(2,i) = diff_types(j)
                  noP2 = .false.
              end if
              if( stype(i).eq.'C' ) then
                  if( diff_types(j).eq.'C6'.and.vals(j).ne.0  ) then
                     p2o_diff(i) = vals(j)
                     saved_dtype(2,i) = diff_types(j)
                     noP2 = .false.
                  endif
              endif  

          end do
      end do

      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         do i = 1,num_data
            write(*,150) ep, i,  stype(i), chan_diff(i),
     .                   P1o_diff(i), saved_dtype(1,i), 
     .                   P2o_diff(i), saved_dtype(2,i)
 150        Format('RAWDIFF EP ',i3,' C ',I2,1x,a1,1x,I3,1x,
     .             2(F14.3,1x,a,1x))
         enddo
      endif


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
*         See if we are processing these data and if we are
*         set the PRN 
          OK = .false.
          if( index(gnsstouse, stype(i)).gt.0 ) then
*             This type is in list so update the PRN
*             number
              chan_diff(i) = chan_diff(i)+conoff(stype(i))
              OK = .true.
          endif          
          if( .not. OK ) then
!             write(*,210) date, k, stype(i)
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
          else if( abs(p1o_diff(i)-p2o_diff(i)).gt.200.d0 .and.
     .                 p2o_diff(i).ne.0.d0 .and. rep_dels ) then
              write(*,220) date, sectag, chan_diff(i), 
     .                     abs(p1o_diff(i)-p2o_diff(i))
 220          format('DELETING REF  DATA AT ',i4,4(1x,i2),1x,F6.2,
     .               ' PRN ',i3,' P1-P2= ',f12.1,' m different')
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
 
      include 'svsp3.h'
 
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
     .                 (data_types(str+i), i=1,rem)
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
          stop 'svsp3: Wrong type data file'
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
 
      include 'svsp3.h'
 
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
              write(*,120) trim(data_file), diff_rxver
 120          format('Data File ',a,' Rinex Version ',F5.2)
              if( diff_rxver.ge.3.0 ) then
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
      include 'svsp3.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
 
 
      integer*4 len_run, i, rcpar
      integer*4 inend, inan  !  Index at end of name
 
*  runstring   - Elements of runstring
 
      character*256 runstring
      character*256 ucrun 
      character*32  opt_str ! Combination of GNSS type and
             ! analsys type. For example, add +GR_P1 to file name

 
*     Get the first runstring parameter

      gnsstouse = 'G'
      anal_typ = 'PC'

*     Initialize
      wn = 0.d0
      wn(4) = 1.d6
      data_noise = 100.d0
      site_xyz = 0.d0
      diff_xyz = 0.d0
      out_spacing = 3600   ! 1-hr spaceing
      out_type = 'XYZ'
      debug_start = 0   ! No debug 
      debug_end   = 0
      proc_start  = 0   ! Process all data
      proc_end    = 0
      rep_dels    = .false.   ! Turn on with -rep_dels options
      atm_apr_sig = 0   ! Set zero so values not estimated (-atm 
                        ! option can be used to set)

      len_run = rcpar(1,nav_file)
      if( len_run.le.0 ) then
          call proper_runstring('svsp3.hlp','svsp3/nav file',1)
      end if
 
*     See if data difle names passed
      len_run = rcpar(2,data_file)
      if( len_run.le.0 ) then
          call proper_runstring('svsp3.hlp','svsp3/data file',1)
      end if

*     Start decoding the GNSS and analysis type: String the follows
*     the either the O or RMX at the end of the name
      ucrun = data_file
      call casefold(ucrun)
      inend = index(ucrun,'O+') 
      if( inend.eq.0 ) then   ! Try long name
          inend = index(ucrun,'X+') 
      end if
      if( inend.gt.0 ) then
          opt_str = ucrun(inend+2:)
          data_file = data_file(1:inend)
          inan = index(opt_str,'_')
          if( inan.gt.0 ) then
             anal_typ = opt_str(inan+1:)
             if( inan.gt.2 ) then
                gnsstouse = opt_str(1:inan-1)
             endif
          else
             gnsstouse = opt_str
          endif
      end if

*     See if diff difle names passed
      len_run = rcpar(3,diff_file)
      if( len_run.le.0 ) then
          diff_file = 'NONE'
          write(*,'(a)') 'No difference RX file, Point Positioning'
          Nodiff = .true.
      else
          Nodiff = .false.
*         call proper_runstring('svsp3.hlp','svsp3/diff file',1)
      end if

****  See if set runstring or new optional string

*     Get process noise for the position
      len_run = rcpar(4, runstring)
      if( runstring(1:1).ne.'-' ) then   ! Old fixed runstring
         if( len_run.gt.0 ) then
             read(runstring,*) wn(1)
             wn(2) = wn(1)
             wn(3) = wn(1)
         else
             do i = 1,3
                wn(i) = 0.d0
             end do
         end if

*        Get clock process noise
         len_run = rcpar(5, runstring)
         if( len_run.gt.0 ) then
             read(runstring,*) wn(4)
         else
             wn(4) = 1.d6
         end if
         wn(5) = 0.000d0
         wn(6) = 0.000d0

*        Get data noise
         len_run = rcpar(6, runstring)
         if( len_run.gt.0 ) then
             read(runstring,*) data_noise 
         else
             data_noise = 100.d0
         end if

*        Now see if positions passed 
         do i = 1,3
            len_run = rcpar(6+i,runstring)
            if(  len_run.gt.0 ) then
                read(runstring,*) diff_xyz(i)
            else
                diff_xyz(i) = 0.d0
            end if
         end do

*        Get data spacing (results output when modulo this value
*        past the hourt
         len_run = rcpar(10, runstring)
         if( len_run.gt.0 ) then
              read(runstring,*) out_spacing
         else
              out_spacing = 3600.0d0  ! 1-hr
         end if

*        Get output coordinate type
         len_run = rcpar(11, out_type )
         if( len_run.eq.0 ) out_type = 'XYZ'

* MOD    SCM Added apriori site coordinates on command line. 6/1/00
*        Now see if new site positions passed 
         do i = 1,3
            len_run = rcpar(11+i,runstring)
            if(  len_run.gt.0 ) then
                read(runstring,*) site_xyz(i)
            else
                site_xyz(i) = 0.d0
            end if
         end do

*        See if debug wanted
         len_run= rcpar(15, runstring)
         if( len_run.gt.0 ) then
             read(runstring,*) debug_start
         else
             debug_start = 0   
         end if

*        See if debug or processing start (+ve if debug)
         if( debug_start.lt.0 ) then
             proc_start = -debug_start 
             debug_start = 0
         else
             proc_start = 0 
         endif

*        See if debug wanted
         len_run= rcpar(16, runstring)
         if( len_run.gt.0 ) then
             read(runstring,*) debug_end   
         else
*            See if end of processing or number of epochs
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
         write(*,320) trim(data_file), gnsstouse, anal_typ
 320     format('For Data File ',a,' GNSS ',a,' Data ',a)  
     
         if( proc_start.gt.0 ) then
            write(*,340) proc_start, proc_end
 340        format('Processing data between epoch ',i6, 
     .             ' and ',i6) 
         endif

          if( debug_start.gt.0 ) then
            write(*,360) debug_start, debug_end
 360        format('Debug output between epoch ',i6, 
     .             ' and ',i6) 
         endif
      ELSE     ! New option driven runstring.

****     Loop over runrstring looking for options.
         call opt_runstr
      ENDIF



****  Thats all
      return
      end


 
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'svsp3.h'
 
*   t       - Time for computation (day number)
 
 
      real*8 t
 
*   i       - Satellite number
 
 
      integer*4 i
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j, bin, ein
      real*8 err, rot_mat(3,3), xa(100), ya(100), za(100)
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
      bin = j - 5
      if( bin.le.0 ) bin = 1
      ein = bin + 9
      if( ein.gt.num_sp3 ) then
          ein = num_sp3
          bin = ein - 9
      end if

***** Now do interporlation
      do j = 0, 9
          xa(j+1) = sp3_xyz(1,i,bin+j)
          ya(j+1) = sp3_xyz(2,i,bin+j)
          za(j+1) = sp3_xyz(3,i,bin+j)
      end do 
      call polint(sp3_time(bin), xa, 10, t, 
     .            svs_xyz(1,i), err)
C     if( abs(err).gt.0.01 ) write(*,500) i,t, 'X', err
      call polint(sp3_time(bin), ya, 10, t, 
     .            svs_xyz(2,i), err)
C     if( abs(err).gt.0.01 ) write(*,500) i,t, 'Y', err

      call polint(sp3_time(bin), za, 10, t, 
     .            svs_xyz(3,i), err)
C     if( abs(err).gt.0.01 ) write(*,500) i,t, 'Z', err

  500 format('*Interpolation error PRN ',i3,' Epoch ',f12.4,
     .       ' Component ',a1,' Magnitude ', f8.4,' m')
 
*     Now convert to latitude and longitude
      if( abs(svs_xyz(1,i)+svs_xyz(2,i)+svs_xyz(3,i)).gt.1.) 
     .   call xyz_to_geod( rot_mat, svs_xyz(1,i), svs_loc(1,i))
 
****  Thats all
      return
      end
 
CTITLE SP3_to_CLK
 
      subroutine sp3_to_clk(t, i)

      implicit none
 
*     Routine to compute clock error based on time by interpolating
*     clock times.  Clock adjustment is returned in meters.
 
      include 'svsp3.h'
      include '../includes/const_param.h'
 
*   t       - Time for computation (day number)
 
      real*8 t
 
*   i       - Satellite number
      integer*4 i
 
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j, bin, ein
      real*8 err, xa(100)
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
      bin = j - 5
      if( bin.le.0 ) bin = 1
      ein = bin + 9
      if( ein.gt.num_sp3 ) then
          ein = num_sp3
          bin = ein - 9
      end if

***** Now do interporlation
      do j = 0, 9
          xa(j+1) = sp3_clk(i,bin+j)
      end do 
      call polint(sp3_time(bin), xa, 10, t, 
     .            svs_dt(i), err)

C      svs_dt(i) = svs_dt(i)*1.d-6*vel_light
C MOD TAH 18032: SP3_LIB converts time to seconds, so use 
C     directly.
      svs_dt(i) = svs_dt(i)*vel_light
      err = err*1.d6

      if( abs(err).gt.0.1 ) write(*,500) prn_sp3(i),t, err

  500 format('*Interpolation error PRN ',i3,' Epoch ',f12.4,
     .       ' Clock.  Error ', f8.4,' us')
  
****  Thats all
      return
      end
      
      subroutine polint(xa,ya,n,x,y,dy)

      implicit none

      integer*4 nmax, n, i,m, ns
      parameter (nmax=21) 
      real*8 xa(n), ya(n), c(nmax),d(nmax), dif, dift, 
     .       den, w, hp, ho, x,y, dy

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
C         if(den.eq.0.)pause
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

      subroutine get_elev( pos, svs, elev)

      implicit none

      include '../includes/const_param.h'

      real*8 pos(3), svs(3), elev

      real*8 dpos(3), dot, mag1, mag2
      integer*4 i

      do i = 1,3
          dpos(i) = svs(i) - pos(i)
      end do
      mag1 = 0.d0
      mag2 = 0.d0
      dot = 0.d0
      do i = 1,3
          dot = dot + pos(i)*dpos(i)
          mag1 = mag1 + pos(i)**2
          mag2 = mag2 + dpos(i)**2
      end do
      dot = dot / sqrt(mag1*mag2)
      elev = 90.d0 - acos(dot)*180/pi
      end 
 
      subroutine mit_dry( elev, part )

      implicit none

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


CTITLE OPT_RUNSTR
 
      subroutine opt_runstr
 
      implicit none

*     Routine to decode the new option driven runstring for svsp3
*
 
      include '../includes/const_param.h'
      include 'svsp3.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
*  nr          - Current runstring number.
 
 
      integer*4 len_run, i, rcpar, nr
 
      logical done    ! Set true when empty runstring.

*  runstring   - Elements of runstring
 
      character*256 runstring
 

* SVSP3 command line options
* -proc_noise <position m^2/ep> <clock m^2/ep> <range noise (m)>
* -ref_apr <X> <Y. <Z> -- Apriori coordinates for Reference (2nd) site
* -dat_apr <X> <Y. <Z> -- Apriori coordinates for first site.
* -out  <out spacing> <out type (XYZ/NEU)> 
* -debug <start epoch> <end epoch>
* -span <start epoch> <end epoch>
* -gnss <systems>  G-GPS, R-Glonass, E-Galileo, C-Beidou
* -anal_type <Type> From PC, P1, P2, P5
* -rep_dels  -- Report deletes

 

****  Start decoding
      nr = 3    ! 1 back from current position in runstring
      done = .false.
      do while ( .not. done ) 
         nr = nr + 1
         len_run = rcpar(nr, runstring)
         if ( len_run.eq. 0 ) then
            done = .true.
         else
*           Start decoding
            call casefold(runstring)
            if( runstring(1:2).eq.'-P' ) then
*              -proc_noise <position m^2/ep> <clock m^2/ep> <range noise (m)>
               nr = nr + 1
               len_run = rcpar(nr, runstring)
               read(runstring,*) wn(1)
               wn(2) = wn(1)
               wn(3) = wn(1)

*              Get clock process noise
               nr = nr + 1
               len_run = rcpar(nr, runstring)
               read(runstring,*) wn(4)
               wn(5) = 0.000d0
               wn(6) = 0.000d0

*              Get data noise
               nr = nr + 1
               len_run = rcpar(nr, runstring)
               if( len_run.gt.0 ) then
                   read(runstring,*) data_noise 
               end if


            elseif ( runstring(1:4).eq.'-REF' ) then
*              -ref_apr <X> <Y. <Z> -- Apriori coordinates for Reference (2nd) site
               do i = 1,3
                  nr = nr +1 
                  len_run = rcpar(nr,runstring)
                  read(runstring,*) diff_xyz(i)
               end do
 
            elseif ( runstring(1:3).eq.'-DA' ) then
*              -dat_apr <X> <Y. <Z> -- Apriori coordinates for first site.
               do i = 1,3
                  nr = nr +1 
                  len_run = rcpar(nr,runstring)
                  read(runstring,*) site_xyz(i)
               end do
 
            elseif ( runstring(1:2).eq.'-O' ) then
*               -out  <out spacing> <out type (XYZ/NEU)> 
                nr = nr + 1
                len_run = rcpar(nr, runstring)
                read(runstring,*) out_spacing
                nr = nr + 1
                len_run = rcpar(nr, out_type )
                call casefold(out_type)

            elseif ( runstring(1:3).eq.'-DE' ) then
*               -debug <start epoch> <end epoch>
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) debug_start
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) debug_end

            elseif ( runstring(1:2).eq.'-S' ) then
*               -span <start epoch> <end epoch>
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) proc_start
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) proc_end

            elseif ( runstring(1:2).eq.'-G' ) then
*               -gnss <systems>  G-GPS, R-Glonass, E-Galileo, C-Beidou
                nr = nr + 1
                len_run = rcpar(nr, gnsstouse) 
                call casefold(gnsstouse)

            elseif ( runstring(1:3).eq.'-AN' ) then
*               -anal_type <Type> From PC, P1, P2, P5
                nr = nr + 1
                len_run = rcpar(nr, anal_typ )
                call casefold(anal_typ)

            elseif ( runstring(1:3).eq.'-AT' ) then
*               -atm_apr_sig <Site 1> <Site 2>
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) atm_apr_sig(1)
                nr = nr + 1
                len_run = rcpar(nr, runstring) 
                read(runstring,*) atm_apr_sig(2)

            elseif ( runstring(1:4).eq.'-REP' ) then
*               -rep_dels
                rep_dels = .true.
            else
                 write(*,'("UNKNOWN OPTION ",a)') trim(runstring)
            endif
         endif
      end do

*     Report setup for run
      write(*,320) trim(data_file), gnsstouse, anal_typ
 320  format('For Data File ',a,' GNSS ',a,' Data ',a)       

      if( proc_start.gt.0 ) then
         write(*,340) proc_start, proc_end
 340     format('Processing data between epoch ',i6, 
     .          ' and ',i6) 
      endif

      if( debug_start.gt.0 ) then
         write(*,360) debug_start, debug_end
 360     format('Debug output between epoch ',i6, 
     .          ' and ',i6) 
      endif

      write(*,420) wn(1), wn(4), data_noise
 420  format('* PROCESS NOISE Site ',F10.2,' m^2; CLOCK ',
     .       F10.2,' m^2; Data Noise ',F6.2,' m')
      write(*,440) atm_apr_sig
 440  format('* ATM APRIORI SIGMAS ',2(F10.2,1x),' m')


***** Thata all
      return
      end



CTITLE PCOBS

      real*8 function pcobs(p1o, p2o, sn )

      implicit none

*     Function to return the requested linear combination for
*     processing.

      include 'svsp3.h'

* PASSED Values
      real*8 p1o  ! Observed range at L1 
     .,      p2o  ! Observed range at L2

      integer*4 sn     ! Satellite number in list of satellites

* LOCAL
      real*8 pcf1, pcf2   ! Factors to form ion-free combination
                       ! based on the frequencies for the specific
                       ! satellite

****  See what has been requested but only compute PC is L2 range is 
*     not zero 

      if( p2o.gt.0 .and. anal_typ.eq.'PC' ) then
*        Form coefficients specific to the frequnency of 
*        this satellite.
         pcf1 =  fL1(sn)**2/(fL1(sn)**2-fL2(sn)**2) 
         pcf2 = -fL2(sn)**2/(fL1(sn)**2-fL2(sn)**2) 
         pcobs = p1o*pcf1 + p2o*pcf2
      elseif( anal_typ.eq.'P1' ) then
         pcobs = p1o
      elseif( anal_typ.eq.'P2' ) then
         pcobs = p2o
      else     ! Default to P1 data.
         pcobs = p1o
      end if

****  Thats all
      return
      end

CTITLE ORB_DROT

      subroutine orb_drot( svs_xyz, svs_xyz_i, pc )

      implicit none

*     Routine to apply rotation to satellite positions
*     to account for Earth rotation during light propagation
*     time (pc expressed as range).

      include '../includes/const_param.h'

* PASSED
      real*8  svs_xyz(3)   ! Earth fixed XYZ at trsnsmit time
     .,       svs_xyz_i(3) ! Earth fixed XYZ rotated to account
                           ! light progragation time delay on the
                           ! rotation of the Earth. (OUTPUT)
     .,       pc           ! Range to satellite in meter (light
                           ! propagation range.

* LOCAL
      real*8 dtheta        ! Amount of Earth rotation in light
                           ! progation time (rad)

     

****  Now compute the Earth's rotation contributions
      dtheta = -(pc/vel_light/86400.d0)*2*pi

*     Rotate about the Z-axis (small angle approximation)
      svs_xyz_i(1) = svs_xyz(1) - dtheta*svs_xyz(2)                 
      svs_xyz_i(2) = dtheta*svs_xyz(1) + svs_xyz(2)
      svs_xyz_i(3) = svs_xyz(3)

***** Thats all
      return
      end


CTITLE CLEAN_DATA

      subroutine clean_data(pomc, pomc_OK, num, scale, type, ep)

      implicit none

*     Routine to look at the omc vector and delete any large
*     outliers from both omc and the original data.  The scale
*     factor gives the multiplier on the data noise to assess
*     what is deleted. 
*     The type string says if this is 'data' or 'diff' sites.

      include 'svsp3.h'

* PASSED
      integer*4 num       ! Number of values in OMC vector
     .,         ep        ! Epoch number for reporting edits if ask

      real*8 pomc(num)    ! O-minus-C values in meters passed array
     .,      scale        ! Mulitplier on data noise

      logical pomc_OK(num)   ! Logically set to true if good; passed array

      character*(*) type    ! Either 'data' or 'diff' depending on
                            ! which site this is.
    
* LOCAL

      integer*4 nu   ! Number used after editing
     .,         i, j    ! Loop counters
     .,         ib      ! Index of biggest residual

      real*8 av         ! Average value of residuals
     .,      big        ! Largest residual (compared with scale*data_noise)

      logical OK        ! Set true when either only 2 data points or no
                        ! more outliers are found.

***** Keep iterating until no more edits or we run out of data
      OK = .false.
      pomc_OK = .true.

      do while ( .not. OK )
         av = 0.d0
         nu = 0
         do i = 1, num
            if( pomc_OK(i) ) then  ! Good, so add to average
               av  = av + pomc(i)
               nu  = nu + 1
            end if
         end do
         if( nu.le.2 ) OK = .true.

****     Compute average and RMS
         av = av/nu

*        Find largest residual
         big = 0
         do i = 1, num
            if( pomc_OK(i) ) then    
               if( abs(pomc(i)-av).gt.big ) then
                  big = abs(pomc(i)-av)
                  ib = i
               endif
            endif
          end do
          if( big.gt. scale*data_noise ) then
*             Outlier so remove
              pomc_OK(ib) = .false.
              if( rep_dels ) then    ! Report the delete
                 write(*,220) type, ep, chan_data(ib), big
 220             format('DELETING ',a,' DATA EPOCH ',I5,' PRN ',i3,
     .                  ' MAGNITUDE ',F20.2,' m')
              end if 
          else
*             All done
              OK = .true.
          end if
      end do

****  If we deleted data, remove from list of data
      if( nu.lt.num ) then
*        Delete data: BASE case
         if( type(1:2).eq.'BA' ) then    ! Data array
*           Run backwards so only delete later data
            do i = num,1,-1
               if( .not.pomc_OK(i) ) then   ! Bad so delete
                   do j = i+1, num_data
                      chan_data(j-1) = chan_data(j)
                      p1o_data(j-1) = p1o_data(j)
                      p2o_data(j-1) = p2o_data(j)
                   end do
                   num_data = num_data - 1
                endif 
            end do
         else      ! Work with the diff data arrays 
*           Run backwards so only delete later data
            do i = num,1,-1
               if( .not.pomc_OK(i) ) then   ! Bad so delete
                   do j = i+1, num_diff
                      chan_diff(j-1) = chan_diff(j)
                      p1o_diff(j-1) = p1o_diff(j)
                      p2o_diff(j-1) = p2o_diff(j)
                   end do
                   num_diff = num_diff - 1
                endif 
            end do
         end if
      end if

****  Now align the omc and omc_OK values with the editing
      i = 0
      do while ( i.lt. nu )
         i = i + 1
         if(  .not.pomc_OK(i) ) then
*            Move data down to replace this value
             do j = i+1, num
                pomc_OK(j-1) = pomc_OK(j)
                pomc(j-1)    = pomc(j)
             end do
             i = i - 1
         endif
      end do

***** Thats all
      return
      end

CTITLE ATM_DEL

      real*8 function atm_del( mjd, geod, elev )

      implicit none

      include '../includes/const_param.h' 

*     Function to return a GPT/GMF model for the atmospheric delay
*

* PASSED  
      real*8 mjd       ! MJD at which delay is to be computed
     .,      geod(3)   ! Geodetic coordinates (colat, long, ht; rads/m) 
     .,      elev      ! Elevation angle (degs): NOTE UNITS

* LOCAL
      real*8 press, temp_C  ! Pressure and temperature
     .,      undu      ! Geoid height (not used).
     .,      dry_zen, wet_zen   ! Dry and Wet zenith delay
     .,      zend      ! Zenith angle (rad)
     .,      dry_map, wet_map   ! Dry and wet mapping function
     .,      lat       ! Latitude (rad)

***** Call GPT to get pressute and temperatue

      lat = pi/2-geod(1)

      call gpt(mjd,lat,geod(2),geod(3), press, temp_C, undu) 

*     Compute dry zenith delay
      call dry_saas_zen( temp_C,lat ,geod(3), press,   dry_zen)

*     Get mapping function
      call gmf(mjd,lat,geod(2),geod(3), (90-elev)*pi/180, 
     .        dry_map, wet_map)

*     Compute the dry delay (set wet to 0.1 meters) 
      wet_zen = 0.1       
      atm_del = dry_zen*dry_map + wet_zen*wet_map

***** That all
      return 
      end

CTITLE SET_GEOD

      subroutine set_geod

      implicit none

*     Routine to compute geodetic positions based on apriori XYZ coords

      include 'svsp3.h'

*     Use XYZ_to_GEOD to get Geodetic co-lat, longitude (rads) and height from
*     XYZ coordinates


* LOCAL
      real*8 rot_mat(3,3)    ! Rotation to topocentric frame

      call XYZ_to_GEOD(rot_mat, site_xyz, site_geod)
      if( diff_xyz(1).ne.0 ) then
          call XYZ_to_GEOD(rot_mat, diff_xyz, diff_geod)
      else
          diff_geod = 0
      endif

****  Thats all
      return
      end


