      program ke_anal
 
*     This program will carry out group and phase delay analysis for
*     the list of experiments passed in the runstring.
*
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include 'ke_anal.h'
 
*     Local variables
 
*   ierr    - IOSTAT error
*   rcpar   - Read runstring rouitne
*   trimlen - Length of string
*   ir      - Counter for position in runstring
*   i       - Loop counter
*   luo     - Unit number for output
*   date(5) - Epoch of observation
*   unw     - Local unweight flag (Frnge_qual 8 or 9 )
*   lunw    - Line unweight indicator (1 if good/ 2 if bad)
 
      integer*4 ierr, rcpar, trimlen, ir, i, luo, date(5), unw, lunw
 
*   sectag  - Seconds tag of observation
*   dgrp, sgrp  - Group - phase (mm) and group delay sigma (mm)
 
      real*8 sectag, dgrp, sgrp
 
*   out_data    - NAme of the output file (or unit)
 
 
      character*128 out_data
 
*     Get the data output file name
 
      ir = rcpar(1, out_data)
      if( ir.gt.0 ) then
          luo = 200
          call open_lu(luo, out_data, ierr, 'append')
      else
          call proper_runstring('ke_anal.hlp', 6,1)
      end if
 
 
*     Initialize the summation variables
 
      call init_sums( sig_sums, elev_sums, rate_sums, frng_sums,
     .            overal_sums, max_bins )
 
*     Now start loop over the names of the files in runstring
      ir = 2
      ne = 1
      do while ( rcpar(ir, KO_names(ne)).gt.0 )
 
          call open_KalObs( KO_names(ne), ierr )
          call report_error('FmpOpen',ierr,'open',KO_names(ne), 0,
     .            'ke_anal')
 
*                                     ! Continue processing
          if( ierr.eq.0 ) then
              call rw_KalObs_header('R', ierr)
 
*             Write header to output
              write(luo, 100) ko_names(ne)(1:trimlen(ko_names(ne)))
  100         format('* Group phase comparision from ',a,/,
     .               '*   Date     Group-phase (mm)    +-  source unw',
     .               ' lunw Qual Code    Rate (ps/s)  El (dg)',
     .               ' Az (dg) Temp (C)' )
 
              call init_lsum
 
*             Now loop over data
              do i = 1, num_obs
                  call rw_KalObs_block('R', 'DATA',site, ierr, i)
 
*                 Get the group delay difference and sigma
                  dgrp = (observable(1) + grp_amb*num_grp_amb -
     .                   observable(2) - phs_amb*num_phs_amb +
     .                   (feed_rot_cont(2,1)-feed_rot_cont(1,1)))*0.3d0 
                  sgrp = sqrt(db_sigma(1)**2+db_sigma(2)**2)*0.3d0
 
                  unw = 1
                  if( frnge_qual(2:2).eq.'9' .or.
     .                frnge_qual(2:2).eq.'8' ) unw = 0
                  unw = or( data_flag, unw)

*                 set the line unweight indicator (2 for bad/ 1 forgood)
                  lunw = 2
                  if( unw.eq.0 ) lunw = 1
 
*                 Write out the observation
                  call jd_to_ymdhms( epoch, date, sectag)
                  write(luo,120) date, sectag, dgrp, sgrp, source,
     .                unw, lunw, frnge_qual, observable(4)/1.d3,
     .                elev(1,1)*180.d0/pi, azimuth(1,1)*180.d0/pi,
     .                temp_C(1,1), feed_rot_cont(1,1)*0.3, 
     .                feed_rot_cont(2,1)*0.3
 
 120              format(i5,4i3,1x,f6.2, 2f12.1, i3, i8, i2, a2, f12.3,
     .                f8.2,1x,f8.2, 1x,f6.1, 2f7.2)
 
****              Now start accumulating statistics on these data.
 
                  call accum_ke_stat( dgrp, sgrp, unw, frnge_qual,
     .                elev(1,1)*180.d0/pi, observable(4)/1.d3,
     .                epoch)

                  call accum_soln( dgrp, sgrp, unw, elev(1,1), 
     .                             azimuth(1,1), epoch)
 
*                         ! Looping over data
              end do

*             Now get the local mean (use the overall statistics)
              if( overal_lsum(1).gt.1 ) then
                  exper_mean(ne) = overal_lsum(3)/overal_lsum(2)
                  exper_chi(ne) = (overal_lsum(4) -
     .                    overal_lsum(3)*exper_mean(ne))/
     .                    (overal_lsum(1)-1)
                  exper_sig(ne) = sqrt(exper_chi(ne)/
     .                     overal_lsum(2))
                  exper_epoch(ne) = mid_epoch
                  exper_num(ne) = overal_lsum(1)

****              Now remove mean from other summations and save
*                 total accumualtion
                  call rem_mean( exper_mean(ne))

*                 Get the station position solution
                  call solve_pos
              end if

              ne = ne + 1
              call close_KalObs(ierr)
*                         ! Open OK
          end if
          ir = ir + 1
*                         ! Looping over experiments
      end do

      ne = ne - 1
 
****  Now finishing the Statistics
 
      call fin_ke_stats
      call out_ke_stats
 
****  Thats all
      end
 
CTITLE INIT_SUMS
 
      subroutine init_sums
 
*     Routine to clear all of the sumation variables
 
      include 'ke_anal.h'
 
*   max_elev        - Maxumium number of elevation bins
 
 
 
      integer*4 max_elev
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
*     Loop over all of the bins
 
      do i = 1, 5
          overal_sums(i) = 0.d0
          do j = 1, max_bins
              sig_sums(i,j)  = 0.d0
              elev_sums(i,j) = 0.d0
              rate_sums(i,j) = 0.d0
              frng_sums(i,j) = 0.d0
              time_sums(i,j) = 0.d0
          end do
      end do
 
****  Thats al
      return
      end
 
CTITLE INIT_LSUM
 
      subroutine init_lsum
 
*     Routine to clear all of the local experiment sumation variables
 
      include 'ke_anal.h'
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
*     Loop over all of the bins
 
      do i = 1, 5
          overal_lsum(i) = 0.d0
          do j = 1, max_bins
              sig_lsum(i,j)  = 0.d0
              elev_lsum(i,j) = 0.d0
              rate_lsum(i,j) = 0.d0
              frng_lsum(i,j) = 0.d0
              time_lsum(i,j) = 0.d0
          end do
      end do

*     Clear the normal equations
      do i = 1,4
         bvec(i) = 0.d0
         do j = 1,4
            norm_eq(i,j) = 0.d0
         end do
      end do
 
****  Thats al
      return
      end
 
      block data KE_anal_bd
 
*     Initialization of the bins of KE_anal
 
      include 'ke_anal.h'
 
      data elev_bins / 0.d0,  5.d0,  7.5d0, 10.d0, 15.d0, 20.d0,
     .                30.d0, 40.d0, 50.d0,  60.d0, 70.d0, 80.d0,
     .                90.d0 /
 
      end
 
CTITLE ACCUM_KE_STAT
 
      subroutine accum_ke_stat( dgrp, sgrp, unw, fq, el, rate, epoch )
 
*     This routine accumulates the statistics for this observation
 
      include 'ke_anal.h'
 
*   unw     - Edit flag, 0-good, <>0-bad
 
      integer*4 unw
 
*   dgrp, sgrp  - group delay difference and sigma (mm)
*   el          - Elvation (degs)
*   rate        - Rate (us/s)
*   epoch       - JD of obervation
 
      real*8 dgrp, sgrp, el, rate, epoch
 
*   fq          - Fringe quality code
 
 
      character*2 fq
 
* LOCAL VARIABLES
 
*   ib          - Bin number
*   decimaltoint    - Ascii code to number value
*   ierr        - IOSTAT error on decimaltoint
 
      integer*4 ib, decimaltoint, ierr
 
*   hour        - UTC Hour of day
 
      real*8 hour
 
****  Start with overal statistics
      if( unw.eq.0 ) then
          call sum_ke( dgrp, sgrp, 1.d0, overal_lsum )
 
****      Now get the elevation bin
          ib = 1
          do while ( ib.lt. max_elev .and.
     .               el.gt. elev_bins(ib) )
              ib = ib + 1
          end do
 
*         IB now points to the top of the range for the bins
          call sum_ke( dgrp, sgrp, el, elev_lsum(1,ib-1))
 
****      Sum into time bins
          hour = mod(epoch + 0.5d0,1.d0)*24.d0
          ib = hour + 1
          call sum_ke( dgrp, sgrp, hour, time_lsum(1,ib))
 
****      Sum into rate bins of 25 ps/sec
          ib = abs(rate/25.d0) + 1
          call sum_ke( dgrp, sgrp, abs(rate), rate_lsum(1,ib))
 
****      Sum into sigma bins ( 0.5 mm wide)
          if( sgrp.lt.5 ) then
              ib = sgrp/0.5d0  + 1
          else
              ib = sgrp/3.d0 + 11
          end if
          if( ib.gt.30 ) ib = 30
          call sum_ke( dgrp, sgrp, sgrp, sig_lsum(1,ib))
 
      end if
 
      if( and(unw,-2).eq.0 ) then
 
****      Now get the frnge quality bin
          if( fq(2:2).ge.'0' .and. fq(2:2).le.'9') then
              ib = decimaltoint(fq, ierr) + 1
              if( ierr.ne.0 ) ib = 1
              call sum_ke( dgrp, sgrp, 1.d0*(ib-1), frng_lsum(1,ib))
          end if
      end if
 
****  Thats all
      return
      end
 
CTITLE SUM_KE
 
      subroutine sum_ke( dgrp, sgrp, bin, sums)
 
*     Routine to actually accumulate values
 
*   dgrp, sgrp  - group delay difference and sigma (mm)
*   bin         - Binning information
*   sums(5)     - Summatation bins
 
      real*8 dgrp, sgrp, bin, sums(5)
 
*   wgh         - Weight in sum
 
      real*8 wgh
 
      if( sgrp.gt.0 ) then
          wgh = 1.d0 / sgrp**2
 
          sums(1) = sums(1) + 1
          sums(2) = sums(2) + wgh
          sums(3) = sums(3) + dgrp*wgh
          sums(4) = sums(4) + dgrp**2*wgh
          sums(5) = sums(5) + bin*wgh
 
      else
          write(*,'('' ** ERROR ** Sigma <= 0, Value '',d13.4)') sgrp
      end if
 
****  Thats all
      return
      end
 
CTITLE REM_MEAN
 
      subroutine rem_mean (mean)
 
*     This routine will remove the means from the local sums and
*     sum them into the total bins
 
      include 'ke_anal.h'
 
*   mean        - Mean for this experiment
 
      real*8 mean
 
*   ib          - Bin number
 
      integer*4 ib
 
****  Loop over all of the summations
      call rm( overal_lsum, overal_sums, mean )
 
      do ib = 1, max_elev
          call rm( elev_lsum(1,ib), elev_sums(1,ib), mean )
      end do
 
      do ib = 1, max_bins
          call rm( time_lsum(1,ib), time_sums(1,ib), mean )
          call rm( frng_lsum(1,ib), frng_sums(1,ib), mean )
          call rm( rate_lsum(1,ib), rate_sums(1,ib), mean )
          call rm( sig_lsum(1,ib),  sig_sums(1,ib),  mean )
      end do
 
****  Thats all
      return
      end
 
CTITLE RM
 
      subroutine rm( lsum, sums, mean)
 
*     This routine will remove the mean for local sum and acumulate
*     demeaned total sums
 
*   lsum(5), sums(5)    - Local and total summation quantities
*   mean            - Mean to be removed
 
 
      real*8 lsum(5), sums(5), mean
 
*     See if we need anything
      if( lsum(1).gt.0 ) then
 
*         Adjust mean
          sums(1) = sums(1) + lsum(1)
          sums(2) = sums(2) + lsum(2)
          sums(3) = sums(3) + lsum(3) - mean*lsum(2)
          sums(4) = sums(4) + lsum(4) - 2*mean*(lsum(3)-mean*lsum(2))
     .                    - mean**2* lsum(2)
          sums(5) = sums(5) + lsum(5)
 
      end if
 
****  Thats all
      return
      end
 
CTITLE FIN_KE_STATS
 
      subroutine fin_ke_stats
 
*     This routine will finish up the statistic calculation for
*     the total set of experiments analysed.
 
      include 'ke_anal.h'
 
*   ib      - Loop counter
 
      integer*4 ib
 
***** Loop over all the types
 
      call fin( overal_sums, overal_fin )
 
      do ib = 1, max_elev
          call fin(  elev_sums(1,ib), elev_fin(1,ib))
      end do
 
      do ib = 1, max_bins
          call fin( time_sums(1,ib), time_fin(1,ib))
          call fin( frng_sums(1,ib), frng_fin(1,ib))
          call fin( rate_sums(1,ib), rate_fin(1,ib))
          call fin( sig_sums(1,ib),  sig_fin(1,ib) )
      end do
 
****  Thats all
      return
      end
 
CTITLE FIN
 
      subroutine fin( sums, finals)
 
*     This routine will compute the final statsics for each of the
*     distributions
 
*   sums(5)     - The sumation values
*   finals(6)   - Final statistics
 
      real*8 sums(5), finals(6)
 
***** See if we need do anything
*                             ! Compute vlaues
      if( sums(1).gt.0 ) then
 
*****     Compute weighted mean
          finals(1) = sums(3)/ sums(2)
 
*         chi**2 (about 0)
          finals(4) = sums(4) / sums(1)
 
*         Sigma of mean
          finals(2) = sqrt(finals(4)/sums(2))
 
*         RMS about 0
          finals(3) = finals(2)*sqrt(sums(1))
 
*         Number of values
          finals(5) = sums(1)
 
*         Average bin value
          finals(6) = sums(5) / sums(2)
 
      end if
 
****  Thats all
      return
      end
 
CTITLE OUT_KE_STATS
 
      subroutine out_ke_stats
 
****  This routine will output the final results from the analysis
 
      include 'ke_anal.h'
 
* LOCAL VARIABLES
 
*   ib      - Bin counter
*   i       - Loop counter
*   trimlen - Length of string
*   date(5) - Date of experimet
 
      integer*4 ib, i,j, trimlen, date(5)
 
*   sectag  - Seconds tag
 
      real*8 sectag
 
*     Write out information about each experiment
      write(*,100)
 100  format('* Summary of statistics from following data',/,
     .        '* Experiment Date    Mean (mm)  +-     Chi**2 ',
     .        ' #   Name')
 
      do i = 1, ne
          call jd_to_ymdhms( exper_epoch(i), date, sectag )
          write(*,120) date, exper_mean(i), exper_sig(i),
     .        exper_chi(i), exper_num(i),
     .        ko_names(i)(1:trimlen(ko_names(i)))
 120      format(i5,4i3, f12.2, f10.2, f12.4, 1x,i4, 1x, a)
      end do
 
***** Now put out the statistics
 
      call out('Overall', 1, overal_fin )
 
      do ib = 1, max_elev
          call out( 'Elevation angle', ib, elev_fin(1,ib))
      end do
 
      do ib = 1, max_bins
          call out( 'Time (UTC hrs)', ib, time_fin(1,ib))
      end do
      do ib = 1, max_bins
          call out( 'Fringe quality', ib, frng_fin(1,ib))
      end do
      do ib = 1, max_bins
          call out( 'Rate value (ps/s)', ib, rate_fin(1,ib))
      end do
      do ib = 1, max_bins
          call out( 'Sigma (mm)', ib, sig_fin(1,ib) )
      end do

****  Now output the position estimates

      write(*,200) 
 200  format(/'* Summary of Position estimates from following data',/,
     .        '* Experiment Date    Offset   +-',7x,'N (mm) +-',
     .        7x,'E (mm)    +-',7x,'U (mm)    +-',1x,'#   Name')

      do i = 1, ne
          call jd_to_ymdhms( exper_epoch(i), date, sectag )
          write(*,220) date, (exper_soln(j,i), exper_snsg(j,i),
     .        j=1,4), exper_num(i),
     .        ko_names(i)(1:trimlen(ko_names(i)))
 220      format(i5,4i3, 4(f8.2,1x, f8.2,1x),i4, 1x, a)
      end do
 
****  Thats all
      return
      end
 
 
CTITLE OUT
 
      subroutine out( label, ib, finals )
 
*     Routine to output the final answers
 
*   ib      - Bin being output (used to see if header
*           - should be writen)
 
      integer*4 ib
 
*   finals(6)   - Final statsics
 
      real*8 finals(6)
 
 
      character*(*) label
 
***** See if we should write header
 
      if( ib.eq.1 ) then
          write(*,100) label
 100      format(/,'* Final statistics for ',a,' binning',/,
     .        '*  Mean (mm)     +-           WRMS (mm)   chi**2/n',
     .        '   #    Average Bin')
      end if
 
*                                 ! Output
      if( finals(5).ne.0 ) then
          write(*,200) finals
  200     format(f12.2, 1x, f10.2,1x, f12.2, 1x, f12.4,1x, f6.0,1x,
     .          f12.4)
 
      end if
 
***** Thats all
      return
      end
 
CTITLE ACCUM_SOLN    
 
      subroutine accum_soln( dgrp, sgrp, unw, el, az, epoch )
 
*     This routine accumulates the statistics for this observation
 
      include 'ke_anal.h'
 
*   unw     - Edit flag, 0-good, <>0-bad
 
      integer*4 unw
 
*   dgrp, sgrp  - group delay difference and sigma (mm)
*   el          - Elvation (rads)
*   az          - Azimuth (rads)
*   epoch       - JD of obervation
 
      real*8 dgrp, sgrp, epoch
      real*4 el, az 
 
* LOCAL VARIABLES

* i,j           - Loop counters
* apart(4)      - Offset, NEU partials
* wgh           - Weight for this observation

      integer*4 i,j
      real*8 apart(4), wgh

****  Increment the normal equations if data OK

      if( unw.eq.0 ) then
          apart(1) = 1.d0
          apart(2) = -cos(el)*cos(az)
          apart(3) = -cos(el)*sin(az)
          apart(4) = -sin(el)

          wgh = 1.d0/sgrp**2

*         Increment normal equations
          do i = 1,4
             bvec(i) = bvec(i) + apart(i)*dgrp*wgh
             do j = 1,4
                norm_eq(i,j) = norm_eq(i,j) + apart(i)*apart(j)*wgh
             end do
          end do
      end if

***** Thats all
      return
      end
          
CTITLE SOLVE_POS

      subroutine solve_pos 

*     This routine will solve the normal equations for the 
*     group-phase station position.

      include 'ke_anal.h'

* LOCAL VARIABLES

* scale    -- Used to scale the matrux

      real*8 scale(4)

* ipivot(4) - Pivot elements
* i,j       - Loop counters

      integer*4 ipivot(4), i, j

***** Solve the system

      call invert_vis( norm_eq, bvec, scale, ipivot, 4,4,1)

*     Now save the results

      do i = 1,4
         exper_soln(i,ne) = bvec(i)
         exper_snsg(i,ne)  = sqrt(norm_eq(i,i))
      end do

***** Thats all
      return
      end


