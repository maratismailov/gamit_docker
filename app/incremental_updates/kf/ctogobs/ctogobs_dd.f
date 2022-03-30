CTITLE CLEAN_DD
 
      subroutine clean_dd(pass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )

      implicit none
 
*     This routine will clean the phase data in double difference mode
*     The procedure used is loop over all one-way data and when a gap
*     or bias flag is found, a double differece combination is formed
*     that will allow the cycle slip to be fixed.  If the patch is
*     reliable then the bias flag is removed.
*
*     The order in which the data is cleaned is based on the range
*     data rms for sites and clock allan standard deviation for
*     satellites.  As the data is clean these list are used for the
*     first guess at finding the double difference combination.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   pass            - Indicates which pass this is (for iteration of 
*                     cleaning.  On second pass gaps are ignored
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), pass
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*   i,j,k   - Counters of stations, satellites and epochs
*   ns, lv  - Current station and satellite list number
*   js, kv  - Double difference site and SV if used.
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
*   first_ep    - Epoch number of first data point after the first
*           - bias flag
*   step_ep - Number of epochs between each data point for this
*           - site.
*   ct      - Temporary channel number (used for flagging)
*   trimlen - Length of string
*   ierr    - IOSTAT error on opening the dd report file
*   dd_clean_iter - Number of iteratations through the dd cleaning
*             for a single station/satellite.  Each time data is
*             flagged or a bias flag added we iterate the clean loop.
*   latest_dd_edit - Epoch number of the highest dd edit.  Used to see
*             if we should push a bias flag to this point.
 
      integer*4 i,j,k, l, ns, lv, ch, ltoc, first_ep, step_ep, ct,
     .          js, kv, trimlen, ierr, dd_clean_iter, latest_dd_edit
 
*   data_OK - function returns true if data OK.  phs_mask does
*           - no check bias flags
*   kbit    - Checks if bit is set.
*   first_bias  - Indicates that first bias flag has been
*           - found
*   in_gap  - Indicates that we are in gap
*   bias_last - Indicates that we added a bias flag in check_continuity
*   checked   - Indicates to check_continuity that these data have not
*               been checked (not really used here)
 
 
      logical data_OK, kbit, first_bias, in_gap, bias_last, checked,
     .        good_bf

***** See if we should open the double difference report file.
      if( trimlen(dd_outfile).gt.0 ) then
          open(202, file=dd_outfile, status='unknown', iostat=ierr)
          call report_error('IOSTAT',ierr,'open',dd_outfile,0,
     .                      'clean_dd')
      end if

***** Check the status of flagging gaps.  If we the user has said
*     to flag gaps then set the ingore_gaps true, otherwize do what
*     use said
      if( gaps_flagged .or. pass.gt.1 ) then
          ignore_gaps = .true.
      else
          ignore_gaps = usr_ignore_gaps
      end if

****  Set up the scaling on slip detector so with each pass it
*     gets larger and large.

      tol_scale = (1.0 + (pass-1)/4.0)
 
***** First order the stations and satellites for cleaning
 
      call get_clean_order('CLEANING')
 
*     Now loop over the one-way data
      do i = 1, num_cfiles
          ns = dd_site_list(i)
          dd_clean_iter = 1

*         Loop over the satellites as a do-while loop
*         because we may need to repeat a satellite if there
*         is edited data or bias flags added.
          j = 0
          do while (j.lt.num_sat)
              j = j + 1
              lv = dd_svs_list(j)
              first_bias = .false.
              in_gap     = .false.
 
              curr_L1_slip = 0.d0
              curr_L2_slip = 0.d0
              dd_bias_added = .false.

****          Set up the default search lists
              call create_dd_search( ns, lv )
 
****          Scan up to the first bias flag (always marked) and
*             step at the data sampling interval for this site
              k = 0
              do while( k.lt.num_ep .and. .not.first_bias )
                  k = k + 1
                  ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! SV observed at this time
                  if( ch.gt.0 ) then
 
c                      if( kbit(data_flag_cse(ch,ns,k),31).or.
c     .                    kbit(data_flag_cse(ch,ns,k),32) ) then
                       if( good_bf(data_flag_cse(ch,ns,k),0, 
     .                             phs_mask) ) then
                          first_bias = .true.
                      end if
                  end if
              end do
 
*             Based on the sampling interval at this site compute
*             number of epochs bewteen data points (this is so that
*             we know when there are gaps).
              if( orig_sampling(ns).gt.sampling ) then
                  step_ep = orig_sampling(ns)/sampling
              else
                  step_ep = 1
              end if
 
*             Start at the first point after the first bias flag
C             first_ep = k + step_ep
*             start at the first point so that form_data will 
*             check the quality even though it will not be able
*             to patch the data
              first_ep = k
 
              write(*,120) cf_codes(ns), prn_list(lv), first_ep,
     .                    step_ep, dd_clean_iter
  120         format(' Processing ',a4,' PRN ',i2.2,' Start ep ',i4,
     .                ' Step ',i3,' epochs: Iteration ',i2)
 
              do k = first_ep, num_ep, step_ep
 
*                 see if we have the end of a gap or a bias flag
                  ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! SV observed at this time
                  if( ch.gt.0 ) then
 
*                     See if we have found a good data point after a
*                     gap or with a bias flag.
                      if( data_OK(data_flag_cse(ch,ns,k),0,phs_mask)
     .                    .and. ( in_gap .or.
     .                       kbit(data_flag_cse(ch,ns,k),31).or.
     .                       kbit(data_flag_cse(ch,ns,k),32))    ) then
 
*                         Set the number of unflagged cycle slips to
*                         zero
                          unflg_num = 0

*                         Clear the double difference second site and
*                         satellite so thatwe will know double diffs were
*                         used
                          js = 0
                          kv = 0
                          call scan_and_fix(ns, lv, js, kv, k, step_ep,
     .                        L1r_phs_cse, L2r_phs_cse,
     .                        L1r_rng_cse, L2r_rng_cse,
     .                        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                        data_flag_cse, bf_type_cse, params_cse,
     .                        par_flag_cse, pass )

*                         See if we edited any data because there was
*                         not enough of it
                          do l = 1, dd_edit_num

*                             Set the iteration needed flag
                              dd_bias_added = .true.
                              write(*,210) dd_edit_ep(l), cf_codes(ns),
     .                                prn_list(lv), k
 210                          format(' FLAGGING EP ',i4,' Site ',A4,
     .                               ' PRN ',i2.2,' Too little data,',
     .                               ' Patching epoch ',i5)
                              ct = ltoc(ctol_cse(1,ns,dd_edit_ep(l)),
     .                                  lv, actual_max_chan)
                              if( ct.gt.0 ) then

*                                 MOD TAH 980806: If this is a good
*                                 data point that we are edting and it
*                                 has a bias flag, make sure that we 
*                                 push it to the next good point
                                  call check_pb( data_flag_cse,
     .                                 ctol_cse, ns, lv,                                
     .                                 dd_edit_ep(l) )                                                           
                                  call sbit(data_flag_cse(ct,ns,
     .                                      dd_edit_ep(l)), 6,1)
                              else
                                   write(*,*) 'ERROR ct = ',ct
                              end if
                          end do

*                         Flag a cycle slip at this epoch in case it was
*                         a gap thatwe originally found
                          if( dd_edit_num.gt.0 ) then
                              call sbit(data_flag_cse(ch,ns,k),31,1)
*                             If this is a forward edit (i.e., the edited epochs
*                             are after the epoch being processed) then push the
*                             bias flag forard to the next good point.  Look
*                             forward to latest epoch editted
                              latest_dd_edit = dd_edit_ep(1)
                              do l = 2, dd_edit_num
                                  if( dd_edit_ep(l).gt.latest_dd_edit)
     .                                latest_dd_edit = dd_edit_ep(l)
                              end do
                              if( latest_dd_edit.ge.k ) then
                                  call push_ddbf(latest_dd_edit, 
     .                                  ns, lv, data_flag_cse,
     .                                  bf_type_cse,ctol_cse)
                              end if
                          end if

*                         Check to see if there were unflagged cycle
*                         slips.
                          do l = 1,unflg_num

*                             Set the iteration needed flag
                              ct = ltoc(ctol_cse(1,ns,unflg_ep(1,l)),
     .                                  lv, actual_max_chan)

*                             Resolve where the slip is:
                              call set_scan_list(unflg_ep(1,l) )
                              bias_last = .false.
                              if( num_scan.gt.0 .and. 
     .                            unflg_type(l).eq.'LCDD' ) then
                                  dd_bias_added = .true.
                                  checked = .false.
                                  call find_dd_slip(ns,lv, 
     .                                unflg_sv(1,l),unflg_sv(2,l),
     .                                checked, bias_last, L1r_phs_cse,
     .                                L2r_phs_cse, L1_cyc_cse, 
     .                                L2_cyc_cse,  ctol_cse,
     .                                data_flag_cse, bf_type_cse)
                                   bias_last = .true.
                              end if

*                             If we did not flag anthing, mark the
*                             current one-way.  Will happen with WL 
*                             on a long search and for points flagged
*                             imediately after a bias flag.  If this
*                             second pass and the slip is WL then
*                             don't flag.

                              if( ct.gt.0 .and. .not.bias_last .and.
     .                            (unflg_type(l).eq.'WLOW'.and.
     .                             pass.eq.1 ) ) then
                                  dd_bias_added = .true.
* MOD TAH 200617: Updated ep to I5 format
                                  write(*,220) unflg_ep(1,l), 
     .                                     cf_codes(ns),prn_list(lv),
     .                                     unflg_mag(l), unflg_type(l)
  220                             format(' UNFLAGGED slip at epoch ',I5,
     .                               ' Site ',a4,' PRN ',i2.2,
     .                               ': Magnitude ',F10.2,' cyc,',
     .                               ' type ',a4)
                                  call sbit(data_flag_cse(ct,ns,
     .                                    unflg_ep(1,l)),31,1)

*                                 Set the type of bias flag
                                  if( unflg_type(l).eq.'WLOW' ) 
     .                            call sbit(bf_type_cse(ct,ns,
     .                                      unflg_ep(1,l)),6,1)
                                  if( unflg_type(l).eq.'LCDD' ) 
     .                            call sbit(bf_type_cse(ct,ns,
     .                                      unflg_ep(1,l)),7,1)
                              else if( ct.le.0 ) then
                                  write(*,*) 'ERROR ct = ',ct
                              end if

                          end do

 
                          in_gap = .false.

                      end if
 
*                     See if this ia a bad data point.  Treat as gap
                      if(.not.
     .                    data_OK(data_flag_cse(ch,ns,k),
     .                            0,phs_mask) ) then
                          in_gap = .true.

*                         if we are ignoring gaps then set in_gap to .false. 
*                         even though we are in a gap
                          if( ignore_gaps ) in_gap = .false.
                      end if
 
*                     Apply the current number of cycle slips to
*                     this point
 
                      L1_cyc_cse(ch,ns,k) = L1_cyc_cse(ch,ns,k) +
     .                                    curr_L1_slip
                      L2_cyc_cse(ch,ns,k) = L2_cyc_cse(ch,ns,k) +
     .                                    curr_L2_slip
 
                  else  
 
****                  No data, so say we in a gap
                      in_gap = .true.

*                     if we are ignoring gaps then set in_gap to .false. 
*                     even though we are in a gap
                      if( ignore_gaps ) in_gap = .false.
                  end if
*                         ! Looping over epochs
              end do

****          See if we need to iterate the double difference cleaning
*             loop.
              if( dd_bias_added .and. dd_clean_iter.lt.5 ) then
                  dd_clean_iter = dd_clean_iter + 1
                  j = j - 1
              else
                  dd_clean_iter = 1
              end if
*                         ! Looping over satellites
          end do
*                         ! Looping over stations.
      end do
 
***** Thats all
      return
      end
 
CTITLE GET_CLEAN_ORDER
 
      subroutine get_clean_order(type)

      implicit none
 
*     This routine will make the list of the order in which to
*     clean the double differences.  For stations the order is based
*     on the rms scatter of the range residuals, for the satellites on
*     rms of the clock errors.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES 

* type  - Reason for calling routine (for output only)

      character*(*) type
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
*   rms(max_cfiles+max_gprn)   - RMS quantity used for
*                   - sorting (rearranged in sort
*                   - and therefore a copy is used
*  WARNING: Dimension here should be the larger of max_gprn
*           and max_cfiles. 
*  MOD TAH 950818: Problem avioded by using sum.
 
      real*8 rms(max_cfiles+max_gprn)
 
***** First sort the sites by the rms of the delay residuals.
 
*     Make a copy of the rms list since this will be changed in the
*     sort routine
 
      do i = 1, num_cfiles
          rms(i) = rng_noise(i)
      end do
 
*     Now get the sorted list
      call sort_dd_list( rms, num_cfiles, dd_site_list)
 
      write(*,120) type, (i, cf_codes(dd_site_list(i)), i=1,num_cfiles)
  120 format(/,' DOUBLE DIFFERENCE ',a,/,
     .        ' Cleaning data with following site order:',/,
     .        30(5(I3,1x,a4,3x):,/))
 
****  Now do the same for the satellites
      do i = 1, num_sat
          rms(i) = apr_clk_var(num_cfiles+i)
      end do
      call sort_dd_list( rms, num_sat, dd_svs_list)
      write(*,140) (i, prn_list(dd_svs_list(i)), i=1,num_sat)
  140 format(' Cleaning data with following satellite order:',/,
     .        30(5(I3,' PRN ',i2.2,1x):,/))
 
****  Thats all
      return
      end
 
CTITLE SORT_DD_LIST
 
      subroutine sort_dd_list( rms, num, dd_list )
 
*     Routine to return the list of sites or satellites sorted in
*     decending order by the rms quantity given.
 
* PASSED VARIABLES
 
*   num     - Number of values to be sorted
*   dd_list(num)  - list of values sorted in ascending rms order.
 
 
      integer*4 num, dd_list(num)
 
*  rms(num)  - Rms of values (used for sorting).  The order of this list
*     is changed during the sort.
 
 
      real*8 rms(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - Smallest integer in current pass.
*   swap   - Value used to swap real*8
*   iswap  - the index to be swapped.
 
      integer*4 i,j, smallest_one , iswap
 
 
      real*8 swap
 
****  set up list for indexing
      do i = 1, num
         dd_list(i) = i
      end do
 
****  Start loop using exchange sort
 
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( rms(j).lt. rms(smallest_one) ) then
                  smallest_one = j
              end if
          end do
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              swap = rms(smallest_one)
              iswap = dd_list(smallest_one)
              rms(smallest_one) = rms(i)
              dd_list(smallest_one) = dd_list(i)
              rms(i) = swap
              dd_list(i) = iswap
          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end
 
CTITLE SCAN_AND_FIX
 
      subroutine scan_and_fix(ns, lv, js, kv, ep, step_ep,
     .        L1r_phs_cse, L2r_phs_cse, L1r_rng_cse, L2r_rng_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse, data_flag_cse,
     .        bf_type_cse, params_cse, par_flag_cse, pass )

      implicit none
 
*     This routine will scan about gap or a bias flag and fix the
*     cycle slip and remove the bias flag if it is deeemed OK.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number
*   js, kv   - Double difference pair if used. (set in the routine below)
*   ep       - Epoch at which the gap or bias flag occurrs
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   pass - Number of iteratations through the dd cleaning
*             for a single station/satellite.  Each time data is
*             flagged or a bias flag added we iterate the clean loop.
 
 
      integer*4 ns, lv, ep, step_ep, js, kv,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), pass
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*       ch      - Channel number for current PRN
*     ltoc      - Function to get channel number of a given
*                 satellite list number
*     i         - Loop countr
*     gbi       - Indicates if Gap or Bias that is being
*                 patched.
      integer*4 ch, ltoc, i, gbi

*   reliable - Logical to say if patch is reliable
*   force_one_way - Logical to say that their is a slip
*              on all channels at one site and that we
*              should force the one-way to have bias
*              flag removed.
*   kbit     - Check status of bit

      logical reliable, force_one_way, kbit

 
*   norm_cyc(2,2)   - Normal equations for L1 and L2 cycle
*               - slip estimate
*   b_cyc(2)        - B-vector for estimation of number of cycle
*               - slips 
*   dL1_slip, dL2_slip - number of additional cycles needed
*                 accross a patch
*   chiqual     - Quantitative estimate of quality of cycle
*                 slip fixing.
*   norm_trl(2,2), b_trl(2) - Trial values for the nomral equations
*                 when we want to try one-way ionosphere fix.
*   norm_onebg(2,2), b_onebg(2) - Trial values when we want to see if
*                 allowing one bias or gap will allow us to fix cycle
*                 slips.
*   chi_cyc, chi_trl, chi_onebg  - Prefit chi**2 for trial and actual 
*                 estimate
 
      real*8 norm_cyc(2,2), b_cyc(2), dL1_slip, dL2_slip,
     .       chiqual, norm_trl(2,2), b_trl(2), 
     .       norm_onebg(2,2), b_onebg(2), chi_cyc, chi_trl, chi_onebg

*    gap_or_bias(2) - String to denoted if a bias or gap being
*                     patched.  (Used for output only)
*    small_one_bg   - Logical that is true if the gap is small
*                     enough for allow one bias/gap to be used.

      character*4 gap_or_bias(2)
      logical     small_one_bg

      data gap_or_bias / 'GAP ', 'BFLG' /
      
 
***** Start: Initialize the estimator for the number of cycle slips
 
      call init_cyc_est( norm_cyc, b_cyc, chi_cyc )
      reliable = .false.
      force_one_way = .false.
      dd_edit_num = 0
      dL1_slip = 0
      dL2_slip = 0

 
*     collect the one-way for this site.  If we have L2 ranges first
*     collect a lot and so that we can widelane, narrow lane and LG
*     patch.
 
      ch = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)

*     If we have an L2 range measurement then try to use the
*     wide-lane and narrow lane observables.
      if( lambda(lv,4,ns).eq.1 ) then
 
*         Collect the widelane data.  The data are returned through
*         common in data_lft and data_rgh
          call get_one_way(ns, lv, ep, step_ep, max_wl_ret,
     .            L1r_phs_cse, L2r_phs_cse,
     .            L1r_rng_cse, L2r_rng_cse,
     .            L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .            data_flag_cse, params_cse)
 
*         Add this contribution to the estimate of the number of
*         cycles in slip. Copy the one-way data to the work arrays
          call copy_dtow
          call inc_cyc_est(ns, lv, 0, 0, 'WL', 'OW', norm_cyc, b_cyc,
     .                     chi_cyc )
      end if

****  Now collect the double difference data with this site and
*     satellite.  Again the data are returned in the data_lft and
*     data_rgh arrays.
      call get_one_way(ns, lv, ep, step_ep, 2*max_dd_ret,
     .        L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, params_cse)

****  See how much data is availiable.  If it is too small then edit it
*     out.
      if( num_rgh.lt. min_good_bias .or. num_lft.lt.min_good_bias ) then

*         not enough data so save the epochs and flag data
          if( num_rgh.lt.min_good_bias ) then
              do i = 1, num_rgh
                 dd_edit_ep(i) = ep_rgh(i)
              end do
              dd_edit_num = num_rgh
          end if
          if( num_lft.lt.min_good_bias ) then
              do i = 1, num_lft
                 dd_edit_ep(i+dd_edit_num) = ep_lft(i)
              end do
              dd_edit_num = dd_edit_num + num_lft
          end if

* MOD TAH 970110: Only exit the loop if there is finite amount of data.
*         If zero in left segment, we are at start and want to use
*         the absolute fix (most appropriate when postfit residuals
*         are used).
          if( num_lft.ne.0 ) RETURN
          RETURN
      end if

*     If the gap in the data is less than tolerance for one-way fixing
*     try to fix now.  If it is deemed reliable then don't bother with
*     double difference fixing.  Only try this for data with L2 range.

* MOD TAH 970110: added check on the number of data being greater than
*     zero in both segments.
 
      if( (ep_rgh(1)-ep_lft(1))*sampling.le. tol_one_way_fix .and.
     .     lambda(lv,4,ns).eq.1 .and. 
     .     num_lft.gt.0 .and.num_rgh.gt.0 ) then
 
*         Copy one-ways to work arrays and Increment the estimator
          call copy_dtow

*         Copy the normal equations for the test.  This is so that
*         we will not add the ionopshere twice if have to form
*         double differences
          call copy_norm(norm_cyc, b_cyc, chi_cyc, 
     .                   norm_trl, b_trl, chi_trl)
          call inc_cyc_est(ns, lv, 0, 0,'LG', 'OW', norm_trl, b_trl,
     .                   chi_trl )
 
*         solve for the number of cycles in L1 and L2.  Since this is
*         a one-way trial pass 0 as the second station and satellite.
          call est_dd_cyc( ep, ns, lv, 0, 0, norm_trl, b_trl, chi_trl,
     .                reliable, chiqual, dL1_slip, dL2_slip )
      end if
 
*     Go to double difference fixing if we do not have a relaible
*     fix yet.  Only try the double difference cleaning if we
*     have more than one cfiles and the gap in the one ways is
*     less than twice the duration over which we will allow the
*     patch.
      if( .not. reliable ) then

*         Initially do not allow any bias flags or gaps in second
*         satellite in single difference.  If we don't find any
*         data; then we will allow one bias or gap.
          js = 0

*         Clear the normal equations so that double differences
*         will be done only with double differnece
* MOD TAH 970110: Added num_lft check since we can get here with
*         no data in the left segment.
          call init_cyc_est( norm_cyc, b_cyc, chi_cyc )
          if( num_cfiles.gt.1 .and. num_lft.gt.0 .and.
     .        (ep_rgh(1)-ep_lft(1))*sampling.le.2*dchi2_max_sep ) then
              allow_one_bg = .false.
              call copy_dtow
              call get_dd_data( ns, lv, js, kv, ep, step_ep, 
     .                    norm_cyc, b_cyc, chi_cyc,
     .                    L1r_phs_cse, L2r_phs_cse,
     .                    L1r_rng_cse, L2r_rng_cse,
     .                    L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                    data_flag_cse, reliable )

*             If no doble differnces; try with one bias or gap allowed
*             if the user has said it is OK.
              if( js.eq.0 .and. kv.eq.0 .and. do_one_bg ) then
                  allow_one_bg = .true.
                  call copy_dtow
                  call copy_norm(norm_cyc, b_cyc, chi_cyc, 
     .                           norm_onebg, b_onebg, chi_onebg)
                  call get_dd_data( ns, lv, js, kv, ep, step_ep, 
     .                    norm_onebg, b_onebg, chi_onebg,
     .                    L1r_phs_cse, L2r_phs_cse,
     .                    L1r_rng_cse, L2r_rng_cse,
     .                    L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                    data_flag_cse, reliable )

*                 See if we could now reliably fix biases.  If we could
*                 then set force_one_way to true and set the combination
*                 back to zero so that we wil force LG in the one-ways.
*                 Check to see if we edted data which occurrs when there
*                 are some data between bias flags but not enough for
*                 patch (i.e. <5 data points)
* MOD TAH 940113: Put more conditions on allowing one bias or gap:
*                 New conditions are:
*                 (1) Not allowed on the first iteration;
*                 (2) Not allowed if the gap in the one-way sequence
*                     is more than 5% of maximum gap over which patching
*                     will be done.
*                 See if we match condition for small enough gap for
*                 one-bias-gap to be used.
* MOD TAH 940406: Set the small_one_bg condition so that it grows with
*                 passes through the iteration loop.
                  small_one_bg = .false.
                  if( (ep_rgh(1)-ep_lft(1))*sampling
     .                         .le.(pass-1)*0.05*dchi2_max_sep .or.
     .                 ep_rgh(1)-ep_lft(1).le.step_ep*gap_size(ns) 
     .                                                          ) then
* MOD TAH 990518: Allowed gap to encompass gap_size range.
C    .                 ep_rgh(1)-ep_lft(1).le.step_ep ) then
                       small_one_bg = .true.
                  end if

                  if( js.gt.0          .and. reliable .and. 
     .                dd_edit_num.eq.0 .and. 
     .                pass.gt. 1 .and. small_one_bg ) then

                      js = 0
                      force_one_way = .true.

                  end if
              end if
          end if

*         If still no double difference were found then js (the second
*         site number should return = 0, if this is the case then
*         apply the 1-one way ion-constraint
          if( js.eq.0 ) then
              call copy_dtow
              if( lambda(lv,4,ns).eq.1 ) then
*                 Just use LG and WL.  If no widelane then use LC as well
*                 MOD TEST: Changed fix to LC instead of ion (LG).  Next
*                 MOD used new LA feature.  Introduced WA option for
*                 absolute check on wide-lane.  Since special code now
*                 to process one-way postfits, go back to LGWL
                  if( force_one_way ) then
                     call inc_cyc_est(ns,lv, 0, 0,'LGWL','OW',norm_cyc, 
     .                                b_cyc, chi_cyc )
                  else
                     call inc_cyc_est(ns,lv, 0, 0,'LAWA','OW', 
     .                                norm_cyc, b_cyc, chi_cyc )                        
                  end if
              else
* MOD TAH 990517: Added check to see if we have L2 phase data.              
                  if( lambda(lv,2,ns).ne.0 ) then
                      call inc_cyc_est(ns, lv, 0, 0,'LCLG', 'OW', 
     .                                 norm_cyc, b_cyc, chi_cyc )
                  else
                      call inc_cyc_est(ns,lv, 0, 0,'L1','OW', 
     .                                 norm_cyc, b_cyc, chi_cyc )
                  endif
              end if
          end if

*         If good double differences were found then the normal
*         equations will have been updated so Just estimate
*         number of cycles 
*         Now solve for the number of cycles in L1 and L2

         call est_dd_cyc( ep, ns, lv, js, kv, norm_cyc,b_cyc, chi_cyc, 
     .                reliable, chiqual, dL1_slip, dL2_slip )

*         See if we are forcing one_way data
          if( force_one_way ) then
              reliable = .true.
          else
              if( js.eq.0 ) reliable = .false.
          end if

*         To be sure set the estimate unreliable if there was
*         no data in the left or write segments
          if( reliable ) then
              if( num_wl.lt.3 .or. num_wr.lt.3 ) then
                  write(*,310) num_wl, num_wr
 310              format(' ** WARNING ** Num WL and WR are ',2i3,
     .                   ' yet est_dd_cyc marked reliable ')
                  reliable = .false.
              end if
          end if
      end if

*     See if we just patched a bias or a gap:  This just for output to
*     let user know if there are slips across gaps in the data that
*     have not been flagged.
      if( kbit(data_flag_cse(ch,ns,ep),31).or.
     .    kbit(data_flag_cse(ch,ns,ep),32)   ) then
          gbi = 2
      else
          gbi = 1
      end if
 
*     Update the number of cycles to be applied to the data
      write(*,200) ep, cf_codes(ns), prn_list(lv), curr_L1_slip ,
     .            dL1_slip, curr_L2_slip , dL2_slip, reliable,
     .            chiqual, gap_or_bias(gbi), allow_one_bg, force_one_way
 200  format(' Epoch ',i4,' Site ',a4,' PRN ',i2.2, ' L1 from ',f8.2,
     .        ' to ',F8.2,' L2 from ',f8.2,' to ', f8.2,' Reliable? ',
     .        L2, F8.2,1x,a4,' OneBG ',L1,' Force ',L1)
 
      curr_L1_slip = dL1_slip
      curr_L2_slip = dL2_slip
 
*     If the estimate of the number of cycles is considered reliable
*     then remove the bias flag.  If it is not reliable then add a
*     flag (in case this is a gap)
 
      if( reliable ) then
          call sbit( data_flag_cse(ch,ns, ep), 31, 0 )
          call sbit( data_flag_cse(ch,ns, ep), 32, 0 )
      else
          call sbit( data_flag_cse(ch,ns, ep), 31, 1 )
*         If we are setting the bias flag here, record as a 
*         gap.
          if( .not. kbit(data_flag_cse(ch,ns, ep), 31) .and.
     .        .not. kbit(data_flag_cse(ch,ns, ep), 32) ) then
              call sbit( data_flag_cse(ch,ns, ep), 31, 1 )
              call sbit( bf_type_cse(ch,ns, ep), 4,1 )
          end if

      end if
 
****  Thats all
      return
      end
 
CTITLE INIT_CYC_EST
 
      subroutine init_cyc_est(  norm_cyc, b_cyc, chi_cyc )

      implicit none
 
*     Routine to initialize the normal equations and bvector for
*     estimating the number of cycle slips.
 
* PASSED VARIABLES
 
*   norm_cyc(2,2)   - Normal equations for L1 and L2 cycle
*               - slip estimate
*   b_cyc(2)        - B-vector for estimation of number of cycle
*               - slips
*   chi_cyc     - Prefit chi**2 value
 
      real*8 norm_cyc(2,2), b_cyc(2), chi_cyc
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
****  Clear the matrix and array
      do i = 1,2
          b_cyc(i) = 0.d0
          do j = 1,2
              norm_cyc(i,j) = 0.d0
          end do
      end do

      chi_cyc = 0.d0
 
****  Thats all
      return
      end
 
CTITLE COPY_DTOW
 
      subroutine copy_dtow
 
      implicit none

*     Routine to copy the one-way to the work arrays where it will
*     be used.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES - None
 
* LOCAL VARIABLES
 
*   i,j     - Loop parameters
 
      integer*4 i,j
 
****  Copy over the left and right segments
      num_wl = num_lft
      do i = 1, num_lft
          ep_wl(i) = ep_lft(i)
          do j = 1,4
              data_wl(j,i) = data_lft(j,i)
          end do
      end do
 
      num_wr = num_rgh
      do i = 1, num_rgh
          ep_wr(i) = ep_rgh(i)
          do j = 1,4
              data_wr(j,i) = data_rgh(j,i)
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE GET_ONE_WAY
 
      subroutine get_one_way(ns, lv, ep, step_ep, max_ret,
     .            L1r_phs_cse, L2r_phs_cse,
     .            L1r_rng_cse, L2r_rng_cse,
     .            L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .            data_flag_cse, params_cse)

      implicit none
 
*     This routine will scan forward and backwards in time about
*     epoch ep getting as much data as possible upto max_ret
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number
*   ep       - Epoch at which the gap or bias flag occurrs
*   step_ep  - Number of epochs bewteen data points at this site.
*   max_ret  - Maximum number of values to return.
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
      integer*4 ns, lv, ep, step_ep, max_ret,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep) - Clock estimates at the L1 frequency
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep), 
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*   k       - Counter for epochs
*   nl,nr   - Number of data in left and right segments
*   ltoc    - Function to return channel number given list
*           - number
*   ch      - Channel number (-1 if satellte not observed)
 
      integer*4 k, nl,nr, ltoc, ch
 
*   data_OK - Indicates that data is OK
*   kbit    - Checks the status of a bit in a word.
*   finished    - Indicates that search is finished (either
*           - beginning or end of data; enough data ; or bias
*           - flag found.
 
 
      logical data_OK, kbit, finished

*   dclk   - Contribution to one-ways from clock in reciever and
*            satellite (L1 cycles)

      real*8 dclk
	
 
***** Start by working backwards from current ep until we have reached
*     max_ret, or the beginning of the data, or a bias flag
 
      k = ep
      finished = .false.
      nl = 0
      do while ( .not. finished )
          k = k - step_ep
          if( k.le.0 ) finished = .true.
*                                         ! See if this point OK
          if( .not.finished ) then
 
*             Get the channel number at this epoch
              ch = ltoc( ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! We have data at this time
              if( ch.gt.0 ) then
 
*                 See if data is OK
                  if(data_OK(data_flag_cse(ch,ns,k),0, phs_mask)) then
 
*                     Add these data to the left segment
                      nl = nl + 1
                      dclk = params_cse(ns,k) - 
     .                       params_cse(num_cfiles+lv,k)
                      data_lft(1,nl) = L1r_phs_cse(ch,ns,k) +
     .                            L1_cyc_cse(ch,ns,k) - dclk
                      if( lambda(lv,2,ns).ne.0 ) then
                         data_lft(2,nl) = L2r_phs_cse(ch,ns,k) +
     .                       L2_cyc_cse(ch,ns,k) - fL2(lv)*dclk/fL1(lv)
                      else
                         data_lft(2,nl) = 0.d0
                      endif
                      data_lft(3,nl) = L1r_rng_cse(ch,ns,k) - dclk
                      data_lft(4,nl) = L2r_rng_cse(ch,ns,k) - 
     .                                          fL2(lv)*dclk/fL1(lv)
                      ep_lft(nl) = k
 
*                     See if we have enough data
                      if( nl.eq.max_ret ) finished = .true.
 
*                     See if this data point has a bias flag
                      if( kbit(data_flag_cse(ch,ns,k),31) .or.
     .                    kbit(data_flag_cse(ch,ns,k),32)  )
     .                        finished = .true.
                  end if
              end if
          end if
      end do
 
      num_lft = nl
 
****  Now do the right segment of data.  Same conditions apply except
*     then we check for bias flag before adding data to list
      k = ep - step_ep
      finished = .false.
      nr = 0
      do while ( .not. finished )
          k = k + step_ep
          if( k.gt.num_ep ) finished = .true.
*                                         ! See if this point OK
          if( .not.finished ) then
 
*             Get the channel number at this epoch
              ch = ltoc( ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! We have data at this time
              if( ch.gt.0 ) then
 
*                 See if this data point has a bias flag (but ignore
*                 the bias flag on the data point we are trying to
*                 patch.
                  if( (kbit(data_flag_cse(ch,ns,k),31) .or.
     .                kbit(data_flag_cse(ch,ns,k),32)) .and.
     .                k.ne.ep   )
     .                    finished = .true.
              else
                  finished = .true.
              end if
 
*             See if data OK
              if( ch.gt.0 .and. .not.finished ) then
 
*                 See if data is OK
                  if(data_OK(data_flag_cse(ch,ns,k),0, phs_mask)) then
 
*                     Add these data to the left segment
                      nr = nr + 1
                      dclk = params_cse(ns,k) - 
     .                       params_cse(num_cfiles+lv,k)
                      data_rgh(1,nr) = L1r_phs_cse(ch,ns,k) +
     .                            L1_cyc_cse(ch,ns,k) - dclk
                      if( lambda(lv,2,ns).ne.0 ) then
                         data_rgh(2,nr) = L2r_phs_cse(ch,ns,k) +
     .                       L2_cyc_cse(ch,ns,k) - fL2(lv)*dclk/fL1(lv)
                      else
                         data_rgh(2,nr) = 0.d0
                      end if
                      data_rgh(3,nr) = L1r_rng_cse(ch,ns,k) - dclk
                      data_rgh(4,nr) = L2r_rng_cse(ch,ns,k) - 
     .                                            fL2(lv)*dclk/fL1(lv)
                      ep_rgh(nr) = k
 
*                     See if we have enough data
                      if( nr.eq.max_ret ) finished = .true.
 
                  end if
              end if
          end if
      end do
 
      num_rgh = nr
 
***** Thats all
      return
      end
 
CTITLE INC_CYC_EST
 
      subroutine inc_cyc_est( ns, lv, js, kv, dat_type, obs_type, 
     .                        norm_cyc, b_cyc, chi_cyc )

      implicit none
 
*     This routine will increment the estimates of the L1 and L2
*     slip using the current data in the left and right segments.
*     The data are checked first to get the rms for the type of patch
*     to be made, and to check quality.  The estimator is then
*     incremented.  There is the possibility that the quality checking
*     may add a bias flag to the right segment data.  (A left segment
*     problem will be issued as a warning.)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ns,lv        - Site number and satellite list at which cycle 
*                  is being fixed
*   js,kv        - Second site and satellite if double differences

      integer*4 ns, lv, js, kv
 
*   norm_cyc(2,2)   - Normal equations for L1 and L2 cycle
*               - slip estimate
*   b_cyc(2)        - B-vector for estimation of number of cycle
*               - slips
*   chi_cyc     - Prefit chi squared
 
      real*8 norm_cyc(2,2), b_cyc(2), chi_cyc
 
*   dat_type    - Gives the type of information to be used in
*               - incrementing the estimator:
*               - WL - Wide lane based on phase and range
*               - NL - Narrow lane based on phase and range
*               - LC - LC phase observable
*               - LG - LG ionospheric delay.  If the cap is
*               -      large this obervable will be down
*               -      weighted relative to the others.
*               - LA - Force right segment data to have 0 average residual
*               - WA - Force right segment data to have 0 average WL.
*               - NA - Force right segmant data to have 0 average NL computed
*                      phase and range data.
*  obs_type     - Observable type:
*               - OW -- One way
*               - DD -- Double difference
 
      character*(*) dat_type, obs_type
 
 
* LOCAL VARIABLES
 
*   ep_ref      - Reference epoch for patching at.  For the
*                 one-way data this is the center of the one-way gap;
*                 for double differences it is the center of the
*                 double difference gap.

*   apart(4)        - Partials of phase and range with respect to
*               - the observable being pactched.
*   dobs, sobs  - Difference of obervable accross gap or bias
*               - flag and it sigma based on the rms of the
*               - polynomial fit
 
      real*8 ep_ref, apart(4), dobs, sobs
 
***** Scan the type and increment the estimate
      call casefold( dat_type )
      call casefold( obs_type )
      if( obs_type(1:2).eq.'OW' ) then
          ep_ref = (ep_lft(1)+ep_rgh(1))/2.0
      else
          ep_ref = (dd_ep(1,1)+dd_ep(1,2))/2.0
      end if

*     See if widelane WL
      if( index(dat_type,'WL').gt.0 ) then
          apart(1) = -1.d0
          apart(2) = +1.d0
          apart(3) = dfsf(lv)
          apart(4) = dfsf(lv)
 
          call form_data(ns,lv, js,kv, ep_ref,apart, 'WL', obs_type, 
     .                   dobs, sobs, 0)
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
* MOD TAH 041202: debug
c          if( ns.eq.1 .and.lv.eq.16 .and. js.eq.0 .and. kv.eq.0 ) then
c              print *,'INC_CYC_EST WL',dobs, sobs
c          endif 
	  
      end if
 
*     See if Narrow lane NL
      if( index(dat_type,'NL').gt.0 ) then
          apart(1) = 1.d0
          apart(2) = 1.d0
          apart(3) = -sfdf(lv)
          apart(4) = +sfdf(lv)
 
          call form_data(ns,lv,js,kv, ep_ref,apart, 'NL', obs_type,
     .                   dobs, sobs, 0)
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
      end if
 
*     See if LC observable
      if( index(dat_type,'LC').gt.0 ) then
          apart(1) = lcf1(lv)
          apart(2) = lcf2(lv)
          apart(3) = 0.d0
          apart(4) = 0.d0
          
* MOD TAH 970110: If we are patching a wide-gap, use 0 order poly
*         nomial fit.          
          if((ep_wr(1)-ep_wl(1))*sampling.le.2*dchi2_max_sep ) then
              call form_data(ns, lv,js,kv, ep_ref,apart, 'LC', 
     .                  obs_type, dobs, sobs, 1)
          else
              call form_data(ns, lv,js,kv, ep_ref,apart, 'LC', 
     .                  obs_type, dobs, sobs, 0)
          end if
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
* MOD TAH 041202: debug
c          if( ns.eq.1 .and.lv.eq.16 .and. js.eq.0 .and. kv.eq.0 ) then
c              print *,'INC_CYC_EST LC',dobs, sobs
c          endif 

      end if
      
* MOD TAH 990517: See if L1 observable only
      if( index(dat_type,'L1').gt.0 ) then
          apart(1) = 1.d0
          apart(2) = 0.d0
          apart(3) = 0.d0
          apart(4) = 0.d0
          
* MOD TAH 970110: If we are patching a wide-gap, use 0 order poly
*         nomial fit.          
          if((ep_wr(1)-ep_wl(1))*sampling.le.2*dchi2_max_sep ) then
              call form_data(ns, lv,js,kv, ep_ref,apart, 'L1', 
     .                  obs_type, dobs, sobs, 1)
          else
              call form_data(ns, lv,js,kv, ep_ref,apart, 'L1', 
     .                  obs_type, dobs, sobs, 0)
          end if
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
      end if
      
* MOD TAH 970110: Introduced new option which is LA -- Absolute fit
*     using only Right segment of data.
      if( index(dat_type,'LA').gt.0 ) then
          apart(1) = lcf1(lv)
          apart(2) = lcf2(lv)
          apart(3) = 0.d0
          apart(4) = 0.d0
          
          call form_rhda(ns, lv,js,kv, ep_ref,apart, 'LC', 
     .                   obs_type, dobs, sobs, 0)

*         If this is LA, downweight relative to WL to preserve the
*         WL ambiquity.
          if( sobs.lt.0.00001d0  ) then
              sobs = 10
          end if
          sobs = sobs*4
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
      end if      
 
*     See if widelane WL
      if( index(dat_type,'WA').gt.0 ) then
          apart(1) = -1.d0
          apart(2) = +1.d0
          apart(3) = dfsf(lv)
          apart(4) = dfsf(lv)
 
          call form_rhda(ns,lv, js,kv, ep_ref,apart, 'WL', obs_type, 
     .                   dobs, sobs, 0)
*         Remove the mean one way bias in the wide lane
*         Add bias because dobs has been changed in sign.
          dobs = dobs + WL_bias_site(ns) + WL_bias_svs(lv)
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
      end if

*     See if Narrow lane NA to applied to right segment.
      if( index(dat_type,'NA').gt.0 ) then
          apart(1) = 1.d0
          apart(2) = 1.d0
          apart(3) = -sfdf(lv)
          apart(4) = +sfdf(lv)
 
          call form_rhda(ns,lv,js,kv, ep_ref,apart, 'NL', obs_type,
     .                   dobs, sobs, 0)
 
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
      end if
 
 
*     See if LG observable
      if( index(dat_type,'LG').gt.0 ) then
          apart(1) = lgf1(lv)
          apart(2) = lgf2(lv)
          apart(3) = 0.d0
          apart(4) = 0.d0

*         if we are patching accross a gap larger than the amount
*         overw which we will be allowed to remove a bias, then
*         fit 0 order polynomial rather than a first order one. 
          if((ep_wr(1)-ep_wl(1))*sampling.le.2*dchi2_max_sep ) then
              call form_data(ns,lv,js,kv,ep_ref,apart, 'LG', obs_type, 
     .                       dobs, sobs, 1)
          else  
              call form_data(ns,lv,js,kv,ep_ref,apart, 'LG', obs_type, 
     .                       dobs, sobs, 0)
          end if
          call inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)
* MOD TAH 041202: debug
c          if( ns.eq.1 .and.lv.eq.16 .and. js.eq.0 .and. kv.eq.0 ) then
c              print *,'INC_CYC_EST LG',dobs, sobs
c         endif 

      end if
 
****  Thats all
      return
      end
 
CTITLE INC_NORM
 
      subroutine inc_norm( apart, dobs, sobs, norm_cyc, b_cyc, chi_cyc)

      implicit none
 
*     This routine will the normal equations for the estimates of the
*     number of cycle slips.
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   apart(2)        - L1 and L2 apartial with respect to the
*               - observable dobs
*   dobs        - Observable difference across boundary
*   sobs        - Sigma of the observable based on rms in
*               - form_data
*   norm_cyc(2,2)   - Normal equations
*   b_cyc(2)        - Solution vector
*   chi_cyc     - Prefit chi**2
 
 
      real*8 apart(2), dobs, sobs, norm_cyc(2,2), b_cyc(2), chi_cyc
 
* LOCAL PARAMETERS
 
*   i,j     - Loop counters
 
      integer*2 i,j
 
*   wgh     - Weight of observation in LSQ estimator
 
      real*8 wgh
 
****  Increment the normal equations and solution vector
      if( sobs.le.0 ) then
          write(*,*) 'SOBS BAD', sobs 
          RETURN 
      end if 
      wgh = 1.d0 / sobs**2
      do i = 1,2
          b_cyc(i) = b_cyc(i) + apart(i)*dobs*wgh
          do j = 1,2
              norm_cyc(i,j) = norm_cyc(i,j) + apart(i)*apart(j)*wgh
          end do
      end do

      chi_cyc = chi_cyc + dobs**2*wgh
 
****  Thats all
      return
      end
 
CTITLE FORM_DATA
 
      subroutine form_data( ns, lv,js,kv, ep_ref, apart, dat_type, 
     .                      obs_type, dobs, sobs, poly_order )

      implicit none
 
*     This routine will form the obervable given by the partials in
*     apart from the left and right segments of data. A poly nominal of
*     order poly_order is fit to segments and the polynomial is
*     evaluated at the midpointbetween the left and right segments.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ns,lv          - Site and satellite list number at which data is 
*                    being formed 
*   js,kv          - Second site and satellite for double diffs.
*   poly_order  - Order of polynomial to fit (0 is mean,
*               - 1 is slope etc.)
 
      integer*4 ns, lv, js,kv, poly_order
 
*   ep_ref    - Refernce epoch for computations (may be
*           - fractional)
*   apart(4)        - Partials of obervable with respect to L1,L2
*               - P1 and P2
*   dobs        - O-C computed at the midpoint of the spans of
*               - data.
*   sobs        - Sigma of dobs based on the quality of the fit
*               - to data in each segment.
 
      real*8 ep_ref, apart(4), dobs, sobs
 
*   dat_type    - Gives the type of information to be used in
*               - incrementing the estimator:
*               - WL - Wide lane based on phase and range
*               - NL - Narrow lane based on phase and range
*               - LC - LC phase observable
*               - LG - LG ionospheric delay.  If the cap is
*               -      large this obervable will be down
*               -      weighted relative to the others.
*  obs_type     - Observable type:
*               - OW -- One way
*               - DD -- Double difference
 
      character*(*) dat_type, obs_type

* LOCAL PARAMETERS
 
*   i,j     - loop counters
*   num_ul, num_ur  - Number of values to use in the left and right
*             segments (number reduced for the LG observables)
 
      integer*4 i,j, num_ul, num_ur
 
*   obs(max_max_wl_ret)     - Observable for left or right
*           - segment of data.
*   obs_refr, obs_refl  - Value of polynomial fit to the right
*           - and left segments of data at the reference epoch
*   sig_refr, sig_refl  - Sigmas of the obs_refr and obs_refl
*           - values
*   apl(max_patch_poly,max_patch_poly), bpl(max_patch_poly) -
*             Covariance matrix and solution from left data segment
*   apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly) -
*             Covariance matrix and solution from right data segment

 
      real*8 obs(max_max_wl_ret), obs_refr, obs_refl,
     .    sig_refr, sig_refl, 
     .    apl(max_patch_poly,max_patch_poly), bpl(max_patch_poly),
     .    apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly)
 
***** Start by forming the observable in the left and right segments
*     fiting the polynomial
 
*     Get the number of valuess to use.  Use all of not LG dat_type
      if( dat_type(1:2).ne.'LG' ) then
          num_ul = num_wl
          num_ur = num_wr
      else
          num_ul = min(max_lg_use, num_wl)
          num_ur = min(max_lg_use, num_wr)
      end if
 
      do i = 1, num_ul
          obs(i) = 0
          do j = 1, 4
              obs(i) =  obs(i) + apart(j)*data_wl(j,i)
          end do
      end do
 
****  Now fit to the left segment and check the quality
      call fit_obs( ns, lv,js,kv, obs, ep_wl, num_ul, dat_type, 
     .              obs_type, poly_order, ep_ref, obs_refl, sig_refl,
     .              apl, bpl )

*     Update the working values for the number of data in the segment
*     in case we flagged data in the fit obs routine (don't do for
*     LG since we don't flag for this data type).
      if( dat_type(1:2).ne.'LG' ) then
          num_wl = num_ul
      end if

****  Now fit to the right segment of data
      do i = 1, num_ur
          obs(i) = 0
          do j = 1, 4
              obs(i) =  obs(i) + apart(j)*data_wr(j,i)
          end do
      end do
 
****  Now fit to the left segment and check the quality
      call fit_obs( ns, lv, js, kv, obs, ep_wr, num_ur, dat_type,
     .              obs_type,poly_order, ep_ref, obs_refr, sig_refr,
     .              apr, bpr )

*     Update the working values for the number of data in the segment
*     in case we flagged data in the fit obs routine (don't do for
*     LG since we don't flag for this data type).
      if( dat_type(1:2).ne.'LG' ) then
          num_wr = num_ur
      end if

****  Now get the jump values for this data types.  Change the sign of
*     dobs so we estimate the number of cycles to remove (rather than
*     the number accross the gap).
      if( dat_type(1:2).ne.'LC' .and. dat_type(1:2).ne.'L1' ) then
*         Get a simple estimate (do not force polynomials to be the same)
          dobs = -(obs_refr - obs_refl)
          sobs = sqrt(sig_refr**2 + sig_refl**2)
      else
          call force_lc_poly(poly_order+1,apl,bpl,apr,bpr, dobs, sobs)
      end if
      if( sobs.le.0 ) then
          write(*,*) 'SOBS FDD', sobs, sig_refr, sig_refl, dat_type
          write(*,*) apl, bpl, poly_order, apr, bpr 
      end if
 
***** Thats all
      return
      end
 
CTITLE FORM_RHDA
 
      subroutine form_rhda( ns, lv,js,kv, ep_ref, apart, dat_type, 
     .                      obs_type, dobs, sobs, poly_order )

      implicit none
 
*     This routine will form the obervable given by the partials in
*     apart for only the right segment of data. Only a mean value
*     is fit to the data (ie., poly_order=0)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ns,lv          - Site and satellite list number at which data is 
*                    being formed 
*   js,kv          - Second site and satellite for double diffs.
*   poly_order  - Order of polynomial to fit (0 is mean,
*               - 1 is slope etc.)
 
      integer*4 ns, lv, js,kv, poly_order
 
*   ep_ref    - Refernce epoch for computations (may be
*           - fractional)
*   apart(4)        - Partials of obervable with respect to L1,L2
*               - P1 and P2
*   dobs        - O-C computed at the midpoint of the spans of
*               - data.
*   sobs        - Sigma of dobs based on the quality of the fit
*               - to data in each segment.
 
      real*8 ep_ref, apart(4), dobs, sobs
 
*   dat_type    - Gives the type of information to be used in
*               - incrementing the estimator:
*               - WL - Wide lane based on phase and range
*               - NL - Narrow lane based on phase and range
*               - LC - LC phase observable
*               - LG - LG ionospheric delay.  If the cap is
*               -      large this obervable will be down
*               -      weighted relative to the others.
*  obs_type     - Observable type:
*               - OW -- One way
*               - DD -- Double difference
 
      character*(*) dat_type, obs_type

* LOCAL PARAMETERS
 
*   i,j     - loop counters
*   num_ul, num_ur  - Number of values to use in the left and right
*             segments (number reduced for the LG observables)
 
      integer*4 i,j, num_ur
 
*   obs(max_max_wl_ret)     - Observable for left or right
*           - segment of data.
*   obs_refr, obs_refl  - Value of polynomial fit to the right
*           - and left segments of data at the reference epoch
*   sig_refr, sig_refl  - Sigmas of the obs_refr and obs_refl
*           - values
*   apl(max_patch_poly,max_patch_poly), bpl(max_patch_poly) -
*             Covariance matrix and solution from left data segment
*   apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly) -
*             Covariance matrix and solution from right data segment

 
      real*8 obs(max_max_wl_ret), obs_refr,  sig_refr,
     .    apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly)
 
***** Start by forming the observable in the left and right segments
*     fiting the polynomial
 
      num_ur = num_wr
 
****  Now fit to the right segment of data
      do i = 1, num_ur
          obs(i) = 0
          do j = 1, 4
              obs(i) =  obs(i) + apart(j)*data_wr(j,i)
          end do
      end do
 
****  Now fit to the left segment and check the quality
      call fit_obs( ns, lv, js, kv, obs, ep_wr, num_ur, dat_type,
     .              obs_type,poly_order, ep_ref, obs_refr, sig_refr,
     .              apr, bpr )

****  Now get the jump values for this data types.  Change the sign of
*     dobs so we estimate the number of cycles to remove (rather than
*     the number accross the gap).
      dobs = -obs_refr
      sobs = sig_refr
 
***** Thats all
      return
      end
 
CTITLE FIT_OBS
 
      subroutine fit_obs( ns, lv, js, kv, obs, ep, num, dat_type, 
     .                    obs_type, po, ep_ref, obs_ref, sig_ref, 
     .                    ap, bp )

      implicit none
 
*     This routine will fit a polynomial of order po to the obs data
*     and return the obs at the reference epoch.  This routine will
*     also check the quality of the fit.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ns,lv,js,kv      - Site number and satellite numbers in dd.
*   num     - Number of data in segment
*   ep(num) - Epoch numbers of data in the segment
*   po      - Polynomial order to fit to the data
 
      integer*4 ns, lv,js,kv, num, ep(num), po
 
*   obs(num)    - Observed values
*   ep_ref  - Reference epoch to the fit to the data
*   obs_ref - polynomial value at the reference epoch
*   sig_ref - Sigma at the reference epoch
 
      real*8 obs(num), ep_ref, obs_ref, sig_ref

*   ap(max_patch_poly, max_patch_poly)  - Normal equations for
*               - fitting polynomial
*   bp(max_patch_poly)  - Solution vector for poly fit

      real*8 ap(max_patch_poly, max_patch_poly), bp(max_patch_poly)
 
*   dat_type    - Gives the type of information to be used in
*               - incrementing the estimator:
*               - WL - Wide lane based on phase and range
*               - NL - Narrow lane based on phase and range
*               - LC - LC phase observable
*               - LG - LG ionospheric delay.  If the cap is
*               -      large this obervable will be down
*               -      weighted relative to the others.
*  obs_type     - Observable type:
*               - OW -- One way
*               - DD -- Double difference
 
      character*(*) dat_type, obs_type

* LOCAL PARAMETERS
 
*   i,j,k       - Loop counters
*   nus         - Number of data to use (as we iterate we may
*               - discard some data due to jumps)
*   ipivot(max_patch_poly)  - Pivot values in inversion
*   imax        - Point number with the largest jump
*   ep_jump     - The epoch that we should the bias flag at
*                 when a jump is found
*   ep_comp     - Epoch next to the jump for comparison with
 
      integer*4 i,j,k, nus, ipivot(max_patch_poly), imax,
     .          ep_jump, ep_comp
 
*   scale(max_patch_poly)   - Scale factors in inversion
*   res, prev_res       - Residual to poly fit at current and
*                   - previous epoch (latter used for rms of
*                   - differences)
*   res_save(max_max_wl_ret) - Array of residuals for further
*                     checking.
*   ressq, diffsq   - Sum of residual and residual differences
*               - squares
*   rms_res, rms_diff   - RMS of residuals and their differences
*   max_jump        - Max jump in residuals while looking for jumps
*   used_tol        - Tolerance actually used in looking for jumps
*   dep             - difference in epoch numbers for computing
*                     polynomials.
 
      real*8 scale(max_patch_poly), res, prev_res, ressq, diffsq,
     .    rms_res, rms_diff, max_jump, res_save(max_max_wl_ret),
     .    used_tol, dep
 
*   finished        - Indicates that we have finished looking for
*               - jumps (based on total rms being greater than
*               - rms of differences)
*   flg_jump    - Logical to indicate that we should flag data
*                 associated wuth a jump
*   jump_already_found - True if we have already recorded the jump.
 
 
      logical finished, flg_jump, jump_already_found

* PT960925: add a tmp char variable to get around compile problems with
*           string concatenation using g77

      character*4 tmp_var
 
****  Keep repeating the fit until all looks Ok (or we run out of data)
      finished = .false.
* MOD TAH 041202: debug
c          if( ns.eq.1 .and.lv.eq.16 .and. js.eq.0 .and. kv.eq.0 ) then
c              print *,'FITOBS',dat_type,num,(obs(j),j=1,num)
c          endif 

 
*     set the used number of data to the full set processed
      nus = num
      imax = 1
      do while ( .not.finished )
 
*         Make sure we have enough data
          if( nus.lt.po+3 ) then
*             We don;t have enough data; set the obs_ref to zero
*             give it a very large sigma
              if( nus.gt.0 ) then
                  obs_ref = obs(1)
              else
                  obs_ref = 0
              end if
              sig_ref = 1.d2
              num = nus

* Reset the polynomial estimation coefficients
              do i = 1, po+1
                 bp(i) = 0.0d0
                 do j = 1, po+1
                    ap(i,j) = 0.d0
                    if( i.eq.j ) ap(i,j) = sig_ref**2
                 end do
              end do 
              RETURN
          end if
 
*         Clear the polyomial fit arrays
          do i = 1, po+1
              bp(i) = 0.d0
              do j = 1, po+1
                  ap(i,j) = 0.d0
              end do
          end do
 
*         Accumulate current data into fit
* MOD TAH 980513: Use non-zero dep to avoild 0**0 Dec problems.
          do i = 1, nus
              dep = ep_ref-ep(i) + 0.001d0
              do j = 1,po+1
                  bp(j) = bp(j) + obs(i)*(dep)**(j-1)
                  do k = 1,po+1
                      ap(j,k) = ap(j,k) + (dep)**(j+k-2)
                  end do
              end do
          end do
 
*****     Now solve for the polynomial
          call invert_vis(ap,bp, scale, ipivot, po+1, max_patch_poly,1)

*         Clear the rms square variables
          ressq = 0.d0
          diffsq = 0.d0
 
*****     Now check the quality of the data.  Here we compute the rms
*         of the residual and their time differences
          do i = 1, nus
              res = obs(i)
              dep = ep_ref-ep(i) + 0.001d0
              do j = 1, po+1
                  res = res - bp(j)*(dep)**(j-1)
              end do
              res_save(i) = res
 
****          Accumulate stats and if this is after first point
*             accumulate difference stats
              ressq = ressq + res**2
              if( i.gt.1 ) then
                  diffsq = diffsq + (res-prev_res)**2
              end if
              prev_res = res
          end do
 
****      Get the rms scatter
          rms_res = sqrt(ressq/(nus-po-1))
*         For the difference rms we divide by the same degrees of
*         freedom because we lost an order in the polyomial when we
*         differenced.
          rms_diff = sqrt(diffsq/(nus-po-1))
 
*         Now the rms_diff should be greater than the rms_res by
*         square root of 2 (1.414).  If it is less than the rms_res
*         then we may have a jump in the data.  If this is the case
*         find the largest jump and save the epoch for later recording.
*         Then remove the data after the "bad" point (largest jump) and
*         try again.  Only use this test for Double difference LC
*         and one-way WL.
          tmp_var = dat_type//obs_type
          if( (rms_diff.lt.rms_res .or. 
     .        (rms_res.gt.0.1 .and. obs_type.eq.'DD')) .and.
     .        (tmp_var.eq.'WLOW' .or. tmp_var.eq.'LCDD')     ) then
c     .        (dat_type//obs_type.eq.'WLOW' .or.
c     .         dat_type//obs_type.eq.'LCDD')     ) then

*             Search for the largest jump and do not use this data
              max_jump = 0
              prev_res = res_save(1)
              do i = 2, nus
                  res = res_save(i) 
 
*                 If this is after the first point check the
*                 size of the jump
                  if( abs(max_jump).lt.abs(res-prev_res)) then
                      imax = i
                      max_jump = res-prev_res
                  end if
                  prev_res = res
              end do
 
****          Now save the epoch of the jump and its value, and
*             reset the number of data to use.

*             See if we should fix jump.
              flg_jump = .false.

*             First check if these are double differnce, LC. The
*             tolerance is actually computed in the check_jump_bad
*             routine.
              if( obs_type.eq.'DD' ) then
                  call check_jump_bad( ns, nus, imax, ep, res_save, 
     .                 dd_lc_tol, used_tol, flg_jump, dat_type  )
              end if

              if( obs_type.eq.'OW' ) then
                  call check_jump_bad( ns, nus, imax, ep, res_save, 
     .                 dd_wl_tol, used_tol, flg_jump, dat_type  )
              end if

*             See if we really have a jump to process.
              if( flg_jump ) then 

*                 Now the epoch of the jump depends on if we are in
*                 left or right segment
                  if( ep(1).lt.ep_ref ) then
                      ep_jump = ep(imax-1)
                      ep_comp = ep(imax)
                  else
                      ep_jump = ep(imax)
                      ep_comp = ep(imax-1)
                  end if 

                  nus = imax -1

*****             See if we have already found this one.
                  jump_already_found = .false.
                  do i = 1, unflg_num
                      if( unflg_ep(1,i).eq. ep_jump ) then
                          jump_already_found = .true.
                      end if
                  end do


*                 If we have already found the jump then ignore.
                  if( .not. jump_already_found ) then 
                      unflg_num = unflg_num + 1
                      if( unflg_num.gt. max_unflg ) then
                          write(*,300) (i, unflg_ep(1,i), unflg_mag(i), 
     .                                unflg_type(i), i = 1,unflg_num-1)
 300                      format(' **WARNING** Too many unflaged jumps',
     .                        ' replacing last one: ',/,
     .                        '  #    Epoch   Magnitude (cycles) type',/,
     .                        200(I4,1x,i4,1x,F20.2,2x,a4:/))
                          unflg_num = unflg_num -1
                      end if

                      unflg_ep(1,unflg_num) = ep_jump
                      unflg_ep(2,unflg_num) = ep_comp
                      unflg_sv(1,unflg_num) = js
                      unflg_sv(2,unflg_num) = kv
                      unflg_type(unflg_num) = dat_type // obs_type
                      unflg_mag(unflg_num) = max_jump
                      if( obs_type.eq.'OW' )
     .                write(*,320) max_jump, cf_codes(ns), prn_list(lv),
     .                             unflg_ep(1,unflg_num),
     .                             dat_type,obs_type, used_tol
 320                  format(' JUMP of ',F9.2,' cycles found at site ',
     .                     a4,' PRN ',i2.2,'  Epoch ',I4,
     .                        ' with data type ',2a2,
     .                        ': tol ',f6.2)
                      if( obs_type.eq.'DD' )
     .                write(*,330) max_jump, cf_codes(ns), prn_list(lv),
     .                             unflg_ep(1,unflg_num),
     .                             cf_codes(js), prn_list(kv),
     .                             dat_type,obs_type, used_tol
 330                  format(' JUMP of ',F9.2,' cycles found at site ',
     .                     a4,' PRN ',i2.2,' Epoch ',I4,' with ',a4,
     .                        ' PRN ',i2.2,' and data type ',2a2,
     .                       ': tol ',f6.2)
                   end if
              else
                  finished = .true.
              end if
          else
 
****          All looks OK, so say we are finished
              finished = .true.
          end if
      end do

***** Scale the covaraince matrix by the rms
      do i = 1, po+1
         do j = 1, po+1
            ap(i,j) = ap(i,j)*(ressq/(nus-po-1))
         end do
      end do
 
****  Return the observable and its sigma
      obs_ref = bp(1)
      sig_ref = sqrt(ap(1,1))
      num = nus
 
***** Thats all
      return
      end
 
CTITLE est_dd_cyc
 
      subroutine est_dd_cyc( ep, ns, lv, js, kv, norm_cyc, b_cyc, 
     .                chi_cyc, reliable, chiqual, dL1_slip, dL2_slip )

      implicit none
 
*     This routine will take the accumulation information on the
*     number of cycle slips accross a gap or bias flag and form the
*     integer estimate
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   ep          - Epoch at which bias is being tested. 
*   ns, lv      - Station number and satellite list number
*   js, kv      - Second station and satellite if double difference
*                 (js=0 if one-way)
 
      integer*4 ep, ns, lv, js, kv
 
*   norm_cyc(2,2), b_cyc(2) - Normal equations and solution
*               - vector for the estimate
*   chi_cyc     - Prefit chi**2
*   dL1_slip, dL2_slip      - Estimated cycle slips at L1 and
*               - L2 modulo the wavelength factor for the site
*               - satellites
*   chiqual     - Quantitiave estimate of quality of fix.

      real*8 norm_cyc(2,2), b_cyc(2), dL1_slip, dL2_slip,
     .       chiqual, chi_cyc
 
*   reliable        - Logical to indicate the slip is fixed
*               - reliably.
*   kbit        - Checks to see if bit is set
 
      logical reliable, kbit
 
* LOCAL PARAMETERS
 
*   det     - Determinate of the normal equations
*   cyc_est(2)  - Floating point estimates of the number of
*               - cycles at L1 and L2
*   cov_est(2,2)    - Covariance matrix of the estimated number
*               - of cycles
*   rho_12      - Correlation between the estimates
*   dL1_0, dL2_0    - Estimates of integer number of cycles with
*               - dL1_0 obtained directly and dL2_0 based on
*               - change imposed at L1
*   L1, L2      - Trial values of L1 and L2 cycles which we
*               - use to constrast the chi**2 values
*   L2_est      - Floating point estimate of L2 cycles based on
*               - the integer number of L1 cycles
*   dL1, dL2        - Differences between the L1 and L2 integer
*               - values and the floating point estimates of
*               - their value
*   dchi2min(2) - Two smallest changes in chi**2 with 1 being
*               - the smallest
*   L1_min(2), L2_min(2)    - Number of cycles corresponding to the
*               - the two minium chi**2
*   dchi2       - Change in chi**2 as value taken to integer
*   arg         - Angular argument used for atan of ratio of gap to
*                 minumum data.  Used to scale the chiqual for a large
*                 gap with little data
*   chiscale    - scaling on chiqual.
*   chi_fin     - Final chi**2 for fitting to data.

 
      real*8 det, cyc_est(2), cov_est(2,2), rho_12, dL1_0, dL2_0,
     .    L1, L2, L2_est, dL1, dL2, dchi2min(2), L1_min(2), L2_min(2),
     .    dchi2, arg, chiscale, chi_fin

      real*8 wf1, wf2  ! Wavelength factor (for codeless).  Introcuded to make
                 ! code easier to read for Glonass. L1 and L2 value
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

* MOD TAH 200505: Added to remove floating do loops.
      integer*4 CR,   ! cycle ranges in terms of wave lengths
     .          ic    ! loop variable for testing cycles.

*   dont_make_reliable - Logical to indicate that the det is
*      near singular on the estimate so don't set the reliable flag.
*      We do this to "flattenen" out LG even when there is no 
*      range data ans we only have a one-way LG fix.

      logical dont_make_reliable
 
***** Solve the system of equations to get the floating point
*     estimates of the number of cycles.  Check the det. first to
*     make sure not singular
 
      reliable = .false.
      dont_make_reliable = .false.

*     Set the default chiqual to a very small value (but not zero)
*     this way we know that some double difference was formed.
      chiqual = 0.01d0
 
      det = norm_cyc(1,1)*norm_cyc(2,2) - norm_cyc(1,2)**2
*                                         ! Estimate is singular
      if( det.lt.max(1.d-4,1.d-8*norm_cyc(1,1)) ) then

*         Tell routine not to make reliable fix
          dont_make_reliable = .true.
 
*         Force the number of L1 cycles to zero and solve for L2
*         cycles
          cyc_est(1) = 0.d0

*         See if any data.  Changed to L1 check. 990517:  Added
*         L2 only check to handle L1 only data
          if( norm_cyc(1,1).eq.0.d0 ) then
* TEST MOD: Do Not reset
c             dL1_slip = 0
c             dL2_slip = 0
              RETURN
          end if
          
* MOD TAH 990517: Check to see if we can fix only L1 cycles (i.e., no L2 data).
          if( norm_cyc(2,2).eq.0.d0 ) then 
*             This is a much simplier decision function so complete
*             all the calculations here:
              dont_make_reliable = .false.
              cyc_est(2) = 0.d0
              cyc_est(1) = b_cyc(1)/norm_cyc(1,1) 
* MOD TAH 180320: Introduced re-mapping factors for Glonass
C             dL1_0 = float(nint(cyc_est(1)*abs(lambda(lv,1,ns))))/
C     .               abs(lambda(lv,1,ns))
              wf1 = abs(lambda(lv,1,ns))
              fr1 = fL1u(lv)/fL1(lv)
* MOD TAH 200617: Added fr1 factor to convert back to nominal freq.
              dL1_0 = float(nint(cyc_est(1)*wf1*fr1))/(wf1*fr1)
              dchi2min(1) = 1.d10
              dchi2min(2) = 1.d10
* MOD TAH 180320: L1 only code (not updated for Glonass)
* MOD TAH 200505: Removed floating point loop.
*             do L1 = dL1_0-1, dL1_0+1, 1/wf1
              CR = nint(1/wf1)
              do ic = -CR, CR
* MOD TAH 200520: Fixed line below.
                 L1 = dL1_0 + ic*wf1
                 dL1 = L1/fr1-cyc_est(1)
                 dchi2 = norm_cyc(1,1)*dL1**2
                 if( dchi2.lt.dchi2min(1) ) then
*                    Move the current samllest to the next smallest
                     dchi2min(2) = dchi2min(1)
                     L1_min(2) = L1_min(1)
                     L2_min(2) = 0.d0
 
*                    Save the smallest values
                     dchi2min(1) = dchi2
                     L1_min(1) = L1
                     L2_min(1) = 0.d0
*                See if smaller than the next smallest
                 else if( dchi2.lt.dchi2min(2) ) then
 
*                    Save values
                     dchi2min(2) = dchi2
                     L1_min(2) = L1
                     L2_min(2) = 0.d0
                end if
             end do
             dL1_slip = L1_min(1)/fr1
             dL2_slip = 0.d0
             
*            Now quantify the quality of the patch.  This uses the ratio
*            next-best to a modified best chi**2 value.  The modification
*            based on setting a min_chi**2 to allow.  Limit the exp value
*            so that we won't get underflows.

             if( dchi2min(1)/dchi2_min_val.lt.10 ) then
                 chiqual = dchi2min(2)/
     .             (dchi2min(1)+dchi2_min_val*
     .                          exp(-dchi2min(1)/dchi2_min_val))
             else
                 chiqual = dchi2min(2)/dchi2min(1)
             end if

*            Get the scaling for chiqual
*            MOD TAH 980513: Make sure that arg does not go to infinity.
             if( min(num_wl,num_wr).gt.0 ) then
                 arg = float(ep_wr(1)-ep_wl(1))/
     .                 float(min(num_wl,num_wr))
             else
                 arg = 1.d4
             end if
             chiscale = 1.d0 + dchi2_gap_fact*atan(arg)
             chiqual = chiqual/chiscale

*            Now based on the chi**2 values see if reliable.
             if( kbit(status_rep,6) )
     .           write(*,200) ep, ns, lv, js, kv, 
     .               dL1_slip, dL2_slip,
     .               num_wl, num_wr, ep_wl(1), ep_wr(1), dchi2min,
     .               chiqual
 200             format(' Ep ',i5,' S1/C1 ',2i3,' S2/C2 ',2i3,
     .                  ' dL1   SLP ',2F9.2,' cyc ',
     .                  ' NumLR ',2i3,' EpLR ',2i6,' dchi, ChiQ ',
     .                    3(F6.1,1x))

*            See if actual estimates are wanted:
             if( kbit(status_rep,7) )
     .            write(*,210) ep, ns, lv, js, kv, 
     .                 cyc_est, sqrt(cov_est(1,1)), sqrt(cov_est(2,2)),
     .                 ep_wl(1), ep_wr(1), chiqual
 210              format(' Ep ',i5,' S1/C2 ',2i3,' S2/C2 ',2i3,
     .                   ' dL1   ests ',2F8.2,' +- ',2F6.2,
     .                   ' EpLR ',2i6,' Chiqual ',F7.1,1x)

             if( chiqual.ge.dchi2_ratio .and. 
     .          (ep_wr(1)-ep_wl(1))*sampling.lt.dchi2_max_sep ) then

*               This where we set the bias as being reliably resolved.  
*               Make sure the don't make reliable flag is not set.
                if( .not. dont_make_reliable ) then
                   reliable = .true.
                end if
             end if

*            If this is a reliable fix and the user wants ALL or
*            FIXED biases output, then write this one (even though
*            it is not a double difference)
             if( reliable .and. (dd_out_opts(1:1).eq.'A' .or.
     .                    dd_out_opts(1:1).eq.'F')  ) then
                 call write_dd_scan(202, ep, lv, kv, ns, js, dL1_slip,
     .                       dL2_slip, chiqual, reliable)
             end if
             if( .not.reliable .and. js.gt.0 .and.
     .            (dd_out_opts(1:1).eq.'A' .or. 
     .             dd_out_opts(1:1).eq.'N')) then
                 call write_dd_scan(202, ep, lv, kv, ns, js, dL1_slip,
     .                       dL2_slip, chiqual, reliable)
             end if
 
****         Thats all
             RETURN
              
          else     
                          
             cyc_est(2) = b_cyc(2)/norm_cyc(2,2)
 
             cov_est(1,1) = 1.d0
             cov_est(2,2) = 1.d0/norm_cyc(2,2)
             cov_est(1,2) = 0.d0
             cov_est(2,1) = 0.d0
          end if
 
      else
 
*         Solve the system
          cyc_est(1) = ( norm_cyc(2,2)*b_cyc(1) -
     .                norm_cyc(1,2)*b_cyc(2) )/det
          cyc_est(2) = ( norm_cyc(1,1)*b_cyc(2) -
     .                norm_cyc(1,2)*b_cyc(1) )/det
 
          cov_est(1,1) = norm_cyc(2,2)/det
          cov_est(2,2) = norm_cyc(1,1)/det
          cov_est(1,2) = - norm_cyc(1,2)/det
          cov_est(2,1) = - norm_cyc(2,1)/det
 
      end if

* MOD TAH 950822: Compute the increment in chi**2 when the offsets
*     are forced to be the same from all the observables.  NOTE:
*     when on LC and LG or LG and WL are used this should be zero.
      chi_fin = chi_cyc - (cyc_est(1)*b_cyc(1)+cyc_est(2)*b_cyc(2))
 
*     Now try to resolve to the integral number of wavelengths
*     estimate the reliablity.
 
      rho_12 = cov_est(1,2)/sqrt(cov_est(1,1)*cov_est(2,2))
 
****  Now try to fix to integer numbers of cycles
      dchi2min(1) = 1.d10
      dchi2min(2) = 1.d10
* MOD TAH 180320: Update to allow for frequency ratios when 
*     Gkonass re-mapped frequencies are used.
C     dL1_0 = float(nint(cyc_est(1)*abs(lambda(lv,1,ns))))/
C    .                abs(lambda(lv,1,ns))
      wf1 = abs(lambda(lv,1,ns))
      fr1 = fL1u(lv)/fL1(lv)
*     Keep as integer for looping.(just apply the wavelength factor)
* MOD TAH 200617: Added fr1 factor back to nominal freq.
      dL1_0 = float(nint(cyc_est(1)*wf1*fr1))/(wf1*fr1)
 
*     Now loop around these values getting the smallest increment
*     in chi**2
* MOD TAH 200505: Replaced floating point loop
C     do L1 = dL1_0-3,dL1_0+3 , 1/wf1
      CR = nint(1/wf1)
      do ic = -3*CR, 3*CR
* MOD TAH 200520: Fixed line below. L1 value at nominal frequency
*         but need to step by fLlu cycles (divide by fr1)
          L1 = dL1_0 + ic*wf1/fr1
 
*         Make estimate of L2 based on L1 value
* MOD TAH 180320: Map the integer value back to original frequency
* MOD TAH 200617: Need to keep at nominal frequency (cyc_est freq).
*         Removed the fr1 factor (included above)
          dL1 = L1 - cyc_est(1)
          L2_est = cyc_est(2) + rho_12*dL1*
     .                        sqrt(cov_est(2,2))/
     .                        sqrt(cov_est(1,1))
* MOD TAH 180320: Update for mapped frequencies
C         dL2_0 = float(nint(L2_est*abs(lambda(lv,2,ns))))/
C    .                    abs(lambda(lv,2,ns))
          wf2 = abs(lambda(lv,2,ns))
          fr2 = fL2u(lv)/fL2(lv)
*         Keep the dL2_0 value as integer (don't /fr2
* MOD TAH 200617: Same logic as above dL2_0 is in nominal
*         frequency unit.
          dL2_0 = float(nint(L2_est*wf2*fr2))/(wf2*fr2)

*****     Compute chi**2 change:  No need to search L2
*         since value computed above is the least squares
*         result. 
          L2 = dL2_0
*         Map back now. 
* MOD TAH 200617: Scale factor included above so not needed.
          dL2 = L2 - cyc_est(2)
          dchi2 =  norm_cyc(1,1)*dL1**2 +
     .             norm_cyc(2,2)*dL2**2 +
     .           2*norm_cyc(1,2)*dL1*dL2 + abs(chi_fin)
 
****      Save the top two dchi**2. See if it smaller than the
*         smallest
          if( dchi2.lt.dchi2min(1) ) then
 
*             Move the current samllest to the next smallest
              dchi2min(2) = dchi2min(1)
              L1_min(2) = L1_min(1)
              L2_min(2) = L2_min(1)
 
*             Save the smallest values
              dchi2min(1) = dchi2
              L1_min(1) = L1
              L2_min(1) = L2
*             See if smaller than the next smallest
          else if( dchi2.lt.dchi2min(2) ) then
 
*             Save values
              dchi2min(2) = dchi2
              L1_min(2) = L1
              L2_min(2) = L2
          end if
      end do
 
****  Save the number of cycles in slips for the min chi**2 change
* MOD TAH 180320: Map the dL1/2_slip values back to original frequency
* MOD TAH 200617: Not needed, values are computed in nominal freq above.
C     dL1_slip = L1_min(1)/fr1
C     dL2_slip = L2_min(1)/fr2
      dL1_slip = L1_min(1)
      dL2_slip = L2_min(1)

*     Now quantify the quality of the patch.  This uses the ratio
*     next-best to a modified best chi**2 value.  The modification
*     based on setting a min_chi**2 to allow.  Limit the exp value
*     so that we won't get underflows.

      if( dchi2min(1)/dchi2_min_val.lt.10 ) then
          chiqual = dchi2min(2)/
     .             (dchi2min(1)+dchi2_min_val*
     .                          exp(-dchi2min(1)/dchi2_min_val))
      else
          chiqual = dchi2min(2)/dchi2min(1)
      end if

*     Get the scaling for chiqual
* MOD TAH 980513: Make sure that arg does not go to infinity.
      if( min(num_wl,num_wr).gt.0 ) then
          arg = float(ep_wr(1)-ep_wl(1))/float(min(num_wl,num_wr))
      else
          arg = 1.d4
      end if
      chiscale = 1.d0 + dchi2_gap_fact*atan(arg)
      chiqual = chiqual/chiscale

*     Now based on the chi**2 values see if reliable.
      if( kbit(status_rep,6) )
     .    write(*,300) ep, ns, lv, js, kv, 
     .            dL1_slip, dL2_slip,
     .            num_wl, num_wr, ep_wl(1), ep_wr(1), dchi2min,
     .            chiqual
 300      format(' Ep ',i5,' S1/C1 ',2i3,' S2/C2 ',2i3,
     .       ' dL1/2 SLP ',2F9.2,' cyc ',
     .       ' NumLR ',2i3,' EpLR ',2i6,' dchi, ChiF ',3(F6.1,1x))

*     See if actual estimates are wanted:
      if( kbit(status_rep,7) )
     .    write(*,310) ep, ns, lv, js, kv, 
     .            cyc_est, sqrt(cov_est(1,1)), sqrt(cov_est(2,2)),
     .            ep_wl(1), ep_wr(1), chiqual
 310      format(' Ep ',i5,' S1/C2 ',2i3,' S2/C2 ',2i3,
     .       ' dL1/2 ests ',2F8.2,' +- ',2F6.2,
     .       ' EpLR ',2i6,' Chiqual ',F7.1,1x)

      if( chiqual.ge.dchi2_ratio .and. 
     .    (ep_wr(1)-ep_wl(1))*sampling.lt.dchi2_max_sep ) then

*         This where we set the bias as being reliably resolved.  Make
*         sure the don't make reliable flag is not set.
          if( .not. dont_make_reliable ) then
              reliable = .true.
          end if
      end if

*     If this is a reliable fix and the user wants ALL or
*     FIXED biases output, then write this one (even though
*     it is not a double difference)
      if( reliable .and. (dd_out_opts(1:1).eq.'A' .or.
     .                    dd_out_opts(1:1).eq.'F')  ) then
          call write_dd_scan(202, ep, lv, kv, ns, js, dL1_slip,
     .                       dL2_slip, chiqual, reliable)
      end if
      if( .not.reliable .and. js.gt.0 .and.
     .    (dd_out_opts(1:1).eq.'A' .or. dd_out_opts(1:1).eq.'N')) then
          call write_dd_scan(202, ep, lv, kv, ns, js, dL1_slip,
     .                       dL2_slip, chiqual, reliable)
      end if
 
****  Thats all
      return
      end

CTITLE CHECK_JUMP_BAD

      subroutine check_jump_bad( ns, nus, imax, ep, res_save, 
     .                           tol_par, tol, flg_jump, dat_type  )

      implicit none

*     This routine will check to see if the maxium flagged jump is
*     just bad point or truly a jump in the residuals.  It checks
*     the values stradling the bad point to see if they are OK.  If
*     the bad point is at the end of sequence it flags it as a jump.

* INCLUDES

      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   ns  - Site number (for output)
*   nus   - Number of data points in sequence
*   imax  -  Site with maximum jump
*   ep(nus) - List of epochs 

      integer*4 ns, nus, imax, ep(nus)

*   res_save(nus)  - Residual to polynomial fit
*   tol_par(3) - (1) scale, (2) minimum and (3) maximum elements of the
*                tolerances for this data
*   tol      - Computed locally based on min bound and scale on rms

      real*8 res_save(nus), tol_par(3), tol

*   flg_jump  - Changed to true if this is a jump and and not a bad
*               data point.  If bad data point then added to edit list

      logical flg_jump

*   dat_type  - Type of data (LC/WL)

      character*(*) dat_type

* LOCAL VARIABLES

*  loc_rms  - Sum of the difference of residuals around the point
*             squared (used for getting local RMS)

      real*8 loc_rms

*  start_loc, end_loc - Start of end point numbers in computing rms
*  num_loc            - Number of local points
*  i                  - Loop counter

      integer*4 start_loc, end_loc, num_loc, i

***** Compute the local RMS

      start_loc = max(1,imax-4)
      end_loc   = min(nus,start_loc+9)

      loc_rms = 0.d0
      num_loc = 0
      do i = start_loc, end_loc-1
         if( i.ne. imax .and. i+1.ne.imax ) then
             loc_rms = loc_rms + (res_save(i)-res_save(i+1))**2
             num_loc = num_loc + 1
         end if
      end do

      if( num_loc.gt.0 ) then
          loc_rms = sqrt(loc_rms/num_loc)
      else
          loc_rms = tol_par(2)
      end if

      tol = min(max(tol_par(2), tol_par(1)*loc_rms), tol_par(3))*
     .      tol_scale

*     See if the jump is large enough to trip the tolerance
      if ( imax.gt.1 ) then
          if( abs(res_save(imax)-res_save(imax-1)).gt.tol ) then 

*            Now see if bad point or a jump.  Since we handle bias flags
*            and finding the correct one-way correctly in version 1.11;
*            change the editing so that it get flagged as a bad point.
*            No editing option; tell program to flag as a bias parameter.
             flg_jump = .true.
          end if 
      end if

****  Thats all
      return
      end
      
CTITLE CHECK_BP

      subroutine check_pb( data_flag_cse, ctol_cse, ns, lv, ep )

      implicit none
      
*     This routine will check if a data point being edited
*     is good and has a bias flag.  If this is the case then
*     bias flag will be pushed to the next good point.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   ns, lv - Cfile number and satellite list entry
*   ep    - Epoch of point edited.

      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), ns, lv, ep
     
* LOCAL VARIABLES
 
*   k     - counter for epochs
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.

      integer*4 k, ch, ltoc
      
*   data_OK, good_bf  -- Logical function to indicate data_OK and
*        good bias flag.
*   pushed -- Logical to indicate that bias flag pushed.

      logical data_OK, good_bf, pushed
      

***** See if this data point contains a good bias flag
      ch = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)      
      if( ch.gt.0 ) then
         if( good_bf(data_flag_cse(ch,ns,ep),0, phs_mask ) ) then

*            This is good bias flag point.  Push the bias forward
             pushed = .false.
             k = ep
             do while ( .not.pushed .and. k.lt.num_ep ) 
                k = k + 1
                ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
                if( ch.gt.0 ) then
*                   See if OK while ignoring the bias flags                 
                    if( data_OK(data_flag_cse(ch,ns,k),0,
     .                          phs_mask) ) then

*                       See if bias flag already present
                        if( good_bf(data_flag_cse(ch,ns,k),0,
     .                          phs_mask) ) then 
                            pushed = .true.
                        else
                            call sbit( data_flag_cse(ch,ns,k),31,1)
C                           write(*,120) cf_codes(ns), prn_list(lv),
C    .                           ep, k
C120                        format(' Pushed Bias Flag for ',a4,' PRN ',
C    .                             I2.2,' from Epoch ',i5,' to ',I5)                         
                            pushed = .true.
                        end if
                    end if    
                 end if
             end do
          end if
       else
          write(*,140) cf_codes(ns), prn_list(lv), ep
 140      format(' ERROR not channel for ',a4,' PRN ',
     .            I2.2,' at Epoch ',i5) 
       end if
       
****   Thats all 
       return
       end 
       
