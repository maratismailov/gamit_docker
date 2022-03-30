CTITLE CLEAN_OW
 
      subroutine clean_ow(pass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )

      implicit none
 
*     This routine will estimate the cycles slips at the remianing
*     bias flags to produce "flat" oww-way residuals.  This will only
*     work well when postfit residuals are used.
 
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
*   trimlen - Length of string
*   ierr    - IOSTAT error on opening the dd report file
*   dd_clean_iter - Number of iteratations through the dd cleaning
*             for a single station/satellite.  Each time data is
*             flagged or a bias flag added we iterate the clean loop.
 
      integer*4 i,j,k, ns, lv, ch, ltoc, first_ep, step_ep,
     .          dd_clean_iter
 
*   data_OK - function returns true if data OK.  phs_mask does
*           - no check bias flags
*   kbit    - Checks if bit is set.
*   first_bias  - Indicates that first bias flag has been
*           - found
*   in_gap  - Indicates that we are in gap (not used)
 
      logical data_OK, kbit, first_bias
c     logical in_gap

***** First order the stations and satellites for cleaning
 
      call get_clean_order('CLEANING')
 
*     Now loop over the one-way data
      do i = 1, num_cfiles
          ns = dd_site_list(i)

*         Loop over the satellites as a do-while loop
*         because we may need to repeat a satellite if there
*         is edited data or bias flags added.
          j = 0
          do while (j.lt.num_sat)
              j = j + 1
              lv = dd_svs_list(j)
              first_bias = .false.
 
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
 
                      if( kbit(data_flag_cse(ch,ns,k),31).or.
     .                    kbit(data_flag_cse(ch,ns,k),32) ) then
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
 
*             Start at the first point so that form_data will 
*             check the quality even though it will not be able
*             to patch the data
              first_ep = k
 
              write(*,120) cf_codes(ns), prn_list(lv), first_ep,
     .                    step_ep, dd_clean_iter
  120         format(' Oneway Set ',a4,' PRN ',i2.2,' Start ep ',i4,
     .                ' Step ',i3,' epochs: Iteration ',i2)
 
              do k = first_ep, num_ep, step_ep
 
*                 see if we have the end of a gap or a bias flag
                  ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! SV observed at this time
                  if( ch.gt.0 ) then
 
*                     See if we have found a good data point after a
*                     bias flag.
                      if( data_OK(data_flag_cse(ch,ns,k),0,phs_mask)
     .                    .and. 
     .                      (kbit(data_flag_cse(ch,ns,k),31).or.
     .                       kbit(data_flag_cse(ch,ns,k),32))    ) then
 
*                         Set the number of unflagged cycle slips to
*                         zero
                          unflg_num = 0

*                         Clear the double difference second site and
*                         satellite so thatwe will know double diffs were
*                         used
                          call set_oneway(ns, lv, k, step_ep,
     .                        L1r_phs_cse, L2r_phs_cse,
     .                        L1r_rng_cse, L2r_rng_cse,
     .                        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                        data_flag_cse, bf_type_cse, params_cse,
     .                        par_flag_cse, pass )

                      end if
 
*                     Apply the current number of cycle slips to
*                     this point
 
                      L1_cyc_cse(ch,ns,k) = L1_cyc_cse(ch,ns,k) +
     .                                    curr_L1_slip
                      L2_cyc_cse(ch,ns,k) = L2_cyc_cse(ch,ns,k) +
     .                                    curr_L2_slip
 
                  end if
              end do

*                         ! Looping over satellites
          end do
*                         ! Looping over stations.
      end do
 
***** Thats all
      return
      end
 
CTITLE SET_ONEWAY  
 
      subroutine set_oneway(ns, lv, ep, step_ep,
     .        L1r_phs_cse, L2r_phs_cse, L1r_rng_cse, L2r_rng_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse, data_flag_cse,
     .        bf_type_cse, params_cse, par_flag_cse, pass )

      implicit none
 
*     This routine will scan ethe one-way data and set the cycle
*     slips so that we residuals and WL are close to zero.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number
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
 
 
      integer*4 ns, lv, ep, step_ep, 
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
      integer*4 ch, ltoc

*   norm_cyc(2,2)   - Normal equations for L1 and L2 cycle
*               - slip estimate
*   b_cyc(2)        - B-vector for estimation of number of cycle
*               - slips 
*   dL1_slip, dL2_slip - number of additional cycles needed
*                 accross a patch
*   chiqual     - Quantitative estimate of quality of cycle
*                 slip fixing.
*   chi_cyc     - Prefit chi**2 for trial and actual 
*                 estimate
 
      real*8 norm_cyc(2,2), b_cyc(2), dL1_slip, dL2_slip,
     .       chiqual, chi_cyc

*    gap_or_bias(2) - String to denoted if a bias or gap being
*                     patched.  (Used for output only)
*    reliable       - Not used in this routine.

      character*4 gap_or_bias(2)
      logical     reliable

      data gap_or_bias / 'GAP ', 'BFLG' /
 
***** Start: Initialize the estimator for the number of cycle slips
 
      call init_cyc_est( norm_cyc, b_cyc, chi_cyc )
 
*     collect the one-way for this site.  If we have L2 ranges first
*     collect a lot and so that we can widelane, narrow lane and LG
*     patch.
 
      ch = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)

 
*     Collect the widelane data.  The data are returned through
*     common in data_lft and data_rgh
      call get_one_way(ns, lv, ep, step_ep, max_wl_ret,
     .        L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, params_cse)
 
*     Add this contribution to the estimate of the number of
*     cycles in slip. Copy the one-way data to the work arrays
*     If we have an L2 range measurement then try to use the
*     wide-lane and narrow lane observables.
      call copy_dtow
      if( lambda(lv,4,ns).eq.1 ) then
          call inc_cyc_est(ns, lv, 0, 0, 'WALA', 'OW', norm_cyc, 
     .                     b_cyc, chi_cyc )
      else
* MOD TAH 990518: Check to see if we have L2
          if( lambda(lv,2,ns).ne.0 ) then
             call inc_cyc_est(ns, lv, 0, 0, 'LALG', 'OW', norm_cyc, 
     .                     b_cyc,  chi_cyc )
          else
             call inc_cyc_est(ns, lv, 0, 0, 'L1', 'OW', norm_cyc, 
     .                     b_cyc, chi_cyc )
          end if
      end if

*     Now solve for the number of cycles in L1 and L2
      call est_dd_cyc( ep, ns, lv, 0, 0, norm_cyc,b_cyc, chi_cyc, 
     .            reliable, chiqual, dL1_slip, dL2_slip )

*     Update the number of cycles to be applied to the data
      if( curr_L1_slip.ne.dL1_slip .or. curr_L2_slip.ne.dL2_slip )
     .write(*,200) ep, cf_codes(ns), prn_list(lv), curr_L1_slip ,
     .            dL1_slip, curr_L2_slip , dL2_slip, chiqual
 200  format(' EPOW  ',i4,' Site ',a4,' PRN ',i2.2, ' L1 from ',f7.1,
     .        ' to ',F7.1,' L2 from ',f7.1,' to ', f7.1,' ChiQual ',
     .        F8.2)
 
      curr_L1_slip = dL1_slip
      curr_L2_slip = dL2_slip
 
****  Thats all
      return
      end
 
