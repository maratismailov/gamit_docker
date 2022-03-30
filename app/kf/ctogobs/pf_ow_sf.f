CTITLE pf_ow_sf 
 
      subroutine pf_ow_sf(pass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )

      implicit none
 
*     This routine is the one-way scan and fix routine to be used after
*     phase clocks have been used.  Strictly should be with postfit
*     residuals but this is not forced on the user.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   pass            - Indicates which pass this is (for iteration of 
*                     cleaning.  Not used at this point.
 
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
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
*   first_ep    - Epoch number of first data point after the first
*           - bias flag
*   step_ep - Number of epochs between each data point for this
*           - site.
 
      integer*4 i,j,k,  ns, lv, ch, ltoc, first_ep, step_ep
 
*   kbit    - Checks if bit is set.
*   first_bias  - Indicates that first bias flag has been
*           - found
*   good_bf   - Function to indicate bias flag on good data.
 
      logical kbit, first_bias, good_bf

*     Now loop over the one-way data.  This pass is only done on one-way
*     data so just use the station order.
      do i = 1, num_cfiles
          ns = i

*         Loop over the satellites 
          do j = 1, num_sat
              lv = j

              curr_L1_slip = 0.d0
              curr_L2_slip = 0.d0

****          Scan up to the first bias flag (always marked) and
*             step at the data sampling interval for this site
              k = 0
              first_bias = .false.
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
 
*             start at the first point so that form_data will 
*             check the quality even though it will not be able
*             to patch the data
              first_ep = k
 
              write(*,120) cf_codes(ns), prn_list(lv), first_ep,
     .                    step_ep
  120         format(' One-way Scan ',a4,' PRN ',i2.2,' Start ep ',i4,
     .                ' Step ',i3,' epochs')

*             Start at point after first bias flag. 
              do k = first_ep+step_ep, num_ep, step_ep
 
*                 see if we have the end of a gap or a bias flag
                  ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                                     ! SV observed at this time
                  if( ch.gt.0 ) then
 
*                     See if we have found a good data point after a
*                     gap or with a bias flag.
                      if( good_bf(data_flag_cse(ch,ns,k),
     .                           0,phs_mask)           ) then
 
                          call ow_scan_fix(ns, lv, k, step_ep,
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
*                         ! Looping over epochs
              end do
*                         ! Looping over satellites
          end do
*                         ! Looping over stations.
      end do
 
***** Thats all
      return
      end
 
CTITLE OW_SCAN_FIX   
 
      subroutine ow_scan_fix(ns, lv, ep, step_ep,
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
*     i         - Loop countr
*     gbi       - Indicates if Gap or Bias that is being
*                 patched.
      integer*4 ch, ltoc, gbi

*   reliable - Logical to say if patch is reliable

      logical reliable

 
*   norm_cyc(2,2)   - Normal equations for L1 and L2 cycle
*               - slip estimate
*   b_cyc(2)        - B-vector for estimation of number of cycle
*               - slips 
*   dL1_slip, dL2_slip - number of additional cycles needed
*                 accross a patch
*   chiqual     - Quantitative estimate of quality of cycle
*                 slip fixing.
*   chi_cyc, chi_trl  - Prefit chi**2 for trial and actual 
*                 estimate
 
      real*8 norm_cyc(2,2), b_cyc(2), dL1_slip, dL2_slip,
     .       chiqual,  chi_cyc

*    gap_or_bias(2) - String to denoted if a bias or gap being
*                     patched.  (Used for output only)
*    small_one_bg   - Logical that is true if the gap is small
*                     enough for allow one bias/gap to be used.

      character*4 gap_or_bias(2)

      data gap_or_bias / 'GAP ', 'BFLG' /
 
***** Start: Initialize the estimator for the number of cycle slips
 
      call init_cyc_est( norm_cyc, b_cyc, chi_cyc )
      reliable = .false.
 
*     collect the one-way for this site.  If we have L2 ranges first
*     collect a lot and so that we can widelane, narrow lane and LG
*     patch.
 
      ch = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)

****  Now collect the double difference data with this site and
*     satellite.  Again the data are returned in the data_lft and
*     data_rgh arrays.
      call get_one_way(ns, lv, ep, step_ep, 2*max_dd_ret,
     .        L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, params_cse)


*     Try to fix the one-way cycle slip.  Use WL if it is available.  Check the range
*     first and if too large do not try to patch (maximum separation given by
*     dchi2_max_sep in seconds)
      if(  ep_rgh(1)-ep_lft(1).gt.0 .and. 
     .    (ep_rgh(1)-ep_lft(1))*sampling.le. dchi2_max_sep  ) then
          call copy_dtow
          if( lambda(lv,4,ns).eq.1 ) then
               call inc_cyc_est(ns,lv, 0, 0,'LCLGWL','OW',norm_cyc, 
     .                          b_cyc, chi_cyc )
          else if( lambda(lv,2,ns).ne.0 ) then
              call inc_cyc_est(ns, lv, 0, 0,'LCLG', 'OW', norm_cyc, 
     .                         b_cyc, chi_cyc )
          else
              call inc_cyc_est(ns, lv, 0, 0,'L1', 'OW', norm_cyc, 
     .                         b_cyc, chi_cyc )
          end if

*         Now solve for the number of cycles in L1 and L2

          call est_dd_cyc( ep, ns, lv, 0, 0, norm_cyc,b_cyc, chi_cyc, 
     .                reliable, chiqual, dL1_slip, dL2_slip )
* MOD TAH 041202: Debug
c          if( ns.eq.1 .and. lv.eq.16 .and.   
c     .       (ep.eq.1137 .or. ep.eq.1186 .or. ep.eq.1258 .or.
c     .        ep.eq.1312 ) ) then
c              print *,'PS_OW_SF ',ns,lv, ep, chi_cyc, dl1_slip, 
c     .                dl2_slip
c              if( abs(dl1_slip).gt.1.d3 ) dL1_slip = 0
c              if( abs(dl2_slip).gt.1.d3 ) dL2_slip = 0
c          end if
          if( abs(dl1_slip).gt.1.d3 .or. abs(dl2_slip).gt.1.d3 ) then
               write(*,150) ep,ns,lv,  dL1_slip, dL2_slip
 150           format('PF_OW_SF: Slip correction to big ',
     .             'EP: ',i5,2i4,2e20.7)
               dL1_slip = 0
               dL2_slip = 0
          end if


*         All these patches are Bias flags so set, gbi = 2 (=1 is for
*         gaps).
          gbi = 2
     
*         Update the number of cycles to be applied to the data
          write(*,200) ep, cf_codes(ns), prn_list(lv), curr_L1_slip ,
     .                dL1_slip, curr_L2_slip , dL2_slip, reliable,
     .                chiqual, gap_or_bias(gbi)
 200      format(' Epoch ',i4,' Site ',a4,' PRN ',i2.2, ' L1 from ',f8.2
     .          ,' to ',F8.2,' L2 from ',f8.2,' to ', f8.2,' Reliable? '
     .          , L2, F8.2,1x,a4,' PF_OW_SF')
 
          curr_L1_slip = dL1_slip
          curr_L2_slip = dL2_slip
 
*         If the estimate of the number of cycles is considered reliable
*         then remove the bias flag.  If it is not reliable then add a
*         flag (in case this is a gap)
 
          if( reliable ) then
              call sbit( data_flag_cse(ch,ns, ep), 31, 0 )
              call sbit( data_flag_cse(ch,ns, ep), 32, 0 )
*             Reset the bias flag type (since we have now removed it).
              bf_type_cse(ch,ns, ep) = 0
          end if
      end if
      
****  Thats all
      return
      end
 
