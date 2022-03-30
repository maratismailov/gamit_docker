CTITLE PROC_PHSFIN
  
      subroutine proc_phsfin(opass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )
 
*     This routine will use the clock estimated from the range data
*     to set the initial values of the L1 and L2 cycle parameters.
*     These values are reset for each bias flag encountered.
                              
      implicit none

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   opass    -- Cleaning pass number.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass
 
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
 
*   i,j        - Epoch loop counter
  
      integer*4 i, j
 
*   data_ok     - Logical function returns true if data is good
*   kbit        - Function to check status of bit
*   cyc_updated - Return from est_clk_phs indicating that cycle values
*                 have been updated.  If this is done then we re-estimate
*                 clocks with the new cycles included.
 
      logical cyc_updated 

*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch. 

      real*8 curr_cyc(2,max_gsvs, max_cfiles)

*  Arrays for getting the estimates of phase offsets
*  sums_lc(2, max_gsvs, max_cfiles)  - Current LC one offsets for mean and rms
*  nums_lc(max_gsvs, max_cfiles) - Number of values in computation
*  epoc_lc(max_gsvs, max_cfiles) - Start epoch for current accumulations.  When
*             cycles updated, we need to go back to this epoch.

      real*8 sums_lc(2, max_gsvs, max_cfiles) 
      integer*4 nums_lc(max_gsvs, max_cfiles), 
     .          epoc_lc(max_gsvs, max_cfiles)

*  cyc_set(max_gsvs, max_cfiles) - Indicates that data has already been
*                 found on this one-way
*  site_used(max_cfiles)  - Set true when site used.

      logical cyc_set(max_gsvs, max_cfiles), site_used(max_cfiles)

 
****  set the data mask to be used for editing the phase data.  Two masks
*     are used.  The phs_bias_mask check the bias flags as well as
*     editing, the simple phs_mask just checks to see that the data
*     are OK.
      call set_phs_mask( phs_mask, phs_bias_mask )

****  Check the maxium size of the RMS residuals and remove any station
*     that is too large
      if( opass.gt.pc_start_edit*2 .and. edit_postfit ) then
          call pf_check_rms(opass,data_flag_cse)
      end if

***** Set up the initial use of data flags
      do i = 1, num_cfiles
          site_used(i) = .false.
          do j = 1, num_sat

*             On odd numbered iterations we fit not using
*             data after bias flags until we have enough
*             data to get a "resonably" good estimate of the
*             number of cycles. 
              curr_cyc(1,j,i) = 0.d0
              curr_cyc(2,j,i) = 0.d0
*             Test of odd (-1) or even pass (0)
              if( int(opass/2)*2-opass.lt.  0 ) then
                  cyc_set(j,i) = .false.
              else
                  cyc_set(j,i) = .true. 
              end if

              nums_lc(j,i) = 0
              sums_lc(1,j,i)  = 0.d0
              sums_lc(2,j,i)  = 0.d0
              epoc_lc(j,i) = 1
          end do
      end do
 
*     Set the initial epoch for the clocks
      do i = 1, num_param
         apr_clk_epoch(i) = 0
      end do
 
*     Start looping over all epochs of data
      do i = 1, num_ep
 
*         clear the normal equations and set the apriori variance
*         for the clocks.  (Derived from the range solution)
          call init_norm(i, .true. )
 
*         Get the apriori values of the clocks for this epoch from
*         range solution values
          do j = 1, num_param
              apr_clk_val(j) = params_cse(j,i)
          end do


*         Now get the new estimates at this epoch and sum the
*         range residual rms estimates.  Update the cycles for
*         those sites and channels which we do not know yet.
*         First call accounting for cycles.  
          call est_clk_fin(i, phs_mask, cyc_updated, 
     .                    L1r_phs_cse(1,1,i), L2r_phs_cse(1,1,i),
     .                    L1_cyc_cse(1,1,1), L2_cyc_cse(1,1,1),
     .                    L1r_rng_cse(1,1,i),  L2r_rng_cse(1,1,i),
     .                    ctol_cse(1,1,1), data_flag_cse(1,1,1),
     .                    bf_type_cse(1,1,i), params_cse(1,i), 
     .                    par_flag_cse(1,i), 
     .                    cyc_set, curr_cyc, sums_lc, nums_lc,
     .                    epoc_lc, opass)

      end do

*     Now finish up any biases that have been left unresolved. 
*     Only do this on odd passes where we tried to update the
*     bias parameters.
      if( int(opass/2)*2-opass.ne.0 ) then
         call end_resolve( L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                   data_flag_cse, 
     .                   cyc_set, curr_cyc, sums_lc, nums_lc,
     .                   epoc_lc, opass )
      end if

****  Thats all
      return
      end
 
CTITLE EST_CLK_FIN
 
      subroutine est_clk_fin(ep, mask, cyc_updated,
     .                L1r_phs_cs, L2r_phs_cs, L1_cyc_cse, L2_cyc_cse,
     .                L1r_rng_cs, L2r_rng_cs, ctol_cse, data_flag_cse, 
     .                bf_type_cs, params_cs, par_flag_cs, 
     .                cyc_set, curr_cyc, sums_lc, nums_lc,
     .                epoc_lc, opass )
 
*     Subroutine to estimate the clocks at this epoch using the
*     phase data.  Only data for which the biases are known are used
*     in the estimation.  The number of cycles associated with
*     the other phases for which the number of cycles is not known
*     are resolved after the clocks have been estimated.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cse(num_chan, num_cfiles,num_ep)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cse(num_chan, num_cfiles,num_ep) - Data flags (to see if any
*                 - good.)
*   bf_type_cs(num_chan, num_cfiles) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   par_flag_cs(num_param)         - parameter estimate flag
*   ep      - Epoch number being processed.
*   mask    - Mask to use in deciding if data should be included
*             On first pass it will not use data with bias flag
*             set, on second pass these are ignored.
*   opass   - Indicates which pass through the clock estimation this is:
*             On odd numbered passes, biases are estimated.  On even
*             numbers these are not estimated.  Routine should be called
*             an even number of times.
 
      integer*4 ctol_cse(num_chan, num_cfiles,num_ep),
     .    data_flag_cse(num_chan, num_cfiles,num_ep),
     .    bf_type_cs(num_chan, num_cfiles),
     .    par_flag_cs(num_param), ep, mask, opass

*   cyc_updated - Logical to indicate that cycles values have been
*                 estimated

      logical cyc_updated
 
*   L1_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L1 cycles to be
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L2 cycles to be
*                     added (size as above).
*   L1r_phs_cs(num_chan, num_cfiles)    - L1 phase measurements (L1 cylces)
*   L2r_phs_cs(num_chan, num_cfiles)    - L2 phase measurements (L2 cylces)
*   L1r_rng_cs(num_chan, num_cfiles)    - L1 range measurements (L1 cylces)
*   L2r_rng_cs(num_chan, num_cfiles)    - L2 range measurements (L2 cylces)
*   params_cs(num_param)    - Clock parameter estimates (L1 cycles)
 
      real*8 L1_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L1r_phs_cs(num_chan, num_cfiles),
     .    L2r_phs_cs(num_chan, num_cfiles),
     .    L1r_rng_cs(num_chan, num_cfiles),
     .    L2r_rng_cs(num_chan, num_cfiles), params_cs(num_param)

*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch. 

      real*8 curr_cyc(2,max_gsvs, max_cfiles)

*  Arrays for getting the estimates of phase offsets
*  sums_lc(2, max_gsvs, max_cfiles)  - Current LC one offsets for mean and rms
*  nums_lc(max_gsvs, max_cfiles) - Number of values in computation
*  epoc_lc(max_gsvs, max_cfiles) - Start epoch for current accumulations.  When
*             cycles updated, we need to go back to this epoch.
*  dL1c, dL2c  - currrent change in the number of cycles
*  mean, rms   - mean and rms change

      real*8 sums_lc(2, max_gsvs, max_cfiles), dL1c, dL2c, mean, rms

* MOD TAH 200617 Added frequency ratio (only fr1 needed)
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

      integer*4 nums_lc(max_gsvs, max_cfiles), 
     .          epoc_lc(max_gsvs, max_cfiles)

*  cyc_set(max_gsvs, max_cfiles) - Indicates that data has already been
*                 found on this one-way

      logical cyc_set(max_gsvs, max_cfiles)

* LOCAL VARIABLES
 
 
*   i,j,k       - Loop counters
*   ns          - parameter number of satellite
*   lv          - Satellite number
*   pivot(max_cfiles_max_gsvs)  - Pivot elmenent for inversion
*   par_loc(max_neq)  - array to indicate which parameters
*                 can be estimated.  (Initially set to zero and
*                 set to one when normal equations are incremented).
*  loc_mask     - Set the local mask to ignore bias flags if this
*                 is after the first pass. 
*  ch           - channel number at site
*  ltoc         - Function to return channel given satellite number
*  ss_used(max_cfiles+max_gsvs, max_cfiles+max_gsvs) -- sets which site
*                 and satellites have been used.
 
      integer*4 i,j, ns, lv, pivot(max_neq),
     .    par_loc(max_cfiles+max_gsvs), loc_mask, ch, ltoc,
     .    ss_used(max_cfiles+max_gsvs, max_cfiles+max_gsvs)
 
*   data_ok     - Logical function which is true if data OK
*   good_bf     - Returns true if bias flag on good data point.
 
      logical data_ok, good_bf
 
*   res         - Generic residual values
*   phs_omc     - Function to return range residual
*   data_var    - Data varaince at one site
*   scale(max_neq)  - Scale factors for invert_vis
*   L1_cyc_sav, L2_cyc_sav  - Values of the cycles before we
*                 estimate new values.  These are used to
*                 see if we have updated the values.

*   lc_est       - estimate of LC error expected given WL resolved
*   lc_err       - Difference between mean and lc_est
*   lc_sig       - Sigma of mean difference based on rms and number
*                  of values.
 
      real*8 res, phs_omc, data_var, scale(max_neq),
     .       lc_est, lc_err, lc_sig

      logical kbit

* Push bias flags
      integer*4 tep, tmask, cho
      logical pushed, nodd_bfgood
 
***** Initialize the parameter local set array
       
      do i = 1, num_param
         par_loc(i) = 0
         do j = 1, num_param
            ss_used(i,j) = 0
         end do
      end do
      
      dL1c = 0.d0
      dL2c = 0.d0

*     Make sure that we have the mask set to check the bias flags
*     on the first past and ignore them afterwards.
      loc_mask = mask
*     If below is -1 for odd passes and 0 for even passes.  On
*     odd passes bias flags are checked for quality of data.
      if( int(opass/2)*2-opass.lt.  0 ) then
         call sbit(loc_mask,31,1)
         call sbit(loc_mask,32,1)
      else
         call sbit(loc_mask,31,0)
         call sbit(loc_mask,32,0)
*        Set ignore noDD if past first few iterations
         tmask = loc_mask
         if( opass.ge.8 ) then
            call sbit(loc_mask,23,0)
         end if

      end if

***** Loop over the data incrementing the normal equations
      do i = 1,num_cfiles

          data_var = phs_noise(i)**2

*         Now compute the next iteration on the clock values
          do j = 1, actual_max_chan

****          if pass one then update number of cycles
              lv = ctol_cse(j,i,ep)
	      IF( lv.gt.0 ) THEN
              L1_cyc_cse(j,i,ep) = L1_cyc_cse(j,i,ep) +
     .                             curr_cyc(1,lv,i)
              L2_cyc_cse(j,i,ep) = L2_cyc_cse(j,i,ep) +
     .                             curr_cyc(2,lv,i)

*             Use only data that has no error and for which
*             we already know the bias parameter values
              if( data_ok(data_flag_cse(j,i,ep),0, loc_mask).and.
     .            cyc_set(lv,i)    ) then
 
*                 Compute OminusC
* MOD TAH 040703: Removed resetting the NoDD bit.  (This might be needed
*                 if data is deleted and then restored)?
C                 call sbit(data_flag_cse(j,i,ep),23,0)
                  ns = num_cfiles + lv 
                  res = phs_omc(L1r_phs_cs(j,i),L2r_phs_cs(j,i),
     .                    L1_cyc_cse(j,i,ep), L2_cyc_cse(j,i,ep),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(lv), fL2(lv) )
 
*                 Increment normal equations based on the data quality
                  norm_eq(i,i)   = norm_eq(i,i)   + 1.d0/data_var
                  norm_eq(ns,ns) = norm_eq(ns,ns) + 1.d0/data_var
                  norm_eq(ns,i)  = norm_eq(ns,i)  - 1.d0/data_var
                  norm_eq(i,ns)  = norm_eq(i,ns)  - 1.d0/data_var
                  sol_eq(i)      = sol_eq(i)      + res/data_var
                  sol_eq(ns)     = sol_eq(ns)     - res/data_var
                  par_loc(i) = 1
                  par_loc(ns) = 1
                  ss_used(i,i)   = ss_used(i,i)   + 1
                  ss_used(ns,ns) = ss_used(ns,ns) + 1
                  ss_used(ns,i)  = ss_used(ns,i)  + 1
                  ss_used(i,ns)  = ss_used(i,ns)  + 1
              end if
	      END IF
          end do
      end do

****  Now check to see if we elements that have been used only once
* MOD TAH 040703: Only check the double differences on the 
*     even numbered passed
      if( int(opass/2)*2-opass.eq. 0 ) then

         do i = 1, num_param             
            if( ss_used(i,i).eq.1 ) then
*               Find the other site or satellite with only one observation
                if( i.le.num_cfiles ) then
*                   Must be a site with only one satellite observed.  
*                   Find SV
                    do j = num_cfiles+1, num_param
                       if( ss_used(i,j).eq.1 .and. cyc_set(lv,i) ) then
                           lv = j - num_cfiles
                           ch = ltoc(ctol_cse(1,i,ep),lv,
     .                                       actual_max_chan)
                           cho = ch
                           if( ch.gt.0 ) then
                              nodd_bfgood = good_bf(
     .                                         data_flag_cse(ch,i,ep),
     .                                         0,tmask)
                              call sbit(data_flag_cse(ch,i,ep),23,1)
* MOD TAH 040703: See if need to push a bias flag forward
                              if( nodd_bfgood ) then
                                 pushed = .false.
                                 tep = ep
                                 do while (.not.pushed) 
                                    tep = tep + 1
                                    if( tep.le.num_ep ) then
                                        ch = ltoc(ctol_cse(1,i,tep),lv,
     .                                       actual_max_chan)
                                        if( ch.gt.0 ) then
                                        if ( data_OK(
     .                                       data_flag_cse(ch,i,tep),
     .                                       0,tmask)          ) then
                                             call sbit(
     .                                          data_flag_cse(ch,i,tep),
     .                                          31,1)
                                             pushed = .true.
                                        endif
                                        endif
                                    else
                                        pushed = .true.
                                    endif
                                 end do
                              endif
c                             write(*,996)'PusghNoDD Site ',i,lv,ep,tep,
c    .                            data_flag_cse(cho,i,ep)
c996                          format(a,1x,4i5,1x,o16) 
                           end if          
                       end if
                    end do
                else
*                   More likely case of satellite visible from only one
*                   station.  See which site
                    lv = i - num_cfiles
                    do j = 1, num_cfiles
                       if( ss_used(i,j).eq.1 .and. cyc_set(lv,j) ) then
                           ch = ltoc(ctol_cse(1,j,ep),lv,
     .                                    actual_max_chan)
                           cho = ch
                           if( ch.gt.0 ) then
                              nodd_bfgood = good_bf(
     .                                         data_flag_cse(ch,j,ep),
     .                                         0,tmask)
                              call sbit(data_flag_cse(ch,j,ep),23,1)
* MOD TAH 040703: See if need to push a bias flag forward
                              if( nodd_bfgood ) then
                                 pushed = .false.
                                 tep = ep
                                 do while (.not.pushed) 
                                    tep = tep + 1
                                    if( tep.le.num_ep ) then
                                        ch = ltoc(ctol_cse(1,j,tep),lv,
     .                                       actual_max_chan)
                                        if( ch.gt.0 ) then
                                        if ( data_OK(
     .                                       data_flag_cse(ch,j,tep),
     .                                       0,tmask)          ) then
                                             call sbit(
     .                                          data_flag_cse(ch,j,tep),
     .                                          31,1)
                                             pushed = .true.
                                        endif
                                        endif
                                    else
                                        pushed = .true.
                                    endif
                                 end do
                              endif 
c                             write(*,996)'PusghNoDD SVS  ',j,lv,ep,tep,
c    .                            data_flag_cse(cho,j,ep)

                           endif
                       end if
                    end do
                end if
            end if
         end do
      end if
 
****  Check the parameter flag to see if we have data on a s
*     parameter.  If we no not set its diagonal to one.
      do i = 1, num_param
*                                         ! No data
          if( par_loc(i).eq.0 )  then
              norm_eq(i,i) = 1.d0
*             No data at this time so mark clock as not determined
              call sbit(params_cs(i),1,1)
          else
              call sbit(params_cs(i),1,0)
          end if
          norm_eq(i,i) = norm_eq(i,i)*(1.d0+1.d-6)
      end do
 
*     Now invert the system
      call invert_vis(norm_eq, sol_eq, scale, pivot, num_param,
     .                max_neq,1)
 
*     Save as apriori and save values in parameter array
      do i = 1, num_param
          params_cs(i) = apr_clk_val(i) + sol_eq(i)
          apr_clk_val(i) = params_cs(i)
 
*         if we have data, save this epochs value as the
*         apriori for the next epoch
          if( par_loc(i).eq.1 )  then
              apr_clk_epoch(i) = ep
          end if
      end do      

*     If this is an even pass, we are not estimating biases so
*     just return.
      if( int(opass/2)*2-opass.eq.0 ) RETURN

****  Now if this is the first pass then check to see if
*     have reached a bias flag but not resolved the bias.

      do i = 1,num_cfiles

*         Now check on the values of the bias parameters.     
          do j = 1, actual_max_chan

*             Use only data that has no error and for which
*             we already know the bias parameter values
              lv = ctol_cse(j,i,ep)
	      IF( LV.gt.0 ) THEN
* MOD TAH 200617: Added frequency ratio
              fr1 = fL1u(lv)/fL1(lv)
              if( good_bf(data_flag_cse(j,i,ep),0, mask) ) then

*                 we have reached a bias flag, but have not 
*                 yet set the previous bias flags values.  Force
*                 if to be set now, and apply the numbers of cycles
*                 back through earlier data.
                  if( nums_lc(lv,i).gt.0 .and. .not.cyc_set(lv,i) ) then
                      mean = sums_lc(1,lv,i)/nums_lc(lv,i)
                      rms  = sqrt(abs(sums_lc(2,lv,i)/nums_lc(lv,i) -
     .                            mean**2))
                      if( kbit(status_rep,12) )
     .                write(*,300) cf_codes(i), prn_list(lv),
     .                       epoc_lc(lv,i), mean, rms,
     .                       nums_lc(lv,i)
 300                  format(' For site ',a4,' PRN ',i2.2,
     .                       ' Unresolved bias from epoch ',i4,
     .                       ' Mean, RMS ',2f6.2,' cyles, # ',i6)


*                     Resolve to nearest integer assuming that
*                     the widelane is set correctly.  Use negative
*                     so that mean will be removed.
                      if( opass.le.pc_non_int ) then
                          if( fdd_L2_fact.ne.0 ) then
* MOD TAH 200617: Account for frequency differences
                              dL1c = -nint((mean/(lcf1(lv)+lcf2(lv)))
     .                                                      *fr1)/fr1
                          else
                              dL1c = -nint(mean*fr1)/fr1
                          end if
                      else
                          if( fdd_L2_fact.ne.0 ) then
                              dL1c = -mean/(lcf1(lv)+lcf2(lv))
                          else
                              dL1c = -mean
                          end if
                      end if
* MOD TAH 990518: Check for L1 only data.
                      if( lambda(lv,2,i).ne.0 ) then
                         dL2c = dL1c
                      else
                         dL2c = 0.0d0
                      end if

                      if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Added passing of edit flag
                          call update_cyc(dL1c, dL2c, i, lv, 
     .                         epoc_lc(lv,i), ep, opass,
     .                         L1_cyc_cse(1,1,1),
     .                         L2_cyc_cse(1,1,1),
     .                         ctol_cse(1,1,1),data_flag_cse,'N' )
                          curr_cyc(1,lv,i) = curr_cyc(1,lv,i) + dL1c
                          curr_cyc(2,lv,i) = curr_cyc(2,lv,i) + dL2c
                      end if
                  end if

*                 Clear the summation arrays and set the next bias
*                 flag epoch.
                  sums_lc(1,lv,i) = 0.d0
                  sums_lc(2,lv,i) = 0.d0
                  nums_lc(lv,i)   = 0
                  epoc_lc(lv,i)   = ep
                  cyc_set(lv,i)   = .false.
              end if
	      END IF
          end do
      end do

****  Now check the entries for which we have not yet set the number
*     of cycles.  We accumulate their statistics and once we have more
*     data than the mimimum to fix a bias flag we starting checking to
*     see if we can resolve the number of biases.

      do i = 1,num_cfiles
*         Now check on the values of the bias parameters.     
          do j = 1, actual_max_chan

*             Use only data that has no error and for which
*             we already know the bias parameter values
* ** Use mask passed since this should be set to ignore bias flags.
              lv = ctol_cse(j,i,ep)
	      IF( LV.gt.0 ) THEN
              fr1 = fL1u(lv)/fL1(lv)
              if( data_OK(data_flag_cse(j,i,ep),0, mask).and.
     .            .not. cyc_set(lv,i)    ) then

*                 Compute OminusC
                  ns = num_cfiles + ctol_cse(j,i,ep)
                  res = phs_omc(L1r_phs_cs(j,i),L2r_phs_cs(j,i),
     .                    L1_cyc_cse(j,i,ep), L2_cyc_cse(j,i,ep),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(lv), fL2(lv) )

*                 Add this residual to current statistics
* MOD TAH 980522: Use only data which is part of double difference.
                  if( .not.kbit(data_flag_cse(j,i,ep),23) ) then
                     sums_lc(1,lv,i) = sums_lc(1,lv,i) + res
                     sums_lc(2,lv,i) = sums_lc(2,lv,i) + res**2
                     nums_lc(lv,i)   = nums_lc(lv,i)   + 1   
                  end if

*                 Now see if we can fix
                  if( nums_lc(lv,i).gt.min_good_bias*10 ) then
                      mean = sums_lc(1,lv,i)/nums_lc(lv,i)
                      rms  = sqrt(abs(sums_lc(2,lv,i)/nums_lc(lv,i) -
     .                            mean**2))

****                  See if OK.
                      if( lambda(lv,2,i).ne.0 ) then
                         lc_est = nint((mean/(lcf1(lv)+lcf2(lv)))*
     .                       (lcf1(lv)+lcf2(lv))*fr1)/fr1
                      else
                         lc_est = nint(mean*fr1)/fr1
                      endif
                      lc_err = mean-lc_est
                      lc_sig = sqrt(rms**2/nums_lc(lv,i))
                      if( abs(lc_err) .lt.0.1d0 .and.
     .                    lc_sig      .lt.0.05d0 ) then
                          if( kbit(status_rep,12) )
     .                    write(*,400) cf_codes(i), prn_list(lv),
     .                           epoc_lc(lv,i), mean  , lc_sig,
     .                           nums_lc(lv,i)
 400                      format(' For site ',a4,' PRN ',i2.2,
     .                           ' Resolved bias from epoch ',i4,
     .                           ' Diff +- ',2f6.2,' cyles, # ',i6)

*                         Resolve to nearest integer assuming that
*                         the widelane is set correctly.  Use negative
*                         so that mean will be removed.
                          if( opass.le.pc_non_int ) then
                              dL1c = -nint((mean/(lcf1(lv)+lcf2(lv)))
     .                                                      *fr1)/fr1
                          else   
                              dL1c = -(mean/(lcf1(lv)+lcf2(lv)))
                          end if
* MOD TAH 990518: Check if L2 available.
                          if( lambda(lv,2,i).ne.0 ) then
                             dL2c = dL1c
                          else
                             dL2c = 0.d0
                          end if
                          cyc_set(lv,i) = .true.

*                         Only apply the corrections if they are
*                         non-zero.
                          if( dL1c.ne.0 .or. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data
                              call update_cyc(dL1c, dL2c, i, lv, 
     .                             epoc_lc(lv,i), ep, opass,
     .                             L1_cyc_cse(1,1,1),
     .                             L2_cyc_cse(1,1,1),
     .                             ctol_cse(1,1,1),data_flag_cse, 'N' )
                              curr_cyc(1,lv,i) = curr_cyc(1,lv,i) + dL1c
                              curr_cyc(2,lv,i) = curr_cyc(2,lv,i) + dL2c
                          end if

*                         Clear the summation arrays
                          sums_lc(1,lv,i) = 0.d0
                          sums_lc(2,lv,i) = 0.d0
                          nums_lc(lv,i)   = 0
*                         Set the epoch for applying the next set of cycle
*                         slip changes to -1 so that they will not be applied
*                         until nexy bias flag found.
                          epoc_lc(lv,i)  = -1
                      end if
                  end if
              end if
	      END IF
          end do
      end do

****  Thats all
      return
      end

CTITLE UPDATE_CYC
      
      subroutine update_cyc(dL1c, dL2c, ns, lv, 
     .           ep_start, ep_end, opass, L1_cyc_cse, L2_cyc_cse, 
     .           ctol_cse, data_flag_cse, edit_opt )

*     Routine to apply the current estimate of the number
*     of cycle slips back to the epoch when the bias flag
*     was first encountered.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   ns, lv   - Site and satellite list number
*   ep_start, ep_end - Start and stop epoch numbers.
*   opass     - Cleaning pass number.
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

      integer*4 ctol_cse(num_chan, num_cfiles, num_ep), 
     .          ns, lv, ep_start, ep_end, opass,
     .          data_flag_cse(num_chan, num_cfiles, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   dL1c, dL2c      - Change to the number of cycles.
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    dL1c, dL2c

* edit_opt  -- Editing option:  If set to E then data is edited at the
*     same time the cycles are applied.  This is used in allign_phs to
*     same pieces of data early in the processing

      character*1 edit_opt

* LOCAL VARIABLES

*  i  -- loop counter
*  ch -- Channel number
*  ltoc - Function to return channel number for specific satellite

      integer*4 i, ch, ltoc

*  kbit - Check bit setting

      logical kbit

***** Loop over epochs, finding the entries to change
*     Make sure ep_start is greater than 0.  This case should
*     always be.
      if( ep_start.le.0 ) then
          write(*,100) cf_codes(ns), prn_list(lv), ep_end
 100      format('**WARNING** Attempt to fix non-existant bias at ',
     .           a4,' PRN ',i2.2,' Epoch ',i6)
          RETURN
      end if
      if( kbit(status_rep,12) )
     .write(*,120) cf_codes(ns), prn_list(lv), ep_start, ep_end,
     .             dL1c, dL2c, opass
 120  format('Updating cycles for ',a4,' PRN ',i2.2,' from Epoch ',
     .       i5,' to ',i5,' by dL1/2 ',2f6.2,' cycles, pass ',i2)

      do i = ep_start, ep_end
         ch = ltoc( ctol_cse(1,ns,i),lv, actual_max_chan)  
         if( ch.gt.0 ) then
             L1_cyc_cse(ch,ns,i) = L1_cyc_cse(ch,ns,i) + dL1c
             L2_cyc_cse(ch,ns,i) = L2_cyc_cse(ch,ns,i) + dL2c
* MOD TAH 000310: If the E editing option passed, then remove the data
*            at this time.
             if( edit_opt(1:1).eq.'E' ) then
                 call sbit(data_flag_cse(ch,ns,i),6,1)
             end if
         end if
      end do

***** Thats all
      return
      end

CTITLE GET_MEAN_OWWL
 
      subroutine get_mean_owwl(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )
 
*     This routine will compute the mean oneway widelane 
*     observable by station.  These means are later used to
*     allign the widelanes during postfit residual generation.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
 
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
     .    par_flag_cse(num_param, num_ep)
 
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

* i,j,k,l   -- Loop counters over epochs, sites, satellites and
*              to find channel.
* ch        -- Channel for current satellite.
* nums_wl(2) -- Accumulation on number of points in current span
*              and total
* pv        -- Parameter number for satellite means
* pivot(max_gsvs+max_cfiles) - Pivot elements
* nums_var  -- Sum for total at station.

      integer*4 i,j,k, ch, nums_wl(2), pv,
     .          pivot(max_gsvs+max_cfiles), nums_var

* sums_wl(2,2) -- Accumulation of mean and RMS for current span
*              and total
* sums_var     -- Sum for getting total variance
* WL           -- Estimate of Widelane for current measurement
* mean_wl(max_gsvs, max_cfiles) - Mean wide lane by satellite and
*                 site.
* scale(max_gsvs+max_cfiles)   - Scale for inversion

      real*8 sums_wl(2,2), WL, mean_wl(max_gsvs, max_cfiles),
     .       scale(max_gsvs+max_cfiles) , sums_var

*  min_rms  -  Miniumum WL RMS
*  ewl, dwl -  Estimated mean wide (from double difference combination)
*              and integer version.

      real*8 min_rms, ewl, dwl

*  ref_wl, ref_sv - Reference site for WL double differences, and reference
*              satellite (for each site).
*  ltoc   - Function to return channel for a given satellite entry.

      integer*4 ref_wl, ref_sv, ltoc

*  data_OK  - Logical function that returns true if not bits in
*     data mask match those in the data_flag
*  all_svs  - Indicates that all satellites are seen at site

      logical data_OK, all_svs

****  Start, loop over each station, first by satellite and then
*     average all the satellite results together.


      do j = 1, num_cfiles

*         Clear RMS statistics
          sums_var = 0.d0
          nums_var = 0
      
          do k = 1, num_sat

****         Clear the accumulation arrays.  Each sattelite
*            is accumulated separately but we keep the
*            rms for all satellites at a station.
             call biased_clear( sums_wl, nums_wl )
          

*            Now loop down the epochs
             do i = 1, num_ep
             
*               If the data point is OK, and does not have a bias
*               flag then form the WL and add to accumulation arrays
                ch = ltoc( ctol_cse(1,j,i),k, actual_max_chan)  

*               If the satellite is observed this epoch and data OK
*               then process.  (The two if's below are really one, but
*               due to bounds checking are separated.
                if( ch.gt.0 .and. lambda(k,4,j).ne.0 ) then
                if( data_ok(data_flag_cse(ch,j,i),0,phs_bias_mask)) then      
          

*                   Data is good and does not have a bias flag
                    WL = -(L1r_phs_cse(ch,j,i) +
     .                     L1_cyc_cse(ch,j,i))   +
     .                    (L2r_phs_cse(ch,j,i) + 
     .                     L2_cyc_cse(ch,j,i))    + 
     .                     dfsf(k)*(L1r_rng_cse(ch,j,i)+
     .                           L2r_rng_cse(ch,j,i))

*                   Now accumulate
                    sums_wl(1,1) = sums_wl(1,1) + WL
                    sums_wl(2,1) = sums_wl(2,1) + WL**2
                    nums_wl(1)   = nums_wl(1)   + 1
                else

****                See if there is a bias flag
                    if(data_ok(data_flag_cse(ch,j,i),0, phs_mask)) then

*                      data is good, it just has bias flag on it.  If
*                      we have already accumulated some mean WL estimates
*                      solve and save these values. 
                       call biased_stats( sums_wl, nums_wl, 1.d0)
*                      Now save the current data point in the accumulation
*                      arrays
                       WL = -(L1r_phs_cse(ch,j,i) +
     .                        L1_cyc_cse(ch,j,i))   +
     .                       (L2r_phs_cse(ch,j,i) + 
     .                        L2_cyc_cse(ch,j,i))    + 
     .                        dfsf(k)*(L1r_rng_cse(ch,j,i)+
     .                              L2r_rng_cse(ch,j,i))
                       sums_wl(1,1) = WL
                       sums_wl(2,1) = WL**2
                       nums_wl(1)   = 1
                   end if
                end if
                end if

*                             ! Looping over epochs
            end do

*           Call the biased_stats routine to make sure that we add the
*           last block of data
            call biased_stats( sums_wl, nums_wl, 1.d0) 
                              ! End do looping over satellites
            if( nums_wl(2).gt.0 ) then
                mean_wl(k,j) = sums_wl(1,2)/nums_wl(2)
                sums_var = sums_var + sums_wl(2,2)
                nums_var = nums_var + nums_wl(2)
             else
                mean_wl(k,j) = 0.d0
             end if
         end do  

*        Save the site statistics
         if( nums_var.gt.0 ) then
             WL_rms(j) = sqrt(sums_var/nums_var)
         else
             WL_rms(j) = 0.d0
         end if
         WL_num(j) = nums_var
      end do


****  Now remove the "ambiquiuties" from the mean wide-lanes.  We do
*     this by setting the double differneces to zero.  First find "base"
*     station.
      min_rms = 1.d10
      ref_wl = 0
      do i = 1, num_cfiles
         all_svs = .true.
         do j = 1, num_sat
            if( mean_wl(j,i).eq.0 ) all_svs = .false.
         end do
         if( all_svs .and. WL_rms(i).lt.min_rms ) then

*            OK, sees all satellites and has lowest rms. Save
             min_rms = WL_rms(i)
             ref_wl  = i
         end if
      end do

      if( ref_wl.gt.0 ) then
          write(*,120) cf_codes(ref_wl), WL_rms(ref_wl)
 120      format(' WL Using ',a4,' as Reference site, RMS ',f6.2,
     .           ' cycles')
      else
          do i = 1,num_cfiles
             if( WL_rms(i).lt.min_rms ) then
                 min_rms = WL_rms(i)
                 ref_wl  = i
             end if
          end do
          write(*,140) cf_codes(ref_wl), WL_rms(ref_wl) 
 140      format(' WL No complete reference, using ',a4,
     .            ' RMS ',f6.2, 'cycles')
      end if

****  Now fix up any ambiquities: Based on ref_site, set
*     double differneces to be zero
      do i = 1, num_cfiles
         if( i.ne. ref_wl ) then
             ref_sv = num_sat
             do j = num_sat,1,-1
                if( mean_wl(j,i).ne.0 ) ref_sv = j
             end do
             do j = ref_sv + 1, num_sat
                if( mean_wl(j,i).ne.0 ) then
* MOD TAH 20617: Need to think about mixed satellites here.  Ratio
*                not implemented.
                    ewl = mean_wl(ref_sv,i) -
     .                    (mean_wl(ref_sv,ref_wl)-mean_wl(j,ref_wl))
                    dwl = nint(mean_wl(j,i)-ewl)
                    mean_wl(j,i) = mean_wl(j,i) - dwl
                    if( dwl.ne.0 ) then
                        write(*,160) cf_codes(i), prn_list(j),
     .                               mean_wl(j,i), ewl
 160                    format(' WL Update ',a4,' PRN ',i2.2,
     .                         ' Mean and Est widelane ',2f6.2)
                    end if
                end if
             end do
          end if
      end do

****  Now compute the individual site and satellite contributions
      do i = 1, num_cfiles+num_sat
         sol_eq(i) = 0.d0
         do j = 1, num_cfiles+num_sat
            norm_eq(i,j) = 0.0d0
            if( i.eq.j ) norm_eq(i,j) = 1.d-6
         end do
      end do

*     Now loop over sites and satellites
      do i = 1, num_cfiles
         do j = 1, num_sat
            pv = num_cfiles + j
            if( mean_wl(j,i).ne.0.d0 ) then
                sol_eq(i)  = sol_eq(i) + mean_wl(j,i)
                sol_eq(pv) = sol_eq(pv) + mean_wl(j,i)
                norm_eq(i,i) = norm_eq(i,i) + 1
                norm_eq(pv,pv) = norm_eq(pv,pv) + 1
                norm_eq( i,pv) = norm_eq( i,pv) + 1
                norm_eq(pv, i) = norm_eq(pv, i) + 1
            end if
         end do
      end do

      call invert_vis(norm_eq, sol_eq, scale, pivot, 
     .                num_cfiles+num_sat, max_neq,1)

*     Now save the  results
      do i = 1, num_cfiles
         WL_bias_site(i) = sol_eq(i)
         write(*,220) cf_codes(i), WL_bias_site(i), WL_rms(i), WL_num(i)
 220     format(' WL Stats for ',a4,'  OW Bias ',f6.3,' OW RMS ',
     .                           f6.3,' cycles, Num ',i6)
      end do
      do i = 1, num_sat
         pv = num_cfiles+i
         WL_bias_svs(i) = sol_eq(pv)
         write(*,230) prn_list(i), WL_bias_svs(i)
 230     format(' WL Stats for PRN ',i2.2,' OW Bias ',f6.3)
      end do

****  Now write out the comparison.
      do i = 1, num_cfiles
         do j = 1, num_sat
            write(*,250) cf_codes(i), j, mean_wl(j,i), 
     .                   mean_wl(j,i)-(WL_bias_site(i)+WL_bias_svs(j))
 250        format(' WL Differences ',a4,' PRN ',i2.2,2f9.3)
         end do
      end do


****  Thats all
      return
      end


CTITLE BIASED_CLEAR

      subroutine biased_clear( sums, nums )

*     Rouitne to clear the biased statistics arrays

* PASSED VARIABLES

*  sums(2,2)  -- Summ arrays for mean and rms
*  nums(2)    -- number of data in current and total snpan

      integer*4 nums(2)
      real*8 sums(2,2)

****  Set all values to zero
      nums(1) = 0
      nums(2) = 0
      sums(1,1) = 0.d0
      sums(2,1) = 0.d0
      sums(1,2) = 0.d0
      sums(2,2) = 0.d0

***** Thats all
      return
      end 
     
CTITLE Biased_stats

      subroutine biased_stats( sums, nums, ambi)

*     Routine to accumulate statistics for groups of data
*     that may have biases between the means.  If ambi is
*     passed then the means accumulated are assumed to have
*     this width ambiguity,

* PASSED VARIABLES

*  sums(2,2)  -  Accumulation sums for mean and variance for
*                current accumulation (*,1) and the total (*,1)
*  nums(2)    -  Number of values in current and total accumulation
*  ambi       -  Possible ambiquity in mean.  If this value is
*                zero then means are separatley removed from each
*                segment.

      integer*4 nums(2)
      real*8 sums(2,2), ambi

* LOCAL VARAIBLES
*  mean       - Mean from current data group
*  vars       - Sum residuals squares about the mean for 
*               current group of data.
*  prev_mean  - Mean from previous segments of data
*  dmean      - Change in mean due to ambiquities.

      real*8 mean, vars, prev_mean, dmean

****  See if we have already saved values in the totals part.
      if( nums(2).eq.0 ) then

*         Nothing saved in total so just move sums accross,
*         demeaning the squared sum.
          if( nums(1).gt.0 ) then
*             Make sure we have data.
              mean = sums(1,1)/nums(1)
              vars = abs(sums(2,1)/nums(1) - mean**2) 
              if( ambi.ne.0 ) then
                  dmean = nint(mean/ambi)*ambi
              else
                  dmean = 0.d0
              end if
              sums(1,2) = sums(1,1) - dmean*nums(1)
              sums(2,2) = vars*nums(1)
              nums(2) = nums(1)
          end if
       else

*         We alreasy have data, so add this group in
          if( nums(1).gt.0 ) then 
*             Make sure we have data.
              mean = sums(1,1)/nums(1)
*             see if we need to adjust the mean
              if( ambi.ne.0.d0 ) then
                  prev_mean = sums(1,2)/nums(2)
                  dmean = nint((mean-prev_mean)/ambi)*ambi
              else
                  dmean = 0.d0
              end if
              vars = abs(sums(2,1)/nums(1) - mean**2) 
              sums(1,2) = sums(1,2) + (mean-dmean)*nums(1)
              sums(2,2) = sums(2,2) + vars*nums(1)
              nums(2) = nums(2) + nums(1)

*             Now clear the lead accumulation arrays
              nums(1) = 0
              sums(1,1) = 0.d0
              sums(2,1) = 0.d0
          end if
      end if

****  Thats all
      return
      end

CTITLE END_RESOLVE 

      subroutine end_resolve( L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                data_flag_cse,  
     .                cyc_set, curr_cyc, sums_lc, nums_lc,
     .                epoc_lc, opass )

*     Rotuine to finish up the one-way bias resolution at the end of
*     looping over all epochs.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cs(num_chan, num_cfiles,num_ep)  - Conversion from channel
*                 - number to local list of PRN's

*   opass    - Indicates which pass through the clock estimation this is:
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
 
      integer*4 ctol_cse(num_chan, num_cfiles,num_ep), opass,
     .          data_flag_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L2 cycles to be
*                     added (size as above).
 
      real*8 L1_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles,num_ep)

*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch. 

      real*8 curr_cyc(2,max_gsvs, max_cfiles)

*  Arrays for getting the estimates of phase offsets
*  sums_lc(2, max_gsvs, max_cfiles)  - Current LC one offsets for mean and rms
*  nums_lc(max_gsvs, max_cfiles) - Number of values in computation
*  epoc_lc(max_gsvs, max_cfiles) - Start epoch for current accumulations.  When
*             cycles updated, we need to go back to this epoch.
*  dL1c, dL2c  - currrent change in the number of cycles
*  mean, rms   - mean and rms change

      real*8 sums_lc(2, max_gsvs, max_cfiles), dL1c, dL2c, mean, rms
      integer*4 nums_lc(max_gsvs, max_cfiles), 
     .          epoc_lc(max_gsvs, max_cfiles)
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

*  cyc_set(max_gsvs, max_cfiles) - Indicates that data has already been
*                 found on this one-way

      logical cyc_set(max_gsvs, max_cfiles)
 
* LOCAL VARIABLES
 
 
*   i,j,k       - Loop counters
*   lv          - Satellite number 
*   ch          - Channel number
*   ltoc        - Functio to convert sv list to channel 
 
      integer*4 i, lv, ch, ltoc 
      logical kbit 
 
*     Loop over all stations
      do i = 1,num_cfiles

*         Now scan over all satellites to see if any still have
*         unresolved biases.
          do lv = 1, num_sat
                               
              ch = ltoc( ctol_cse(1,j,i),lv, actual_max_chan)  
              if( nums_lc(lv,i).gt.0 .and. .not.cyc_set(lv,i) ) then
                  mean = sums_lc(1,lv,i)/nums_lc(lv,i)
                  rms  = sqrt(abs(sums_lc(2,lv,i)/nums_lc(lv,i) -
     .                        mean**2))
                  if( kbit(status_rep,12) ) 
     .            write(*,300) cf_codes(i), prn_list(lv),
     .                   epoc_lc(lv,i), mean, rms,
     .                   nums_lc(lv,i), opass
 300              format(' For site ',a4,' PRN ',i2.2,
     .                   ' Final bias from epoch ',i4,
     .                   ' Mean, RMS ',2f6.2,' cyles, # ',i6, 
     .                   ' Pass ',i2)

*                 Resolve to nearest integer assuming that
*                 the widelane is set correctly.  Use negative
*                 so that mean will be removed.
                  if( opass.le.pc_non_int ) then
* MOD TAH 990518: Check to see if we L2 data.
                      if( lambda(lv,2,i).ne.0 ) then
                         fr1 = fL1u(lv)/fL1(lv)
                         dL1c = -nint(mean/(lcf1(lv)+lcf2(lv))*fr1)/fr1
                      else
                         dL1c = -nint(mean)
                      endif
                  else
* MOD TAH 990518: Check to see if we L2 data.
                      if( lambda(lv,2,i).ne.0 ) then
                         dL1c = -(mean/(lcf1(lv)+lcf2(lv)))
                      else
                         dL1c = -mean
                      end if
                  end if
* MOD TAH 990518: Check to see if we L2 data.
                  if( lambda(lv,2,i).ne.0 ) then
                     dL2c = dL1c
                  else
                     dL2c = 0.d0
                  end if
                  if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data
                      call update_cyc(dL1c, dL2c, i, lv, 
     .                     epoc_lc(lv,i), num_ep, opass,
     .                     L1_cyc_cse(1,1,1),
     .                     L2_cyc_cse(1,1,1),
     .                     ctol_cse(1,1,1),data_flag_cse, 'N' )
                      curr_cyc(1,lv,i) = curr_cyc(1,lv,i) + dL1c
                      curr_cyc(2,lv,i) = curr_cyc(2,lv,i) + dL2c
                  end if
              end if
          end do
      end do

***** Thats all
      return
      end

CTITLE COMP_PF_RMS

      subroutine comp_pf_rms(lun, opass, L1r_phs_cse, L2r_phs_cse, 
     .           L1_cyc_cse, L2_cyc_cse,
     .           ctol_cse, data_flag_cse, 
     .           params_cse, par_flag_cse)
 
*     Subroutine to increment the post-fit rms statics by station
*     and satellite, accounting for the bias flags in the data.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   lun  - Unit number to write results to 
*   opass - Iteration number.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), lun, opass
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)

 
* LOCAL VARIABLES
 
* i,j,k,l   -- Loop counters over epochs, sites, satellites and
*              to find channel.              
* ch        -- Channel for current satellite.
* nums_lc(2) -- Accumulation on number of points in current span
*              and total
* pv        -- Parameter number for satellite means
* pivot(max_gsvs+max_cfiles) - Pivot elements
* nums_var  -- Sum for total at station.
* ns        -- Parameter number of satellite clock.
* ltoc      -- Function to return channel number for satellite.
* num_rav, num_avg -- Running number of values in average, and total number
*              of averaged points

      integer*4 i,j,k, ch, nums_lc(2), 
     .          nums_var, ns, ltoc, num_rav, num_avg

* sums_lc(2,2) -- Accumulation of mean and RMS for current span
*              and total
* sums_var     -- Sum for getting total variance 
* sum_rav, sum_avg -- Sum for running average and total average variance
* mean_lc(max_gsvs, max_cfiles) - Mean wide lane by satellite and
*                 site.
* scale(max_gsvs+max_cfiles)   - Scale for inversion
* res          -- Residual of LC.
* phs_omc      -- Routine to phase obs-comp'd.

      real*8 sums_lc(2,2), sums_var, res, phs_omc, sum_rav, sum_avg

* data_OK      -- Logical to indicate good data
* kbit         -- Test single bit (used for No Double Difference
*                 bit).

      logical data_OK, kbit

***** Save the previous iteration values if this after the
*     first iteration
      if( opass.gt.1 ) then
          do j = 1, num_cfiles
             PR_rms(j) = LC_rms(j)
             PR_num(j) = LC_num(j)
          end do
          All_prv = All_rms
      else

*         Initialize the previous RMS values to something large
          do j = 1, num_cfiles
             PR_rms(j) = 100.d0
             PR_num(j) = 0
          end do
          All_rms = 100.d0
      end if
 
****  Loop over all sites and satellites
      do j = 1, num_cfiles

*         Clear RMS statistics
          sums_var = 0.d0
          nums_var = 0

          do k = 1, num_sat

****         Clear the accumulation arrays.  Each sattelite
*            is accumulated separately but we keep the
*            rms for all satellites at a station.
             call biased_clear( sums_lc, nums_lc )

*            Clear statistics for averages
             sum_avg = 0.d0
             sum_rav = 0.d0
             num_avg = 0
             num_rav = 0
 
*            Now loop down the epochs
             do i = 1, num_ep
             
*               If the data point is OK, and does not have a bias
*               flag then form the WL and add to accumulation arrays
                res = 0.d0  
                ch = ltoc( ctol_cse(1,j,i),k, actual_max_chan)  

*               If the satellite is observed this epoch and data OK
*               then process.  (The two if's below are really one, but
*               due to bounds checking are separated.  Check we really
*               have double differences (residual will be zero if we
*               don't)
                if( ch.gt.0 ) then
		
                if( data_ok(data_flag_cse(ch,j,i),0,phs_bias_mask) .and.
     .              .not.kbit(data_flag_cse(ch,j,i),23)  ) then      

*                   Data is good and does not have a bias flag  
                    ls = ctol_cse(1,j,i)
                    ns = num_cfiles + k
                    res = phs_omc(L1r_phs_cse(ch,j,i),
     .                    L2r_phs_cse(ch,j,i),
     .                    L1_cyc_cse(ch,j,i), L2_cyc_cse(ch,j,i),
     .                    params_cse(j,i), params_cse(ns,i),
     .                    fL1(k),fL2(k) ) 
     
*                   Now accumulate, but only of non-zero
                    sums_lc(1,1) = sums_lc(1,1) + res
                    sums_lc(2,1) = sums_lc(2,1) + res**2
                    nums_lc(1)   = nums_lc(1)   + 1

*                   Now accumulate the average values
                    num_rav = num_rav + 1
                    sum_rav = sum_rav + res
                    if( num_rav.eq.25 ) then
*                       OK, we have enough data to find the average
*                       and accumulate the statistics
                        num_avg = num_avg + 1
                        sum_avg = sum_avg + (sum_rav/num_rav)**2
*                       Now reset the values
                        num_rav = 0
                        sum_rav = 0.d0
                    end if
                else

****                See if there is a bias flag
                    if(data_ok(data_flag_cse(ch,j,i),0, phs_mask).and.
     .                 .not.kbit(data_flag_cse(ch,j,i),23)  ) then      
*                      data is good, it just has bias flag on it.  If
*                      we have already accumulated some mean WL estimates
*                      solve and save these values. 
                       call biased_stats( sums_lc, nums_lc, 100.d0)

*                      Now save the current data point in the accumulation
*                      arrays
                       sums_lc(1,1) = res
                       sums_lc(2,1) = res**2
                       nums_lc(1)   = 1

                   end if
                end if
                end if

*                             ! Looping over epochs
            end do

*           Call the biased_stats routine to make sure that we add the
*           last block of data
            call biased_stats( sums_lc, nums_lc, 100.d0) 
                              ! End do looping over satellites
            if( nums_lc(2).gt.0 ) then
                lc_svs_bias(k,j) = sums_lc(1,2)/nums_lc(2)
                lc_svs_rms(k,j)  = sqrt(sums_lc(2,2)/nums_lc(2))
                lc_svs_num(k,j)  = nums_lc(2) 
                sums_var = sums_var + sums_lc(2,2)
                nums_var = nums_var + nums_lc(2)
             else
                lc_svs_bias(k,j) = 0.d0
                lc_svs_rms(k,j)  = 0.d0
                lc_svs_num(k,j)  = 0
             end if

*            Now save the RMS of the averaged values
             if( num_avg.gt.0 ) then
                 lc_svs_ams(k,j) = sqrt(sum_avg/num_avg)
               else
                 lc_svs_ams(k,j) = 0
             end if
             lc_svs_anm(k,j) = num_avg
         end do  

*        Save the site statistics
         if( nums_var.gt.0 ) then
             LC_rms(j) = sqrt(sums_var/nums_var)
         else
             LC_rms(j) = 0.d0
         end if
         LC_num(j) = nums_var

*        Now do the average for all of the satellites at this station
         num_avg = 0
         sum_avg = 0 
         do k = 1, num_sat
            sum_avg = sum_avg + lc_svs_ams(k,j)**2*lc_svs_anm(k,j)
            num_avg = num_avg + lc_svs_anm(k,j)
         end do
         if( num_avg.gt.0 ) then
             LC_AMS(j) = sqrt(sum_avg/num_avg)
         else
             LC_AMS(j) = 0.d0
         endif
         LC_anm(j) = num_avg

      end do
      
*     Now print the results

      call report_rms( lun, opass )

****  Thats all
      return
      end

CTITLE REPORT_RMS

      subroutine report_rms( lun, opass )

*     Routine to report the RMS scatters but station and
*     satellite at station.

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*  lun  - Unit number for outpit
*  opass - Output pass number

      integer*4 lun, opass

* LOCAL VARIABLES

* i,j,k  - Loop counters
      integer*4 j,k

* cyc_to_mm -- Conversion from cyc to mm for output. 
* all_svs_rms(max_gsvs)   -- RMS of all sites

      real*8 cyc_to_mm, all_svs_rms(max_gsvs)

* all_svs_num(max_gsvs)  -- Number of data in rmss
      integer*4 all_svs_num(max_gsvs), all_num

****  Get conversion from cyc to mm for output
* RWK 150202: Don't distinguish among GLONASS frequencies for this statistic
      cyc_to_mm = (vel_light/fClk)*1000.d0

***** Start output
      write(lun, 300) opass
 300  format(/,'ONE-WAY POSTFIT RESIDUAL STATISTICS: Pass ',i3,/,
     .         '---------------------------------------------')  
* MOD TAH 200618: Updated 32I to 50I to allow for 45 Beidou satellites
      write(lun,320) opass,(prn_list(k),k=1,num_sat)
 320  format('RMS by site and satellite (mm): Pass ',i3,/,
     .       'RMS  IT Site   All',50I5.2)
                                                            
      do j = 1, num_cfiles
         write(lun,330) opass, cf_codes(j), LC_rms(j)*cyc_to_mm, 
     .               (nint(LC_svs_RMS(k,j)*cyc_to_mm), k=1,num_sat)
 330     format('RMS ',i3,1x,a4,f6.1,50I5)

* MOD TAH 980522: Update the phase noise based on RMS
         if( opass.gt.2 ) then
            phs_noise(j) = ( LC_rms(j)/all_rms + 2)*all_rms/2
         end if
      end do

*     Accumulate up overall statistics
      all_rms = 0.d0
      all_num = 0
      do k = 1, num_sat
         all_svs_rms(k) = 0.d0
         all_svs_num(k) = 0.d0
      end do

*     Now loop over all sites and satellites
      do j = 1, num_cfiles
         do k = 1, num_sat
            all_svs_rms(k) = all_svs_rms(k) + 
     .                       LC_svs_RMS(k,j)**2*LC_svs_NUM(k,j)
            all_svs_num(k) = all_svs_num(k) + LC_svs_NUM(k,j)
            all_rms        = all_rms + 
     .                       LC_svs_RMS(k,j)**2*LC_svs_NUM(k,j)
            all_num        = all_num + LC_svs_NUM(k,j)
         end do
      end do

*     Now complete the statistics computations
      do k = 1, num_sat
         if( all_svs_num(k).gt.0 ) then 
             all_svs_rms(k) = sqrt(all_svs_rms(k)/all_svs_num(k))
         end if
      end do
      if( all_num.gt.0 ) then
         all_rms = sqrt(all_rms/all_num)
      end if

****  Now write the over all statitstics line
* MOD TAH 030317: Increase all output to 0.1 mm.
* MOD MAF 210629: Updated 32(1x,I4) to 50(1x,I4) to allow for 45 Beidou satellites
      write(lun,335) opass, 'ALL ', all_rms*cyc_to_mm, 
     .              (nint(all_svs_RMS(k)*cyc_to_mm*10), k=1,num_sat)
 335     format('RMS ',i3,1x,a4,f6.1,50(1x,I4))

*     Write out number of values in estimates.
* MOD TAH 200618: Updated 32I to 50I to allow for 45 Beidou satellites
      write(lun,340) opass,(prn_list(k),k=1,num_sat)
 340  format(/,'Number of data by site and satellite: Pass ',i3,/,
     .       'NUM  IT Site   All',50I5.2)
 
      do j = 1, num_cfiles
         write(lun,350) opass, cf_codes(j), LC_num(j), 
     .               (LC_svs_num(k,j), k=1,num_sat)
 350     format('NUM ',i3,1x,a4,i6,50I5)
      end do

*     Now write biases.
C     write(lun,360) opass,(prn_list(k),k=1,num_sat)
C360  format(/,'Average OW residuals by site and satellite (mm):',
C    .       ' Pass ',i3,/,
C    .       'AVG  IT Site   RMS',32I4.2)
C
C     do j = 1, num_cfiles
C        write(lun,370) opass, cf_codes(j), LC_rms(j)*cyc_to_mm, 
C    .               (nint(LC_svs_bias(k,j)*cyc_to_mm), k=1,num_sat)
C370     format('AVG ',i3,1x,a4,f6.1,32I4)
C     end do

****  Now write the Average value RMS scatters
* MOD TAH 200618: Updated 32I to 50I to allow for 45 Beidou satellites
      write(lun,420) opass,(prn_list(k),k=1,num_sat)
 420  format(/,'RMS of 25-point averages by site and satellite (mm):',
     .       ' Pass ',i3,/,
     .       'AMS  IT Site   All Ratio ',50I5.2)
 
      do j = 1, num_cfiles
         write(lun,430) opass, cf_codes(j), LC_ams(j)*cyc_to_mm, 
     .               (LC_rms(j)/(LC_ams(j)+1.d-5)),
     .               (nint(LC_svs_AMS(k,j)*cyc_to_mm), k=1,num_sat)
 430     format('AMS ',i3,1x,a4,f6.1,1x,f6.2, 50I5)
      end do

****  Thats all
      return
      end

CTITLE RESET_CYC
      
      subroutine reset_cyc(  L1_cyc_cse, L2_cyc_cse)
 
*     Subroutine to reset the cycles offset in the one way
*     data back to the nearest integer ambiquitity value.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)

 
* LOCAL VARIABLES
      
*  i,j,k  -- Loop variables
      integer*4 i,j,k
      
*  dL1, dL2  -- Number of L1 and L2 cycles
*  dL1_msw, dL1_lsw -- Large and small parts of number of cycles
*               (Only needed for some very old Astech data).     
*  dL2_msw, dL2_lsw -- Large and small parts of number of cycles
*               (Only needed for some very old Astech data). 

      real*8 dL1, dL2,  dL1_msw, dL1_lsw,  dL2_msw, dL2_lsw 
          
      do i = 1, num_ep
         do j = 1, num_cfiles
             do k = 1, actual_max_chan
             
*               Set ambiquity back to integer number of lambda'a
                dL1 =   L1_cyc_cse(k,j,i)                           
                if( abs(dL1).gt.2.d0**31) then
                   dL1_msw = nint(dL1/1.d4)*1.d4 
                   dL1_lsw = dL1 - dL1_msw
                   L1_cyc_cse(k,j,i)  = dL1_msw + 
     .                        (1.d0*nint(dL1_lsw*lambda(k,1,j)))/
     .                                           lambda(k,1,j)
                else
*                  Do computations with standard nint function
                   L1_cyc_cse(k,j,i) = 
     .                        (1.d0*nint(dL1*abs(lambda(k,1,j))))/
     .                                       abs(lambda(k,1,j))
                end if
                dL2 =   L2_cyc_cse(k,j,i) 

*               Make sure we use the L2 factor if the cycles are
*               very large.
                if( abs(dL2).gt.2.d0**30) then
                   dL2_msw = nint(dL2/1.d4)*1.d4 
                   dL2_lsw = dL2 - dL2_msw
                   L2_cyc_cse(k,j,i)  = 
     .                         (1.d0*nint(dL2_lsw*
     .                               lambda(k,2,j)))/
     .                               lambda(k,2,j) 
     .                                     + dL2_msw
                else
*                  Do computations with standard nint function.  
*                  Make sure to multiply by 1.d0 so that comp
*                  is double precision floating point.
                   L2_cyc_cse(k,j,i) = 
     .                       (1.d0*nint(dL2*lambda(k,2,j)))/
     .                                    lambda(k,2,j)  
 
                end if
            end do
         end do
      end do 
 
*     Thats all
      return
      end
                    
CTITLE RM_OW_MEAN  
 
      subroutine rm_ow_mean(opass,
     .                L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                ctol_cse, data_flag_cse, 
     .                params_cse, par_flag_cse )
 
*     Subroutine to remove the mean offsets in each section of
*     data after the phase clocks have been estimated.

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cse(num_chan, num_cfiles,num_ep)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   opass   - Indicates which pass through the clock estimation this is:
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   opass   - number of iterations so far
 
      integer*4 ctol_cse(num_chan, num_cfiles,num_ep),
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass

*   L1_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L1 cycles to be
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L2 cycles to be
*                     added (size as above).
*   L1r_phs_cse(num_chan, num_cfiles)    - L1 phase measurements (L1 cylces)
*   L2r_phs_cse(num_chan, num_cfiles)    - L2 phase measurements (L2 cylces)
*   params_cse(num_param)    - Clock parameter estimates (L1 cycles)
 
      real*8 L1_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)


* LOCAL VARIABLES
 
 
*   i,j,k       - Loop counters
*   ns          - parameter number of satellite
*   ch          - Channel number of sattelite (0 if not observed)
*   ep          - epoch number counter
 
      integer*4 i,j, ep, ch, ns 
 
*   data_ok     - Logical function which is true if data OK
*   good_bf     - Returns true if bias flag on good data point.
 
      logical data_ok, good_bf
 
*  Arrays for getting the estimates of phase offsets
*  sum_lc  - Sum of lc residuals
*  num_lc  - Number of values in sime
*  epoc_lc - Epoch from which the sum started
*  dL1c, dL2c  - currrent change in the number of cycles
*  mean, rms   - mean and rms change
*  curr_cyc(2) - Current offset being applied to cycles.
*  ltoc        - Function to return channel number

      real*8 sum_lc, dL1c, dL2c, mean, curr_cyc(2)
      integer*4 num_lc,  epoc_lc, ltoc

*   res         - Generic residual values
*   phs_omc     - Function to return range residual
*   data_var    - Data varaince at one site

      real*8 res, phs_omc

 
***** Loop over the data computing the mean and removing it
*     from the number of cycles added to the phase.
      do i = 1,num_cfiles
         do j = 1, num_sat

              sum_lc = 0.d0
              num_lc = 0
              epoc_lc = 1
              curr_cyc(1) = 0.d0
              curr_cyc(2) = 0.d0

              do ep = 1, num_ep

*                Now check on the values of the bias parameters.     
                 ch = ltoc( ctol_cse(1,i,ep),j, actual_max_chan)  
 
*                Do we have this satellite?
                 if( ch.gt.0 ) then

*                    See if we have hit a bias flag
                     if( good_bf(data_flag_cse(ch,i,ep),
     .                                        0, phs_mask) ) then

*                         we have reached a bias flag, but have not 
*                         yet set the previous bias flags values.  Force
*                         if to be set now, and apply the numbers of cycles
*                         back through earlier data.
                          if( num_lc.gt.0 ) then
                              mean = sum_lc/num_lc

*                             Resolve to nearest integer assuming that
*                             the widelane is set correctly.  Use negative
*                             so that mean will be removed.
                              if( opass.le.pc_non_int ) then
                                  dL1c = -nint(mean/(lcf1(j)+lcf2(j)))
                              else
                                  dL1c = -mean/(lcf1(j)+lcf2(j))*
     .                              pc_over_shoot
                              end if
                              dL2c = dL1c
                              if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data
                                  call update_cyc(dL1c, dL2c, i, j, 
     .                                 epoc_lc, ep-1, opass,
     .                                 L1_cyc_cse(1,1,1),
     .                                 L2_cyc_cse(1,1,1),
     .                                 ctol_cse,data_flag_cse, 'N' )
                              end if
                          end if

*                         Clear the summation arrays and set the next bias
*                         flag epoch.
                          sum_lc = 0.d0
                          num_lc = 0
                          epoc_lc = ep
                      end if

*                     Apply the current phase offset and accumlate the
*                     residuals
                      ns = num_cfiles + j
                          
                      res = phs_omc(L1r_phs_cse(ch,i,ep),
     .                              L2r_phs_cse(ch,i,ep),
     .                              L1_cyc_cse(ch,i,ep),
     .                              L2_cyc_cse(ch,i,ep),
     .                              params_cse(i,ep),
     .                              params_cse(ns,ep),
     .                              fL1(j), fL2(j) )

                      if( data_OK(data_flag_cse(ch,i,ep),0,
     .                                          phs_mask) ) then

*                         Add this residual to current statistics
                          sum_lc = sum_lc + res
                          num_lc = num_lc + 1   
                      end if
*                                ! Satellite observed this epoch
                  end if
              end do

*             Now we have looped over all epochs, apply any mean offset
*             for last segment of data
              if( num_lc.gt.0 ) then
                  mean = sum_lc/num_lc

*                 Resolve to nearest integer assuming that
*                 the widelane is set correctly.  Use negative
*                 so that mean will be removed.
                  ls = ctol_cse(ch,i,1)
                  if( opass.le.pc_non_int ) then
                      dL1c = -nint(mean/(lcf1(j)+lcf2(j)))
                  else
                      dL1c = -mean/(lcf1(j)+lcf2(j))*pc_over_shoot
                  end if
                  dL2c = dL1c
                  if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data
                      call update_cyc(dL1c, dL2c, i, j, 
     .                     epoc_lc, num_ep, opass,
     .                     L1_cyc_cse(1,1,1),
     .                     L2_cyc_cse(1,1,1),
     .                     ctol_cse(1,1,1),data_flag_cse,'N' )
                  end if
              end if
*                        ! Looping over satellites
          end do
*                        ! Looping over sites
      end do

****  Thats all
      return
      end

CTITLE CHECK_CONVD

      subroutine check_convd(opass, converged ) 

*     Routine to check if the postfit RMS one-way residuals
*     have converged.

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*  opass - Current pass number
      integer*4 opass

*  converged -  Set true if RMS changes have converged

      logical converged

* LOCAL VARIABLES

*  i    - Loop counter over stations
*  ns_max  - Site with maximum change
*  ns_neg  - Site with negative change

      integer*4 i, ns_max, ns_neg

*  drms - Ratio of change in rms (old/new-1).  Compard with
*     limit pc_convd.
*  max_drms - Largest change in rms
*  neg_drms - Largest negative change is rms

      real*8 drms, max_drms, neg_drms

***   Loop over all stations
      converged = .true.
      max_drms  = -1.d6
      neg_drms  =  1.d6
      do i = 1, num_cfiles

*        Only check those sites that seem to be present.
         if( LC_rms(i).gt.0.d0 ) then
             drms = (PR_rms(i)/LC_rms(i) - 1.d0)
*            If rms increased (drms<0) or change larger than
*            tolerance then we have not converged.
             if( drms .lt. -pc_convd .or. drms.gt.pc_convd ) then
                 converged = .false.
             end if
             if( drms.gt.max_drms ) then
                 max_drms = drms
                 ns_max   = i
             end if
             if( drms.lt.neg_drms ) then
                 neg_drms = drms
                 ns_neg   = i
             end if
         end if
      end do

*     Make sure that we have reached editing loop
      if( edit_postfit .and. opass.le. pc_start_edit+2 ) 
     .                                      converged = .false.

*     Check to see if total RMS is diverging
C     if( ALL_rms.gt. All_prv .and. opass.gt.2 ) converged = .true.
* MOD TAH 980907: Changed the divergence test.  If we start diverging
*     change the reduce the overshoot.
      if( ALL_rms.gt. All_prv .and. opass.gt.2 ) then
         pc_over_shoot = sqrt(pc_over_shoot)
      end if
      
*     Write results
      write(*,120) opass, cf_codes(ns_max), max_drms*100, 
     .             cf_codes(ns_neg), neg_drms*100, pc_over_shoot,
     .             converged
 120  format('+CHECK converged: Iteratation ',i3,' Max DRMS at ',a4,
     .       ' Change ',f5.1,'%. Min DRMS at ',a4,' Change ',
     .       f5.1,'%. Over-shoot ',f5.2,' Converged? ',L1)

****  Thats all
      return
      end
                 
CTITLE PF_EDIT  
 
      subroutine pf_edit(opass,
     .                L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                ctol_cse, data_flag_cse, 
     .                params_cse, par_flag_cse, azel_cse )
 
*     Subroutine to do postfit editing of one-way phase data.
*     The procedure used is:
*       (a) Compare the rms of each satellite with the average
*           for the station. All those satellites whose rms is
*           greater than <pf_svs_ratio> times the average for the
*           station are looked at.  Data can be restored on other
*           satellites.
*       (b) Any data that is <pf_nsig> outlier based on the satation
*           rms is deleted, any previously edited point is checked
*           to see if now satisfies the condition. 

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cse(num_chan, num_cfiles,num_ep)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   opass   - Indicates which pass through the clock estimation this is:
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   opass   - number of iterations so far
 
      integer*4 ctol_cse(num_chan, num_cfiles,num_ep),
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass

*   L1_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L1 cycles to be
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L2 cycles to be
*                     added (size as above).
*   L1r_phs_cse(num_chan, num_cfiles)    - L1 phase measurements (L1 cylces)
*   L2r_phs_cse(num_chan, num_cfiles)    - L2 phase measurements (L2 cylces)
*   params_cse(num_param)    - Clock parameter estimates (L1 cycles)
 
      real*8 L1_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
*                      angles to be written to phs_res_root files.
*                      (Real*4, radians)

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)  


* LOCAL VARIABLES
 
 
*   i,j,k       - Loop counters
*   ns          - parameter number of satellite
*   ch          - Channel number of sattelite (0 if not observed)
*   ep          - epoch number counter
*   res_mask    - Mask to be used in checking if we can we
*                 restore data.  These are:
*                 Bit  6  -- Bias flags were set too close
*                 Bit 13  -- Postfit edit (from previous iteration)
*   last_good_ep -  Epoch of the last good data point.  Used to see if
*                 can remove a bias flag from a restored point which
*                 has a bias flag set.

      integer*4 i,j, ep, ch, ns, res_mask, last_good_ep
 
*   data_ok     - Logical function which is true if data OK
*   good_bf     - Returns true if bias flag on good data point.
 
      logical data_ok, good_bf
 
*  ltoc        - Function to return channel number

      integer*4 ltoc

*   res         - Generic residual values
*   phs_omc     - Function to return range residual
*   sig         - Expected sigma based on elevation angle dependent
*                 noise model

      real*8 res, phs_omc, sig

*   push_bias   - Set true if we remove a data point with a 
*                 biad flag and we need to push it forward.

      logical push_bias

***** set the restoring mask
      res_mask = phs_mask
* MOD TAH 060109: Added bit 5, too close to end as well.
* MDO TAH 061017: Removed non-checking bit 5.  Current code can generate
*     lots of bias flags at end of strings
C      call sbit(res_mask,  5, 0)
      call sbit(res_mask,  6, 0)
      call sbit(res_mask, 13, 0)
      call sbit(res_mask, 20, 0)
*     Set to check for no double difference residual
      call sbit(res_mask, 23, 1)
      
***** Loop over the data computing the mean and removing it
*     from the number of cycles added to the phase.
      do i = 1,num_cfiles
* MOD TAH 020731: Only do postfit editing if the lc_rms is greater than
*        zero (suggesting that there is data)
         IF( LC_rms(i).gt.0 ) THEN
         do j = 1, num_sat

*           Scan down the epochs for this satellite
            push_bias = .false.
****        Set last_good_ep well back in time so that first good data
*           will be well after it.
            last_good_ep = -2*num_ep 
            do ep = 1, num_ep

*              Now check on the values of the bias parameters.     
               ch = ltoc( ctol_cse(1,i,ep),j, actual_max_chan)  
 
*              Do we have this satellite?
               if( ch.gt.0 ) then

*                  See if we are pushing a bias flag and this 
*                  data point is OK.
                   if( push_bias .and.
     .                 data_OK(data_flag_cse(ch,i,ep),
     .                              0,phs_mask)) then                      
                       call sbit(data_flag_cse(ch,i,ep),31,1)
                       push_bias = .false.
                   end if

*                  Compute the residual and see if we should
*                  restore or remove.
                   ns = num_cfiles + j
                   res = phs_omc(L1r_phs_cse(ch,i,ep),
     .                           L2r_phs_cse(ch,i,ep),
     .                           L1_cyc_cse(ch,i,ep),
     .                           L2_cyc_cse(ch,i,ep),
     .                           params_cse(i,ep),
     .                           params_cse(ns,ep),
     .                           fL1(j), fL2(j) )
                   sig = sqrt(pf_zdep(1,i)+pf_zdep(2,i)/
     .                        sin(azel_cse(2,ch,i,ep))**2)

*****              Now see what we should do.
C                  if( abs(res/lc_rms(i)).lt. pf_nsig .and.
                   if( abs(res/sig).lt. pf_nsig .and.
     .                 abs(res).lt. 2*pf_maxres .and.
     .                 LC_rms(i).lt. pf_maxres   ) then
                   
*                      Data looks good, see we are not using it
*                      and that it is restorable. 
                       if( .not.data_OK(data_flag_cse(ch,i,ep),
     .                                  0,phs_mask) .and.
     .                      data_OK(data_flag_cse(ch,i,ep),
     .                                  0,res_mask) ) then

*                           If its been a long time since a good data
*                           point, say we should push a bias flag
* MOD TAH 060109: Removed because of better check below
C                           if( (ep-last_good_ep)*sampling .ge.
C    .                           dchi2_max_sep ) push_bias = .true.    
C
****                        Yes, we can restore it.  Turn off any
*                           bias flag that it may have had
* MOD TAH 060109:           Added bit 5 as well.
* MOD TAH 061017: Removed resetting bit 5.  We can get as string of
*                 bias flags very close together with current code
*                 (Bit 5 is bias too close to end of data).
C                            call sbit(data_flag_cse(ch,i,ep), 5,0)
                            call sbit(data_flag_cse(ch,i,ep), 6,0)
* MOD TAH 080528: Do not reset range error bit if we are using lc_autcln
                            if( .not. resolve_wl )
     .                      call sbit(data_flag_cse(ch,i,ep),13,0)
                            call sbit(data_flag_cse(ch,i,ep),20,0)
* MOD TAH 060109: See if we need to add a bias flag at this time.
*                           Can be need if data deleted and later
*                           restored.
                            call check_pushbf(ep,i,j,
     .                           L1_cyc_cse, L2_cyc_cse,
     .                           ctol_cse, data_flag_cse, push_bias)
 
                            if( .not.push_bias ) then
                               call sbit(data_flag_cse(ch,i,ep),31,0)
                               call sbit(data_flag_cse(ch,i,ep),32,0)
                            else
                               call sbit(data_flag_cse(ch,i,ep),31,1)
                               push_bias = .false.
                            end if
                            if( 1.eq.2 )
     .                      write(*,220) ep, cf_codes(i), 
     .                                   prn_list(j), res, sig,
     .                                   opass, push_bias
 220                        format('PF_EDIT: Restoring epoch ',i5,
     .                             1x,a4,' PRN ',i2.2,' LC Res. ',
     .                             f8.3,' RMS ',f6.3,' Cycles,',
     .                             ' Pass ',i3,' PBF ',L1)
                       end if

*                      Save this as the current last good epoch.                       
                       if( data_OK(data_flag_cse(ch,i,ep),
     .                                       0,phs_mask) ) then
                           last_good_ep = ep
                       end if
C                  else if( abs(res/lc_rms(i)).gt. pf_nsig ) then
                   else if( abs(res/sig).gt. pf_nsig ) then
*                      Data is bad.  Flag the data.  Check to see
*                      if it has a bias flag.
                       if( good_bf(data_flag_cse(ch,i,ep),
     .                             0,phs_mask)          ) then
                           push_bias = .true.
                       end if

*                      Restart the counter on the number of good 
*                      observation                       

*                      Only report editing the point if it
*                      is not already flaged.
                       if( data_OK(data_flag_cse(ch,i,ep),
     .                             0,phs_mask)          ) then
                          call sbit(data_flag_cse(ch,i,ep),20,1)
                          if( 1.eq.2 )
     .                    write(*,240) ep, cf_codes(i), 
     .                                 prn_list(j), res, sig,
     .                                 opass, push_bias
 240                      format('PF_EDIT: Deleting  epoch ',i5,
     .                           1x,a4,' PRN ',i2.2,' LC Res. ',
     .                           f8.3,' RMS ',f6.3,' Cycles,',
     .                           ' Pass ',i3,' PBF ',L1)
                      end if
                   end if
C                  write(*,9999) ep, i, j, push_bias, 
C    .                         data_flag_cse(ch,i,ep)
C9999              format('PF: ',3i5,2x,L1,2x,O14)
*                     ! Satellite observed this epoch
               end if
*                       ! Looping over epochs
            end do
*                       ! Looping over satellites.
         end do
         ENDIF          ! If LC_RMS is greater than zero.
*                       ! Looping over sites
      end do

****  Thats all
      return
      end

CTITLE COMP_ELEV_RMS

      subroutine comp_elev_rms(lun, opass, L1r_phs_cse, L2r_phs_cse, 
     .           L1_cyc_cse, L2_cyc_cse,
     .           ctol_cse, data_flag_cse, 
     .           params_cse, par_flag_cse, azel_cse, out)
 
*     Subroutine to fit an elevation angle dependent model to the
*     one-way station residuals.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
*   lun  - Unit number to write results to 
*   opass - Iteration number.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), lun, opass
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)

 
*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
*                      angles to be written to phs_res_root files.
*                      (Real*4, radians)

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep) 

      character*(*) out  ! Output option (REP is write the results to the 
                         ! summary file). 

* LOCAL VARIABLES
 
* i,j,k,l   -- Loop counters over epochs, sites, satellites and
*              to find channel.
* ch        -- Channel for current satellite.
* nums_lc(18) -- Accumulation on number of points in each elevation
*              bin           
* lv        -- Satellite number
* ns        -- Parameter number of satellite clock.
* bin       -- Bin for current elevation

      integer*4 i,j,k, ch, nums_lc(91), ns, bin

* sums_lc(91) -- Accumulation of RMS for each elevation
*               bin (Note: we do not remove mean in RMS comps).
* res          -- Residual of LC.
* phs_omc      -- Routine to phase obs-comp'd.
* eldeg        -- Elevation in degrees
* cyc_to_mm    -- Conversion from L1 cycles to mm
* norm_atm(3)  -- Normal eqn for sig**2 = a**2 + b**2/sin(el)**2
* batm(2)      -- Solution vector
* det          -- Determinate of norm_atm
* zpart        -- Partial for 1/sin(el)**2
* zdep(2)      -- A and B coefficients in above model 

      real*8 sums_lc(91), res, phs_omc, eldeg, cyc_to_mm, 
     .       norm_atm(3), batm(2), det, zpart, zdep(2)

* azdeg, nadeg -- Azimuth and nadir angle
      real*8 azdeg, nadeg

* data_OK      -- Logical to indicate good data
* kbit         -- Test single bit (used for No Double Difference
*                 bit).

      logical data_OK, kbit

****  Loop over all sites and satellites
* MOD RWK: Ignore differences in GLONASS frequencies in the cyc_to_mm concersion
      cyc_to_mm = (vel_light/fClk)*1000.d0
      if( out(1:3).eq.'REP' ) then
         write(lun,150) (k*5,(k+1)*5, k = 0,17)
 150     format(/,'Elevation angle dependent RMS statistics.',
     .         'MODEL: RMS^2 = A^2 + B^2/(sin(elv))^2',/,
     .         'ATELV Site    A     B   ',18(I2,'-',I2.2,1x))
      else
         write( * ,150) (k*5,(k+1)*5, k = 0,17)
      endif

      do j = 1, num_cfiles

*         Clear RMS statistics
          do k = 1, 18
             sums_lc(k) = 0.d0
             nums_lc(k) = 0
          end do
      
*         Now loop down the epochs
          do i = 1, num_ep
             
             do ch = 1, actual_max_chan
                               
                lv = ctol_cse(ch,j,i) 
*               If the satellite is observed this epoch and data OK
*               then process.  (The two if's below are really one, but
*               due to bounds checking are separated.  Check we really
*               have double differences (residual will be zero if we
*               don't)
                if( data_ok(data_flag_cse(ch,j,i),0,phs_mask) .and.
     .              .not.kbit(data_flag_cse(ch,j,i),23)  ) then      
          
*                   Data is good and does not have a bias flag
                    ns = num_cfiles + ctol_cse(ch,j,i)
                    res = phs_omc(L1r_phs_cse(ch,j,i),
     .                    L2r_phs_cse(ch,j,i),
     .                    L1_cyc_cse(ch,j,i), L2_cyc_cse(ch,j,i),
     .                    params_cse(j,i), params_cse(ns,i),
     .                    fL1(lv), fL2(lv) )
     
*                   Now accumulate, but only of non-zero
                    eldeg = azel_cse(2,ch,j,i)*180/pi
                    if( eldeg.gt.0.d0 .and. eldeg.lt.90.d0 ) then
                        bin = eldeg/5.d0 + 1
                        sums_lc(bin) = sums_lc(bin) + res**2
                        nums_lc(bin)   = nums_lc(bin)   + 1
                    end if

                end if

*                             ! Looping over channels 
             end do
*                             ! Looping over epochs
          end do     

*         No finish the computation of the statistics
          do k = 1,3
             norm_atm(k) = 0.d0
          end do
          do k = 1,2
             batm(k) = 0.d0
          end do

          do k = 1, 18
             if( nums_lc(k).gt.0 ) then
                sums_lc(k) = sqrt(sums_lc(k)/nums_lc(k))*cyc_to_mm
             end if

             zpart = 1.d0/sin((k*5-2.5d0)*pi/180.d0)**2
*            Accumulate the normals weighted by the number
*            of data points.  Do not weight by number of data
             if( nums_lc(k).gt.0 ) then
                 norm_atm(1) = norm_atm(1) + 1        ! nums_lc(k)
                 norm_atm(2) = norm_atm(2) + zpart    ! *nums_lc(k)
                 norm_atm(3) = norm_atm(3) + zpart**2 ! *nums_lc(k)
                 batm(1) = batm(1) + sums_lc(k)**2    ! *nums_lc(k)
                 batm(2) = batm(2) + zpart*sums_lc(k)**2 ! *nums_lc(k)
             end if
          end do

*         Now compute deterimate and solve the equations accounting
*         for that both zdep(1) and (2) need to be positive and that
*         we have little or no data.
          det = norm_atm(1)*norm_atm(3) - norm_atm(2)**2
          if( det.gt.0 ) then
             zdep(1) = (batm(1)*norm_atm(3) - batm(2)*norm_atm(2))/det
             zdep(2) = (batm(2)*norm_atm(1) - batm(1)*norm_atm(2))/det
*            If the mean is less than zero; set it to 1 mm and use
*            elevation angle dependence.
             if( zdep(1).lt.0.d0 ) then
* MOD TAH 030830: Make the model match at zenith and 12.5 degrees
                 zdep(1) = (zdep(1)+zdep(2))/2
                 batm(2) = batm(2) - norm_atm(2)*zdep(1)
                 zdep(2) = batm(2)/norm_atm(3)
             end if
*            If elevation term is zero, then just use a constant value.
             if ( zdep(2).lt.0.d0 ) then
                 zdep(1) = batm(1)/norm_atm(1)
                 zdep(2) = 0.d0
             end if
          else
             if( norm_atm(1).gt.0 ) then
                 zdep(1) = batm(1)/norm_atm(1)
                 zdep(2) = 0.d0
             else
                 zdep(1) = 10.d0
                 zdep(2) = 0.d0
             end if
          end if

*         Final check to make sure a non-zero value is given in the
*         n-file
          if( zdep(1).le.0.01d0 ) zdep(1) = 10.d0 

*         Save the coefficinets
          pf_zdep(1,j) = zdep(1)/cyc_to_mm**2
          pf_zdep(2,j) = zdep(2)/cyc_to_mm**2

*         Now write out the values
          if( out(1:3).eq.'REP' ) then
             write(lun, 300 ) cf_codes(j), sqrt(zdep(1)), sqrt(zdep(2)),
     .                    (sums_lc(k), k = 1,18)
 300         format('ATELV ',a4,1x,f5.1,1x,f5.1,2x,18(F5.1,1x))
          else
             write( * , 300 ) cf_codes(j), sqrt(zdep(1)), sqrt(zdep(2)),
     .                    (sums_lc(k), k = 1,18)
          endif
      end do

***** Now loop over data with a 1-deg mean value calculation
* MOD RWK 150202: Ignore differences in GLONASS frequencies in the cyc_to_mm conversion
      cyc_to_mm = (vel_light/fClk)*1000.d0
      if( out(1:3).eq.'REP' ) then
         write(lun,450) (k,k=0,90)
 450     format(/,'Elevation angle dependent mean phase residual (mm)',
     .         /,'ELMEAN Site ',91(2x,i2.2,2x))
      else
         write( * ,450) (k,k=0,90)
      endif

      do j = 1, num_cfiles

*         Clear RMS statistics
          do k = 1, 91
             sums_lc(k) = 0.d0
             nums_lc(k) = 0
          end do
      
*         Now loop down the epochs
          do i = 1, num_ep
             
             do ch = 1, actual_max_chan
                                         
                lv = ctol_cse(ch,j,i)

*               If the satellite is observed this epoch and data OK
*               then process.  (The two if's below are really one, but
*               due to bounds checking are separated.  Check we really
*               have double differences (residual will be zero if we
*               don't)
                if( data_ok(data_flag_cse(ch,j,i),0,phs_mask) .and.
     .              .not.kbit(data_flag_cse(ch,j,i),23)  ) then      
          
*                   Data is good and does not have a bias flag
                    ns = num_cfiles + ctol_cse(ch,j,i)
                    res = phs_omc(L1r_phs_cse(ch,j,i),
     .                    L2r_phs_cse(ch,j,i),
     .                    L1_cyc_cse(ch,j,i), L2_cyc_cse(ch,j,i),
     .                    params_cse(j,i), params_cse(ns,i),
     .                    fL1(lv), fL2(lv) )
     
*                   Now accumulate, but only of non-zero
                    eldeg = azel_cse(2,ch,j,i)*180/pi
                    if( eldeg.gt.0.d0 .and. eldeg.lt.90.d0 ) then
                        bin = int(eldeg)+1
                        sums_lc(bin) = sums_lc(bin) + res
                        nums_lc(bin)   = nums_lc(bin)   + 1
                    end if

                end if

*                             ! Looping over channels 
             end do
*                             ! Looping over epochs
          end do     

*****     Finish up mean value calculation
          do k = 1, 91
             if( nums_lc(k).gt.10 ) then
                sums_lc(k) = (sums_lc(k)/nums_lc(k))*cyc_to_mm
             else
                sums_lc(k) = 99.9d0
             end if
          end do

*         Now write out the values
          if( out(1:3).eq.'REP' ) then
             write(lun, 480 ) cf_codes(j), (sums_lc(k), k = 1,91)
 480         format('ELMEAN ',a4,91(1x,F5.1))
          else
             write( * , 480 )  cf_codes(j), (sums_lc(k), k = 1,91)
          endif
      end do


***** If ALLMEAN status has been requested then 
*     Now loop over data with a 1-deg mean value calcualtion
      if( kbit(status_rep,15) .or. 1.eq.1 ) then
          print *,'ALLMEAN STAT ',kbit(status_rep,15)
          if( out(1:3).eq.'REP' ) then
             write(lun,520) (k,k=0,355,5)
 520         format(/,'Azimmuth dependent mean phase residual (mm)',
     .            /,'AZMEAN Site ',72(2x,i3.3,1x))
          else
            write( * ,520) (k,k=0,355,5)
          endif

          do j = 1, num_cfiles

*            Clear RMS statistics
             do k = 1, 72
                sums_lc(k) = 0.d0
                nums_lc(k) = 0
             end do
      
*            Now loop down the epochs
             do i = 1, num_ep
                
                do ch = 1, actual_max_chan
                                         
                   lv = ctol_cse(ch,j,i) 

*                  If the satellite is observed this epoch and data OK
*                  then process.  (The two if's below are really one, but
*                  due to bounds checking are separated.  Check we really
*                  have double differences (residual will be zero if we
*                  don't)
                   if( data_ok(data_flag_cse(ch,j,i),0,phs_mask) .and.
     .                 .not.kbit(data_flag_cse(ch,j,i),23)  ) then      
             
*                      Data is good and does not have a bias flag
                       ns = num_cfiles + ctol_cse(ch,j,i)
                       res = phs_omc(L1r_phs_cse(ch,j,i),
     .                       L2r_phs_cse(ch,j,i),
     .                       L1_cyc_cse(ch,j,i), L2_cyc_cse(ch,j,i),
     .                       params_cse(j,i), params_cse(ns,i),
     .                       fL1(lv), fL2(lv) )
     
*                      Now accumulate, but only of non-zero
                       azdeg = azel_cse(1,ch,j,i)*180/pi
                       eldeg = azel_cse(2,ch,j,i)*180/pi
                       if( azdeg.gt.0.d0 .and. azdeg.lt.360.0 .and.
     .                     eldeg.lt.35.d0 ) then
                           bin = int(azdeg/5)+1
                           sums_lc(bin) = sums_lc(bin) + res
                           nums_lc(bin)   = nums_lc(bin)   + 1
                       end if

                   end if

*                                ! Looping over channels 
                end do
*                                ! Looping over epochs
             end do     

*****        Finish up mean value calculation
             do k = 1, 72
                if( nums_lc(k).gt.10 ) then
                   sums_lc(k) = (sums_lc(k)/nums_lc(k))*cyc_to_mm
                else
                   sums_lc(k) = 99.9d0
                end if
             end do

*            Now write out the values
             if( out(1:3).eq.'REP' ) then
                write(lun, 530 ) cf_codes(j), (sums_lc(k), k = 1,72)
 530            format('AZMEAN ',a4,72(1x,F5.1))
             else
                write( * , 530 )  cf_codes(j), (sums_lc(k), k = 1,72)
             endif
          end do

*******************************************************************
****      No accumulate averaged NADIR angle residuals by satellite
*         Output values in 0.2 degree bins (70 in all)
          if( out(1:3).eq.'REP' ) then
             write(lun,550) (k*0.2,k=0,69)
 550         format(/,'Nadir dependent mean phase residual by',
     .                ' satellite (mm)',
     .            /,'NAMEAN SVS  ',70(1x,f4.1,1x))
          else
            write( * ,550) (k*0.2,k=0,69)
          endif

          do k = 1, num_sat

*            Clear RMS statistics
             do i = 1, 70
                sums_lc(i) = 0.d0
                nums_lc(i) = 0
             end do
      
*            Now loop down the epochs
             do i = 1, num_ep
                
                do j = 1, num_cfiles

*                  If the satellite is observed this epoch and data OK
*                  then process.  (The two if's below are really one, but
*                  due to bounds checking are separated.  Check we really
*                  have double differences (residual will be zero if we
*                  don't)
                   ch = ltoc( ctol_cse(1,j,i),k, actual_max_chan)  
                   lv = ctol_cse(1,j,i) 
                   IF ( CH.gt.0 ) THEN

                   if( data_ok(data_flag_cse(ch,j,i),0,phs_mask) .and.
     .                 .not.kbit(data_flag_cse(ch,j,i),23)  ) then      
             
*                      Data is good and does not have a bias flag
                       ns = num_cfiles + ctol_cse(ch,j,i)
                       res = phs_omc(L1r_phs_cse(ch,j,i),
     .                       L2r_phs_cse(ch,j,i),
     .                       L1_cyc_cse(ch,j,i), L2_cyc_cse(ch,j,i),
     .                       params_cse(j,i), params_cse(ns,i),
     .                       fL1(lv), fL2(lv) )
     
*                      Now accumulate, but only of non-zero
                       nadeg = asin(6378.d0*cos(azel_cse(2,ch,j,i))/
     .                             26378.d0)*180/pi
                       if( nadeg.gt.0.d0 .and. nadeg.lt.14.0d0 ) then
                           bin = int(nadeg*5.0)+1
                           sums_lc(bin) = sums_lc(bin) + res
                           nums_lc(bin)   = nums_lc(bin)   + 1
                       end if

                   end if
                   ENDIF

*                                ! Looping over channels 
                end do
*                                ! Looping over epochs
             end do     

*****        Finish up mean value calculation
             do i = 1, 70
                if( nums_lc(i).gt.10 ) then
                   sums_lc(i) = (sums_lc(i)/nums_lc(i))*cyc_to_mm
                else
                   sums_lc(i) = 99.9d0
                end if
             end do

*            Now write out the values
             if( out(1:3).eq.'REP' ) then
* MOD TAH 200628: Add GNSS type
                write(lun, 560 ) sv_gnss, prn_list(k), 
     .                                   (sums_lc(i), i = 1,70)
 560            format('NAMEAN ',a1,i2.2,1x,70(1x,F5.1))
             else
                write( * , 560 ) sv_gnss, prn_list(k), 
     .                                   (sums_lc(i), i = 1,70)
             endif
          end do
      ENDIF

****  Thats all
      return
      end

CTITLE PF_CHECK_RMS

      subroutine pf_check_rms(opass,data_flag_cse)

*     Routine to check the RMS scatters from each station and
*     remove any station for which the value is too large

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*   opass -- Current iteration 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .          opass

* LOCAL VARIABLES
* i, k -- Loop counters
* ep   -- epoch number
* trimlen -- Length of string

      integer*4 i,k, ep, trimlen

*  worst_rms  -- Worst RMS among the sites 

      real*8  worst_rms 

*  worst_site -- Site with worst RMS
      integer*4 worst_site

* kbit  -- Check bit status

      logical kbit

* message -- Report_stat message

      character*80 message

****  Loop over the RMS seeing which ones are too large
      worst_rms = 0.d0
      worst_site = 0

      do i = 1, num_cfiles
         if( LC_RMS(i).gt.pf_max_rms .and.
     .       LC_num(i).gt.0 ) then 
             if( LC_RMS(i).gt.worst_rms ) then
                worst_rms = LC_RMS(i)
                worst_site = i
             end if
         endif
      end do

*     Now remove the worst_rms site if one is selected
      if( worst_site.gt.0 ) then
         i = worst_site

*         OK, too large kill all the data so that this site does
*         effect the others in the solution
          do ep = 1, num_ep
             do k = 1, actual_max_chan
                if( .not.kbit(data_flag_cse(k,i,ep),30) ) then
                    call sbit(data_flag_cse(k,i,ep),20,1)
                end if
             end do
          end do
                                           
* MOD RWK 150202: Ignore differences in GLONASS frequencies in this calculation
          write(message,100) cf_codes(i), 
     .                LC_RMS(i)*(vel_light/fClk)*1000.d0,
     .                LC_num(i),  pf_max_rms* (vel_light/fClk)*1000.d0
 100      format('Removing ',a4,' Postfit RMS ',F12.1,
     .           ' mm too large. Num ',i6,' Limit ',F5.1,' mm')
          write(*,'(a)') message(1:trimlen(message))
          write(uns,'(a)') message(1:trimlen(message))
          call report_stat('status','autcln','pf_check_rms',' ',
     .              message,0)

      end if

****  That all
      return
      end

CTITLE CHECK_PUSHBF

      subroutine check_pushbf(ep,ns,js,
     .            L1_cyc_cse, L2_cyc_cse,
     .            ctol_cse, data_flag_cse, push_bias)

*     Rouitine to check if we need to add a bias flag
*     when a data point is restored.  The check is done
*     by scanning back throug L1_cyc_cse, L2_cyc_cse to
*     see if values change before the last valid bias flag.
*     If there is a change, then a new bias flag is neeed.

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cse(num_chan, num_cfiles,num_ep)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   ep  -- Epocj number
*   ns  -- Site number
*   js  -- satellite number
 
      integer*4 ctol_cse(num_chan, num_cfiles,num_ep),
     .    data_flag_cse(num_chan, num_cfiles, num_ep), ep, ns, js

*   L1_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L1 cycles to be
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cse(num_chan, num_cfiles,num_ep) - Number of L2 cycles to be
*                     added (size as above).
 
      real*8 L1_cyc_cse(num_chan, num_cfiles,num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles,num_ep)

*   push_bias   - Set true if we need to push the bais flag to this point.

      logical push_bias

* LOCAL VARIABLES
*   ch          - Channel number of sattelite (0 if not observed)
*   k     - epoch loop counter
*   ltoc  - Function to return channel number for salellite number

      integer*4 ch, k, ltoc

*   L1_cyc_ep, L2_cyc_ep -- Number of cycles at current ep

      real*8 L1_cyc_ep, L2_cyc_ep

*   good_bf     - Returns true if bias flag on good data point.
*   done  - Logical set true when valid bias flag found or we
*           run out of epochs. 
      logical good_bf, done

*
*     Get the number of cycles for current data point
      ch = ltoc( ctol_cse(1,ns,ep),js, actual_max_chan)  
* 
      L1_cyc_ep = L1_cyc_cse(ch,ns,ep)
      L2_cyc_ep = L1_cyc_cse(ch,ns,ep)

*     Now start ruuning backwards to find the last valid bias flags
      k = ep
      done = .false.
      push_bias = .true.
      do while ( .not.done )
          k = k - 1
          if( k.gt.0 ) then
*             Get channel
              ch = ltoc( ctol_cse(1,ns,k),js, actual_max_chan) 
              if( ch.gt.0 ) then
*                See if we have a valid bias flag.  If we do
*                see if number of cycles are the same.  If they 
*                are then we don't need to push bias flag, but if
*                they differ then we do.
                 if( good_bf(data_flag_cse(ch,ns,k),
     .                             0,phs_mask)  ) then
*                    See if cycles are the same
                     if( L1_cyc_ep.eq.L1_cyc_cse(ch,ns,k) .and.
     .                   L2_cyc_ep.eq.L2_cyc_cse(ch,ns,k) ) then
                         push_bias = .false.
                         done = .true.
                     else
                         push_bias = .true. 
                         done = .true.
                     endif
                  end if
              endif
          else
              done = .true.
          endif
      enddo

****  Thats all
      return
      end

CTITLE REPORT_WLSTAT

      subroutine report_wlstat( un )

*     Report quality of MW-WL

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

      integer*4 un   ! Output unit

* LOCAL

      integer*4 i

*     Now save the  results
      write(uns, 100 ) 
 100  format(/,'MW-WL Statistics by site and satellite (cycles)',/,
     .         'MW-WL',1x, 4('Site   Mean   RMS    Num   '))
      write(un,140) (cf_codes(i), WL_bias_site(i), WL_rms(i),
     .                WL_num(i), i=1,num_cfiles)
 140  format(50('MW-WL',1x,4(a4,1x,f6.3,1x,f6.3,1x,I6,2x),/))
      write(un,220)
 220  format('MW-WL',1x,8('PRN    Mean   '))
      write(un,240) (prn_list(i), WL_bias_svs(i),i=1,num_sat)
 240  format(20('MW-WL',1x,8('PRN',i2.2,1x,f6.3,2x),/))
      write(un,'(1x)')

      return
      end
