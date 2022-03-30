CTITLE PROCESS_PHS
 
      subroutine process_phs(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse, azel_cse, cpass  )

      implicit none
 
*     This routine will use the clock estimated from the range data
*     to set the initial values of the L1 and L2 cycle parameters.
*     These values are reset for each bias flag encountered.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
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
*   cpass   - Cleaning pass number.  On pass 1 we do not update
*             the clock estimates.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), cpass
 
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

*  azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
*                      angles to be written to phs_res_root files.
*                      (Real*4, radians)

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)  
 

* LOCAL VARIABLES
 
*   i,j        - Epoch loop counter
 
      integer*4 i, j, k
 
*   data_ok     - Logical function returns true if data is good
*   cyc_set(max_gsvs, max_cfiles) - Indicates that data has already been
*                 found on this one-way
*   kbit        - Function to check status of bit
*   cyc_updated - Return from est_clk_phs indicating that cycle values
*                 have been updated.  If this is done then we re-estimate
*                 clocks with the new cycles included.
 
      logical data_ok, cyc_set(max_gsvs, max_cfiles), kbit, 
     .        cyc_updated 

*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch. 

      real*8 curr_cyc(2,max_gsvs, max_cfiles)

 
****  set the data mask to be used for edting the phase data.  Two masks
*     are used.  The phs_bias_mask check the bias flags as well as
*     editing, the simple phs_mask just checks to see that the data
*     are OK.
      call set_phs_mask( phs_mask, phs_bias_mask )
 
***** Set up the initial use of data flags
      do i = 1, num_cfiles
          do j = 1, num_sat
              cyc_set(j,i) = .false.
              curr_cyc(1,j,i) = 0.d0
              curr_cyc(2,j,i) = 0.d0
              curr_ion(j,i)   = 0.d0
*             Set initial value large so that range ion will
*             be used if both ranges are available
              curr_ion_ep(j,i) = -1000
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

*         Carry the current set biases foward to this epoch of data
          do j = 1, num_cfiles
              do k = 1, actual_max_chan

*                Make sure we actually have data before saving
*                values
* MOD TAH 000309: Only reset the number of cycles on the first pass
                 if( .not.kbit(data_flag_cse(k,j,i),30) .and. 
     .                cpass.lt.2 ) then
                      L1_cyc_cse(k,j,i) = curr_cyc(1,ctol_cse(k,j,i),j)
                      L2_cyc_cse(k,j,i) = curr_cyc(2,ctol_cse(k,j,i),j)
                 end if

*                Check to see of we need to set a bias flag
*                because there is a new channel/site observed
*                at the next epoch.
*                                      ! Data OK at next epoch
                 if( data_ok(data_flag_cse(k,j,i),0, phs_mask)
*                                      ! Bias not yet set for this
*                                      ! satellite/site combination
     .      	     .and. .not.cyc_set(ctol_cse(k,j,i),j) ) then
 
*                     Introduce the ctogobs bias flag and
*                     indicates that the initial bias has been
*                     set for this satellite/site combination
                      call sbit(data_flag_cse(k,j,i),31,1)
                      cyc_set(ctol_cse(k,j,i),j) = .true.
                  end if
              end do
          end do
 
*         Now get the new estimates at this epoch and sum the
*         range residual rms estimates.  Update the cycles for
*         those sites and channels which we do not know yet.
*         First call accounting for cycles.
          call est_clk_phs(i, phs_bias_mask, cyc_updated, 
     .                    L1r_phs_cse(1,1,i), L2r_phs_cse(1,1,i),
     .                    L1_cyc_cse(1,1,i), L2_cyc_cse(1,1,i),
     .                    L1r_rng_cse(1,1,i),  L2r_rng_cse(1,1,i),
     .                    ctol_cse(1,1,i), data_flag_cse(1,1,i),
     .                    bf_type_cse(1,1,i), params_cse(1,i), 
     .                    par_flag_cse(1,i), 1, cpass)

*         If we updated cycles, then reestimate the clocks accounting
*         for the new bias values.  (Use just phs_mask, so that bias
*         flags will be ignored.)
          if( cyc_updated ) then
             call init_norm(i, .true. )
             call est_clk_phs(i, phs_mask, cyc_updated, 
     .                       L1r_phs_cse(1,1,i), L2r_phs_cse(1,1,i),
     .                       L1_cyc_cse(1,1,i), L2_cyc_cse(1,1,i),
     .                       L1r_rng_cse(1,1,i),  L2r_rng_cse(1,1,i),
     .                       ctol_cse(1,1,i), data_flag_cse(1,1,i),
     .                       bf_type_cse(1,1,i), params_cse(1,i), 
     .                       par_flag_cse(1,i), 2, cpass)

          end if

*         See if we need a third pass because of additional updates.         
          if( cyc_updated ) then
             call init_norm(i, .true. )
             call est_clk_phs(i, phs_mask, cyc_updated, 
     .                       L1r_phs_cse(1,1,i), L2r_phs_cse(1,1,i),
     .                       L1_cyc_cse(1,1,i), L2_cyc_cse(1,1,i),
     .                       L1r_rng_cse(1,1,i),  L2r_rng_cse(1,1,i),
     .                       ctol_cse(1,1,i), data_flag_cse(1,1,i),
     .                       bf_type_cse(1,1,i), params_cse(1,i), 
     .                       par_flag_cse(1,i), 3, cpass)

          end if

*         Now save the current cycles so that we can propagate forwards 
          do j = 1, num_cfiles
              do k = 1, actual_max_chan
 
                  if(  .not.kbit(data_flag_cse(k,j,i),30) ) then
                      curr_cyc(1,ctol_cse(k,j,i),j) = L1_cyc_cse(k,j,i) 
                      curr_cyc(2,ctol_cse(k,j,i),j) = L2_cyc_cse(k,j,i) 
                  end if
              end do
          end do
 
*                                 ! Looping over epochs
      end do
 
 
****  Thats all
      return
      end
 
CTITLE EST_CLK_PHS
 
      subroutine est_clk_phs(ep, mask, cyc_updated,
     .                L1r_phs_cs, L2r_phs_cs, L1_cyc_cs, L2_cyc_cs,
     .                L1r_rng_cs, L2r_rng_cs, ctol_cs, data_flag_cs, 
     .                bf_type_cs, params_cs, par_flag_cs, pass, cpass )

      implicit none
 
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
 
*   ctol_cs(num_chan, num_cfiles)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   bf_type_cs(num_chan, num_cfiles) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   par_flag_cs(num_param)         - parameter estimate flag
*   ep      - Epoch number being processed.
*   mask    - Mask to use in deciding if data should be included
*             On first pass it will not use data with bias flag
*             set, on second pass these are ignored.
*   pass    - Indicates which pass through the clock estimation this is:
*             pass = 1 uses the psuedorange clocks
*             pass = 2 uses the phase clocks from the previous pass as the
*                      apririo values
*   cpass   - Cleaning level pass.  One first pass (before double difference
*             scanning, we do not update clock estimates).
 
 
      integer*4 ctol_cs(num_chan, num_cfiles),
     .    data_flag_cs(num_chan, num_cfiles),
     .    bf_type_cs(num_chan, num_cfiles),
     .    par_flag_cs(num_param), ep, mask, pass, cpass

*   cyc_updated - Logical to indicate that cycles values have been
*                 estimated

      logical cyc_updated
 
*   L1_cyc_cs(num_chan, num_cfiles) - Number of L1 cycles to be
*                    added to the phase values (May be fractional
*                    (here we try only to resolve to the nearest
*                    full wavelength)
*   L2_cyc_cs(num_chan, num_cfiles) - Number of L2 cycles to be
*                     added (size as above).
*   L1r_phs_cs(num_chan, num_cfiles)    - L1 phase measurements (L1 cylces)
*   L2r_phs_cs(num_chan, num_cfiles)    - L2 phase measurements (L2 cylces)
*   L1r_rng_cs(num_chan, num_cfiles)    - L1 range measurements (L1 cylces)
*   L2r_rng_cs(num_chan, num_cfiles)    - L2 range measurements (L2 cylces)
*   params_cs(num_param)    - Clock parameter estimates (L1 cycles)
 
      real*8 L1_cyc_cs(num_chan, num_cfiles),
     .    L2_cyc_cs(num_chan, num_cfiles),
     .    L1r_phs_cs(num_chan, num_cfiles),
     .    L2r_phs_cs(num_chan, num_cfiles),
     .    L1r_rng_cs(num_chan, num_cfiles),
     .    L2r_rng_cs(num_chan, num_cfiles), params_cs(num_param)
 
* LOCAL VARIABLES
 
 
*   i,j,k       - Loop counters
*   ls          - Satellite number
*   ns          - parameter number of satellite
*   ipivot(max_neq)  - Pivot elmenent for inversion
*   par_loc(max_cfiles+max_gsvs)  - array to indicate which parameters
*                 can be estimated.  (Initially set to zero and
*                 set to one when normal equations are incremented).
 
 
      integer*4 i,j, ls,ns, ipivot(max_neq),
     .    par_loc(max_cfiles+max_gsvs)
 
*   data_ok     - Logical function which is true if data OK
*   kbit        - Logical function to see if bit is set.
*   biases      - Logical to indicate that bias flags have been set.
 
      logical data_ok, kbit, biases, dump
 
*   res         - Generic residual values
*   phs_omc     - Function to return range residual
*   data_var    - Data varaince at one site
*   scale(max_neq)  - Scale factors for invert_vis
*   L1_cyc_sav, L2_cyc_sav  - Values of the cycles before we
*                 estimate new values.  These are used to
*                 see if we have updated the values.

*   fit_res(max_gsvs) - Prefit residuals of phase to range 
*                 estimates of the clock.  If these are large then
*                 the largest ones are assigned bias flags.
*   mean_res     - Mean value of the prefit residuals (cycles)
*   max_res      - Biggest residual when the fit is not good.
*                  if this differrs from the mean by too much
*                  then a bias flag is set. (cycles)
*   mean_tol, max_tol - Tolerances for match of prefit residuals
*                  The values used depend on whether this is pass 1 or
*                  pass 2. (cycles)
*   dt_ion       - Time differnece in secinds bewteen last estimate of
*                  ionopsheric delay and the current one.
*   trial_ion    - Trial estimate of the ionospheric delay to see
*                  if it matches previuos epoch. This is mainly for
*                  those receivvers that seem to drift by cycles.
*   trial_dion   - Trial value for the change in the ionospheric delay
*   dion_tol     - Tolerance for change in ion delay:  Max of min allowed
*                  or last change multiplied by the user multiplier.
*   dcyc(2)      - Chang in the number of cycles at L1 and L2
 
      real*8 res, phs_omc, data_var, scale(max_neq),
     .    L1_cyc_sav, L2_cyc_sav, fit_res(max_gsvs), mean_res, max_res,
     .    mean_tol, max_tol, dt_ion, dcyc(2), trial_ion, trial_dion,
     .    dion_tol

*   num_res     - Number of residuals fit residual
*   ms_off      - Offset in phase due to millisecond jump.

      integer*4 num_res, ms_off 
 
***** Initialize the parameter local set array
 
      do i = 1, num_param
         par_loc(i) = 0
      end do

***** Initialize the tolerance for the prefit resiudals.  If pass 1 ranges
*     are being used so set larges; if pass 2 then phases so set low.
      if( pass.eq.1 .or. cpass.eq.1 ) then
         mean_tol = phs_fit_tol(1)
         max_tol  = phs_fit_tol(2)
      else
         mean_tol = phs_fit_tol(3)
         max_tol  = phs_fit_tol(4)
      end if
 
***** Loop over the data incrementing the normal equations
      do i = 1,num_cfiles
          data_var = phs_noise(i)**2

*         Initialize values for computing the mean residual at
*         this site.
          num_res = 0
          mean_res = 0.d0
          dump = .false.
          biases = .false. 

          do j = 1, actual_max_chan

*             clear the fit residual (so that we can tell later
*             if there is data here)
              fit_res(j) = 0
 
*             Use only data that has no error and for which
*             there is no bias flag.
* MOD TAH 200616: Added ctol check to make sure data present
              if( ctol_cs(j,i).gt.0 .and. 
     .            data_ok(data_flag_cs(j,i),0, mask) ) then
 
*                 Compute OminusC
                  ls = ctol_cs(j,i)
                  ns = num_cfiles + ctol_cs(j,i)
                  res = phs_omc(L1r_phs_cs(j,i),L2r_phs_cs(j,i),
     .                    L1_cyc_cs(j,i), L2_cyc_cs(j,i),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(ls), fL2(ls) )

* MOD TAH 050321: Check for millisec effect residuals.  Problem here
*                 is some files have MS jumps in both range and phase
*                 while others just have range jumps.  Code below
*                 catches case of range only jumps (and add jump to
*                 phase) 
                  ms_off = nint(res/(fL1(ls)*0.001d0))*
     .                         (fL1(ls)*0.001d0)
                  if( abs(ms_off).gt.0 ) then
*                    Possible ms jump; see if close
                     if( abs(res-ms_off).lt.100.d0 ) then
*                      Correct for ms_jmp
                       res = res - ms_off
* MOD TAH 180320: No change made here for Glonaas -- should not have
*                 millisecond offsets (unique to GPS).
                       L1_cyc_cs(j,i) = L1_cyc_cs(j,i)-ms_off
                       L2_cyc_cs(j,i) = L2_cyc_cs(j,i)-ms_off*
     .                       fL2(ls)/fL1(ls)
                       write(*,90) ep,cf_codes(i),
     .                         prn_list(ctol_cs(j,i)), res, 
     .                         nint(ms_off/(fL1(ctol_cs(j,i))*0.001d0))
 90                      format('MS PHS JUMP Ep ',i4,' Site ',a4,
     .                         ' PRN ',i2.2,' Res ',f8.2,
     .                         ' cylces, Jump ',i4,' ms') 
                     end if
                  endif                               
                  fit_res(j) = res
               
                      
                  if( abs(res).gt.mean_tol ) then 
                      biases = .true.
                  endif
              end if
          end do

*****     If biases status was set then, see if all of the residuals are
*         about the same size.
          do while ( biases )
              call prescan_phs( fit_res, actual_max_chan, 
     .                          data_flag_cs(1,i), bf_type_cs(1,i),
     .                          biases, max_res, max_tol, mean_tol, j )
              if( biases ) then
                   if( kbit(status_rep,3) ) 
     .            write(*,100) ep,cf_codes(i),prn_list(ctol_cs(j,i)), 
     .                         max_res, max_tol
 100              format('JMP BIAS flag added at ',I4,' Site ',a4,
     .                   ' PRN ',i2.2, 1x,f20.2,1x,f10.1)
              end if
          end do

*         Now compute the next iteration on the clock values
          do j = 1, actual_max_chan

*             Use only data that has no error and for which
*             there is no bias flag.
              if( data_ok(data_flag_cs(j,i),0, mask) ) then
 
*                 Compute OminusC
                  ns = num_cfiles + ctol_cs(j,i)
                  ls = ctol_cs(j,i)
                  res = phs_omc(L1r_phs_cs(j,i),L2r_phs_cs(j,i),
     .                    L1_cyc_cs(j,i), L2_cyc_cs(j,i),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(ls), fL2(ls) )
 
*                 Increment normal equations based on the data quality
                  norm_eq(i,i)   = norm_eq(i,i)   + 1.d0/data_var
                  norm_eq(ns,ns) = norm_eq(ns,ns) + 1.d0/data_var
                  norm_eq(ns,i)  = norm_eq(ns,i)  - 1.d0/data_var
                  norm_eq(i,ns)  = norm_eq(i,ns)  - 1.d0/data_var
                  sol_eq(i)      = sol_eq(i)      + res/data_var
                  sol_eq(ns)     = sol_eq(ns)     - res/data_var
                  par_loc(i) = 1
                  par_loc(ns) = 1
              end if
          end do
      end do
 
****  Check the parameter flag to see if we have data on a s
*     parameter.  If we no not set its diagonal to one.
      do i = 1, num_param
*                                         ! No data
          if( par_loc(i).eq.0 )  then
              norm_eq(i,i) = 1.d0
          else 
* MOD TAH 2006016: Add small amount to diagonals in case 
*             norm_eq is singular for small data sets. 
              norm_eq(i,i) = norm_eq(i,i) + 1.d-6
          end if
      end do
 
*     Now invert the system
      call invert_vis(norm_eq, sol_eq, scale, ipivot, num_param,
     .                max_neq,1)
 
*     Save as apriori and save values in parameter array if this is
*     not the first pass.
      if( cpass.gt. 1 ) then
          do i = 1, num_param
              params_cs(i) = apr_clk_val(i) + sol_eq(i)
              apr_clk_val(i) = params_cs(i)
 
*             if we have data, save this epochs value as the
*             apriori for the next epoch
              if( par_loc(i).eq.1 )  then
                  apr_clk_epoch(i) = ep
              end if
          end do
      end if

***** Now loop over data and set the numbers of cycles at those
*     sites with bias flags (These are just approximate calculations
*     aimed at getting the residual error in the number of cycles
*     to < 100 cycles.)
      cyc_updated = .false.
      do i = 1,num_cfiles
          do j = 1, actual_max_chan
                                        
              ls = ctol_cs(j,i)
*             see if we actually have data.
              if( .not.kbit(data_flag_cs(j,i),30) ) then
 
*                  We have data, see if there is a bias flag here
                   if( kbit(data_flag_cs(j,i),31) .or.
     .                 kbit(data_flag_cs(j,i),32)      ) then
 
*                      Estimate L1 and L2 numbers of cycles (NOTE here
*                      that we do not pass in the wavelength factors
*                      because this is an approximate resolution of the
*                      cycles only. i.e., we resolve to nearest full
*                      wavelength
                       ns = num_cfiles + ctol_cs(j,i)

*                      Save current values
                       L1_cyc_sav = L1_cyc_cs(j,i)
                       L2_cyc_sav = L2_cyc_cs(j,i)

*                      Comput time from last ionospheric delay
                       dt_ion = (ep-curr_ion_ep(ctol_cs(j,i),i))*
     .                           sampling

* MOD TAH 950713: Added curr_ion to call to subroutine.  (Apparently
*                 missed in original coding.  Should yeild more 
*                 consistent allignment of phase and range measurements.)
* MOD TAH 180320: Pass in original frequencies as well as remapped ones
*                 (Currently only differ for Glonass)
                       call est_cyc( L1r_phs_cs(j,i), L2r_phs_cs(j,i),
     .                               L1_cyc_cs(j,i),  L2_cyc_cs(j,i),
     .                               L1r_rng_cs(j,i), L2r_rng_cs(j,i),
     .                               curr_ion(ctol_cs(j,i),i),
     .                               dt_ion, apr_clk_val(i),
     .                               apr_clk_val(ns),ep,
     .                               fL1(ls),fL2(ls), fL1u(ls),fL2u(ls))

*****                  Get the number of cycles change
                       dcyc(1) = L1_cyc_sav-L1_cyc_cs(j,i)
                       dcyc(2) = L2_cyc_sav-L2_cyc_cs(j,i)

*****                  Save that we have updated the number of cycles

                       if( dcyc(1).ne.0.d0 .or.
     .                     dcyc(2).ne.0.d0  ) then
                           cyc_updated = .true.
                       end if

*                      Now test to see if values changed provided that
*                      saved values are not zero (indicating that we
*                      have not used this satellites before.)
                       if(  L1_cyc_sav.ne.0.d0 .and.
     .                      L2_cyc_sav.ne.0.d0  ) then

*                          Only print out the message if requested to do
*                          so
                           if( (kbit(status_rep,4) .and. pass.eq.1) .or.
     .                         (kbit(status_rep,5) .and. pass.ge.2) )
     .                     write(*,200) ep, cf_codes(i),
     .                           prn_list(ctol_cs(j,i)), dcyc,pass,
     .                           curr_ion(ctol_cs(j,i),i)
 200                       format('Updating at ',i5,' site ',a4,' PRN ',
     .                           i2.2,' cycles by ', 2F13.2,' Pass ',i2,
     .                           ' Curr Ion ',f13.2)

*****                      Now as a final check, if the number of
*                          cycles to be shift is small, then set it
*                          back to zero. [All changes < 100 are not
*                          added in]
                           if( abs(dcyc(1)).lt.10 ) 
     .                         L1_cyc_cs(j,i) = L1_cyc_sav
                           if( abs(dcyc(2)).lt.10 ) 
     .                         L2_cyc_cs(j,i) = L2_cyc_sav
                               
                       end if
                 end if

*                If this data is OK, then estimate the current value
*                of the ionospheric delay.  We will use this later
*                for getting number of cycles when the there is a bias
*                flag later
                 if( data_OK(data_flag_cs(j,i),0,mask) .and.
     .               L2r_phs_cs(j,i).ne. 0.d0 ) then

*                    compute trial ion delay first.  If there is no
*                    gap or the gap is with-in the user tolerance
*                    then check to see if ion falls with in the
*                    random walk tolerance.  Only do this on the first
*                    pass since since for the second pass tje ion delay
*                    will already be updated.
                     dt_ion = (ep-curr_ion_ep(ctol_cs(j,i),i))*
     .                         sampling
                     trial_ion =
     .                ((L1r_phs_cs(j,i)+L1_cyc_cs(j,i)) -
     .                (fL1(ctol_cs(j,i))/fL2(ctol_cs(j,i)))*
     .                    (L2r_phs_cs(j,i)+L2_cyc_cs(j,i)))/
     .                        (1-(fL1(ls)/fL2(ls))**2)

*                    Since we will pass here on the second pass; only
*                    update the trial_dion value on pass 1; on pass 2
*                    use the first pass value unless the bias flags
*                    are set (in which case, this was not computed on
*                    first pass).	
                     if( pass.eq.1 ) then
                         trial_dion = trial_ion - 
     .                                curr_ion(ctol_cs(j,i),i)
                     else
                         trial_dion = curr_dion(ctol_cs(j,i),i)
                     end if

****                 Now if this is first pass check to see within
*                    tolerance.
                     if( pass.eq.1 .and.
     .                   (dt_ion.le. dt_ion_tol .or.
     .                    dt_ion.eq.orig_sampling(i)) ) then

*                        check the ion delay continuity.  If we have
*                        no history yet, use the maximum value
*                        otherwize compute based on last change.
                         if( curr_dion(ctol_cs(j,i),i).ne.0 ) then
                             dion_tol = min(ion_rw_tol(3,i),
     .                            max(ion_rw_tol(2,i),ion_rw_tol(1,i)*
     .                                abs(curr_dion(ctol_cs(j,i),i))))
                         else
                             dion_tol = ion_rw_tol(3,i)
                         end if
*                        See if we pass th continuity condition  
                         if( abs(trial_dion-curr_dion(ctol_cs(j,i),i))
     .                       .gt.dion_tol .and.
     .                       abs(trial_dion).gt.dion_tol ) then

*                            Ion is discontuous: Add bias flag and tell
*                            user:
                             write(*,250) ep, cf_codes(i), 
     .                           prn_list(ctol_cs(j,i)),
     .                           trial_dion-curr_dion(ctol_cs(j,i),i),
     .                           dt_ion, dion_tol, trial_dion, 
     .                           curr_dion(ctol_cs(j,i),i)
 250                         format('ION BIAS flag EP ',I5,' Site ',a4,
     .                              ' PRN ',i2.2,1x,' Dion ',f7.2,
     .                              ' dT ',f6.1,' sec; tol ',F5.2,2f7.2)
                             call sbit(data_flag_cs(j,i),31,1)
                             call sbit(bf_type_cs(j,i), 3,1)
                         end if
                     end if

*                    say current ion delay.
                     curr_ion(ctol_cs(j,i),i)  = trial_ion
*                    Only update the curr_dion value if it is not
*                    large (so that the next epoch will not be flagged)
                     if( abs(trial_dion).lt.ion_rw_tol(3,i) ) then
                         curr_dion(ctol_cs(j,i),i) = trial_dion
                     else
                         curr_dion(ctol_cs(j,i),i) = 0.0d0
                     end if
                     curr_ion_ep(ctol_cs(j,i),i) = ep
                 end if
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE EST_CYC
 
      subroutine est_cyc( L1r_phs, L2r_phs, L1_cyc, L2_cyc,
     .           L1r_rng, L2r_rng, curr_ion, dt_ion, 
     .           rcv_clk, svs_clk, ep, f1, f2, f1u, f2u )

      implicit none
 
*     This routine will get estimates of the number of cycles for L1
*     and L2 but matching to the apriori clock model.  If dual frequency
*     ranges are available then these will be used to account for the
*     ionospheric delay.
 
* INCLUDES
 
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
*   L1_cyc    - Number of L1 cycles to be
*               added to the phase values (May be fractional
*               (here we try only to resolve to the nearest
*                full wavelength)
*   L2_cyc    - Number of L2 cycles to be
*               added (size as above).
*   L1r_phs   - L1 phase measurements (L1 cylces)
*   L2r_phs   - L2 phase measurements (L2 cylces)
*   L1r_rng   - L1 range measurements (L1 cylces)
*   L2r_rng   - L2 range measurements (L2 cylces)
*   dt_ion    - Time difference from current eppoch and previous 
*               estimate of the ionospheric delay (secs)
*   curr_ion  - Current ionospheric delay (cyc)
*   rcv_clk   - The clock offset at the reciever (L1 cycles)
*   svs_clk   - The clock offset at the satellite (L1 cycles)
*   ep        - Epoch number (just for output).  Value not used for any
*               computations.
*   f1,f2     - L1 and L2 frequencies for this channel/site/epoch 
*   f1u, f2u  - Original L1 and L2 frequencies used for remapped phase values.

      integer*4 ep
 
      real*8 L1_cyc, L2_cyc, L1r_phs, L2r_phs, L1r_rng, L2r_rng,
     .    curr_ion, dt_ion, rcv_clk, svs_clk, f1, f2, f1u, f2u
 
* LOCAL VARIABLES
 
*  ion_L1    - Ionospheric delay for L1 Phase measurements (set to
*            - 50 cycles (typical ionodelay) when not estimate is
*            - availiable.
*  dL1, dL2, dLC - Differernces in L1 L2 and LC.  We select the
*              L2 number of cycles to match dLC as closely as
*              possible
*  dL1_msw, dL1_lsw - The dL1 value splitt into two parts so that
*              we do not have integer*4 overflow.
*  dL2_msw, dL2_lsw - The dL2 value split into two parts so that
*              we do not have integer*4 overflow.
 
      real*8 ion_L1, dL1, dL2, dLC, dL1_msw, dL1_lsw,
     .       dL2_msw, dL2_lsw
      real*8 L2_fract  ! Compuation of L2 fractional.  For non-Glonass
                       ! we could take nint of this value 
 
***** Try to get an estimate of L1 ionosphere
 
      if( L2r_rng.ne.0.d0 ) then

*         use the range estimate of the ion delay if more than
*         one hour has elapsed.
          if( dt_ion.gt. 3600.d0 ) then
              ion_L1 = -(L1r_rng-(f1/f2)*L2r_rng)/(1.d0-(f1/f2)**2)

* MOD TAH 98101: See if this value make sense.
              if( ion_L1.gt.0 .or. ion_L1.lt.-150 ) then
                  ion_L1 = -30 
              end if
C             ion_L1 = -(-L1r_rng*(f2/f1)+L2r_rng)
          else
              ion_L1 = curr_ion
          end if
      else

*         If we don't have range, then use last value unless this
*         if first observation in which case arbitarily set to an
*         "zero ion delay"
          if( curr_ion.ne.0.d0 ) then
              ion_L1 = curr_ion
          else
              ion_L1 =  -30.d0
          end if
      end if
*
*     Now compute the number of cycles which will make the phase match
*     the range estimates of the clocks.
*     FORMULATION:
*     L1 = LC + LG
*     L2 = (f2/f1)*LC + (f1/f2)*LG
*     P1 = LC - LG
*     P2 = (f2/f1)*LC - (f1/f2)*LG
 
      dL1 = -(L1r_phs-(rcv_clk-svs_clk + ion_L1))
      if( L2r_phs.ne. 0.d0 ) then
         dL2 = -(L2r_phs-((f2/f1)*(rcv_clk-svs_clk) +
     .                       (f1/f2)*ion_L1))
      else
         dL2 = 0.0d0
      end if

*     Compute change in LC
      dLC = (dL1 - (f2/f1)*dL2)/(1-(f2/f1)**2)

*     Check to make sure that we do not overflow integer*4 on 
*     nint command
      if( abs(dL1).gt.2.d0**31 .or. abs(dL2).gt.2.d0**31 ) then
          write(*,100) dL1, dL2
 100      format('**WARNING** Integer*4 overflow in set cycles',
     .          /'            Splitting ',2F20.3,' to resolve' )
* MOD TAH 180320: Did not update this code for Glonass (should not be
*         neded for modern data).
          dL1_msw = nint(dL1/1.d4)*1.d4 
          dL1_lsw = dL1 - dL1_msw
          L1_cyc  = nint(dL1_lsw) + dL1_msw
          dL2 = (L1_cyc - dLC*(1-(f2/f1)**2))*(f1/f2)
          dL2_msw = nint(dL2/1.d4)*1.d4 
          dL2_lsw = dL2 - dL2_msw
          L2_cyc  = nint(dL2_lsw) + dL2_msw
          write(*,120) L1_cyc, L2_cyc 
 120      format('            Resolved  ',2F20.3,' cycles')
      else
*         Do computations with standard nint function.
* MOD TAH 180320: Adjust for remapping frequencies.  This will make
*         number of cycles non-integer for Glonass.
          L1_cyc = nint(dL1*f1u/f1)*(f1/f1u)
          if( L2r_phs.ne. 0.d0 ) then
*            L2_cyc = nint((L1_cyc - dLC*(1-(f2/f1)**2))*(f1/f2))
             L2_fract = (L1_cyc - dLC*(1-(f2/f1)**2))*(f1/f2)
             L2_cyc = nint(L2_fract*f2u/f2)*(f2/f2u)
          else
             L2_cyc = 0.d0
          end if
      end if

***** Thats all
      return
      end
 
CTITLE PHS_OMC
 
      real*8 function phs_omc( L1r_phs, L2r_phs, L1_cyc, L2_cyc,
     .    rcv_clk, svs_clk, f1, f2 )

      implicit none
 
*     Function to return the ionospheric free phase residual relative
*     to the current clock values.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   L1_cyc    - Number of L1 cycles to be
*               added to the phase values (May be fractional
*               (here we try only to resolve to the nearest
*                full wavelength)
*   L2_cyc    - Number of L2 cycles to be
*               added (size as above).
*   L1r_phs   - L1 phase measurements (L1 cylces)
*   L2r_phs   - L2 phase measurements (L2 cylces)
*   rcv_clk   - The clock offset at the reciever (L1 cycles)
*   svs_clk   - The clock offset at the satellite (L1 cycles)
*   f1, f2    - L1 and L2 frequencies for the channel (named so as not
*               to conflict with the fL1, fl2 arrays in ctogobs_comb.h
 
 
      real*8 L1_cyc, L2_cyc, L1r_phs, L2r_phs, rcv_clk, svs_clk, f1, f2
 
* LOCAL VARIABLES
 
*   res     - Residual computed accounting for the ionosphere
*           - and the current clock values.
 
 
      real*8 res
 
***** Fitsly compute the ionospheric free residual before adding clock
*     model.
* MOD TAH 990512: See if we have L2 phase measurements.
* MOD TAH 880524: Also check to see if fdd_L2_fact is zero.
      if( L2r_phs.ne. 0.d0 .and. fdd_L2_fact.ne.0 ) then 
         res = ((L1r_phs+L1_cyc)-(f2/f1)*(L2r_phs+L2_cyc))/
     .        (1.d0 - (f2/f1)**2)
      else
         res = (L1r_phs+L1_cyc)
      end if
 
*     Now remove the clock contributions
      res = res - (rcv_clk-svs_clk)
 
*     Set the function return
      phs_omc = res
 
****  Thats all
      return
      end
 
CTITLE WRITE_PHS
 
      subroutine write_phs(root, L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, params_cse, par_flag_cse,
     .    azel_cse  )

      implicit none
 
*     This routine will write out all of the phase residuals
*      using the current estimates of the clocks and cycles.  These
*      residuals will be those used as the prefit residual in the
*      filter. (Unless site positions/ephemeris etc are updated
*      during the filter run.)
 
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
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L! range residuals
*                   - cylces at L1
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

*  azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
*                      angles to be written to phs_res_root files.
*                      (Real*4, radians)

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)  
 
*   root  - Root part of name for files

      character*(*) root
 
* LOCAL VARIABLES
 
*   i,j  - Loop counter over stations and  epochs
*   lv   - Satellite number 
*   num_out  - Records number of values output
*   trimlen  - Length of string
*   ierr     - IOSTAT error
*   plot_on_page - Keeps track of number of plots on the
*                  page (currently 4 per page)
*   plot_flag    - Flag for plotting.  Set to 0 if OK, 1 otherwize
 
      integer*4 i,j, k,l, lv, num_out, trimlen, ierr, plot_on_page,
     .          plot_flag
 
*   data_ok  - Logoical function returns true is data OK
*   kbit     - Checks if bit is set
 
 
      logical data_ok, kbit
 
*   min_sc, max_sc, ds  - Min and max values and difference for
*                         scales
*   min_io, max_io      - Min and Max scales for the ionosphere
*   omin_sc,omax_sc     - Original min and max scales for working
*            where on the plot the labels should go.
*   jump  -  Value needed to adjust plot for jumps in the
*            clock
*   clk   - clock value in usec after accounting for jumps
*   lab_offset - Estimate of offset in label (vertically) to
*           try to stop the text overwriting the data.
*   dL1, dL2, dLC, dLG  - L1 L2 and LC residualsi and LG.  Note our
*             definition of LG here is the ion contribution to L1
*             This is 1.983 times larger than the LG shown in cview
*             plots.
*   dRC     - LC range residual (If CA code range, we use the Phase
*             estimate of the ionsphere.)
*   NL, WL  - Narrow lane and widelane combinations for phase and
*             range measurements.
*   rng_omc - Function to compute reange residual.
 
 
      real*8 min_sc, max_sc, omin_sc, omax_sc, ds, rng_omc,
     .    lab_offset,  dL1, dL2, dLC, dLG, min_io, max_io, dRC,
     .    NL, WL, dP1, dP2 
 
*   outfile - Name of output file
*   plotfile = Name of the plot files to generate
 
      character*128 outfile, plotfile
 
***** Loop over all the sites creating and writing the output
*     file and generating the cplot command file.

      plotfile = root(1:trimlen(root)) // '.plt'  
      open(200, file=plotfile, iostat=ierr, status='unknown')
      plot_on_page = 0
 
****  Start looping over all sites and satellites
      do i = 1, num_cfiles
        if( kbit(phs_res_sites,i) ) then
*                         
          do j = 1, num_sat
 
 
              outfile = root(1:trimlen(root)) // '.' // cf_codes(i)

* MOD TAH 130329: Implemenation of single verus multiple output
*             DPH files (by PRN) 
              if(dph_output(1:3).eq.'SIN'.and.j.eq.1) then 
                 write(outfile(trimlen(outfile)+1:),100) 
 100             format('.PRNAL')
                 open(201,file=outfile,iostat=ierr,status='unknown')

              else if(dph_output(1:3).eq.'MUL') then
                 write(outfile(trimlen(outfile)+1:),110) prn_list(j)
 110             format('.PRN',I2.2)
                open(201,file=outfile,iostat=ierr,status='unknown')

              end if
              write(201,120,iostat=ierr) cf_codes(i),
     .             rcv_types(i), prn_list(j)

* MOD TAH 051101: Added P1 and P2 residuals to output
c 110          format('* Clock information for site ',a,
c     .               ' receiver ',a4,'. PRN ',i2.2,/,
c     .               '* Epoch   L1       L2       LC    LG     RC  ',
c     .               '     WL             NL  LSVS   ',
c     .               ' Azimuth     Elev PF    data_flag')
* MOD TAH 130329: Added PRN to line
 120          format('* Clock information for site ',a,
     .               ' receiver ',a4,'. PRN ',i2.2,/,
     .               '* Epoch  L1 cyc   L2 cyc   P1 cyc   P2 cyc',
     .               '   LC cyc   LG cyc   PC cyc   WL cyc   N cyc',
     .               '  LSV    Azimuth    Elev  PF    data_flag',
     .               '         L1_cycles            L2_cycles',
     .               '    PRN')

*             Set the min and max scales
              min_sc = 100.d9
              max_sc = -100.d9
              min_io = 100.d9
              max_io = -100.d9
              num_out = 0
*             Now loop over all the epochs write out those which
*             are for this site and satellite.
              do k = 1, num_ep
 
*                 Loop over the channels at this epoch and see if we
*                 have this satellite.
                  do l = 1, actual_max_chan
                       lv = ctol_cse(l,i,k) 
                       if ( lv.eq.j ) then
C    .                  data_ok(data_flag_cse(l,i,k),0,phs_mask)) then
*                         satellite match, Compute the L1 L2 and LC residuals
                          dL1 = L1r_phs_cse(l,i,k) + L1_cyc_cse(l,i,k)
     .                         - (params_cse(i,k) -
     .                            params_cse(num_cfiles+j,k))
                          dL2 = L2r_phs_cse(l,i,k) + L2_cyc_cse(l,i,k)
     .                         - (fL2(lv)/fL1(lv))*(params_cse(i,k) -
     .                            params_cse(num_cfiles+j,k))
                          dLC = (dL1-(fL2(lv)/fL1(lv))*dL2) /
     .                          (1.d0 - (fL2(lv)/fL1(lv))**2)
                          dLG = (dL1-(fL1(lv)/fL2(lv))*dL2) /
     .                          (1.d0 - (fL1(lv)/fL2(lv))**2)
C                         dLG = -(fL2(lv)/fL1(lv))*dL1 + dL2
                          dRC = rng_omc(L1r_rng_cse(l,i,k),
     .                                  L2r_rng_cse(l,i,k),
     .                                  params_cse(i,k), 
     .                                  params_cse(num_cfiles+j,k),
     .                                  fL1(lv),fL2(lv) )

* MOD TAH 051102: Added range residuals
                          dP1  = L1r_rng_cse(l,i,k) 
     .                         - (params_cse(i,k) -
     .                            params_cse(num_cfiles+j,k))

                          dP2  = L2r_rng_cse(l,i,k) 
     .                         - (fL2(lv)/fL1(lv))*(params_cse(i,k) -
     .                            params_cse(num_cfiles+j,k))

                          if( lambda(l,4,i).eq.0 ) then
* MOD TAH 990518: See if L2 phase is present
                              if( lambda(l,2,i).ne.0 ) then
* MOD TAH 060105: Set P2 quantiies to zero
C                                dRC = dRC + dLG
                                 dRC = dP1
                                 dP2 = 0
                                 WL = 0
                                 NL = 0
                              else
                                 dL2 = 0.d0
                                 dP2 = 0.d0
                                 dLC = dL1
                                 dRC = dP1
                                 dLG = 0.d0
                                 WL = 0
                                 NL = 0
                              end if
                          else  
                              WL = -(L1r_phs_cse(l,i,k) +
     .                               L1_cyc_cse(l,i,k))   +
     .                              (L2r_phs_cse(l,i,k) + 
     .                               L2_cyc_cse(l,i,k))    + 
     .                               dfsf(lv)*(L1r_rng_cse(l,i,k) +
     .                               L2r_rng_cse(l,i,k))
                              NL = (L1r_phs_cse(l,i,k) +
     .                              L1_cyc_cse(l,i,k))   +
     .                             (L2r_phs_cse(l,i,k) + 
     .                              L2_cyc_cse(l,i,k))   - 
     .                              sfdf(lv)*(L1r_rng_cse(l,i,k)-
     .                              L2r_rng_cse(l,i,k))
                          end if

*                         if the fdd_L2_fact is zero, then we made clocks
*                         with L1 only so save L1 residual as LC
                          if( fdd_L2_fact.eq.0 ) dLC = dL1

*                         Check to see if good or bad for 1 digit flag
                          plot_flag = 1
                          if( data_ok(data_flag_cse(l,i,k),0,
     .                            phs_mask) )  plot_flag = 0

* MOD TAH 990518: If dLC is large then set value to 9999 so that it
*                 written to file (to show editing)
                          if( .not.
     .                    data_ok(data_flag_cse(l,i,k),0,phs_mask)) then
                              if( abs(dL1).gt. 1.d5 ) dL1 =  9999.d0
                              if( abs(dl2).gt. 1.d5 ) dl2 =  9999.d0
                              if( abs(dP1).gt. 1.d5 ) dP1 =  9999.d0
                              if( abs(dP2).gt. 1.d5 ) dP2 =  9999.d0
                              if( abs(dLC).gt. 1.d5 ) dLC =  9999.d0
                              if( abs(dLG).gt. 1.d5 ) dLG =  9999.d0
                              if( abs(WL ).gt. 1.d5 ) WL  =  9999.d0
                              if( abs(NL ).gt. 1.d5 ) NL  =  9999.d0
                              if( abs(dRC).gt. 1.d5 ) dRC =  9999.d0
                          end if

                          if( abs(dLC).lt.10000.d0 .or.
     .                    data_ok(data_flag_cse(l,i,k),0,phs_mask)) then 

* MOD TAH 130329:            See if NONE output option selected 
                             if(dph_output(1:4).ne.'NONE') then
* MOD TAH 980215: Write out cycle offsets if bias flag is set.
                                if( kbit(data_flag_cse(l,i,k),32) .or.
     .                              kbit(data_flag_cse(l,i,k),31) ) then
                                    write(201,300) k, dL1, dL2, dP1,   
     .                                  dP2, dLC, dLG, dRC, WL,NL, j, 
     .                                  azel_cse(1,l,i,k)*180.d0/pi,
     .                                  azel_cse(2,l,i,k)*180.d0/pi, 
     .                                  plot_flag, data_flag_cse(l,i,k),
     .                                  L1_cyc_cse(l,i,k),
     .                                  L2_cyc_cse(l,i,k), prn_list(j)
                                else
                                    write(201,310) k, dL1, dL2, dP1,  
     .                                  dP2, dLC, dLG, dRC, WL,NL, j, 
     .                                  azel_cse(1,l,i,k)*180.d0/pi,
     .                                  azel_cse(2,l,i,k)*180.d0/pi, 
     .                                  plot_flag, data_flag_cse(l,i,k),
     .                                  prn_list(j)
                                end if
                             end if 
 300                         format(I5,1x,4(1x,F8.2), 1x,F8.3,
     .                              4(1x,F8.2),
     .                              1x,i3,1x,2f10.4,1x, i1, 1x, o12,
     .                              2(1x,F20.4),2x,I2.2)
 310                         format(I5,1x,4(1x,F8.2), 1x,F8.3,
     .                              4(1x,F8.2),
     .                              1x,i3,1x,2f10.4,1x, i1, 1x, o12,
     .                              44x, I2.2)
                          end if
                          plot_flag = 1
                          if( data_ok(data_flag_cse(l,i,k),0,
     .                            phs_mask) ) then
                               min_sc = min(min_sc, dLC )
                               max_sc = max(max_sc, dLC )
                               min_io = min(min_io, dLG )
                               max_io = max(max_io, dLG )
                               num_out = num_out + 1
                               plot_flag = 0
                          end if
*                                 ! satellite matches
                      end if
*                                 ! Looping over epoch
                  end do
              end do
 
*             If we output data, make entry in cplot file
              if( num_out.gt.0 ) then
 
                  plot_on_page = plot_on_page + 1
 
*                 Clip to nearest 0.1 cycles for LC and 1 cycle
*                 for LG.  Save original scale
*                 so we know if plot trends down or up
                  omax_sc = max_sc
                  omin_sc = min_sc
                  max_sc = int(max_sc*10)/10.0 + 0.1
                  min_sc = int(min_sc*10)/10.0 - 0.1
                  max_io = int(max_io)    + 1
                  min_io = int(min_io)    - 1
                  ds = max_sc - min_sc

*                 Set whether the labels are at the top or bottom
*                 of the window.  (Since the plot starts at zero
*                 we see if the averge of scales is positive or
*                 negative.
                  if( (omax_sc + omin_sc )/2.lt.0 ) then
                      lab_offset = 0
                  else
                      lab_offset = (max_sc-min_sc) - 5*(ds/20)
                  end if
 
*                 Now do commands
                  write(200,400) 0.05+(4-plot_on_page)*.225,
     .                       0.05+(5-plot_on_page)*.225,
     .                       outfile, num_ep+1, min_sc, max_sc,
     .                       lab_offset+(ds/20)*4,
     .                       lab_offset+(ds/20)*3,
     .                       lab_offset+(ds/20)*2,
     .                       lab_offset+(ds/20)*1,
     .                       num_ep+1, min_io, max_io
 400              format('  head 2 1 ',/,' font 5x7 ',/,
     .                 '  char 2.0',/,
     .                 '  view 0.1 0.90 ',2f8.4,/,
     .                 '  file ',a,/,
     .                 '  x_field 1 1 0',/,
     .                 '  y_field 1 4 0 "LC Residual (cycles)"',/,
     .                 '  p_field 9 12 ',/,
     .                 '  read ',/,
     .                 '  x_scale  0 ',i6,/,
     .                 '  y_scale  ',f10.3,1x,f10.3,/,
     .                 '  point -1',/,
     .                 '  char 1.0',/,
     .                 '  line 0',/,
     .                 '  draw ',/,
     .                 '  char 1.7',/,
     .                 '  xmn -1 0 ',/,
     .                 '  xmx -1 0 ',/,
     .                 '  ymn -1 1 "LC Res. (cycles)"',/,
     .                 '  poly ep cyc  1 1', /,
     .                 '  fit 0',/,
     .                 '  pdr ',/,
     .                 '  label 10 ',f10.3,' 1 0  :p1',/,
     .                 '  label 10 ',f10.3,' 1 0  :p2',/,
     .                 '  label 10 ',f10.3,' 1 0  :p3',/,
     .                 '  label 10 ',f10.3,' 1 0  :h1',/, 
     .                 '  y_field 1 5 0',/,
     .                 '  read ',/,
     .                 '  x_scale 0 ', i6,/,
     .                 '  y_scale  ',f10.3,1x,f10.3,/,
     .                 '  point 1',/,
     .                 '  draw ',/,
     .                 '  ymx -1 1 "LG (L1 Cycles)"')
 
                   if( plot_on_page.eq.1 ) then
                        write(200,'('' xmx -1 1 "Epoch number"'')')
                   end if
                   if( plot_on_page.eq.4 ) then
                       plot_on_page = 0
                       write(200,'('' ident'',/)')
                       write(200,'('' erase'',/)')
                   end if
*                         ! There was something to plot
              end if
*                         ! Looping over satellites
* MOD TAH 130329: Only close file if multiple files
              if(dph_output(1:3).eq.'MUL') close(201)
          end do

* MOD TAH 130329: Only close file if single PRN per file
          if(dph_output(1:3).eq.'SIN') close(201)

*                         ! Site selected.
        end if
*                         ! Looping over sites	
      end do
 
***** Thats all
      close(200)
      return
      end
 
CTITLE prescan_phs

      subroutine prescan_phs( fit_res, num_chan, data_flag_c, 
     .                        bf_type_c, biases,
     .                        max_res, max_tol, mean_tol, 
     .                        chan_with_max_res )

      implicit none

*     This routine will scan the fit_res from the phs_omc calculations
*     and decide if we need to add any bias flags based on the deviations
*     from a constant value.

* PASSED VARIABLES

* num_chan  - number of channels in this data
* data_flag_c(num_chan) - Data flag, bit 31 set to show bias flag
*   bf_type_cse(num_chan) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

* chan_with_max_res - channel with the maxium residual.

      integer*4 num_chan, data_flag_c(num_chan), 
     .          bf_type_c(num_chan), chan_with_max_res

* biases    - Logical to indicate that biases have been found

       logical biases

* fit_res   - Prefit residuals for each channel (if zero then no data or
*             bias flag already set.)
* max_tol - Tolerance for cchecking to see if bias should be added.
* mean_tol - Tolerance on mean value of offset (used if we are
*            reduced to one data point)
* max_res - Size of the largest residual from mean

       real*8 fit_res(num_chan), max_tol, mean_tol, max_res

* LOCAL VARIABLES

* j         - loop counter
* num_res   - number of values in mean compuation.

       integer*4 j, num_res

* kbit      - Tests to see if bit is set

       logical kbit

* mean_res  - Mean value of residuals

       real*8 mean_res


***** Compute the mean value not using any values with fit_res = 0 
*     indicates no data or bias flag already set on this data.
      biases = .false.
      mean_res = 0
      num_res  = 0
      do j = 1, num_chan
         if( fit_res(j).ne.0.0d0 ) then
             mean_res = mean_res + fit_res(j)
             num_res  = num_res  + 1
         end if
      end do

****  If there are no residuals just return
      if( num_res.eq.0 ) RETURN

      mean_res = mean_res / num_res
 
***** See if all of the residuals are about the same size.
      max_res = -1.0
      do j = 1, num_chan
         if( fit_res(j).ne.0.d0 ) then
             if( max_res.lt.abs(fit_res(j)-mean_res) ) then
                 max_res = abs(fit_res(j)-mean_res)
                 chan_with_max_res = j
              end if
         end if
      end do

*     If max_res is too large then set bias flag or if we
*     have just one residual and the mean to too far off.
      if( max_res.gt.max_tol .or.
     .   (abs(mean_res).gt.mean_tol) ) then
C    .   (num_res.eq.1 .and. abs(mean_res).gt.mean_tol) ) then
          j = chan_with_max_res
*         If there is no bias flag currently set, the set one
*         now
          if( .not.kbit(data_flag_c(j),32) .or.
     .        .not.kbit(data_flag_c(j),31) ) then	
              call sbit(data_flag_c(j),31,1)

*             Only set the jump bias flag if we more than data point.
*             (These data will be edited because of bias flags added).
              if( num_res.gt.1 )
     .        call sbit(bf_type_c(j), 2,1)
          end if
          fit_res(j) = 0.d0
          biases = .true.

*         if have just one data point, then return mean_res and
*         set biases false. (no more data to check)
          if( num_res.eq.1 ) then
C             biases  = .false.
              max_res = mean_res
          end if
      end if

****  Thats all
      return
      end

