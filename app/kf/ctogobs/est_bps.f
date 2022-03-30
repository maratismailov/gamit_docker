CTITLE EST_BPS
 
      subroutine est_bps( opass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse,
     .    params_cse, par_flag_cse  )
 
*     This is a new version of the final phase estimation that
*     excplicitly estimates the bias parameters.  The bias parameters
*     are set up for estimation and the clocks are solved impliticlity
*     at each epoch.  Once the bias parameters have been determined
*     the normal clock estimation run is made.
 
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
 
*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch.
 
 
      real*8 curr_cyc(2,max_gsvs, max_cfiles)
 
*   i,j       - Loop counters
*   ep        - Epoch number loop counter
*   most_np   - Maximum number of parameters expected based on all
*               sites observing all channels at the same time
 
      integer*4 i,j, ep, most_np
 
****  Initialize the estimation
      num_allp = num_param
      do i = 1, num_cfiles
          do j = 1, num_sat
              bp_sp(i,j)   = 0
              ep_bp(1,i,j) = 0
              ep_bp(2,i,j) = 0
              curr_cyc(1,j,i) = 0.d0
              curr_cyc(2,j,i) = 0.d0
          end do
      end do
      most_np = num_param + actual_max_chan*num_cfiles
      do i =1, most_np
         sol_eq(i) = 0.d0
         do j = 1, most_np
            norm_eq(i,j) = 0.d0
         end do
      end do 
* MOD TAH 160711: Fixed index from i,j to i,i
      do i = num_param+1, most_np
         norm_eq(i,i) = 1.d0
      end do
 
****  Loop over all of the data
      do ep = 1, num_ep
 
****      Initialize the clock parameter part of the solution
*         (Need to mod init_norm to use larger matrix).
          call init_norm(i, .true. )
 
*         Get the apriori values of the clocks for this epoch from
*         range solution values
          do j = 1, num_param
              apr_clk_val(j) = params_cse(j,ep)
          end do
 
 
****      Now loop over the sites and satellites to add the next
*         epoch.  First check to see if we have any bias parameters
*         that needed to be eliminated
 
          call bp_elim(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .            L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .            ctol_cse, data_flag_cse, bf_type_cse,
     .            params_cse, par_flag_cse, curr_cyc)
 
****      Now increment the normal equations for new clock estimation
*         and increment to the bias parameters
 
          call incr_bp(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .            L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .            ctol_cse, data_flag_cse, bf_type_cse,
     .            params_cse, par_flag_cse, curr_cyc)
 
****      Now do the implicit solution for the clocks
          call implicit_clk( ep )
 
      end do
 
****  Now that we have accumulated all the epochs, resolve any of the
*     remaining biases.
      call fin_bp(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .            L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .            ctol_cse, data_flag_cse, bf_type_cse,
     .            params_cse, par_flag_cse, curr_cyc)
 
****  Thats all
      return
      end
 
CTITLE IMPLICIT_CLK
 
      subroutine implicit_clk( ep )
 
*     Routine a make an implicit solution for the clocks similiar
*     to the way solve handles bias parameters.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
* ep  -- Epoch number (for debug output)

      integer*4 ep
 
* LOCAL VARIABLES
 
* i,j,k -- Loop variabales

      integer*4 i,j,k
 
****  Loop over all the clock parameters implicitly solving each
*     one.  Scale the diagonal to aviod singularity problems since
*     all clocks and biases are estimated.
      do k = 1, num_allp
         norm_eq(k,k) = norm_eq(k,k) + norm_eq(k,k)*1.d-3
      end do
      
*     Now do the implicit solution.
      do k = 1, num_param
 
 
*         Now reomve the effect of this clock from all remaining
*         parameters
          if( norm_eq(k,k).gt. 0.d0 ) then

*             Now remove the clock contributions from the bp_sol
*             array (which saves the effects of a +1 cycle change to
*             solution vector)

*             Now do the implicit clock solution for solution itself.
              do j = k+1, num_allp
                  sol_eq(j) = sol_eq(j) - sol_eq(k)*norm_eq(k,j)/
     .                                    norm_eq(k,k)
                  do i = k+1, num_allp
                      norm_eq(i,j) = norm_eq(i,j) -
     .                    norm_eq(i,k)*norm_eq(k,j)/norm_eq(k,k)
                  enddo
              end do
              

          end if
 
*         Now remove clear the row and column for parameter k
          sol_eq(k) = 0.d0
          do j = 1,num_allp
              norm_eq(j,k) = 0.d0
              norm_eq(k,j) = 0.d0
          end do
      end do

****  Thats all
      return
      end
 
 
CTITLE bp_elim
 
      subroutine bp_elim(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse,
     .        params_cse, par_flag_cse, curr_cyc)
 
*     This routine will look at the bias flags at the current epoch
*     and see if there are any new ones.  If there are it will work
*     out which of the old ones to solve and apply to the phase data
*     The new parameter number slot is computed based on replacing
*     the bias parameter in the same channel or on replacing the oldest
*     one available.
 
*     This routine assumes that the first valid data point is marked
*     with a bias flag.
 
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
*   ep       -- Epoch number
 
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass, ep
 
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
 
*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch.
 
 
 
      real*8 curr_cyc(2,max_gsvs, max_cfiles)
 
* LOCAL VARIABLES
 
*         i,j,k - Loop counters
*       lv      - List number of satellite
*       pn      - Generic parameter number for new parameters
*       oldest_bp   - Epoch number of the oldest bias parameter
*               - being estimated (ie. smallest epoch number).
*               - Used to decide which parameters to eliminate
*       oldest_ch   - Satellite list entry for the the oldest bias
*       pivot(max_neq)  - Pivot elements matrix inversion
*       nc   - Count of number of channels in use
 
      integer*4 i,j,k, lv, pn, oldest_bp, oldest_ch,
     .          nc
 
*   dL1c, dL2c       - Changes to the number of cycles based on
*                   - bias parameter estimation.
*   est             - Estimatof bias from normal equations.
 
      real*8 dL1c, dL2c, est

      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

*  Variables needed for forcing bias parameters

*    found       - Controls loops while searching for entries.
*    kbit  -- Function to check bit setting
 
      logical found, kbit
 
*       good_bf - Logical function to return true is bias flag
*               - on good data point.
 
 
      logical good_bf


****  Loop over the sites and channels at this epoch to see
*     if we need any new bias parameters
      bp_elim_needed = .false.
 
      do i = 1, num_cfiles
 
*         Set all eliminate bias parameters bits to zero (i.e. so
*         that we don't eliminate at end of loop.)  (The 4 value is
*         more than we need.  This would cover 4*32 satellites.)
          do j = 1, 4
              elim_bp(j,i) = 0
          end do
 
          do j = 1, actual_max_chan
 
*             See if we have a bias flag at this epoch for this site
*             and channel.  If we do then mark previous bias parameter
*             on this PRN for elimination or find the oldest bias one.
              if( good_bf(data_flag_cse(j,i,ep),0, phs_mask) ) then
 
*                 OK, see if have an active bias flag for this PRN
                  lv = ctol_cse(j,i,ep)
                  if( bp_sp(i,lv).gt.0 ) then
 
*                     mark bias for elimination.
                      call sbit(elim_bp(1,i),lv,1)
                      pn = bp_sp(i,lv)
                      bp_elim_needed = .true.
 
                  else
 
*                     Current bias parameter is not active so we need
*                     to find the oldest parameter and then mark this
*                     for elimination and set up the new parameter numbers

*                     Check to see if we are using all the channels 
*                     available for this site
                      nc = 0
                      do k = 1, num_sat
                          if( bp_sp(i,k).ne.0 ) nc = nc + 1
                      end do

                      oldest_ch = 0
C                     if( nc.ge.actual_max_chan ) then 
                      if( nc.ge.max_gchannels   ) then 
                          oldest_bp = 100000
                          do k = 1, num_sat

*                             Make sure that we are using the parameter
*                             that is older than the current and that it
*                             is not already marked for elimination.
                              if( ep_bp(2,i,k).gt.0 .and.
     .                            ep_bp(2,i,k).lt.oldest_bp .and.
     .                            .not. kbit(elim_bp(1,i),k) ) then

                                  oldest_bp = ep_bp(2,i,k)
                                  oldest_ch = k
                              end if
                          end do
                      end if 
*                     If we did not find any then we may not have had any
*                     bias parameters at this site yet, so assign a
*                     a parameter number.
                      if( oldest_ch.eq.0 ) then
 
*                         See what channels are available
                          found = .false.
                          k = 0
                          do while ( .not.found .and. k.le.num_sat )
                              k = k + 1
                              if( bp_sp(i,k).eq.0 ) then
                                  oldest_ch = k
                                  found = .true.
                              end if
                          end do
 
*                         Find a parameter number to assign this channel
*                         to.
                          num_allp = num_allp + 1
                          if( num_allp.gt. num_cfiles*
C    .                        (actual_max_chan+1) + num_sat) then
     .                        (max_gchannels+1) + num_sat) then
                              write(*,200) ep, i, j, k
 200                          format('*** DISASTER *** Too many ',
     .                            'bias parameters at epoch ',i4,
     .                            ' Site ',i3,' Channel ',i3,i4)
C                             stop 'Too many bias parameters'
                          end if
                          pn = num_allp
 
                      else
 
*                         We found an old channel that can be used
*                         so mark for elimination
                          call sbit(elim_bp(1,i),oldest_ch,1)

*                         Get the parameter number so that it can be
*                         assigned to the next bias parameter on
*                         a different site and make this one negative
*                         It will be set to zero below to show that
*                         it is not in use (compared to case where
*                         we use same satellite for parameter)
                          pn = bp_sp(i, oldest_ch)
                          bp_sp(i, oldest_ch) = -pn
                          bp_elim_needed = .true.
                      end if
                  end if
 
                  bp_sp(i,lv) = pn
                  sp_bp(1,pn) = i
                  sp_bp(2,pn) = lv
 
*                 We only set if not used for the case when 
*                 we have the same satellite on this parameter estimate.
*                 (i.e., this epoch number will be needed below to
*                 update data).
                  if( ep_bp(1,i,lv).eq.0 ) ep_bp(1,i,lv) = ep

              end if
*         Finish looping over channels for this cfile
          end do
*     Finishing looping over all stations.
      end do
 
****  Now see if we need to eliminate any bias parameters
      if( bp_elim_needed ) then
 
*         Now we need to invert the system to get the bias parameter
*         estimates.  Invert only the bias parameter of the matrix
*         We do this by calling with the first element that referrs
*         to bias parameters and number of parameters equal to only
*         bias parameters.  Num_param is the number of clock parameters
*         So we call with element after this.  This should speed up
*         inverse and clocks have been pre-eliminated so there is no
*         correlation with these parameters.

          do i = 1, num_cfiles
              do j = 1, num_sat
                  if( kbit(elim_bp(1,i),j) ) then
 
*                     This bias parameter marked for elimination.
                      lv = j
                      pn = abs(bp_sp(i,lv))

*                     Get estimate by solving only correlated portion
*                     of the normal equations.                     
                      call get_est ( pn, est )
                      if( opass.le.pc_non_int ) then
* MOD TAH 180320: Remap the frequnecies (for Glonass)
                         fr1 = fL1u(lv)/fL1(lv)
                         dL1c = -nint(est*fr1)/fr1
                      else
                         dL1c = -est
                      end if
                     
                      dL2c = dL1c

*                     Only apply the corrections if they are
*                     non-zero.
                      if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data
                          call update_cyc(dL1c, dL2c, i, lv,
     .                         ep_bp(1,i,lv), ep_bp(2,i,lv), opass,
     .                         L1_cyc_cse(1,1,1),
     .                         L2_cyc_cse(1,1,1),
     .                         ctol_cse(1,1,1),data_flag_cse, 'N')
                          curr_cyc(1,lv,i) = curr_cyc(1,lv,i) + dL1c
                          curr_cyc(2,lv,i) = curr_cyc(2,lv,i) + dL2c
                      end if 

*                     See if we are keeping this bias parameter or
*                     not.
                      if( bp_sp(i,lv).gt.0 ) then 
*                         Set the beginning bias epoch number to which 
*                         this bias estimate applies.
                          ep_bp(1,i,lv) = ep
                      else
*                         These two sets makes the parameter in-active
*                         again.
                          bp_sp(i,lv) = 0
                          ep_bp(1,i,lv) = 0
                      end if

*                     Now set the end parameter for this bias to
*                     be zero so that we know it is not active
                      ep_bp(2,i,lv) = 0
 
*                     Now do a regorous implicit solution for the bias.
                      call implicit_bias(pn)

                  end if
              end do
          end do
 
      end if
 
****  OK, that should be all.  Return to main subroutine and increment
*     in new data
      return
      end
 
CTITLE INCR_BP
 
      subroutine incr_bp(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse,
     .        params_cse, par_flag_cse, curr_cyc)
 
*     This routine will take the phase data at the current epoch
*     and increment the normal equations for clock and bias parameters
*
 
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
*   ep       -- Epoch number
 
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass, ep
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L1 phase residuals
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
 
*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch.
 
 
 
      real*8 curr_cyc(2,max_gsvs, max_cfiles)
 
* LOCAL VARIABLES
 
*         i,j       - Loop coounters
*   lv          - satellite list number
*   ns, pn      - Parameter numbers
 
      integer*4 i,j, lv, ns, pn
 
*      data_var - Phase variance assigned to phase data
*   res         - Phase LC residual
*   phs_omc     - Function to compute Phase o-minus-c value.
 
      real*8 data_var, res, phs_omc

*   data_OK    - Logical function to check data
      logical data_OK
 
*
***** Loop over the data incrementing the normal equations
      do i = 1,num_cfiles
 
          data_var = phs_noise(i)**2
 
*         Now compute the next iteration on the clock values
          do j = 1, actual_max_chan
                        
****          if pass one then update number of cycles
              lv = ctol_cse(j,i,ep)
              L1_cyc_cse(j,i,ep) = L1_cyc_cse(j,i,ep) +
     .                             curr_cyc(1,lv,i)
              L2_cyc_cse(j,i,ep) = L2_cyc_cse(j,i,ep) +
     .                             curr_cyc(2,lv,i)
 
*             Use only data that has no error and for which
*             we already know the bias parameter values
              if( data_ok(data_flag_cse(j,i,ep),0, phs_mask) ) then
 
*                 Compute OminusC
                  ns = num_cfiles + lv
                  res = phs_omc(L1r_phs_cse(j,i,ep),L2r_phs_cse(j,i,ep),
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
 
*                 Now get the bias parameter number and increment these
*                 entries as well
                  pn = bp_sp(i,lv)
                  if( pn.le.0 ) then
                      write(*,300) ep, i, lv, pn, data_flag_cse(j,i,ep)
 300                  format('***DISASTER*** Bias parameter at epoch ',
     .                        I4,' with no parameter number',3i4,
     .                        o15)
                      stop 'No bias parameter number'
                  end if
                  norm_eq(pn,pn) = norm_eq(pn,pn) + 1.d0/data_var
                  norm_eq(pn,i)  = norm_eq(pn,i)  + 1.d0/data_var
                  norm_eq(pn,ns) = norm_eq(pn,ns) - 1.d0/data_var
                  norm_eq(i,pn)  = norm_eq(i,pn)  + 1.d0/data_var
                  norm_eq(ns,pn) = norm_eq(ns,pn) - 1.d0/data_var
                  sol_eq(pn)     = sol_eq(pn)     + res/data_var

*                 Set the end epoch number for this bias parameter
                  ep_bp(2,i,lv) = ep
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE FIN_BP
 
      subroutine fin_bp(opass, ep, L1r_phs_cse, L2r_phs_cse,
     .        L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse,
     .        params_cse, par_flag_cse, curr_cyc)
 
*     This routine finishes up the bias parameter estimation and
*     eliminates the biases for all remaining bias parameters
*     after we have finished going through all the epochs.
*
 
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
*   ep       -- Epoch number
 
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), opass, ep
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L1 phase residuals
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
 
*  curr_cyc(2,max_gsvs, max_cfiles) - Current values for the biases
*                 for L1 and L2 by satellite and site.  These are
*                 propogated forward at each epoch.
 
 
 
      real*8 curr_cyc(2,max_gsvs, max_cfiles)
 
* LOCAL VARIABLES
 
*   i,j       - Loop coounters
*   lv          - satellite list number
*   pn      - Parameter numbers
*   pivot(max_neq)  - Pivot elements matrix inversion
 
      integer*4 i,j, lv, pn, pivot(max_neq)
 
*   scale(max_neq)   - Scale values matrix inversion
*   dL1C, dL2C       - Changes to the number of cycles based on
*                   - bias parameter estimation.
 
      real*8 scale(max_neq), dL1c, dL2c

      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.
 
****  Eliminate the remaining bias parameters.
 
      write(*,300)
 300  format('BPE: Eliminating final bias parameters ')
 
 
*     Now we need to invert the system to get the bias parameter
*     estimates.  Invert only the bias parameter of the matrix
*     We do this by calling with the first element that referrs
*     to bias parameters and number of parameters equal to only
*     bias parameters.  Num_param is the number of clock parameters
*     So we call with element after this.  This should speed up
*     inverse and clocks have been pre-eliminated so there is no
*     correlation with these parameters.
      pn = num_param + 1
      call invert_vis(norm_eq(pn,pn), sol_eq(pn),
     .                 scale, pivot, num_allp-num_param, max_neq,1)
 
*     OK, now find the entries that are marked for elimination.
      do i = 1, num_cfiles
          do j = 1, num_sat
              if( ep_bp(2,i,j).gt.0 ) then
 
*                 This bias parameter marked for elimination.
                  lv = j
                  pn = bp_sp(i,lv)
                  if( opass.le.pc_non_int ) then
                      fr1 = fL1u(lv)/fL1(lv)
                      dL1c = -nint(sol_eq(pn)*fr1)/fr1
                  else
                      dL1c = -sol_eq(pn)
                  end if
                  dL2c = dL1c
 
*                 Only apply the corrections if they are
*                 non-zero.
                  if( dL1c.ne.0 .and. dL2c.ne.0 ) then

* MOD TAH 000310: Passed option not to edit data

                      call update_cyc(dL1c, dL2c, i, lv,
     .                    ep_bp(1,i,lv), ep_bp(2,i,lv), opass,
     .                    L1_cyc_cse(1,1,1),
     .                    L2_cyc_cse(1,1,1),
     .                    ctol_cse(1,1,1),data_flag_cse, 'N')
                      curr_cyc(1,lv,i) = curr_cyc(1,lv,i) + dL1c
                      curr_cyc(2,lv,i) = curr_cyc(2,lv,i) + dL2c
                  end if
 
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE GET_EST

      subroutine get_est( pn, est) 
      
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
      
*   pn    -- Parameter to extract

      integer*4 pn
      real*8 est

*   work  -- Work array to save a row or column in
*   lrho  -- is smallest correlation, only needed if we find too
*            many large correlations
*   rho   -- Is estimate of correlation.

      real*8 bp_sol(2*max_gchannels), 
     .       bp_mat(2*max_gchannels, 2*max_gchannels), 
     .       work(2*max_gchannels), lrho, rho

*   lrin -- Index that points to entry in cpel with the lowest
*           correlation.
     
      integer*4 i,j, cpel(2*max_gchannels), nc, 
     .          pivot(2*max_gchannels), lrin
      
      
      
*     Pull off the appropriate peices of the norm-eq.  Build up the list 
*     entries that are correlated.  We get all entries greater than 0.1
*     unless we run out of array space in which case we take the largest.
*     (Tests on data suggests that 2* the maximum number of channels
*      should not be exceeded by the 0.1 correlation values.)
      nc = 1
      cpel(nc) = pn
      do i = num_param+1, num_allp
         rho = abs(norm_eq(i,pn)/
     .         sqrt(norm_eq(i,i)*norm_eq(pn,pn)))
         if( i.ne.pn .and. rho.gt.pc_max_corr ) then
             nc = nc + 1
             if( nc.le.2*max_gchannels ) then
                cpel(nc) = i
             else
*               Find lowest correlation and replace with this value
                lrho = rho
                lrin = 0
                do j = 2, nc-1
                   if( abs( norm_eq(cpel(j),cpel(1))/
     .                 sqrt(norm_eq(cpel(j),cpel(j))*
     .                      norm_eq(cpel(1),cpel(1)))).lt. lrho ) then
                      lrho =  abs( norm_eq(cpel(j),cpel(1))/
     .                         sqrt(norm_eq(cpel(j),cpel(j))*
     .                              norm_eq(cpel(1),cpel(1))) ) 
                      lrin = j
                   end if
                end do 
                nc = nc - 1
*               Only replace if we find a correlation less than the
*               one for the current element.                
                if( lrin.gt.0 ) cpel(lrin) = i
            end if                     
         end if
      end do
      
****  Now copy the values accross
      do i = 1, nc
         bp_sol(i) = sol_eq(cpel(i))
         do j = 1, nc
            bp_mat(i,j) = norm_eq(cpel(i),cpel(j))
         end do
      end do

*     Get the solution                  
      call invert_vis(bp_mat, bp_sol, work, pivot, nc, 
     .                2*max_gchannels, 1)
     
*     Now save the estimate
      est = bp_sol(1)
      
      return
      end
      
CTITLE IMPLICIT_BIAS
 
      subroutine implicit_bias( pn )
 
*     Routine a make an implicit solution for the bias parameters
*     similar to the way solve handles bias parameters.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
* pn  -- Parameter number to be solved

      integer*4 pn
 
* LOCAL VARIABLES
 
* i,j,k -- Loop variabales

      integer*4 i,j

*     Do the implicit solutions
      do j = num_param+1, num_allp
          if( j.ne.pn ) then 
             sol_eq(j) = sol_eq(j) - sol_eq(pn)*norm_eq(pn,j)/
     .                         (norm_eq(pn,pn)*(1.d0+1.d-3))
             do i = num_param+1, num_allp
                if( i.ne.pn ) then
                    norm_eq(i,j) = norm_eq(i,j) -
     .                  norm_eq(i,pn)*norm_eq(pn,j)/
     .                    (norm_eq(pn,pn)*(1.d0+1.d-3))
                end if
             end do
          end if
      end do              

 
*     Now remove clear the row and column for parameter pn
      sol_eq(pn) = 0.d0
      do j = num_param+1,num_allp
         norm_eq(j,pn) = 0.d0
         norm_eq(pn,j) = 0.d0
      end do
      norm_eq(pn,pn) = 1.d0

****  Thats all
      return
      end
 
