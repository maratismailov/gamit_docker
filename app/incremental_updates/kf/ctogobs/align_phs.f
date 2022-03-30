CTITLE ALIGN_PHS
 
      subroutine align_phs(opass, L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )

      implicit none
 
*     This routine will do "rough" alignment of the phase values
*     assuming the bias flags have already been set.  The algorithm
*     used to compute the wide-lane and the mean difference between
*     phase and range estimates of the clocks.  These values are used 
*     to get over-all estimates of the numbers of L1 and L2 cycles.
 
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
 
*   i,j        - Epoch loop counter: site/satellite
*   ep         - Epocn number
*   ns         - Parameter number for current satellite clock
*   nums_lc, nums_wl -- Number of values in LC and WL estimate
*   start_bf   - First epoch of current bias flag being processed.
*   ch, ltoc   - Channel number and function to return channel number
*                for specific satellite list number
  
      integer*4 i, j, lv, ep, ns, nums_lc, nums_wl, start_bf, ch, ltoc

*   sums_lc, sums_wl -- Sum of LC residual and Widelanes
*   vars_lc, vars_wl -- Sum of LC and WL resiudals squared
*   dL1, dL2  -- Adjustments to L1 and L2 cycles
*   res, WL   -- LC residual and Widelane estimate
*   phs_omc   -- Function to return value of phase omc

      real*8 sums_lc, sums_wl, vars_lc, vars_wl, dL1, dL2, res, WL,
     .       phs_omc

*   data_ok     - Logical function returns true if data is good
*   good_bf     - Logical function returns true if Bias flag is good.
*   first_epoch - True is still looking for first measurement
 
      logical data_OK, good_bf, first_epoch

*   edit_opt  - string set to E if the number of data in the mean is less
*     than the minimum number of values between biases
      character*1 edit_opt 

 
****  set the data mask to be used for editing the phase data.  Two masks
*     are used.  The phs_bias_mask check the bias flags as well as
*     editing, the simple phs_mask just checks to see that the data
*     are OK. 
      call set_phs_mask( phs_mask, phs_bias_mask )

***** Set up the initial use of data flags
      do i = 1, num_cfiles
          do j = 1, num_sat
              nums_lc = 0
              nums_wl = 0
              sums_lc = 0.d0
              sums_wl = 0.d0
              vars_lc = 0.d0
              vars_wl = 0.d0
              first_epoch = .true.
              start_bf = 1 
              ns = num_cfiles + j
 
*             Start looping over all epochs of data
              do ep = 1, num_ep

*                Accumulate the statistics to find means and
*                rms.  When a bias flag is found, compute the
*                results and update the number of cyles
                 ch = ltoc( ctol_cse(1,i,ep), j, actual_max_chan)
                 if( ch.gt.0 ) then

*                   Satellite j is observed at this time, see if
*                   data is OK (including presence of bias flag.
*                   If the latter case, then update the cycles before
*                   proceeding)
                    if( data_OK(data_flag_cse(ch,i,ep),
     .                                        0,phs_mask) ) then

*                       See if this a bias flag
                        if( good_bf(data_flag_cse(ch,i,ep),
     .                                            0,phs_mask) ) then
                            if( first_epoch ) then
                                first_epoch = .false.
                            else
*                               OK, we have our first block of data
*                               to process.  Get the update to number
*                               cycles
                                edit_opt = 'N'
                                if( nums_lc.lt.min_good_bias ) 
     .                                                   edit_opt = 'E'
                                call get_lcwl_cyc(nums_lc, nums_wl,
     .                               sums_lc, sums_wl, vars_lc, vars_wl,
     .                               i, j, ch, ep, dL1, dL2, edit_opt)
                                call update_cyc(dL1, dL2, i, j, 
     .                               start_bf, ep-1, opass,
     .                               L1_cyc_cse(1,1,1),
     .                               L2_cyc_cse(1,1,1),
     .                               ctol_cse,data_flag_cse, edit_opt)
*                               OK, now save the current epoch
                                start_bf = ep
*                               clear the accumulation variables
                                nums_lc = 0
                                nums_wl = 0
                                sums_lc = 0.d0
                                sums_wl = 0.d0
                                vars_lc = 0.d0
                                vars_wl = 0.d0
                            end if
                        end if

*                       OK, data is good and does not have a bias 
*                       flag.  Compute the WL and LC residuals and
*                       accumulate result  
cd                     if( i.eq.2.and.j.eq.9) then         
cd                        print *,'ALIGN_PHS i j ep ch  ns '
cd     .                     ,' l1r_phs_cse  l2r_cyc_cse '
cd     .                     ,  'l1_cyc_cse l2_cyc_cse params_cse fl1 fl2'
cd     .                     ,  i,j,ep,ch,ns
cd     .                     ,  l1r_phs_cse(ch,i,ep),l2r_phs_cse(ch,i,ep)
cd     .                     ,  l1_cyc_cse(ch,i,ep),l2_cyc_cse(ch,i,ep)
cd     .                     ,  params_cse(i,ep),params_cse(ns,ep)
cd     .                     ,  fl1(j),fl2(j)
cd                      endif
                        res = phs_omc(L1r_phs_cse(ch,i,ep),
     .                        L2r_phs_cse(ch,i,ep),
     .                        L1_cyc_cse(ch,i,ep), L2_cyc_cse(ch,i,ep),
     .                        params_cse(i,ep),params_cse(ns,ep),
     .                        fL1(j), fL2(j) )         
cd               if( i.eq.2.and.j.eq.9 ) print *,'ep i j res ',ep,i,j,res
                                
*                       Test to see if we have L2 range for computing
*                       widelane.
                        if( lambda(ch,4,i).ne.0 ) then
                            WL = -(L1r_phs_cse(ch,i,ep) +
     .                             L1_cyc_cse(ch,i,ep))   +
     .                            (L2r_phs_cse(ch,i,ep) + 
     .                             L2_cyc_cse(ch,i,ep))    + 
     .                             dfsf(j)*(L1r_rng_cse(ch,i,ep) +
     .                                   L2r_rng_cse(ch,i,ep))
                        else
                           WL = 0
                        end if

*                       Now accumlate the statics
                        nums_lc = nums_lc + 1
                        nums_wl = nums_wl + 1    
cd        if( i.eq.2.and.j.eq.9 ) print *,'ep res sums_lc ',ep,res,sums_lc
                        sums_lc = sums_lc + res
                        sums_wl = sums_wl + WL
                        vars_lc = vars_lc + res**2
                        vars_wl = vars_wl + WL**2
                    end if
                end if
            end do 

*           Now finish up any biases that have been left unresolved. 
            edit_opt = 'N'
            if( nums_lc.lt.min_good_bias )  edit_opt = 'E'

            call get_lcwl_cyc(nums_lc, nums_wl,
     .                        sums_lc, sums_wl, vars_lc, vars_wl,
     .                        i, j, ch, num_ep, dL1, dL2, edit_opt)
            call update_cyc(dL1, dL2, i, j, 
     .                      start_bf, num_ep, opass,
     .                      L1_cyc_cse(1,1,1),
     .                      L2_cyc_cse(1,1,1),
     .                      ctol_cse(1,1,1),data_flag_cse,edit_opt )


         end do
      end do

****  Thats all
      return
      end

CTITLE GET_LCWL_CYC

      subroutine  get_lcwl_cyc(nums_lc, nums_wl,
     .                         sums_lc, sums_wl, vars_lc, vars_wl,
     .                         i, j, ch, ep, dL1, dL2, edit_opt )

      implicit none

*     Routine to compute the number of cycles needed to align the
*     widelanes and make the LC residuals to the range clock close
*     to zero. 

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED  VARIABLES
 
*   i,j        - Site number and satellite list number
*   ch         - Channel number
*   ep         - Epocn number
*   nums_lc, nums_wl -- Number of values in LC and WL estimate
  
      integer*4 i, j, ch, ep, nums_lc, nums_wl

*   sums_lc, sums_wl -- Sum of LC residual and Widelanes
*   vars_lc, vars_wl -- Sum of LC and WL resiudals squared
*   dL1, dL2  -- Adjustments to L1 and L2 cycles

      real*8 sums_lc, sums_wl, vars_lc, vars_wl, dL1, dL2

* LOCAL VARIABLES
* mean_dlc, mean_wl -- Mean values for phase difference and WL
* rms_dlc, rms_wl   -- RMS scatters

      real*8 mean_dlc, mean_wl, rms_dlc, rms_wl

* edit_opt  -- Editing option:  If set to E then data is edited at the
*     same time the cycles are applied.  This is used in align_phs to
*     same pieces of data early in the processing

      character*1 edit_opt

                          
cd       if( i.eq.2.and.j.eq.9 ) 
cd     .  print *,'GET_LCWL_CYC i j ep code prn sums_lc nums_lc '
cd     .         , i,j,ep,cf_codes(i),prn_list(j),sums_lc,nums_lc  
****  OK, make sure we have data
      if( nums_lc.gt.0 ) then
          mean_dlc = sums_lc/nums_lc
          mean_wl  = sums_wl/nums_wl
          rms_dlc = sqrt(abs(vars_lc/nums_lc - mean_dlc**2))
          rms_wl = sqrt(abs(vars_wl/nums_wl - mean_wl**2))

****      Now compute the change in cycles needed  
C         dL1 = -nint((mean_dlc - lcf2(j)*mean_wl)/(lcf1(j)+lcf2(j)))
C         dL2 = -nint(mean_wl-dL1)
* MOD TAH 200612: Account for frequency effects with Glonass
****      Now compute the change in cycles needed  
          dL1 = -(mean_dlc - lcf2(j)*mean_wl)/(lcf1(j)+lcf2(j))
          dL1 = nint(dL1*fL1u(j)/fL1(j))*(fL1(j)/fL1u(j))
          dL2 = -(mean_wl-dL1)
          dL2 = nint(dL2*fL2u(j)/fL2(j))*(fL2(j)/fL2u(j))
      else
          dL1 = 0
          dL2 = 0
      end if
      
      if( edit_opt.ne.'E' )
     .write(*,200) cf_codes(i), prn_list(j), ep, nums_lc, mean_dlc,
     .             mean_wl, rms_dlc, rms_wl, dL1, dL2, edit_opt
 200  format('PRESCAN ',a4,' PRN',i2.2,' EP ',i5,' # ',I6,' Means ',
     .       2F18.2,' RMS ',2F8.2,' dL1/2 ',2F12.2,1x,a)

****  Thats all
      return
      end

