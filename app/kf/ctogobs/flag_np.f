CTITLE FLAG_NP       
 
      subroutine flag_np(L1r_phs_cse, L2r_phs_cse, 
     .              L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse,
     .              ctol_cse, data_flag_cse,
     .              params_cse, par_flag_cse )

      implicit none
 
*     This routine will scan the data distribution and flag all
*     data that can not be used in a normal point.  (All data
*     for the normal duration must be available for it not to
*     get flagged.
 
* INCLUDES FILES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
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

*  i, j    - Loop counters
*  ltoc    - FUnction to return channels number given satellite list number
*  cf      - Cfile number
*  ls      - satellite number
*  ch      - Channel number
*  ep      - Normal point epoch counter 
*  phs_mask_fel  - Phase editing mask accounting for the data below the
*            final cutoff . Bias flags not checked
*  phs_bias_mask_fel - Final elevation checked plus bias flags checked.
*  num_np_good(max_cfiles)  - Number of good normal points
*  num_np_bad(max_cfiles)   - number of bad normal points
*  first_ep, end_ep, tot_npep - First epoch to use in forming
*     normal, the epoch at the start of the last complete
*     normal point, and the total number of normal points
*     that can be formed.
*  np_ep   - Epoch that the current normal point will refer to
 
      integer*4 i,j,k, ltoc, cf, ls, ch, ep, phs_mask_fel, 
     .          phs_bias_mask_fel,
     .          num_np_good(max_cfiles) , num_np_bad(max_cfiles),
     .          first_ep, end_ep, tot_npep, np_ep

* sum_omc(4,max_gsvs)  - Sum of the omc values corrected for clocks
*     for L1, L2, P1 and P2.
* rcv_clk, svs_clk     - Temporary storage of the rcv and satellite
*     clocks

       real*8 sum_omc(4,max_gsvs), rcv_clk, svs_clk

* push_bias(max_gsvs)   - Logical which indicates that we have to
*     push the bias flag because the data with bias flag could
*     no be used in the normal.
* good_bf               - Logical function that says there is a
*     bias flag on a good data point
* data_OK               - Logical function which indicates that
*     data ppoint is good.
* NP_OK(max_gsvs)       - Logical to indicate that all data in normal
*     point is good.
* NP_possible(max_gsvs) - Logocal to indicate that a normal might be
*     able to be formed if all data is available.
* fnd_this_epoch(max_gsvs) - Logical which checks after the first epoch
*     of data in a normal point, that subsequent epochs are found.  This
*     handles the case of no data being collected on a satellite after
*     first epoch.

      logical push_bias(max_gsvs), good_bf, data_OK, NP_OK(max_gsvs),
     .        NP_possible(max_gsvs), fnd_this_epoch(max_gsvs)

***** Set the phase masks

      call set_phs_mask(phs_mask, phs_bias_mask ) 

*     Now set the mask to account for points below the final cutoff
      phs_mask_fel = phs_mask
      call sbit(phs_mask_fel,27,1)
      phs_bias_mask_fel = phs_bias_mask
      call sbit(phs_bias_mask_fel,27,1)

*     Initial counters for number of data
      do i = 1, num_cfiles
         num_np_good(i) = 0
         num_np_bad(i)  = 0
      end do


*     Compute the epoch ranges for normal points
      first_ep = abs(np_start)
      tot_npep = int((num_ep-first_ep)/np_size)
      end_ep   = first_ep + (tot_npep-1)*np_size

      write(*,120) np_size, first_ep, end_ep, tot_npep
 120  format('Forming ',i3,' epoch normal points from EP ',i4,
     .       ' to ',i4,' Total: ',i4)

****  Loop over the cfiles. 
      do cf = 1, num_cfiles

****      Initialize the push biases.  Since there
*         will always be a bias on the first point, set the push_bias
*         true.  Make the logic in getting the first point to be used
*         in normal pointing somewhat simpler.
          do i = 1, num_sat
             push_bias(i) = .true.
          end do

****      First loop up the last epoch before the point we form normal
*         points
          do i = 1, first_ep-1
*            Scan the channels
             do ch = 1, actual_max_chan

*               See if the data is OK
                ls = ctol_cse(ch,cf,i)
                if( ls.gt.0 ) then

*                   If point is marked good, then set the not
*                   in normal point flag
                    if( data_OK(data_flag_cse(ch,cf,i),0,
     .                          phs_bias_mask_fel) )  then
                        call sbit(data_flag_cse(ch,cf,i),7,1)
                    end if
                end if
             end do
          end do

****      Now loop over the epochs of data (Currently set to start at
*         first epoch
          do ep = first_ep, end_ep, np_size
             
*            Set NP_OK to indicate initially that the normal point is OK
             do j = 1, num_sat
                NP_OK(j) = .true.
                NP_possible(j) = .false.

*               Clear the summation variables for the normal points
                do k = 1,4
                   sum_omc(k,j) = 0.d0
                end do
             end do

*            Loop over the epochs within the normal point
             do j = 1, np_size

*               Set that we have found no data yet.  Set true when
*               data actually found
                do k = 1, num_sat
                   fnd_this_epoch(k) = .false.
                end do
                i = ep + j - 1

*               Scan the channels
                do ch = 1, actual_max_chan

*                  See if the data is OK
                   ls = ctol_cse(ch,cf,i)
                   if( ls.gt.0 ) then   

*                      We have data on satellite ls, see if it is OK
*                      (use the mask accounting for final cutoff)
                       if( .not.data_OK(data_flag_cse(ch,cf,i),0,
     .                        phs_mask_fel) ) then
                           NP_OK(ls) = .false.
                       else
*                          OK, we have some good data so say np is 
*                          possible
*                          Set that we found data
                           fnd_this_epoch(ls) = .true.
                           NP_possible(ls) = .true.
                       end if

*                      if this is after the first point see if
*                      there is a bias flag
                       if( good_bf(data_flag_cse(ch,cf,i),0,
     .                             phs_mask) ) then
                           push_bias(ls) = .true.

*                          If the bias flag is not on the first
*                          first point, then mark the Normal point
*                          as bad.
                           if( j.gt.1 ) NP_OK(ls) = .false.                             
                           
                       end if

*                      Now sum the values for the normal point.  These
*                      are sums of omc - (rcv_clk-svs_clk) 
                       rcv_clk = params_cse(cf,i)
                       svs_clk = params_cse(num_cfiles+ls,i)
                       
                       sum_omc(1,ls) = sum_omc(1,ls) + 
     .                     (L1r_phs_cse(ch,cf,i)+L1_cyc_cse(ch,cf,i)) -
     .                     (rcv_clk-svs_clk)
                       sum_omc(2,ls) = sum_omc(2,ls) + 
     .                     (L2r_phs_cse(ch,cf,i)+L2_cyc_cse(ch,cf,i)) -
     .                     (rcv_clk-svs_clk)*fL2(ls)/fL1(ls)
                       sum_omc(3,ls) = sum_omc(3,ls) + 
     .                      L1r_rng_cse(ch,cf,i)     -
     .                     (rcv_clk-svs_clk)
                       sum_omc(4,ls) = sum_omc(4,ls) + 
     .                      L2r_rng_cse(ch,cf,i)     -
     .                     (rcv_clk-svs_clk)*fL2(ls)/fL1(ls)

                   end if
                end do

*               Now if we did not find data on a satellite, set the 
*               NP_OK status to .false.
                do k = 1, num_sat
                   if( .not.fnd_this_epoch(k)) NP_OK(k) = .false.
                end do
*                       ! Looping over the epochs for this normal point
             end do

****         Now loop over the satelites and see of the normal points 
*            are OK
             do ls = 1, num_sat

                if( .not.np_ok(ls) ) then

*                   We need to flag the data as BAD.
                    do j = 1, np_size
                       i = ep + j - 1

*                      Get the channel number
                       ch = ltoc(ctol_cse(1,cf,i), ls, 
     .                           actual_max_chan)        
*                      Now we may have not observed this satellite
*                      so ch may be zero.
                       if( ch.gt.0 ) then
                           if( data_OK(data_flag_cse(ch,cf,i),0,
     .                                 phs_mask_fel) ) then
                                call sbit(data_flag_cse(ch,cf,i),7,1)
                           end if
                       end if
                    end do

*                   Increment the number of bad normal points
                    if( np_possible(ls) ) then
                        num_np_bad(cf) = num_np_bad(cf) + 1
                    end if
                else
*                   Normal point is fine.  Check to see if we are 
*                   pushing a bias flag and compute epoch number
*                   to which the normal point will refer.
                    np_ep = ep + np_size/2

                    if( push_bias(ls) ) then
*                       Get the channel number
                        if( np_start.gt.0 ) then
                           ch = ltoc(ctol_cse(1,cf,np_ep), ls, 
     .                               actual_max_chan)        
                        else
                           ch = ltoc(ctol_cse(1,cf,ep), ls, 
     .                               actual_max_chan)        
                        end if
*                       Now we may have not observed this satellite
*                       so ch may be zero.
                        if( ch.gt.0 ) then
*                           We really are forming normal points, 
*                           put the bias flag on the normal point
*                           epoch, otherwise put it on the first
*                           data point to be used in the normal
*                           point
                            if( np_start.gt.0 ) then
                                call sbit(data_flag_cse(ch,cf,np_ep),
     .                                                           32,1)
                            else
                                call sbit(data_flag_cse(ch,cf,ep),32,1)
                            end if
                            push_bias(ls) = .false.
                        end if
                     end if

*                    if a normal point was possible, then increment
*                    that we have a good one
                     if( np_possible(ls) ) then
                         num_np_good(cf) = num_np_good(cf) + 1
                     end if

****                 Now compute the average residual values and save
*                    in the residual slot for this normal point epoch.
                     do k = 1, 4
                        sum_omc(k,ls) = sum_omc(k,ls)/np_size
                     end do

*                    Now create the average omc value.  Here we add
*                    back the clocks and remove the cycle slips (the
*                    cycle slips will be added back when the cfile
*                    is written.
*                    Get the channel number
                     if( np_start.gt.0 ) then
                         ch = ltoc(ctol_cse(1,cf,np_ep), ls, 
     .                               actual_max_chan)        
                     else
                         ch = ltoc(ctol_cse(1,cf,ep), ls, 
     .                               actual_max_chan)        
                     end if
                     if( ch.gt.0 .and. np_start.gt.0 ) then
                         rcv_clk = params_cse(cf,np_ep)
                         svs_clk = params_cse(num_cfiles+ls,np_ep)
                         
                         L1r_phs_cse(ch,cf,np_ep) = (sum_omc(1,ls) -
     .                                     L1_cyc_cse(ch,cf,np_ep)) +
     .                               (rcv_clk-svs_clk)
                         L2r_phs_cse(ch,cf,np_ep) = (sum_omc(2,ls) -
     .                                     L2_cyc_cse(ch,cf,np_ep)) +
     .                               (rcv_clk-svs_clk)*fL2(ls)/fL1(ls)
                         L1r_rng_cse(ch,cf,np_ep) =  sum_omc(3,ls) +
     .                               (rcv_clk-svs_clk)
                         L2r_rng_cse(ch,cf,np_ep) =  sum_omc(4,ls) +
     .                               (rcv_clk-svs_clk)*fL2(ls)/fL1(ls)
                     end if
*                          ! Normal point was good
                 end if
*                          ! Looping over satellites
              end do
*                          ! Looping over the normal point epochs
          end do

****      Now downweight the data at the end which will not form a
*         complete normal point
          do i = end_ep+np_size, num_ep        
*            Scan the channels
             do ch = 1, actual_max_chan

*               See if the data is OK
                ls = ctol_cse(ch,cf,i)
                if( ls.gt.0 ) then
                    if( data_OK(data_flag_cse(ch,cf,i),0,
     .                          phs_bias_mask_fel ) ) then
                        call sbit(data_flag_cse(ch,cf,i),7,1)
                    end if
                end if
             end do
          end do

*                          ! Looping over the sites.
      end do

****  Now write out a summary:
      write(uns, 200) np_size, np_start
      if( uns.ne.6 ) write( 6 , 200) np_size, np_start
 200  format(/'Normal point statistics for ',i4,' epoch normal points',
     .        ' Starting at epoch ',i4,/,
     .        6(' SITE  Good Bad'))
      write(uns,220) (cf_codes(cf), num_np_good(cf), num_np_bad(cf),
     .                cf=1,num_cfiles)
      if( uns.ne.6 ) 
     .write( 6 ,220) (cf_codes(cf), num_np_good(cf), num_np_bad(cf),
     .                cf=1,num_cfiles)
 220  format(6(:,1x,a4,1x,i4,1x,i4))   


****  Thats all
      return 
      end


