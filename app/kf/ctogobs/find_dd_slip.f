
CTITLE FIND_DD_SLIP
 
      subroutine find_dd_slip(ns, lv, js, kv, checked, bias_last,
     .        L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse )

      implicit none
 
*     This routine will find the one-way observable which caused
*     a slip on the double difference using ns, lv, js, and kv.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv, js, kv  - Site and satellite list numbers to be checked.
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 ns, lv, js, kv,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
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
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep)
 
*   checked - Indicates whether the first two points have
*           - been previously checked.
*   bias_last   - Indicates that the last point has had a bias
*           - flag added
 
      logical checked, bias_last
 
* LOCAL VARIABLES
 
*   step_ep  - Number of epochs bewteen data points at this site.
*   ch      - Channel number for current PRN
*   i,j     - Pointers to the dd_site_list and dd_svs_list arrays
*   ltoc    - Returns channel number for a given satellite
 
      integer*4 step_ep , ch, i,j, ltoc

*   one_found  - Indicates  that at least one-way was found
*                with a slip.  If not set then all one-ways
*                will be flagged.

      logical one_found

***** First start with the site and satellie thatwe were originally
*     cleaning.  Find out which entries they are in the dd_lists
      call lv_to_list(ns, lv, i, j)
      if( orig_sampling(ns).gt.sampling ) then
          step_ep = orig_sampling(ns)/sampling 
      else
          step_ep = 1
      end if
      one_found = .false.

****  Temporially recude the tolerance
      dd_lc_tol(2) = dd_lc_tol(2) * 0.8d0

      call check_continuity(i,j, step_ep,
     .        checked, bias_last, L1r_phs_cse,
     .        L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse)

      one_found = one_found .or.  bias_last

*     if we did flag any bias with this pair try the other one ways
*     in the list.
      if( .not. bias_last ) then
          write(*,100) cf_codes(ns), prn_list(lv), scan_ep(2)
 100      format(' FIND_DD_SLIP: Slip not found for ',a4,' PRN ',i2.2,
     .           ' Epoch ',i5)
      end if 

*     Now check all the other one-ways to make sure that there
*     are no slips on these
      call lv_to_list(ns, kv, i,j)
      call create_dd_search( ns, kv ) 
      call check_continuity(i,j, step_ep,
     .        checked, bias_last, L1r_phs_cse,
     .        L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse)

      one_found = one_found .or.  bias_last

      if ( .not. bias_last .and. .not.one_found ) then 
          write(*,100) cf_codes(ns), prn_list(kv), scan_ep(2)
      end if

      call lv_to_list(js, lv, i,j)
      call create_dd_search( js, lv ) 
      call check_continuity(i,j, step_ep,
     .        checked, bias_last, L1r_phs_cse,
     .        L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse)

      one_found = one_found .or.  bias_last

      if ( .not. bias_last .and. .not.one_found ) then 
          write(*,100) cf_codes(js), prn_list(lv), scan_ep(2)
      end if
      call lv_to_list(js, kv, i,j)
      call create_dd_search( js, kv ) 
      call check_continuity(i,j, step_ep,
     .        checked, bias_last, L1r_phs_cse,
     .        L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse)

      one_found = one_found .or.  bias_last

      if ( .not. bias_last .and. .not.one_found ) then 
          write(*,100) cf_codes(js), prn_list(kv), scan_ep(2)
      end if

****  If we found no slips at all then flag everyone to be safe
*     (This really should not occurr).

      if( .not.one_found ) then
          write(*,200)
 200      format(' FIND_DD_SLIP: could not find slip.  Flagging ',
     .           ' all combinations')
          write(*,220) scan_ep(2), cf_codes(ns),prn_list(lv)
 220      format(' UNFLAGGED slip at epoch ',I5,' Site ',a4,
     .           ' PRN ',i2.2,' type LCDD')
          ch = ltoc(ctol_cse(1,ns,scan_ep(2)), lv, actual_max_chan)  
          call sbit(data_flag_cse(ch,ns, scan_ep(2)), 31,1)
          call sbit(bf_type_cse(ch,ns,scan_ep(2)), 7,1)

          write(*,220) scan_ep(2), cf_codes(ns),prn_list(kv)
          ch = ltoc(ctol_cse(1,ns,scan_ep(2)), kv, actual_max_chan)  
          call sbit(data_flag_cse(ch,ns, scan_ep(2)), 31,1)
          call sbit(bf_type_cse(ch,ns,scan_ep(2)), 7,1)

          write(*,220) scan_ep(2), cf_codes(js),prn_list(lv)
          ch = ltoc(ctol_cse(1,js,scan_ep(2)), lv, actual_max_chan)  
          call sbit(data_flag_cse(ch,js, scan_ep(2)), 31,1)
          call sbit(bf_type_cse(ch,js,scan_ep(2)), 7,1)

          write(*,220) scan_ep(2), cf_codes(js),prn_list(kv)
          ch = ltoc(ctol_cse(1,js,scan_ep(2)), kv, actual_max_chan)  
          call sbit(data_flag_cse(ch,js, scan_ep(2)), 31,1)
          call sbit(bf_type_cse(ch,js,scan_ep(2)) ,7,1)

      end if

****  Reset the search list so that we can keep on cleaning

      call create_dd_search( ns, lv ) 
      dd_lc_tol(2) = dd_lc_tol(2) / 0.8d0

****  Thats all
      return
      end

CTITLE LV_TO_LIST

      subroutine lv_to_list( ns, lv, i, j)

      implicit none

*     Routine to return the entry in the dd lists for site and satellit
*     ns and lv.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*  ns, lv  - Site and satellite numbers
*   i,  j  - The position of ns and lv in the dd_site_list and dd_sv_list
*            that contain the order in wich data should be cleanned

      integer*4 ns, lv, i,  j

* LOCAL VARIABLES

*   k    - Counters for scanning the list

      integer*4 k

*     Scan down the list to get the values
      i = 0
      do k = 1, num_cfiles
         if( dd_site_list(k).eq.ns ) i = k
      end do
      j = 0
      do k = 1, num_sat
         if( dd_svs_list(k).eq.lv ) j = k
      end do

***** Thats all
      return
      end
 
