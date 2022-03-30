 
CTITLE CREATE_DD_SEARCH
 
      subroutine create_dd_search( ns, lv )

      implicit none
 
*     routine to create the list of satellites and sites to search
*     for double differences. These lists will be modified as the
*     searching for double differences continues.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for current one-
*              way. These will be excluded from the list
 
 
      integer*4 ns, lv
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
***** Make the station list based on dd_site_list
 
      j = 0
      do i = 1, num_cfiles
          if( dd_site_list(i).ne.ns ) then
              j = j + 1
              dd_site_search(j) = dd_site_list(i)
          end if
      end do
 
*     Do the same for the satellites
 
      j = 0
      do i = 1, num_sat
          if( dd_svs_list(i).ne.lv ) then
              j = j + 1
              dd_svs_search(j) = dd_svs_list(i)
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE UPDATE_DD_SEARCH
 
      subroutine update_dd_search( is, jv )

      implicit none
 
*     This routine will move the is'th entry in site list and the
*     jv'th entry in th satellite to the top of the search list.
*     CURRENT VERSION ONLY UPDATES THE SITE LIST SO THAT ALLOW_ONE_BG
*     WILL WORK (order of satellites is not changed).
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   is, jv   - Site number and satellite in the search list which will
*             be moved to the top of the list
 
 
      integer*4 is, jv
 
* LOCAL VARIABLES
 
*   issave, jvsave  - Saved values of the entries to be moved
*   i               - Loop counter
 
 
      integer*4 issave, i
 
***** Push the top list of entries down by one spot and then put the
*     is'th entry at the top of the list
*
 
      if( is.gt.1 ) then
          issave = dd_site_search(is)
 
*         Push the other entries down
          do i = is-1 , 1, -1
              dd_site_search(i+1) = dd_site_search(i)
          end do
          dd_site_search(1) = issave
      end if
 
****  Now do the same for the satellites
C     if( jv.gt.1 ) then
C         jvsave = dd_svs_search(jv)
 
*         Push the other entries down
C         do i = jv-1, 1, -1
C             dd_svs_search(i+1) = dd_svs_search(i)
C         end do
C         dd_svs_search(1) = jvsave
C     end if
 
****  Thats all
      return
      end
 
CTITLE GET_DD_DATA
 
      subroutine  get_dd_data( ns, lv, js, kv, ep, step_ep, 
     .                norm_cyc,b_cyc, chi_cyc,
     .                L1r_phs_cse, L2r_phs_cse,
     .                L1r_rng_cse, L2r_rng_cse,
     .                L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                data_flag_cse, reliable)

      implicit none
 
*     This routine will take the list of one way data in data_lft
*     and data_rgh arrays and find the double differences that match
*     them.  The algorithm is to first scan the satellites at this
*     station, and having found an exceptable group, to scan the
*     stations to find an acceptable one.  As the algorithm works
*     it makes a list of the "best" choices.
 
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number
*   js, kv   - Second site and satellite in DD.
*   ep       - Epoch at which the gap or bias flag occurrs
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 ns, lv, ep, step_ep, js, kv, 
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   reliable       - Logical to indicate that the current patch
*             could be done reliably 
 
      logical reliable

*   norm_cyc(2,2), b_cyc(2) - Normal equations and b-solution
*             for estimating number of cycles (used to see if
*             have a reliable estimate.) 
*   chi_cyc         - Prefit chi**2 value
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
 
      real*8 norm_cyc(2,2), b_cyc(2), chi_cyc,
     .    L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   is, jv  - postitions in the search lists for the sites and
*           - satellites
*   ow_ep(max_max_dd_ret,2) - Epoch numbers for the one ways
*           - in the left and right segments (1-left, 2-right)
*   sd_ep(max_max_dd_ret,2) - Epoch numbers for the satellite
*           - satellite single differences
*   ow_num(2)       - Number of values in the left and right
*                   - one ways
*   sd_num(2)       - Number of values in the left and right
*                   - single differences
*   max_dd_min      - Maximum value for the mininum of the
*                   - left and right segments of data.  (used
*                   - so that if we don't get as much data as
*                   - we want, we can use the most available
*   max_sd_min      - Maximum value for the mininum of the 
*                     for the single differences.
*   best_dd_kv, best_dd_js    - Values of jv, is corresponding to the
*                   - max_dd_min
*   best_sd_kv      - Best satellite single difference found.
*   best_ch_js, best_ch_kv - Best site and satellite combination
*                     based on the chiqual value from est_dd_cyc.
*   kv_save         - Value set to any satellite which we can form
*                     form a single difference to (to check if we
*                     should allow one_bg_gap).
 
      integer*4 i, is, jv, ow_ep(max_max_dd_ret,2),
     .    sd_ep(max_max_dd_ret,2), ow_num(2),
     .    sd_num(2), max_dd_min, best_dd_kv, best_dd_js,
     .    max_sd_min, best_sd_kv, best_ch_js, best_ch_kv,
     .    kv_save  
 
*   svs_OK, site_OK - indicates that data collected for a
*                   - satellite and site are OK
*   fin_svs, fin_site   - Indicates that all conditions are
*                   - satisfied for the satellite and site
*                   - double differences
*   continued       - Used in the station loop if we have
*                     made a trial estimate and decided that
*                     we need to keep searching the double
*                     differences
 
      logical svs_OK, site_OK, fin_svs, fin_site, continued

*   chiqual        - Estimate of the quality of the bias fixing
*                    from est_dd_cyc
*   best_chiqual   - Best best_chiqual value found.

*   norm_trl(2,2), b_trl(2) - Trial solution values for the
*                    estimation of the number of cycle slips.
*   norm_best(2,2), b_best(2) - Best set of normal equations 
*                    and solution vector.
*   dl1_slip, dl2_slip - Estimates of L1 and L2 slip

      real*8 chiqual, best_chiqual, norm_trl(2,2), b_trl(2),
     .       norm_best(2,2), b_best(2) , dl1_slip, dl2_slip,
     .       chi_trl, chi_best

*   no_good_fnd  - Logical to indicate that no good double 
*            differences were found during scanning.
*   obg_nd - Returns true if one_bg_needed to form single
*            difference

      logical no_good_fnd, obg_nd
 
***** First scan over the satellites at this station try to find
*     an acceptable one.  When is found try to find a good station,
*     If no good station can be found then try another satellite and
*     keep going until we run out of data.
 
      fin_svs = .false.
      continued = .false.
      no_good_fnd = .true.
      jv = 0
      max_sd_min = 0
      best_chiqual = 0.d0

      best_sd_kv = 0
      best_dd_kv = 0
      best_dd_js = 0

      kv_save = 0
 
*     Copy the epochs from the one_ways
      do i = 1, num_lft
          ow_ep(i,1) = ep_lft(i)
      end do
      do i = 1, num_rgh
          ow_ep(i,2) = ep_rgh(i)
      end do
      ow_num(1) = num_lft
      ow_num(2) = num_rgh
 
      fin_svs = .false.
      jv = 0
 
      do while ( .not. fin_svs )
 
*         Go to next satellite
          jv = jv + 1
          kv = dd_svs_search(jv)
*         Here we just get the list of epochs which have good data
          call get_sdep( ns, lv, kv, ep, step_ep, 
     .        ow_ep, ow_num, sd_ep, sd_num, ctol_cse, data_flag_cse,
     .        svs_OK, obg_nd)

*         If we were able to form a single difference without invoking
*         one bias gap, then save satellite.
          if( svs_OK .and. .not.obg_nd ) kv_save = kv

 
*****     If the satellite looks good then get the double difference
          if ( svs_OK ) then
*             Loop over the stations with this combination of
              if( .not. continued ) is = 0
              continued = .false.
              fin_site = .false.
              max_dd_min = 0
              do while ( .not.fin_site )
                  is = is + 1
                  js = dd_site_search(is)
*                 Here we just get the list of epochs which have
*                 good data
* MOD TAH 990517: Only form double differences between L1+L2 data
*                 stations if the main station is of this type.
                  site_OK = .true.
                  if( lambda(lv,2,ns).ne.0 .and. 
     .                lambda(lv,2,js).eq.0      ) site_OK = .false.
                  if( site_OK )    
     .            call get_ddep( ns, lv, js, kv, ep, step_ep,
     .                sd_ep, sd_num,
     .                ctol_cse, data_flag_cse, site_OK)

 
****              If this combination is not OK, save the maxiumum
*                 of the data with the minimum of the left and right
*                 segments.  In this way, if we never find anything
*                 really good we can use the best available.
                  if( .not.site_OK ) then
                      if( min(dd_num(1),dd_num(2)).gt.
     .                    max_dd_min ) then
                          max_dd_min = min(dd_num(1),dd_num(2))
                          best_dd_kv = dd_svs_search(jv)
                          best_dd_js = dd_site_search(is)
                      end if
                  else  

*                     Set status to say we are finished.
                      fin_site = .true.
                      fin_svs  = .true.
                      max_sd_min = min(sd_num(1),sd_num(2))
                  end if
 
*****             Check to make sure that we have not run out of
*                 data.
                  if( is.eq.num_cfiles-1 ) fin_site = .true.
              end do
          else            

*             Save the best sd we have in case we need it later.
              if( min(sd_num(1),sd_num(2)).gt. max_sd_min ) then
                  max_sd_min = min(sd_num(1),sd_num(2))
                  best_sd_kv = dd_svs_search(jv)
              end if
          end if
 
*****     Now if we finished with the satellites and we still have
*         not found an acceptable one (probably due to bias flags
*         on all satellites at the this time), if there was never
*         any good data then then fix in one-ways and remove bias flag.
*         else use the best combination we found
          if( jv.eq.num_sat-1  ) fin_svs = .true.

*         See if all is OK and whether we should do a trial estimation
*         of the number of cycle slips
          if( fin_svs ) then

*             Normally we will be here becuase both svs and site are
*             OK; but often we may have run out of stations and in some
*             cases we may have run out of satellites.  We can tell these
*             cases apart by the values of max_sd_min and max_dd_min.

*             If max_sd_min is zero then we never found any single diffs.
*             at all in which case we will just use LG to keep the ion
*             continuous.  If max_sd_min is greater than zero then there
*             are single differences so we will see if thre arre any 
*             double differences
              fin_svs = .false.
* MOD TAH 950828: Changed the tolerance here to be min_good_bias rather than
*             4 data points.
              if( max_sd_min.ge.min_good_bias ) then

*                 There are single differences, see if we found any
*                 double differences.  If svs_OK is not true then we
*                 never searched for double differences becuase none of
*                 single differences were acceptable.  We are desperate
*                 now so try the best single difference we had.
                  if( .not.svs_OK ) then

*                      Get the single difference epochs.
                       kv = best_sd_kv
                       call get_sdep( ns, lv, kv, ep, step_ep, 
     .                      ow_ep, ow_num, sd_ep, sd_num, 
     .                      ctol_cse, data_flag_cse, svs_OK, obg_nd)

*                      now try to find double differences with
*                      set of single differences.

*                      Loop over the stations with this combination of
                       is = 0
                       fin_site = .false.
                       max_dd_min = 0
                       do while ( .not.fin_site )
                           is = is + 1
                           js = dd_site_search(is)
*                          Here we just get the list of epochs which have
*                          good data
* MOD TAH 990517: Only combine dual frequency data.
                           site_OK = .true.
                           if( lambda(lv,2,ns).ne.0 .and. 
     .                         lambda(lv,2,js).eq.0 ) site_OK = .false.
                           if( site_OK )    
     .                     call get_ddep( ns, lv, js, kv, ep, step_ep,
     .                         sd_ep, sd_num, 
     .                         ctol_cse, data_flag_cse, site_OK)
 
****                       If this combination is not OK, save the maxiumum
*                          of the data with the minimum of the left and right
*                          segments.  In this way, if we never find anything
*                          really good we can use the best available.
                           if( .not.site_OK ) then
                               if( min(dd_num(1),dd_num(2)).gt.
     .                             max_dd_min ) then
                                 max_dd_min = min(dd_num(1),dd_num(2))

                                 best_dd_kv = dd_svs_search(jv)
                                 best_dd_js = dd_site_search(is)
                             end if
                         else  
*                            Set status to say we are finished.
                             fin_site = .true.
                             fin_svs  = .true.
                         end if
 
*****                    Check to make sure that we have not run out of
*                        data.
                         if( is.eq.num_cfiles-1 ) fin_site = .true.
                     end do
                  end if

*                 See if we have any double differences.  Max_dd_min
*                 >0 means we found some, and if site_OK we never
*                 set then we will have to use these.
                  if( .not.site_OK .and. max_dd_min.gt.0 ) then
*                      Get the single difference epochs.
                       kv = best_dd_kv
                       js = best_dd_js
                       call get_sdep( ns, lv, kv, ep, step_ep, 
     .                      ow_ep, ow_num, sd_ep, sd_num, 
     .                      ctol_cse, data_flag_cse, svs_OK, obg_nd)
* MOD TAH 990517: Only combine dual frequency data.
                       site_OK = .true.
                       if( lambda(lv,2,ns).ne.0 .and. 
     .                     lambda(lv,2,js).eq.0 ) site_OK = .false.
                       if( site_OK )    
     .                 call get_ddep( ns, lv, js, kv, ep, step_ep,
     .                      sd_ep, sd_num, 
     .                      ctol_cse, data_flag_cse, site_OK)

                  end if

*****             If we have DD data in both the left and right
*                 segments then do a trial estimate
* MOD TAH 950828: Change min values to ge.min_good_bias from gt.4.
                  if( dd_num(1).ge.min_good_bias .and. 
     .                dd_num(2).ge.min_good_bias       ) then

*                     Set logical to indicate that good data has 
*                     been found so that we will not edit later.
*                     [Data gets edited if there are no double 
*                     diffferences to the last single difference
*                     found]
                      no_good_fnd = .false.

                      call form_dd( ns, lv, js,
     .                    kv, ep, 
     .                    L1r_phs_cse, L2r_phs_cse,
     .                    L1r_rng_cse, L2r_rng_cse,
     .                    L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                    data_flag_cse)

*                     increment the normal equations
                      call copy_norm(norm_cyc, b_cyc, chi_cyc, 
     .                               norm_trl, b_trl, chi_trl)
* MOD TAH 990517: Only use LC if L2 data is available     
                      if( lambda(lv,2,ns).ne.0 ) then
                          call inc_cyc_est(ns,lv,js,kv, 'LCLG','DD',
     .                                   norm_trl, b_trl, chi_trl)
                      else
                          call inc_cyc_est(ns,lv,js,kv, 'L1','DD',
     .                                   norm_trl, b_trl, chi_trl)
                      endif
                      

*                     If widelane available include this as well.
                      if( lambda(lv,4,ns).eq.1 .and. 
     .                    lambda(kv,4,ns).eq.1 .and.
     .                    lambda(lv,4,js).eq.1 .and.
     .                    lambda(kv,4,js).eq.1      ) then
                          call inc_cyc_est(ns,lv,js,kv,'WL','DD',
     .                         norm_trl, b_trl, chi_trl) 
                      end if

*                     Estimate the number of cycles
                      call est_dd_cyc(ep, ns, lv, js, kv, norm_trl, 
     .                     b_trl, chi_trl, reliable, chiqual, 
     .                     dL1_slip, dL2_slip)

*                     Save the best values of chi**2 quality estimate
                      if( chiqual.gt.best_chiqual ) then
                          best_chiqual = chiqual
                          best_ch_js = js
                          best_ch_kv = kv
                          call copy_norm(norm_trl, b_trl, chi_trl,
     .                                   norm_best, b_best, chi_best)
                      end if

*                     If the bias fixing was not reliable then
*                     try again.  Drop the satellite back (so that
*                     we will continue where we left off once it
*                     it is incremented.
                      if( .not.reliable ) then

*                         Make sure we are not out of data
                          if( is.ne.num_cfiles-1 .and.
     .                        jv.ne.num_sat-1 ) then

*                             Only drop the satellite number back if
*                             we still have stations to check, otherwisw
*                             just go to next satellite.
                              if( is.ne.num_cfiles-1 ) then
                                  jv = jv - 1   ! Use pointer in search list
                                  continued = .true.
                              end if
                              fin_svs = .false.
                              max_sd_min = 0
                          end if
                      else

*                         Set fin_svs to say that we are done.
                          fin_svs = .true.
                      end if
                  else     
*                     not enough data for patch, see if we should edit
*                     Get the single difference epochs and double
*                     differnece epochs if we had any.
                      if( best_dd_kv.gt.0 .and. best_dd_js.gt.0 .and.
     .                    no_good_fnd ) then
                          kv = best_dd_kv
                          js = best_dd_js
                          call get_sdep( ns, lv, kv, ep, step_ep, 
     .                         ow_ep, ow_num, sd_ep, sd_num, 
     .                         ctol_cse, data_flag_cse, 
     .                         svs_OK, obg_nd)
* MOD TAH 990517: Only combine dual frequency data.
                          site_OK = .true.
                          if( lambda(lv,2,ns).ne.0 .and. 
     .                        lambda(lv,2,js).eq.0 ) site_OK = .false.
                          if( site_OK )    
     .                    call get_ddep( ns, lv, js, kv, ep, step_ep,
     .                         sd_ep, sd_num, 
     .                         ctol_cse, data_flag_cse, site_OK)

*                         Call the routine to form the data so that
*                         any editing of double differences will be
*                         done.  This checks for bad pointts on the
*                         side where we have data.
                          call form_dd( ns, lv, js,
     .                        kv, ep, 
     .                        L1r_phs_cse, L2r_phs_cse,
     .                        L1r_rng_cse, L2r_rng_cse,
     .                        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                        data_flag_cse)

                          call save_toofew(dd_num(1), dd_ep(1,1),
     .                         dd_edit_num, dd_edit_ep, max_dd_edit)
                          call save_toofew(dd_num(2), dd_ep(1,2),
     .                         dd_edit_num, dd_edit_ep, max_dd_edit)
                      end if
                  end if   ! num_dd > 0
             else

*                 Not enough data for patch in the single differences
*                 see if we should edit
*                 Get the single difference epochs.
                  if( best_sd_kv.gt.0 .and. no_good_fnd ) then
                      kv = best_sd_kv
                      call get_sdep( ns, lv, kv, ep, step_ep, 
     .                          ow_ep, ow_num, sd_ep, sd_num, 
     .                          ctol_cse, data_flag_cse, svs_OK, obg_nd)
                      call save_toofew(sd_num(1), sd_ep(1,1),
     .                     dd_edit_num, dd_edit_ep, max_dd_edit)
                      call save_toofew(sd_num(2), sd_ep(1,2),
     .                     dd_edit_num, dd_edit_ep, max_dd_edit)
                  end if
             end if        ! max_sd_min > 0
             if( jv.eq.num_sat-1  ) fin_svs = .true.
          end if           ! fin_svs
      end do

****  If the best_chiqual is greater than zero, then something
*     found so retrun updated normal equations
      if( best_chiqual.gt.0 ) then
          call copy_norm(norm_best, b_best, chi_best,
     .                   norm_cyc, b_cyc, chi_cyc)
          js = best_ch_js
          kv = best_ch_kv
      else
*         We never found anything suitable so set the site number
*         and satellite number to zero (when routine exits we can
*         use this to tell no double differences were formed
          js = 0

*         If we could find no satellites with out invoking one_bg
*         then kv_save will still be zero.
          kv = kv_save
      end if

****  Now see if  bith svs and site were good then update the
*     dd search list 
*     MODIFIED to only update the site list so that allow_one_bg
*     will work.
      if( svs_OK .and. site_OK ) then
          call update_dd_search( is, jv )
      end if
 
****  Thats all
      return
      end
 
CTITLE GET_SDEP
 
      subroutine get_sdep( ns, lv, kv, ep, step_ep, 
     .        ow_ep, ow_num, sd_ep, sd_num, ctol_cse, data_flag_cse,
     .        svs_OK, obg_nd)

      implicit none
 
*     Routine to the get the epochs of data corresponding to the
*     one-way data for satellite kv.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for one way
*   kv       - Satellite number to difference with.
*   ep       - Epoch at which the gap or bias flag occurrs
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   ow_ep(max_max_dd_ret,2) - Epoch numbers for the one ways
*           - in the left and right segments (1-left, 2-right)
*   sd_ep(max_max_dd_ret,2) - Epoch numbers for the satellite
*           - satellite single differences
*   ow_num(2)       - Number of values in the left and right
*                   - one ways
*   sd_num(2)       - Number of values in the left and right
*                   - single differences
 
      integer*4 ow_ep(max_max_dd_ret,2), sd_ep(max_max_dd_ret,2),
     .    ow_num(2), sd_num(2)
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 ns, lv, kv, ep, step_ep,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   svs_OK      - Set OK if we get at least half the amount
*               - of data in the single differences as in the
*               - one ways
*   obg_nd      - Set true if one bias gap needed to get a single
*                 difference
 
      logical svs_OK, obg_nd
 
* LOCAL VARIABLES
 
*   i,j     - Epoch loop counter
*   ch      - Channel number for sd satellite
*   ltoc    - Function to return channel number
*   num_bf_fnd - Number of bias flags found.  When allow_one_gb is
*             set true we allow one gap or bias flag in the data.
*             This lets us try to patch across a gap (or bias) on
*             all satellites at a station.
*   num_gap_fnd - Number of gaps found (one allowed when allow_one_gb
*             is set.
 
      integer*4 i, j, ch, ltoc, num_bf_fnd, num_gap_fnd
 
*   finished - Indicates that we have finished getting the data
*   data_OK - Logical function which indicates that the data is
*           - OK.
*   kbit    - Logical function for checking bits in a word
*   good_bf - Logical function that returns true if there is 
*             a bias flag on a good data point.
*   good_fnd - Logical to indicate some good data found.  Once
*             we have found good data we no longer allow another
*             gap or bias bias flag.
*   still_allow_bg - Logical that indicates that we will still
*             allow a bias or gap
 
      logical finished, data_OK, good_bf, good_fnd,
     .        still_allow_bg
 
****  Start:
      svs_OK = .false.
      obg_nd = .false.
      num_bf_fnd = 0
      num_gap_fnd = 0
      still_allow_bg = allow_one_bg 

*     Loop over the left and right segmenets of data.
      sd_num(1) = 0

*     Since this is at the same station as the one-ways we only
*     need to check the data at the epochs collected at this station
      i = ep
      finished = .false.
      do while ( .not.finished )
 
          i = i  - step_ep
 
*         Get the channel number of the new satellite at this epoch
          ch = ltoc(ctol_cse(1,ns,i), kv, actual_max_chan)
 
          if( ch.gt. 0 ) then
 
*             We have data see if bias flag.  Also check
*             here for bad data or gaps.  
              if( good_bf(data_flag_cse(ch,ns,i),0,phs_mask) ) then
 
*                 There is a bias flag so stop collecting data if 
*                 allow_one_bg is not set.
                  num_bf_fnd = num_bf_fnd + 1
                  if( .not. allow_one_bg ) finished = .true.
                  if( .not. still_allow_bg ) finished = .true.
                  if( allow_one_bg .and. num_bf_fnd.gt.1 ) 
     .                                     finished = .true.
              end if 
 
*             Add this epoch to the single differences, only if there
*             is a one-way to go with it
              if( data_OK(data_flag_cse(ch,ns,i),0,phs_mask) )
     .                                                then
*                 loop over the one-way epochs to make sure 
*                 here
                  do j = 1, ow_num(1) 
                      if( i.eq. ow_ep(j,1) ) then
                          sd_num(1) = sd_num(1) + 1
                          sd_ep(sd_num(1),1) = i
                          good_fnd = .true.
                          still_allow_bg = .false.
                      end if
                  end do

*                 If we are allowing one_gap and the number of gaps 
*                 has been set (says we have started one gap) then
*                 increment the number.
* MOD TAH 950824: If gaps were flagged then we don't need to worry
*                 about the gap if it does not have a bias flag.
                  if( num_gap_fnd.gt.0 .and. 
     .               .not.ignore_gaps ) num_gap_fnd = num_gap_fnd+1
              else  
*                 Treat a bad data point like a missing data point.
*                 (Same as code below)
                  if( .not.ignore_gaps ) finished = .true.

*                 See if we are allowing one gap (see below):
                  if( still_allow_bg ) then
                      finished = .false.
                      if( num_gap_fnd.gt.1 ) finished = .true.
                      num_gap_fnd = 1      
                  end if
              end if
          else
*             Terminate search if there is a gap provided we are not
*             ignoring caps.
              if( .not.ignore_gaps ) finished = .true.

*             See if we are allowing one gap: The way this works is
*             that num_gap_fnd is set to one and will stay in that
*             state until a good data point is found at which time
*             it will be incremented and the next time through the
*             gap code finished will be set true. (indicates being
*             of next gap).
              if( still_allow_bg ) then
                  finished = .false.
                  if( num_gap_fnd.gt.1 ) finished = .true.
                  num_gap_fnd = 1      
              end if
          end if

*         See if this is the last point we need to check
          if( i.le. ow_ep(ow_num(1),1) ) finished = .true.
      end do

****  See if needed to invoke one bias gap
      if( num_bf_fnd.gt.1 .or. num_gap_fnd.gt.1 ) obg_nd = .true.

***** Now do the right segment of data.
      sd_num(2) = 0
      num_gap_fnd = 0
      num_bf_fnd  = 0
      still_allow_bg = allow_one_bg

*     Since this is at the same station as the one-ways we only
*     need to check the data at the epochs collected at this station
      i = ep - step_ep
      finished = .false.
      do while ( .not.finished )
 
          i = i  + step_ep
 
*         Get the channel number of the new satellite at this epoch
          ch = ltoc(ctol_cse(1,ns,i), kv, actual_max_chan)
 
          if( ch.gt. 0 ) then
 
*             We have data see if bias flag.  Also check
*             here for bad data or gaps.  
              if( good_bf(data_flag_cse(ch,ns,i),0,phs_mask) ) then
 
*                 There is a bias flag so stop collecting data if 
*                 allow_one_bg is not set.
                  num_bf_fnd = num_bf_fnd + 1
                  if( .not. still_allow_bg ) finished = .true.
                  if( still_allow_bg .and. num_bf_fnd.gt.1 ) 
     .                                     finished = .true.
              end if

*             If data is good and we are not finished, add data
*             point in. Make  sure we have a one-way to go withit
              if( data_OK(data_flag_cse(ch,ns,i),0,phs_mask) .and.
     .             .not.finished ) then
*                 loop over the one-way epochs to make sure 
*                 here
                  do j = 1, ow_num(2) 
                      if( i.eq. ow_ep(j,2) ) then
                          sd_num(2) = sd_num(2) + 1
                          sd_ep(sd_num(2),2) = i
                          still_allow_bg = .false.
                      end if
                  end do
*                 If we are allowing one_gap and the number of gaps 
*                 has been set (says we have started one gap) then
*                 increment the number.
* MOD TAH 950824: If gaps were flagged then we don't need to worry
*                 about the gap if it does not have a bias flag.
                  if( num_gap_fnd.gt.0 .and. 
     .               .not.ignore_gaps ) num_gap_fnd = num_gap_fnd+1
              else     
*                 Treat bad data like misssing data and add to
*                 cap condition (same code as below)
                  if( .not.ignore_gaps ) finished = .true.

*                 See if we are allowing one gap
                  if( still_allow_bg ) then
                      finished = .false.
                      if( num_gap_fnd.gt.1 ) finished = .true.
                      num_gap_fnd = 1      
                  end if
              end if
          else

*             Terminate search if there is a gap provided we are not
*             ignoring caps.
              if( .not.ignore_gaps ) finished = .true.

*             See if we are allowing one gap: The way this works is
*             that num_gap_fnd is set to one and will stay in that
*             state until a good data point is found at which time
*             it will be incremented and the next time through the
*             gap code finished will be set true. (indicates being
*             of next gap).
              if( still_allow_bg ) then
                  finished = .false.
                  if( num_gap_fnd.gt.1 ) finished = .true.
                  num_gap_fnd = 1      
              end if
          end if

*         See if this is the last point we need to check
          if( i.ge. ow_ep(ow_num(2),2) ) finished = .true.
      end do

****  See if needed to invoke one bias gap
      if( num_bf_fnd.gt.1 .or. num_gap_fnd.gt.1 ) obg_nd = .true.
 
***** Now see if we have enough data
 
      if( sd_num(1).gt.4 .and. sd_num(2).gt.4 ) then
          svs_OK = .true.
      end if
 
***** Thats all
      return
      end
 
CTITLE GET_DDEP
 
      subroutine get_ddep( ns, lv, js, kv, ep, step_ep,
     .                sd_ep, sd_num, 
     .                ctol_cse, data_flag_cse, site_OK)

      implicit none
 
*     Routine to get the double difference epochs.
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for one way
*   js, kv   - Site and Satellite number to difference with.
*   ep       - Epoch at which the gap or bias flag occurrs
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   sd_ep(max_max_dd_ret,2) - Epoch numbers for the sd
*           - in the left and right segments (1-left, 2-right)
*   sd_num(2)       - Number of values in the left and right
*                   - single differences
*                   - double differences
 
      integer*4 sd_ep(max_max_dd_ret,2), sd_num(2)
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 ns, lv, js, kv, ep, step_ep,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   site_OK     - Set OK if we get at least half the amount
*               - of data in the single differences as in the
*               - one ways
 
      logical site_OK
 
* LOCAL VARIABLES
 
*   i,j     - Epoch loop counter
*   c1,c2   - Channel numbers for the satellites at the second
*           - station
*   ltoc    - Function to return channel number
*   js_step_ep -  Number of epoch bewteen data points for second
*             site
*   js_cnt     - Counter used to see if we have found a gap
*             in data (used particularly for data sampled at
*             greater than the sampling interval for this data.)
 
      integer*4 i, j, c1,c2, ltoc, js_step_ep, js_cnt
 
*   finished - Indicates that we have finished getting the data
*   data_OK - Logical function which indicates that the data is
*           - OK.
*   kbit    - Logical function for checking bits in a word
*   good_bf - Logical function that returns true if there is 
*             a bias flag on a good data point.
 
 
      logical finished, data_OK, good_bf
 
****  Start:
      site_OK = .false.

*     Get the sampling step size for this site
      if( orig_sampling(js).gt.sampling ) then
           js_step_ep = orig_sampling(js)/sampling
      else
           js_step_ep = 1
      end if

      js_cnt = 0
 
*     Loop over the left and right segmenets of data.
      dd_num(1) = 0
*     Since the new station may have a differenct spacing we
*     check ever epoch hear
      i = ep
      finished = .false.

*     Make sure we have some sd data
      if( sd_num(1).eq.0 ) finished = .true.

      do while ( .not.finished )
 
          i = i  - 1
*         See if this is the last point we need to check
          if( i.le. sd_ep(sd_num(1),1) ) finished = .true.
 
*         Get the channel number of the satellites at this epoch
          c1 = ltoc(ctol_cse(1,js,i), lv, actual_max_chan)
          c2 = ltoc(ctol_cse(1,js,i), kv, actual_max_chan)
 
          if( c1.gt. 0 .and. c2.gt.0 ) then

*             Found good data so reset the counter for sampling
*             of data
              js_cnt = 0
 
*             We have data see if bias flag.  We could also check
*             here for bad data or gaps.  For the moment ignore these
*             for terminating the search, but only save epoch if data
*             is good.  (Later we could add a bit to data_flag to say
*             that data has been cleaned --- somewhat safer this later
*             way.
C             if( (kbit(data_flag_cse(c1,js,i),31).or.
C    .            kbit(data_flag_cse(c1,js,i),32))      .or.
C    .            (kbit(data_flag_cse(c2,js,i),31).or.
C    .            kbit(data_flag_cse(c2,js,i),32))     )then
              if( good_bf(data_flag_cse(c1,js,i),0,phs_mask) .or.
     .            good_bf(data_flag_cse(c2,js,i),0,phs_mask) ) then
 
*                 There is a bias flag so stop collecting data
                  finished = .true.

              end if
 
*             Add this epoch to the double differences if both
*             data are good
              if( data_OK(data_flag_cse(c1,js,i),0,phs_mask) .and.
     .           data_OK(data_flag_cse(c2,js,i),0,phs_mask) )
     .                                                then

*                 Make sure that we have a single difference
*                 to go with it.
                  do j = 1, sd_num(1)
                      if( i.eq.sd_ep(j,1) ) then 
                          dd_num(1) = dd_num(1) + 1
                          dd_ep(dd_num(1),1) = i
                      end if
                  end do
               else

*                  Treat bad data point as a gap.
                   js_cnt = js_cnt + 1
                   if( js_cnt.ge.js_step_ep ) then
                      finished = .true.
                   end if
               end if
          else
*              Terminate search on gap if we have skipped more points
*              than the sampling allows
               js_cnt = js_cnt + 1
               if( js_cnt.ge.js_step_ep ) then
                  finished = .true.
               end if
          end if
      end do
 
***** Now do the right segment of data.
      dd_num(2) = 0

*     Reset the counter for number of skipped data
      js_cnt = 0

*     Since the new station may have a differenct spacing we
*     check ever epoch hear
      i = ep - 1
      finished = .false.

*     Make sure we have some data
      if( sd_num(2).eq.0 ) finished = .true.

      do while ( .not.finished )
 
          i = i  + 1
*         See if this is the last point we need to check
          if( i.ge. sd_ep(sd_num(2),2) ) finished = .true.
 
*         Get the channel number of the satellites at this epoch
          c1 = ltoc(ctol_cse(1,js,i), lv, actual_max_chan)
          c2 = ltoc(ctol_cse(1,js,i), kv, actual_max_chan)
 
          if( c1.gt. 0 .and. c2.gt.0 ) then

*             Found good data so reset the counter for number of
*             skipped data
              js_cnt = 0
 
*             We have data see if bias flag.  We could also check
*             here for bad data or gaps.  For the moment ignore these
*             for terminating the search, but only save epoch if data
*             is good.  (Later we could add a bit to data_flag to say
*             that data has been cleaned --- somewhat safer this later
*             way.
C             if( ( (kbit(data_flag_cse(c1,js,i),31).or.
C    .               kbit(data_flag_cse(c1,js,i),32))      .or.
C    .            (kbit(data_flag_cse(c2,js,i),31).or.
C    .            kbit(data_flag_cse(c2,js,i),32))    ) )then
              if( good_bf(data_flag_cse(c1,js,i),0,phs_mask) .or.
     .            good_bf(data_flag_cse(c2,js,i),0,phs_mask) ) then
 
*                 There is a bias flag so stop collecting data
                  finished = .true.
              else
 
*                 Add this epoch to the double differences if both
*                 data are good
                  if( data_OK(data_flag_cse(c1,js,i),0,phs_mask) .and.
     .                data_OK(data_flag_cse(c2,js,i),0,phs_mask) )
     .                                                then

*                     Make sure that we have a single difference
                      do j = 1, sd_num(2)
                          if( i.eq. sd_ep(j,2) ) then
                              dd_num(2) = dd_num(2) + 1
                              dd_ep(dd_num(2),2) = i
                          end if
                      end do
                  else
*                     Treat bad data like a gap (missing data)
                      js_cnt = js_cnt + 1
                      if( js_cnt.ge.js_step_ep ) then
                         finished = .true.
                      end if
                  end if
               end if
          else
*              Terminate search on gap if we have skipped more points
*              than the sampling allows
               js_cnt = js_cnt + 1
               if( js_cnt.ge.js_step_ep ) then
                  finished = .true.
               end if
          end if
      end do
 
 
***** Now see if we have enough data
 
      if( (dd_num(1).gt.sd_num(1)/5 .and. dd_num(1).gt.4) .and.
     .    (dd_num(2).gt.sd_num(2)/5 .and. dd_num(2).gt.4) ) then
          site_OK = .true.
      end if
 
***** Thats all
      return
      end
 
CTITLE FORM_DD
 
      subroutine form_dd( ns, lv, js, kv, ep, 
     .                    L1r_phs_cse, L2r_phs_cse,
     .                    L1r_rng_cse, L2r_rng_cse,
     .                    L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                    data_flag_cse)

      implicit none
 
*     Routine to form the double differences and to put them
*     into the work arrays.  The data quality was checked as
*     the dd_ep array was formed and so all of this data should be
*     good.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number
*   js, kv   - Site number and satellite list number for the second pair
*   ep       - Epoch at which the gap or bias flag occurrs
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 ns, lv, js, kv, ep,
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
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
 
*   i       - loop counter
*   c11, c12, c21, c22 - channels number of each station
*             satellite combination in the double differences
*             (station is first index and satellite is second)
*   ltoc    - Function to return channel number
*   ed      - Epoch number loop variabel

 
      integer*4 i,  c11, c12, c21, c22, ltoc, ed

****  Form the double differences for the left segment
      num_wl = dd_num(1)
 
      do i = 1, dd_num(1)
          ed = dd_ep(i,1)
          ep_wl(i) = ed
 
*         Get all the channel numbers for the double differences
          c11 = ltoc(ctol_cse(1,ns,ed), lv, actual_max_chan)
          c12 = ltoc(ctol_cse(1,ns,ed), kv, actual_max_chan)
          c21 = ltoc(ctol_cse(1,js,ed), lv, actual_max_chan)
          c22 = ltoc(ctol_cse(1,js,ed), kv, actual_max_chan)
 
*         Make sure that all the channel numbers are greater than
*         zero.  We have a bug if they are not:
          if( c11.gt.0 .and. c12.gt.0 .and.
     .        c21.gt.0 .and. c22.gt.0 ) then
 
*             L1 Phase difference
              data_wl(1,i) = (L1r_phs_cse(c11,ns,ed)+
     .                        L1_cyc_cse(c11,ns,ed)) -
     .                    (L1r_phs_cse(c12,ns,ed)+
     .                        L1_cyc_cse(c12,ns,ed)) -
     .                    (L1r_phs_cse(c21,js,ed)+
     .                        L1_cyc_cse(c21,js,ed)) +
     .                    (L1r_phs_cse(c22,js,ed)+
     .                        L1_cyc_cse(c22,js,ed))
 
*             L2 Phase double difference
              if( lambda(lv,2,ns).ne.0 ) then
                 data_wl(2,i) = (L2r_phs_cse(c11,ns,ed)+
     .                        L2_cyc_cse(c11,ns,ed)) -
     .                    (L2r_phs_cse(c12,ns,ed)+
     .                        L2_cyc_cse(c12,ns,ed)) -
     .                    (L2r_phs_cse(c21,js,ed)+
     .                        L2_cyc_cse(c21,js,ed)) +
     .                    (L2r_phs_cse(c22,js,ed)+
     .                        L2_cyc_cse(c22,js,ed))
              else
                 data_wl(2,i) = 0.0d0
              end if
 
*             L1 Range residuals
              data_wl(3,i) = L1r_rng_cse(c11,ns,ed) -
     .                    L1r_rng_cse(c12,ns,ed) -
     .                    L1r_rng_cse(c21,js,ed) +
     .                    L1r_rng_cse(c22,js,ed)
 
*             L2 range residuals
              if( lambda(lv,4,ns).ne.0 ) then
                 data_wl(4,i) = L2r_rng_cse(c11,ns,ed) -
     .                    L2r_rng_cse(c12,ns,ed) -
     .                    L2r_rng_cse(c21,js,ed) +
     .                    L2r_rng_cse(c22,js,ed)
              else
                 data_wl(4,i) = 0.d0
              end if
 
          else
 
*             We have a problem
              write(*,400) ed, c11, c12, c21, c22, ns, lv, js, kv
 400          format(' *** ERROR *** Left segment DD',
     .                ' differences NF ',
     .                ' Ep ',i4,' Cs ',4i3,' Pairs ',4i3)
*             Set number of data to 0
              num_wl = 0
              num_wr = 0
              RETURN
          end if
 
      end do
 
***** Now do the right segment of data.
      num_wr = dd_num(2)
 
      do i = 1, dd_num(2)
          ed = dd_ep(i,2)
          ep_wr(i) = ed
 
*         Get all the channel numbers for the double differences
          c11 = ltoc(ctol_cse(1,ns,ed), lv, actual_max_chan)
          c12 = ltoc(ctol_cse(1,ns,ed), kv, actual_max_chan)
          c21 = ltoc(ctol_cse(1,js,ed), lv, actual_max_chan)
          c22 = ltoc(ctol_cse(1,js,ed), kv, actual_max_chan)
 
*         Make sure that all the channel numbers are greater than
*         zero.  We have a bug if they are not:
          if( c11.gt.0 .and. c12.gt.0 .and.
     .        c21.gt.0 .and. c22.gt.0 ) then
 
*             L1 Phase difference
              data_wr(1,i) = (L1r_phs_cse(c11,ns,ed)+
     .                        L1_cyc_cse(c11,ns,ed)) -
     .                    (L1r_phs_cse(c12,ns,ed)+
     .                        L1_cyc_cse(c12,ns,ed)) -
     .                    (L1r_phs_cse(c21,js,ed)+
     .                        L1_cyc_cse(c21,js,ed)) +
     .                    (L1r_phs_cse(c22,js,ed)+
     .                        L1_cyc_cse(c22,js,ed))
 
*             L2 Phase double difference
              if( lambda(lv,2,ns).ne.0 ) then
                 data_wr(2,i) = (L2r_phs_cse(c11,ns,ed)+
     .                        L2_cyc_cse(c11,ns,ed)) -
     .                    (L2r_phs_cse(c12,ns,ed)+
     .                        L2_cyc_cse(c12,ns,ed)) -
     .                    (L2r_phs_cse(c21,js,ed)+
     .                        L2_cyc_cse(c21,js,ed)) +
     .                    (L2r_phs_cse(c22,js,ed)+
     .                        L2_cyc_cse(c22,js,ed))
              else
                 data_wr(2,i) = 0.d0
              end if
 
*             L1 Range residuals
              data_wr(3,i) = L1r_rng_cse(c11,ns,ed) -
     .                    L1r_rng_cse(c12,ns,ed) -
     .                    L1r_rng_cse(c21,js,ed) +
     .                    L1r_rng_cse(c22,js,ed)
 
*             L2 range residuals
              if( lambda(lv,4,ns).ne.0 ) then
                 data_wr(4,i) = L2r_rng_cse(c11,ns,ed) -
     .                    L2r_rng_cse(c12,ns,ed) -
     .                    L2r_rng_cse(c21,js,ed) +
     .                    L2r_rng_cse(c22,js,ed)
              else 
                 data_wr(4,i) = 0.d0
              endif
 
          else
 
*             We have a problem
              write(*,410) ed, c11, c12, c21, c22, ns, lv, js, kv
 410          format(' *** ERROR *** Right segment DD',
     .                ' differences NF ',
     .                ' Ep ',i4,' Cs ',4i3,' Pairs ',4i3)
*             Set number of data to 0
              num_wl = 0
              num_wr = 0
              RETURN
          end if
 
      end do
 
****  Thats all
      return
      end
 
CTITLE COPY_NORM

      subroutine copy_norm(norm_cyc, b_cyc, chi_cyc,
     .                     norm_trl, b_trl, chi_trl )

      implicit none

*     Routine to copy the normal equations to trial values fro
*     testing how well we can the number of cycle slips

* PASSED VALIABLES

* norm_cyc(2,2), b_cyc(2) - Original normal eqautions and solution
*     vector
* norm_trl(2,2), b_trl(2) - Copy of the original for testing
* chi_cyc, chi_trl        - Prefit chi**2 values

      real*8 norm_cyc(2,2), b_cyc(2), norm_trl(2,2), b_trl(2),
     .       chi_cyc, chi_trl

* LOCAL VARIABLES

*  i,j   - Loop counters

      integer*4 i,j

***** Copy the matrices
      do i = 1,2
         b_trl(i) = b_cyc(i)
         do j = 1,2
            norm_trl(i,j) = norm_cyc(i,j)
         end do
      end do

      chi_trl = chi_cyc

***** Thats all
      return
      end

 
