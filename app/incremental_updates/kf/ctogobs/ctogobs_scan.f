CTITLE SCAN_DD
 
      subroutine scan_dd(iter, L1r_phs_cse, L2r_phs_cse,
     .     L1_cyc_cse, L2_cyc_cse, ctol_cse, data_flag_cse,
     .     bf_type_cse )

      implicit none
 
*     This routine will scan the double differences for the stations
*     given in the scan_sites list (and specified by the scan_sites
*     command).  The algorithm for checking the double differences uses
*     groups of three sucessive points. Nominally, if the first to have
*     already been checked, the change between the second and the third
*     is compared with the first to second change.  If the 2-3 diffence
*     is outdie tolerance then 3 is flagged.  In some cases, 1-2 will
*     not have been checked and is this case the bias flag may be put
*     on point 2.  For an outlier point, the double pair of bias flags
*     will cause the point to be down weighted in clean_dd.
 
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

*   iter  - Iteration for scanning.  On iteration 2 we reverse
*           the order the satellites are scanned in.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), iter
 
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
*   svs_switch - Temp storage for switching satellite order
 
      integer*4 i,j,k, ns, lv, ch, ltoc, first_ep, step_ep
 
*   kbit    - Checks if bit is set.
*   first_bias  - Indicates that first bias flag has been
*           - found
 
      logical kbit, first_bias
 
***** First order the stations and satellites for cleaning
 
      call get_clean_order('SCANNING')

* MOD TAH 990519: Comment out this code which was not used 
*     because iter was always passed as one.  We not use 
*     iter to see if we should ignore gaps/
      if( iter.gt.1 ) then
          ignore_gaps = .true.
          write(*,110) iter
 110      format(' For Iter ',i2,' Ignoring gaps')
      else
          ignore_gaps = usr_ignore_gaps
      endif


C     if( iter.eq.2 ) then
C         write(*,'(a)') ' Iteration 2: Reversing order of scan'
C         do i = 1, num_sat/2
C            svs_switch = dd_svs_list(i)
C            dd_svs_list(i) = dd_svs_list(num_sat+1-i)
C            dd_svs_list(num_sat+1-i) = svs_switch
C         end do
C     end if 
 
*     Now loop over the one-way data
      do i = 1, num_cfiles
          ns = dd_site_list(i)
 
*         Only do this site if we have been asked to do so (scan_sites
*         command).
          if( kbit(scan_sites,ns) ) then
              do j = 1, num_sat
                  lv = dd_svs_list(j)
                  first_bias = .false.
 
****              Scan up to the first bias flag (always marked) and
*                 step at the data sampling interval for this site
                  k = 0
                  do while( k.lt.num_ep .and. .not.first_bias )
                      k = k + 1
                      ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                                         ! SV observed at this time
                      if( ch.gt.0 ) then
 
                          if( kbit(data_flag_cse(ch,ns,k),31).or.
     .                        kbit(data_flag_cse(ch,ns,k),32) ) then
                              first_bias = .true.
                          end if
                      end if
                  end do
 
*                 Based on the sampling interval at this site compute
*                 number of epochs bewteen data points.
                  if( orig_sampling(ns).gt.sampling ) then
                      step_ep = orig_sampling(ns)/sampling
                  else
                      step_ep = 1
                  end if
 
*                 Start at the first point after the first bias flag
*                 after saving the first point
 
                  first_ep = k
                  write(*,120) cf_codes(ns), prn_list(lv), first_ep,
     .                        step_ep
 120              format(' Scanning ',a4,' PRN ',i2.2,' Start ep ',i4,
     .                    ' Step ',i3,' epochs')
 
                  call create_dd_search( ns, lv )
 
                  call do_scan_eps(i,j, first_ep, step_ep,
     .                L1r_phs_cse, L2r_phs_cse,
     .                L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                data_flag_cse, bf_type_cse)
*                         ! Looping over satellites
              end do
*                         ! We are checking this site
          end if
*                         ! Looping oversites
      end do
 
***** Thats all
      return
      end
 
CTITLE DO_SCAN_EPS
 
      subroutine do_scan_eps(i,j, first_ep, step_ep, L1r_phs_cse,
     .        L2r_phs_cse,L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, bf_type_cse)

      implicit none
 
*     This routine actually does the scanning along the one ways
*     checking all continuous triplets of data.  It is scanning the
*     i'th station in the dd_site_list,and the j'th satellite in the
*     dd_svs_list.
 
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
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
 
*   i,j     - Site number in dd_site_list and svs number in
*           - dd_svs_list
*   first_ep    - Epoch number of first data point after the first
*           - bias flag
*   step_ep - Number of epochs between each data point for this
*           - site.
 
      integer*4 i,j, first_ep, step_ep
 
* LOCAL VARIABLES
 
*   k       - Counter over epochs
*   ns, lv  - Current station and satellite list number
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
      integer*4 k, ns, lv, ch, ltoc
 
*   data_OK - function returns true if data OK.  phs_mask does
*           - no check bias flags
*   checked - Set true if we have checked two points
*   bias_last   - Indicates that a bias flag was added to the
*           - last point of the triplet being checked.
 
      logical data_OK, checked, bias_last

      logical debug_out

***** Set the site and satellite list numbers
      ns = dd_site_list(i)
      lv = dd_svs_list(j)
 
      if( ns.eq.1 .and. lv.eq.3 ) then
         debug_out = .true.
      else
         debug_out = .false.
      end if 
         debug_out = .false.
*     Initialize the triplet values
      num_scan = 1
      scan_ep(1) = first_ep
 
      do k = first_ep+step_ep, num_ep, step_ep
 
*         see if we have the end of a gap or a bias flag.
*         Following pair of if's is really just one
*         but bounds checking will fail if we put
*         both in same if statement when ch=0.
          ch = ltoc(ctol_cse(1,ns,k), lv, actual_max_chan)
*                             ! SV observed at this time
          if( debug_out ) then
              write(*,*) 'Epoch ',k,' ch ',ch, num_scan,
     .                   scan_ep(1), scan_ep(2), scan_ep(3)
          end if
          if( ch.gt.0 ) then
*             See if we have found a good data point
*             NOTE: By using phs_bias_mask, both bad data and
*             bias flags will be detected.
 
              if( data_OK(data_flag_cse(ch,ns,k), 0,
     .            phs_bias_mask) ) then
 
*                 Increment number of points in scan list
                  num_scan = num_scan + 1
                  if( debug_out ) write(*,*) 'Good data ',num_scan
                  if( num_scan.eq.2 ) then
*                     Add next to list
                      scan_ep(2) = k
 
*                     See if continuous with previous point
C                      if( scan_ep(1)+step_ep.ne.k ) then
* MOD TAH 990518: Allow gaps up to gap size.
                      if( debug_out ) write(*,*) scan_ep(1), 
     .                   gap_size(ns)*step_ep, k

                      if( scan_ep(1)+gap_size(ns)*step_ep.lt.k ) then
*                         Discontinuous, reset the lists:
*                         (NOTE: We really don't need to worry about
*                         single gaps providedi ngore_gaps is set
*                         false.
                          scan_ep(1) = k
                          num_scan = 1
                      end if
                      if( debug_out ) write(*,*) scan_ep(1), 
     .                   gap_size(ns)*step_ep, k
 
*                     If we have just incremented to two this means
*                     that they have not been checked.  For
*                     continuous data we will not pass though this
*                     code.
                      checked = .false.
                  end if
 
*                 See if we have three continuous points, in
*                 which case we should check them
                  if( num_scan.eq.3 ) then
                      scan_ep(3) = k
c                      if( scan_ep(2)+step_ep.ne.k ) then
* MOD TAH 990518: Allow checks over interval upto gap_size
                      if( scan_ep(2)+gap_size(ns)*step_ep.lt.k ) then
*                         data is discontinuous.  Reset the
*                         number of points to check
                          num_scan = 2
                      end if
 
*****                 If the current epochs have not been
*                     checked then call the continuity
*                     checker.
                      if( debug_out ) write(*,*) 'Check_c1 ',i,j
                      call check_continuity(i,j, step_ep,
     .                    checked, bias_last, L1r_phs_cse,
     .                    L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                    ctol_cse, data_flag_cse, bf_type_cse)
 
*                     If we had 3 data points and not bias flagg
*                     was set during the continuity check, ripple
*                     values back by one, and continue looking.  If
*                     this is not the case, then go back to just
*                     current point
                      if( num_scan.eq.3 .and. .not.bias_last ) then
                          num_scan = 2
                          scan_ep(1) = scan_ep(2)
                          scan_ep(2) = scan_ep(3)
                      else
                          num_scan = 1
                          scan_ep(1) = k
                      end if
*                             ! Three values in list
                  end if
 
              else if(data_OK(data_flag_cse(ch,ns,k), 0,
     .                phs_mask) ) then
****              If this point is OK, except thatit has a bias flag, then
*                 if we have two points check them and reset the scan counters
*                 back
                  if( num_scan.eq.2 ) then
                      if( debug_out ) write(*,*) 'Check_c2 ',i,j
                      call check_continuity(i,j, step_ep,
     .                       checked, bias_last, L1r_phs_cse,
     .                       L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                       ctol_cse, data_flag_cse, bf_type_cse)
                  end if
                  num_scan = 1
                  scan_ep(1) = k
*                             ! data OK
              end if
*                             ! Channel was non-zero
          end if
*                             ! Looping over epochs.
      end do
 
****  Finally, we may still have two points hanging at the end which
*     have not been checked.  If so, then check now.
      if( .not.checked .and.num_scan.eq.2 ) then
          if( debug_out ) write(*,*) 'Check_cF ',i,j
          call check_continuity(i,j, step_ep,
     .        checked, bias_last, L1r_phs_cse,
     .        L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse)
      end if
 
****  Thats all
      return
      end
 
 
CTITLE CHECK_CONTINUITY
 
      subroutine check_continuity(i, j, step_ep, checked, bias_last,
     .        L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, bf_type_cse )

      implicit none
 
*     This routine will check the continuity of the data epochs passed
*     through common, and flag the one ways if the discontinuity is
*     found.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   i, j     - Site number in the dd_site_list and satellite in the
*              dd_svs_list number
*   step_ep  - Number of epochs bewteen data points at this site.
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 i, j, step_ep,
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

* MOD TAH 200617: Added additional test.
      logical verified  ! Set true to jump seems to occur at epoch
                    ! where bias flag to be added (TAH 200614).
 
* LOCAL VARIABLES
 
*   ch      - Channel number for current PRN
*   ns, lv  - Site number and satellite list number.
*   is1, is2    - Index in dd_sites_search for 2nd and 3rd station
*   jv1, jv2    - Index in dd_svs_search to the second and third
*               - sateliites in double difference
*   js1, js2    - Second and third station numbers
*   kv1, kv2    - second and third satellite numbers
*   bs, bv      - The bad station and satellite that needs a bias
*               - flag added.
*   bep         - Epoch at which the bad data is.
*   ltoc        - Returns channel number for given satelite/station
*                 epoch
*   best_bs, best_bv  - These are non-zero values returned by the
*                 scanners for bad station and bad satellite 
*                 as they loop over trial satellites.  (used
*                 when the satellite or station can not be uniquely 
*                 determined.
*  last_bs, last_bv  - The last valid other station and satellite
*                 checked.  We use in case we run out of stations
*                 or satellites but still have a double difference
*                 jump.

 
      integer*4 ch, ns, lv, is1, is2, jv1, jv2, js1, kv1,
     .    bs, bv, bep, ltoc, best_bs, best_bv, last_bs, last_bv
 
*   finished        - Indicates that we have finished checking
*                   - current one-way combination
*   all_fnd         - Indicates that all stations and satellites
*                   - have been found.
*   svs_OK          - indicates that we have found a good satellite
*                   - for forming the primary double differences.
*   dd_OK           - Set true when we find a valid double differnece
*                   - with this site/SVS.
*   bias_flag       - Indiates that when finished a scan we
*                   - found a bias flag somewhere in the combinations.
*   qual_bad            - True if double difference looks bad
*   kbit            - Checks if bit set
 
      logical finished, all_fnd, svs_OK, dd_OK, bias_flag, qual_bad, 
     .        kbit
 
*   scan_qual       - Numerical value of the quality checker
*                   - for the initial double difference
*   sq(2)           - quality checker for combinations with
*                   - another station and satellite. Double difference
*                     changes when (1) satellite is changed; (2) when
*                     station changed
 
 
      real*8 scan_qual, sq(2)
 
***** Start trying to find double differences with this one-way
*     data
      ns = dd_site_list(i)
      lv = dd_svs_list(j)

****  Start by scanning over sateliites to get single difference
 
      finished = .false.
      jv1 = 0
      bep = 0
      bs  = 0
      bv  = 0
      best_bs = 0
      best_bv = 0
      last_bs = 0
      last_bv = 0
      bias_last = .false.
 
      do while ( .not. finished )
 
*         Go to next satellite
          jv1 = jv1 + 1
          kv1 = dd_svs_search(jv1)
*         Here we just get the list of epochs which have good data
          call get_scan_sd( ns, lv, kv1, step_ep, ctol_cse,
     .            data_flag_cse, svs_OK)
 
*****     If the satellite looks good then get the double difference
          if ( svs_OK ) then
*             Loop over the stations with this combination of
              is1 = 0
              call get_scan_dd('SITE',ns, lv, is1, jv1, ctol_cse,
     .            data_flag_cse, dd_OK, bias_flag)
             
 
*             If we have found a second site for the double difference
*             check the quality of the dd's.  If there are not good,
*             then search for a third site and satellite.  By
*             comparision of who is good and bad, decide which of the
*             orginal pair.  (There is the possibility that the third
*             satellite is wrong in which case the bias flag will be
*             put in the wong place.  However if some satellites are OK
*             then this third one will be flagged latter, and the
*             correctly flagged point should be patched later in
*             clean_dd.  If all satellites are bad, then the last one
*             will never be checked (becuase by the time it gets there
*             there will no inflagged data to form a double difference
*             with.  Howver, this data will not be used in cleaning
*             other data until it has been cleaned and therefore there
*             should not be a problem.  (The one-ways at this site may
*             not be clean, but the single difference at the station
*             should be clean.)
              if( dd_OK ) then
 
*****             Scan the quality of this double difference
                  js1 = dd_site_search(is1)
                  last_bs = js1
                  last_bv = kv1
                  scan_qual = -2*dd_lc_tol(3)
                  call check_dd_quality( ns, lv, js1, kv1, bep,
     .                L1r_phs_cse, L2r_phs_cse,
     .                L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                data_flag_cse, scan_qual, qual_bad, checked )
                  if( bs.ne.0 .or. bv.ne.0 ) then
cd                     write(*,195) bep, ns, lv, js1, kv1,
cd   .                         scan_qual, qual_bad
cd 195                  format(' CHKDD  EP ',i5,' S1/C1',2i3,' S2/C2 ',
cd     .                       2i3, ' Scan qual ',F7.2, L2) 
                      if( abs(scan_qual).gt.dd_lc_tol(3)) 
     .                                               qual_bad = .true.
                  end if
                   
*                 If the double difference, is not OK then search
*                 for a third station and satellite
                  if( .not. qual_bad ) then
 
*                     Then we are finished, it not then we will
*                     need to look around
                      finished = .true.
                      checked = .true.
                      bep = 0
                  else
 
*                     Start the site scan at the current site and
*                     scan up
                      call scan_chk_2nd_ss(ns,lv,is1,jv1,is2, jv2,
     .                        bs, bv, sq, all_fnd, 
     .                        L1r_phs_cse, L2r_phs_cse,
     .                        L1_cyc_cse, L2_cyc_cse, 
     .                        ctol_cse,data_flag_cse, checked)
                      if( kbit(status_rep,11) ) 
     .                write(*,200) bep, ns, lv, dd_site_search(is1), 
     .                         dd_site_search(is2),dd_svs_search(jv1),
     .                         dd_svs_search(jv2), bs, bv, 
     .                         scan_qual, sq
 200                  format(' DD2SS  EP ',i5,' S1/C1',2i3,' S2/S3 ',
     .                       2i3, ' C2/C3 ',2i3,' BS/BV',2i3,
     .                       ' Scan qual ',3F7.2) 
*                     If all_fnd is true then both an addition satelite
*                     and station were found and therefre we are finished.
                      if( bs.ne.0 ) best_bs = bs
                      if( bv.ne.0 ) best_bv = bv
                      if( all_fnd ) finished = .true.
                  end if
*                             ! double differnce was found for the
              end if
*                             ! original one way data
*                             ! Satellite found
          end if
 
****      See if we have run out of satellites
          if( jv1.eq.num_sat - 1 ) finished = .true.
      end do

***** Now flag the point if a bad one was found
      if( bep.ne.0 ) then

*****     Use the best_bs and best_bv available (NOTE: If we came directly
*         here after scan_chk_2nd_ss, the best values will be set to the 
*         lastest values.
          bs = best_bs
          if( bs.eq.0 ) js1 = last_bs
          bv = best_bv
          if( bv.eq.0 ) kv1 = last_bv

*         Flag the point with the bias.  See if it was uniques.  If the
*         site or satellite could be resolved flag both stations or
*         sites. We need to do this becuase we allow oene-way data fixing
*         later on in clean_dd.

          if( bs.ne.0 ) then
               if( bv.ne.0 ) then
                  ch = ltoc(ctol_cse(1,bs,bep),bv, actual_max_chan)

* MOD TAH 200613: Add another continity check to see if really as slip
*                 or just close to tolerance.  Only do this if slip is
*                 less than dd_lc_tol(3).  Only unique (UQ) types tested.
                  verified = .true.
                  if( abs(scan_qual).lt.dd_lc_tol(3) ) then
*                     Secondary check. 
                      call  verify_scan_dd( bep, ns, lv, js1, kv1, 
     .                   L1r_phs_cse, L2r_phs_cse,
     .                   L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .                   data_flag_cse, verified )
                  endif
                  if( verified ) then 
                     write(*,400) bep, cf_codes(bs), prn_list(bv),
     .                            scan_qual, ch, 'UQ'
 400                 format(' UNFLAGGED slip at epoch ',i5,' Site ',a4,
     .                      ' PRN_',I2.2,' Jump: ',F9.2,' cyc, CH ',i2,
     .                      ' TY ',a2)

*                    Increment number of flags set for this site and svs.
                     num_dd_flags(bs, bv) = num_dd_flags(bs, bv) + 1
 
                     if( kbit(status_rep,10) ) 
     .               write(*,420) bep, ns, lv, js1, dd_site_search(is2),
     .                            kv1, dd_svs_search(jv2), bs, bv, 
     .                            scan_qual, sq
 420                 format(' DDSCAN EP ',i5,' S1/C1',2i3,' S2/S3 ',2i3,
     .                      ' C2/C3 ',2i3,' BS/BV',2i3,
     .                      ' Scan qual ',3F7.2) 
                     if( ch.gt.0 ) then
                         call sbit(data_flag_cse(ch,bs,bep),31,1)
                         call sbit(bf_type_cse(ch,bs,bep),5,1)
                     else
                       write(*,*) 'Channel < 0 Line 578'
                     end if
                  endif 
              else
*                 Flag both satellites.
                  ch = ltoc(ctol_cse(1,bs,bep),lv, actual_max_chan)
                  write(*,400) bep, cf_codes(bs), prn_list(lv),
     .                     scan_qual, ch, 'LL'
*                 Increment number of flags set for this site and svs.
                  num_dd_flags(bs,lv) = num_dd_flags(bs, lv) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,bs,bep),31,1)
                      call sbit(bf_type_cse(ch,bs,bep),5,1)
                  end if
                  ch = ltoc(ctol_cse(1,bs,bep),kv1, actual_max_chan)
                  write(*,400) bep, cf_codes(bs), prn_list(kv1),
     .                     scan_qual, ch, 'LL'
*                 Increment number of flags set for this site and svs.
                  num_dd_flags(bs, kv1) = num_dd_flags(bs, kv1) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,bs,bep),31,1)
                      call sbit(bf_type_cse(ch,bs,bep),5,1)
                  end if
                  if( kbit(status_rep,10) ) 
     .            write(*,420) bep, ns, lv, js1, dd_site_search(is2),
     .                         kv1, dd_svs_search(jv2), bs, bv, 
     .                         scan_qual, sq
              end if
          else
              if( bv.ne.0 ) then
                  ch = ltoc(ctol_cse(1,ns,bep),bv, actual_max_chan)
                  write(*,400) bep, cf_codes(ns), prn_list(bv),
     .                         scan_qual, ch, 'SS'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(ns, bv) = num_dd_flags(ns, bv) + 1

                  if ( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,ns,bep),31,1)
                      call sbit(bf_type_cse(ch,ns,bep), 5,1)
                  end if
                  ch = ltoc(ctol_cse(1,js1,bep),bv, actual_max_chan)
                  write(*,400) bep, cf_codes(js1), prn_list(bv),
     .                         scan_qual, ch,'SS'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(js1, bv) = num_dd_flags(js1, bv) + 1

                  if ( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,js1,bep),31,1)
                      call sbit(bf_type_cse(ch,js1,bep), 5,1)
                  end if
                  if( kbit(status_rep,10) ) 
     .            write(*,420) bep, ns, lv, js1, dd_site_search(is2),
     .                         kv1, dd_svs_search(jv2), bs, bv, 
     .                         scan_qual, sq
              else
*                 Flag both satellites and sites.
                  ch = ltoc(ctol_cse(1,ns,bep),lv, actual_max_chan)
                  write(*,400) bep, cf_codes(ns), prn_list(lv),
     .                         scan_qual, ch, 'SL'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(ns , lv) = num_dd_flags(ns , lv) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,ns,bep),31,1)
                      call sbit(bf_type_cse(ch,ns,bep), 5,1)
                  end if
                  ch = ltoc(ctol_cse(1,js1,bep),lv, actual_max_chan)
                  write(*,400) bep, cf_codes(js1), prn_list(lv),
     .                         scan_qual, ch, 'SL'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(js1, lv) = num_dd_flags(js1, lv) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,js1,bep),31,1)
                      call sbit(bf_type_cse(ch,js1,bep), 5,1)
                  end if
                  ch = ltoc(ctol_cse(1,ns,bep),kv1, actual_max_chan)
                  write(*,400) bep, cf_codes(ns), prn_list(kv1),
     .                     scan_qual, ch, 'SL'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(ns ,kv1) = num_dd_flags(ns , kv1) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,ns,bep),31,1)
                      call sbit(bf_type_cse(ch,ns,bep), 5,1)
                  end if
                  ch = ltoc(ctol_cse(1,js1,bep),kv1, actual_max_chan)
                  write(*,400) bep, cf_codes(js1), prn_list(kv1),
     .                     scan_qual, ch, 'SL'

*                 Increment number of flags set for this site and svs.
                  num_dd_flags(js1, kv1) = num_dd_flags(js1, kv1) + 1

                  if( ch.gt.0 ) then
                      call sbit(data_flag_cse(ch,js1,bep),31,1)
                      call sbit(bf_type_cse(ch,js1,bep), 5,1)
                  end if
                  if( kbit(status_rep,10) ) 
     .            write(*,420) bep, ns, lv, js1, dd_site_search(is2),
     .                         kv1, dd_svs_search(jv2), bs, bv, 
     .                         scan_qual, sq
              end if
          end if
          if( bep.eq.scan_ep(num_scan) ) then
*             Bias flag added to last point scanned so mark this
*             point now are unchecked
              checked = .false.
              bias_last = .true.
          end if
       end if

 
***** Thats all
      return
      end
 
CTITLE CHECK_DD_QUALITY
 
      subroutine check_dd_quality( ns, lv, js, kv, bep,
     .        L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .        ctol_cse, data_flag_cse, scan_qual, qual_bad, checked)

      implicit none
 
*     This routine will check the quality of the double differences
*     in nominally three contiguous points.  At times there may be
*     only two points in which case we default tolerances
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv, js, kv  - First station, first satellite, second station
*                     second statellite in double difference 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   bep             - The epoch at which the bias flag should be added
*                     if one needs to be added.
 
      integer*4 ns, lv, js, kv, bep,
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
*   scan_qual       - Measure of the discontinuity (used to flag as
*                     good of bad)
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep), scan_qual

*   qual_bad        - Logical, true if data does not pass quality
*                     check
*   checked         - Indicates that we have checked the first to points

      logical qual_bad, checked
 
* LOCAL PARAMETERS
 
 
* LOCAL VARIABLES

*   c11, c12, c21, c22 - Channel numbers for each one-way in
*                 double difference
*   ltoc        - Function to return channel number for a 
*                 site, satellite list and epoch.
*   ep          - Epoch for data
*   i           - Loop counter

      integer*4 c11, c12, c21, c22, ltoc, ep,  i
 
 
*   dd_obs(3)   - The double differences at the three points
*               - to be checked.
*   dd_L1, dd_L2 - Double differences for L1 and L2
*   dfirst, dsecond - Changes in double difference between the first
*                  pair and second pair
*   tol         - Value that the tested changes in the double differnces
*                 must be less than or will be flagged bad.
 
      real*8 dd_obs(3), dd_L1, dd_L2, dfirst, dsecond, tol
 
***** Loop over the four one-way data and form the double difference
      if( num_scan.lt.2 ) RETURN 
      if( ns.eq.js .or. lv.eq.kv ) then
          write(*,100) ns, js, lv, kv
 100      format('*** WARNING *** ns=js or lv=kv ',4i4)
          RETURN
      end if
      do i = 1, num_scan
          ep = scan_ep(i)
          c11 = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)
          c12 = ltoc(ctol_cse(1,ns,ep), kv, actual_max_chan)
          c21 = ltoc(ctol_cse(1,js,ep), lv, actual_max_chan)
          c22 = ltoc(ctol_cse(1,js,ep), kv, actual_max_chan)
 
          If( c11.gt.0 .and. c12.gt.0 .and.
     .        c21.gt.0 .and. c22.gt.0          ) then
              dd_L1 = (L1r_phs_cse(c11,ns,ep)+
     .                    L1_cyc_cse(c11,ns,ep)) -
     .                (L1r_phs_cse(c12,ns,ep)+
     .                    L1_cyc_cse(c12,ns,ep)) -
     .                (L1r_phs_cse(c21,js,ep)+
     .                    L1_cyc_cse(c21,js,ep)) +
     .                (L1r_phs_cse(c22,js,ep)+
     .                    L1_cyc_cse(c22,js,ep))
 
*             L2 Phase double difference
* MOD TAH 990512: Make sure we have L2 everywhere:
              if( L2r_phs_cse(c11,ns,ep).ne.0.d0 .and.
     .            L2r_phs_cse(c12,ns,ep).ne.0.d0 .and.
     .            L2r_phs_cse(c21,js,ep).ne.0.d0 .and.
     .            L2r_phs_cse(c22,js,ep).ne.0.d0  ) then
                 dd_L2 = (L2r_phs_cse(c11,ns,ep)+
     .                    L2_cyc_cse(c11,ns,ep)) -
     .                   (L2r_phs_cse(c12,ns,ep)+
     .                    L2_cyc_cse(c12,ns,ep)) -
     .                   (L2r_phs_cse(c21,js,ep)+
     .                    L2_cyc_cse(c21,js,ep)) +
     .                   (L2r_phs_cse(c22,js,ep)+
     .                    L2_cyc_cse(c22,js,ep))
                       
* RWK 150203: This may need to be reformulated to account for frequency 
*             differences  between Glonass SVs. 
                 dd_obs(i) = lcf1(1)*dd_L1 + lcf2(1)*dd_L2
              else
                 dd_obs(i) = dd_L1
              end if
          else
              write(*,*) 'ERROR CHECK_DD_QUALITY ns, js, lv, kv ',
     .                    ns, js, lv, kv, num_scan, scan_ep
          end if
      end do
 
****  Now check the quality
      qual_bad = .false.
      if( num_scan.eq.3 ) then
 
*         Here we compare the change between the first pair and the
*         second pair
C         write(*,900) 'ep, ns, js, lv, kv, dobs',  num_scan, ep, ns,
C    .                  js, lv, kv, dd_obs
C900      format(a,6i4,3F9.2)
          dfirst = dd_obs(2) - dd_obs(1)
          dsecond = dd_obs(3) - dd_obs(2)

*         if we know that dfirst has been checked then the second
          if ( checked ) then
             tol = min(dd_lc_tol(3),
     .                max(dd_lc_tol(2),dd_lc_tol(1)*abs(dfirst)))
             scan_qual = dsecond
             if( abs(dsecond).gt.tol ) then
                 bep = scan_ep(3)
                 qual_bad = .true.
             else
                 bep = 0
             end if

          else
*            We don't know which could be bad, so check both
             tol = min(dd_lc_tol(3),
     .                max(dd_lc_tol(2),
     .                    dd_lc_tol(1)*min(abs(dfirst),abs(dsecond))))
             if( abs(dfirst).gt.tol ) then
                 bep = scan_ep(2)
                 qual_bad = .true.
                 scan_qual = dfirst
             else if( abs(dsecond).gt.tol ) then
                 bep = scan_ep(3)
                 qual_bad = .true.
                 scan_qual = dsecond
             end if
          end if
      else

*         Only two points:  See if there difference is large
C         write(*,900) 'num_scan, ep, ns, js, lv, kv, dobs', num_scan,
C    .            ep, ns, js, lv, kv, (dd_obs(i),i=1,2)
          tol = dd_lc_tol(2)
          dfirst = dd_obs(2) - dd_obs(1)
          scan_qual = dfirst
          if( abs(dfirst).gt.tol ) then
              bep = scan_ep(2)
              checked = .false.
              qual_bad = .true.
          end if
      end if
 
****  Thats all
      return
      end
 
CTITLE GET_SCAN_SD
 
      subroutine get_scan_sd( ns, lv, kv, step_ep, ctol_cse,
     .                data_flag_cse, svs_OK)

      implicit none
 
*     Routine to the get the epochs of data corresponding to the
*     one-way data for satellite kv making sure that these are
*     continuous along with the primary one way data set.
* MOD TAH 990519: Changed to operate in the same way as get_scan_dd
*     because with allowing gaps in the data we now need to check
*     to see if there are bias flags in the sequence for the next
*     satellite at the main site. 

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for one way
*   kv       - Satellite number to difference with.
*   step_ep  - Number of epochs bewteen data points at this site.
 
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 ns, lv, kv, step_ep,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   svs_OK      - Set OK if we get at least half the amount
*               - of data in the single differences as in the
*               - one ways
 
      logical svs_OK
 
* LOCAL VARIABLES
 
*   i       - Epoch loop counter
*   ch      - Channel number for sd satellite
*   ltoc    - Function to return channel number
*   begin_ep - First epoch in scan list.  At this epoch we allow a bias
*             flag.  A bias flag on other will cause svs_OK to be set
*             false.
*   end_ep  - Last epoch in list.
*   ep      - Epoch being processed.
 
      integer*4 i, j, ch, ltoc, begin_ep, ep, end_ep
 
*   data_OK - Logical function which indicates that the data is
*           - OK.
*   kbit    - Tests if a bit is set
 
      logical data_OK, kbit
      logical debug_out 
      if( ns.eq.23 .and. lv.eq.24 ) then
         debug_out = .true.
      else
         debug_out = .false.
      endif
      debug_out = .false.
 
****  Start:
      svs_OK = .true.
      if( num_scan.le.1 ) then
          svs_OK = .false.
          RETURN
      end if

*     Get the first epoch (when this routine is called during cleaning
*     the epochs may be in reverse order.)
      begin_ep = scan_ep(1)
      end_ep = scan_ep(1)
      do i = 2, num_scan
         ep = scan_ep(i)
         if( ep.lt.begin_ep ) begin_ep = ep
         if( ep.gt.end_ep ) end_ep = ep 
      end do
 
*     See if the data is available for all the epochs needed
      do ep = begin_ep, end_ep , step_ep 
 
*         Get the channel number of the new satellite at this epoch
          ch = ltoc(ctol_cse(1,ns,ep), kv, actual_max_chan)
          if( ch.gt. 0 ) then
 
*             We have data see if it is OK.  If it is not then
*             set svs_OK false
              if( .not. data_OK(data_flag_cse(ch,ns,ep),
     .                    0,phs_mask) ) then
*                 Check to see if the bad data is at one of the
*                 epochs we are testing
                  do j = 1, num_scan
                     if( ep.eq. scan_ep(j) ) svs_OK = .false.
                  end do
              end if

*             See if bias_flag (on a good point) and this is
*             not the first epoch
              if( (kbit(data_flag_cse(ch,ns,ep),31) .or.
     .             kbit(data_flag_cse(ch,ns,ep),32) ) .and.
     .             ep .ne.begin_ep ) then
                  svs_OK = .false.
              end if
          else
*             See if we can ignore gaps and if point is one of
*             scanned values.
              if( ignore_gaps ) then
                  do j = 1, num_scan
                     if( j.eq. scan_ep(j) ) svs_OK = .false.
                  end do
              else
                  svs_OK = .false.
              end if
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE GET_SCAN_DD
 
      subroutine get_scan_dd( type, ns, lv, is, jv,
     .                ctol_cse, data_flag_cse, dd_OK,bias_flag)
 
      implicit none

*     Routine to scan for either a double difference either by
*     changing the station or the satellite; which is scanned is
*     defined by type (SITE or SVS or NONE if no serach to be made)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for one way
*   is, jv   - Site and Satellite in the index list dd_site/svs_search
*              It will increment from this value (making sure that it
*              does not repeat on the first site/svs in the double
*              difference..
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
 
      integer*4 ns, lv, is, jv,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   dd_OK       - Set OK if we find all the double differences
*               - need.
*   bias_flag   - Returns true if a bias flag is found between
*               - points in the scan list (if a signle point is
*               - tested then pushed bias on non-equally sampled
*               - data will not be found).
 
      logical dd_OK, bias_flag
 
*   type        - SITE if sites to be searched, SVS if
*               - satellites are to be searched and NONE
*               - If we are just to test this data combination
 
      character*(*) type
 
* LOCAL VARIABLES
 
*   js, kv  - Second site and source number
 
      integer*4 js, kv
 
*   finished - Indicates that we have finished getting the data
 
      logical finished
 
      logical debug_out 
      if( ns.eq.1 .and. lv.eq.3 ) then
         debug_out = .true.
      else
         debug_out = .false.
      endif
      debug_out = .false.

****  Start:
      dd_OK = .false.
      bias_flag = .false.
 
****  Treat each of the case separately.
*                                         ! No searching, just check
      if( type(1:2).eq.'NO' ) then
 
*         Loop from the first to last epoch in the scan_ep array
*         checking each point so that we don't miss any bias flags.
          js = dd_site_search(is)
          kv = dd_svs_search(jv)
          call test_scan_dd( ns, lv, js, kv, ctol_cse,
     .                 data_flag_cse, dd_OK,bias_flag)
*                             ! No searching
      end if
 
****  Scan for site using the current epochs and satellite
 
      if( type(1:2).eq.'SI' ) then
 
*         We should scan up from the current site
          finished = .false.
          if( is.ge.num_cfiles-1 ) finished = .true.
          do while ( .not. finished )
              is = is + 1

*             Make sure that this is not the first site in the
*             double difference
              js = dd_site_search(is)
              kv = dd_svs_search(jv)
              if( debug_out ) write(*,*) 'SCAN_DD site ',js, kv
              call test_scan_dd( ns, lv, js, kv, ctol_cse,
     .                        data_flag_cse, dd_OK,bias_flag)
              if( dd_OK .and. .not.bias_flag ) then
                  finished = .true.
              end if

              if( is.ge.num_cfiles-1 ) finished = .true.
          end do
      end if
 
****  Scan for satellite using the current epochs and stations
      if( type(1:2).eq.'SV' ) then
 
*         We should scan up from the current site
          finished = .false.
          if( jv.ge.num_sat-1 ) finished = .true.
          do while ( .not. finished )
              jv = jv + 1
             
*             Make sure that this is not the first site in the
*             double difference
              js = dd_site_search(is)
              kv = dd_svs_search(jv)
              if( debug_out ) write(*,*) 'SCAN_DD SVS  ',js, kv
              call test_scan_dd( ns, lv, js, kv, ctol_cse,
     .                    data_flag_cse, dd_OK,bias_flag)
              if( dd_OK .and. .not.bias_flag ) then
                  finished = .true.
              end if
              if( jv.ge.num_sat-1 ) finished = .true.
          end do
      end if
 
****  That's all
      return
      end
 
CTITLE TEST_SCAN_DD
 
      subroutine test_scan_dd( ns, lv, js, kv, ctol_cse,
     .                 data_flag_cse, dd_OK,bias_flag)

      implicit none
 
*     Routine to test a specific combination of stations and satellites
*     and return site_OK if all the data is there and no bias, and
*     bias_flag true if the data is OK but there is a bias.  (If
*     bias_flag is true then site_OK is always false).
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns, lv   - Site number and satellite list number for one way
*   js, kv   - Site and Satellite numbers of the second stationa and
*             site
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 ns, lv, js, kv,
     .    data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
*   dd_OK       - Set OK if we  find all the double differences
*               - need.
*   bias_flag   - Returns true if a bias flag is found between
*               - points in the scan list (if a signle point is
*               - tested then pushed bias on non-equally sampled
*               - data will not be found).
 
      logical dd_OK, bias_flag
 
* LOCAL VARIABLES
 
*   ep      - Epoch loop counter
*   i,j     - Loop variable 
*   c11,c12,c21,c22   - Channel numbers for the satellites at the second
*           - station
*   ltoc    - Function to return channel number
*   begin_ep - The earilest epoch in the scan list
*   end_ep   - Last epoch in list
 
      integer*4 ep, c11,c12,c21,c22, ltoc, begin_ep, end_ep, i,j
 
*   data_OK - Logical function which indicates that the data is
*           - OK.
*   kbit    - Logical function for checking bits in a word
 
 
      logical data_OK, kbit

      logical debug_out 
      if( ns.eq.1 .and. lv.eq.3 ) then
         debug_out = .true.
      else
         debug_out = .false.
      endif
         debug_out = .false.
 
****  Start:
      dd_OK = .true.
      bias_flag = .false.
      if( num_scan.eq.0 ) then
          dd_OK = .false.
          RETURN
      end if

*     Find the first epoch first so that we can ignore a bias flag
*     on this data
      begin_ep = scan_ep(1)
      end_ep = scan_ep(1)
      do i = 2, num_scan
         ep = scan_ep(i)
         if( ep.lt.begin_ep ) begin_ep = ep
         if( ep.gt.end_ep ) end_ep = ep
      end do
 
*     Loop from the first to last epoch in the scan_ep array
*     checking each point so that we don't miss any bias flags.
      if( debug_out) write(*,*) 'Scanning epochs ',begin_ep, end_ep,
     .               js, kv 
      do ep = begin_ep, end_ep
 
*         Get the channel numbers for his ep at the second site
*         and satellite
          c11 = ltoc(ctol_cse(1,ns,ep), lv, actual_max_chan)
          c12 = ltoc(ctol_cse(1,ns,ep), kv, actual_max_chan)
          c21 = ltoc(ctol_cse(1,js,ep), lv, actual_max_chan)
          c22 = ltoc(ctol_cse(1,js,ep), kv, actual_max_chan)
 
*         If both channels are present see if data OK
          if( c11.gt.0 .and. c12.gt.0 .and. 
     .        c21.gt.0 .and. c22.gt.0        ) then
 
*             See if the data is BAD (check both with and
*             with out the bias_flag).  If the data good then
*             check to see if the bias flag is set.
              if( data_OK(data_flag_cse(c11,ns,ep),0,phs_mask).and.
     .            data_OK(data_flag_cse(c12,ns,ep),0,phs_mask).and.
     .            data_OK(data_flag_cse(c21,js,ep),0,phs_mask).and.
     .            data_OK(data_flag_cse(c22,js,ep),0,phs_mask)
*                                             ! data is good
     .                                ) then
*                 Check the bias flags, but not on the first point
*                 Here we can have a bias flag.
                  if((kbit(data_flag_cse(c11,ns,ep),31) .or.
     .                kbit(data_flag_cse(c11,ns,ep),32) .or.
     .                kbit(data_flag_cse(c12,ns,ep),31) .or.
     .                kbit(data_flag_cse(c12,ns,ep),32) .or.
     .                kbit(data_flag_cse(c21,js,ep),31) .or.
     .                kbit(data_flag_cse(c21,js,ep),32) .or.
     .                kbit(data_flag_cse(c22,js,ep),31) .or.
     .                kbit(data_flag_cse(c22,js,ep),32)) .and.
     .                ep .ne. begin_ep    ) then
*                     There is a bias flag
                      bias_flag = .true.
                      dd_OK = .false.
                  end if
              else
*                 if we are ignoring gaps and this not one of the points
*                 points in the orginial one-way then don't set false.
                  if( ignore_gaps ) then
                      do j = 1, num_scan
                         if( ep.eq.scan_ep(j) ) dd_OK = .false.
                      end do
                  else
                      dd_OK = .false.
                  end if
              if( debug_out) write(*,*) 'Bad data: dd_OK ', dd_ok
              end if
*                     ! No data at all.
          else
*             if we are ignoring gaps and this not one of the points
*             points in the orginial one-way then don't set false.
              if( ignore_gaps ) then
                  do j = 1, num_scan
                     if( ep.eq.scan_ep(j) ) dd_OK = .false.
                  end do
              else
                  dd_OK = .false.
              end if
              if( debug_out) write(*,*) 'Missdata: dd_OK ', dd_ok
          end if
*                         ! Loop over the data
      end do
 
***** Thats all
      return
      end
 
 
CTITLE SCAN_CHK_2ND_SS
 
      subroutine scan_chk_2nd_ss(ns, lv, is1, jv1, is2, jv2, 
     .        bs, bv, sq, all_fnd, 
     .        L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse, 
     .        ctol_cse, data_flag_cse, checked )

      implicit none
 
*     This routine will scan to a secnd site and satellite to
*     compare with the first double difference to see which on the
*     onely-way data is bad.  If all_fnd is returned false then
*     we could not find another station and a satellite to go with
*     it thats could be compared to the original data.  If a
*     combination is found then bs and bv return the bad station
*     and satellite which then can be flagged.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ns,lv       - Original one-way being scanned
*   is1,jv1     - index in dd_site_search and dd_svs_search to
*               - to the double difference with the problem
*   is2,jv2     - index in dd_site_search and dd_svs_search for
*               - a second satellite and site to double diff with.
*   bs,bv       - The bad site and satellite numbers (not indices)

      integer*4 ns,lv, is1,jv1, is2,jv2, bs,bv
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4  data_flag_cse(num_chan, num_cfiles, num_ep),
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
 
*   sq(2)           - quality checker for combinations with
*                   - another station and satellite. (1) when
*                     satellite switched (2) when station switched.
 
      real*8 sq(2)
 
*   all_fnd     - Indicates that all combinations were found and we
*               - were able to check data
*   checked         - Indicates that the first two points are checked.
 
      logical all_fnd, checked
 
* LOCAL VARIABLES
 
*   js1, js2    - Second and third station numbers
*   kv1, kv2    - second and third satellite numbers
*   bet         - Epoch at which the bias flag should be placed.
 
      integer*4 js1, js2, kv1, kv2, bet
 
*   dd_OK(2)        - Set true when we find a valid double differnece
*                   - with this site/SVS.
*   qb(2)           - Quality bad combinations on a third
*                   - station and satellite (for determining
*                   - which one-way is bad).
*   bias_flag       - Indicates bias flag found in double differences
 
      logical dd_OK(2), qb(2), bias_flag
 
****  Start: Go the next station in the dd_site_search liast and
*     see if we find some data
 
      all_fnd = .true.
      js1 = dd_site_search(is1)
      kv1 = dd_svs_search(jv1)

*     Look for another satellite
      jv2 = jv1
      call get_scan_dd('SVS',ns, lv, is1, jv2,
     .            ctol_cse,data_flag_cse, dd_OK(1), bias_flag)
      kv2 = dd_svs_search(jv2)
      all_fnd = all_fnd .and. dd_OK(1)

*     Look for another station
      is2 = is1
      call get_scan_dd('SITE',ns, lv, is2, jv1,
     .                ctol_cse,data_flag_cse, dd_OK(2),
     .                bias_flag)
      js2 = dd_site_search(is2)
      all_fnd = all_fnd .and. dd_OK(2)
 
 
****  If we have found a second satellite and station then
*     check the quality.
*     Check each of the combinations
      if( dd_OK(1) ) call check_dd_quality( ns, lv, js1, kv2, bet,
     .        L1r_phs_cse, L2r_phs_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, sq(1), qb(1), checked )
      if( dd_OK(2) )  call check_dd_quality( ns, lv, js2, kv1, bet,
     .        L1r_phs_cse, L2r_phs_cse,
     .        L1_cyc_cse, L2_cyc_cse, ctol_cse,
     .        data_flag_cse, sq(2), qb(2), checked )

****  Now based on th which combinates looked bad, isolate the bad
*     station/satelite combination

*     Check sateliites
      if( dd_ok(1) ) then
          if( qb(1) ) then
*              Switching satellites did not help.  Must be first satelite
               bv = lv
          else
*              Looks good, must be second satellite
               bv = kv1
          end if
      else
          bv = 0
      end if

*     Check stations
      if( dd_ok(2) ) then
          if( qb(2) ) then
*              Switching sites did not help.  Must be first site
               bs = ns
          else
*              Looks good, must be second site
               bs = js1
          end if
      else
          bs = 0
      end if
 
****  Thats all
      return
      end
