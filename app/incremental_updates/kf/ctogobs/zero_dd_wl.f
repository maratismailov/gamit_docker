CTITLE zero_dd_wl
 
      subroutine zero_dd_wl(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, bf_type_cse, azel_cse)

      implicit none
 
*     This routine will use the bf_table from flat_dd to sequentially
*     resolve the wide lane ambiguities to integers.  If this can be
*     done reliably then the WL ambiquities can be set to zero in solve
*     rather than using the ionospheric constraint as is currently done.
*
*     During these calculations the mean and quality of L1, L2, L2-L1 and
*     LC phase only values are assessed.  On short baselines, these may
*     provide more reliable estimates than the MW-WL using psuedornage.
*
*     The sequence starts with the shortest baseline and works it way 
*     up from this point.

* MOD TAH 031221: Modified algorithm to make the initial acceptance of bias stringent 
*     and then reduce the criteria with each iteration.
 
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
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep)
     .,    ctol_cse(num_chan, num_cfiles, num_ep)
     .,    bf_type_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
       real*4 azel_cse(2,num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
      integer*4 shortest_pair(2)  ! Two stations with the shortest 
                                  ! baseline length
     .,         llbf_sv     ! Satellite for long last bias flag at shortest
                            ! baseline pair of stations
     .,         llbf_dur(2) ! Duration of long last bias flag for pair
     .,         llbf_bf(2)  ! BF_TABLE entry of long last bias for pair
     .,         overlap     ! Number of epochs of overlap between bias falgs
     .,         overnext    ! Overlap for next biasflag on the same satellite
     .,         eps, epe    ! Start and stop epochs of overlap in double diffs
     .,         i, is, j    ! Loop counter over bf_table entries and sort index
     .,         bfs, bfj, bfv  ! BF_table entries for 2nd satellite and 2nd site to 
                            ! form double differences
     .,         chk_overlap ! Function to check overlap of bias flags
     .,         num_new, num_resd  ! Number of new biases resolved and total number
                            ! resolved.
     .,         last_bias   ! Number of bf_table entry for the last bias flag
                            ! on a satellite at a site.
     .,         ns, lv      ! Site and satellite numbers
     .,         new_ref     ! new reference bias to be fixed
     .,         longest_ref ! Longest sequence of bias that is not fixable
     .,         ep          ! epoch counter
     .,         num_wl_resd, num_nl_resd ! Number of resolved biases
     .,         bf_sort(max_bf_table)  ! Sorted reference list to bf_table
     .,         ierr        ! IOSTAT error flag
     .,         iter        ! Iteration counter (starts at 4 and works down)

      logical all_resolved  ! Set true if all bias flags are resolved. Mainly
                            ! use to indicate too many bias flags
     .,       kbit          ! Bit testing routine
     .,       new           ! Logical set true when ever new biases are resolved
     .,       apply         ! Logical set true is we mark the reference satellite
                            ! as fixed (not marked fixed if the next bias flag is
                            ! on the same satellite and overlaps with the reference)
     .,       wl_resd, nl_resd       ! Set true if widelane or narrowlane is resolved
     .,       some_unconnect  ! Logical set true if there are some un-connected
                            ! bias flags that do not overlap with resolvable
                            ! bias flags
     .,       single_ref    ! Set true if single reference satellite can be found (visible
                            ! at all sites)

      character*256 stat    ! Status message saying which site and satellite to start
                            ! with for wide-lane bias fixing.
      character*4 type      ! Sets type of new reference selected.  CON overlaps
                            ! with a fixed bias, UNC is unconnected and does not
                            ! overlap with a resolved bias.

***** Variables for dynamic memory allocation for estimation of biases
C      integer*4 cma_data(1)
C      integer*8 ineq, ibvec   ! Memory slots for normal eqautions and bvec

****  Start: Compute the baseline lengths betweeen all the sites
*     and return the pair of stations with the shortest length.  
      call get_blens( shortest_pair)

***** Remove data that has no double differences
      call scan_nodd( ctol_cse, data_flag_cse, bf_type_cse) 

****  Rescan the bias flags to make the bf_table (this is in case some
*     bias flags were removed in pf_remove_bias operation.  This operation may
*     move to before trim_oneways in ctogobs main.
      call get_biases( data_flag_cse, ctol_cse, all_resolved, 'REP' )

***** Scan the MW-WL in the oneway bias flags and report on quality. 
*     Here we are looking for WL's that are systematic.  We also set
*     a correlation time which is used to get the sigma of the means
*     Bit 5 of the bf_table(5,*) is set to mark those widelanes that
*     are systematic enough not to be used in the resolution.
      call scan_wl( L1r_phs_cse, L2r_phs_cse,
     .             L1_cyc_cse, L2_cyc_cse, 
     .             L1r_rng_cse, L2r_rng_cse, 
     .             ctol_cse, data_flag_cse, azel_cse)
     
*
*     See which station and satellite we will choose as base.  The
*     measure here is the longest last bias flag.  Satellite and 
*     length is returned
      single_ref = .true.
      call find_long_last_bf( shortest_pair, llbf_sv, 
     .                        llbf_dur,llbf_bf, single_ref)
c     call find_long_last_bf( shortest_pair(2), llbf_sv(2), 
c    .                        llbf_dur(2),llbf_bf(2), single_ref)

***** See if we found a satellite
      if( llbf_sv.eq.0 ) then
          single_ref = .false.
      endif
      if( .not.single_ref ) then
         call find_long_last_bf( shortest_pair, llbf_sv, 
     .                        llbf_dur,llbf_bf, single_ref)
c        call find_long_last_bf( shortest_pair(2), llbf_sv(2), 
c    .                        llbf_dur(2),llbf_bf(2), single_ref)
      end if


*     See which is better (more quality tests could be add here
*     at the moment use just the length)
      wl_ref_svs  = llbf_sv
      if( llbf_dur(1).gt.llbf_dur(2) ) then
         wl_ref_site = shortest_pair(1)
         wl_ref_dur  = llbf_dur(1)
         wl_ref_bf   = llbf_bf(1)
      else
         wl_ref_site = shortest_pair(2)
         wl_ref_dur  = llbf_dur(2)
         wl_ref_bf   = llbf_bf(2)
      endif
      write(stat,120) cf_codes(wl_ref_site), prn_list(wl_ref_svs),
     .                wl_ref_bf, wl_ref_dur
 120  format('WL DD Reference site ',a4,' PRN_',i2.2,' BF ',i4,
     .       ' Duration ',I4,' Epochs')
      call report_stat('status','autcln','Zero_dd_wl',' ',stat,0)

****  Since we have the last bias flag at these stations, there should
*     not be a problem with complete lose of lock at one of these 
*     stations (i.e., an epoch when there is a bias flag to all 
*     satellites).  If this condition exists then there will be
*     biases that can not be resolved with out selecting another
*     reference satellite. Later code might fix this problem if we
*     want to resolve all biases (currently solve only resolves last
*     bias so this is OK).

*     Mark all the biases at reference station as resolved
      do i = 1, tot_bf
         if( bf_table(2,i).eq.  wl_ref_site ) then
*          Code needs to be modified to not fix all since there
*          can be gaps
* MOD TAH 040507:
*          Find the bias flag that over laps with the reference
*          satellite.  If there are more than one, choose the last
*          one
           overlap = chk_overlap(i, wl_ref_bf, 0, 0, eps, epe)
           if( overlap.gt.0 ) then
               apply = .true.
*              Make sure next bias flag is not overlapping and
*              on same satellite
               if( i.lt.tot_bf ) then
                  if( bf_table(3,i+1).eq.bf_table(3,i) ) then
                      overnext = chk_overlap(i+1, wl_ref_bf, 0, 0, 
     .                                         eps, epe)
*                     For the moment, just take the next one.  Really
*                     should take the longest one.  Implement this
*                     later.
                      if( overnext.gt.0 ) apply = .false.
                  endif
               endif
               if( apply ) then     
                  call sbit(bf_table(5,i),1,1)
                  call sbit(bf_table(5,i),2,1)
                  call sbit(bf_table(5,i),3,1)
                  call sbit(bf_table(5,i),4,1)
               endif
           endif
         elseif( bf_table(3,i).eq.wl_ref_svs ) then 
*            See if this site and satellite overlaps with reference
             overlap = chk_overlap(i, wl_ref_bf, 0, 0, eps, epe)
             if( overlap.gt.0 ) then
* MOD TAH 031220: Check the next bias flag and see if it is the same
*                satellite and also has overlap
                 apply = .true.
                 if( i.lt.tot_bf ) then
                     if( bf_table(3,i+1).eq.bf_table(3,i) ) then
                         overnext = chk_overlap(i+1, wl_ref_bf, 0, 0, 
     .                                          eps, epe)
*                        For the moment, just take the next one.  Really
*                        should take the longest one.  Implement this
*                        later.
                         if( overnext.gt.0 ) apply = .false.
                     endif
                 endif
*                OK: If the next bf_table is different satellite
*                or there is no overlap on the next bias flag with this
*                satellite, apply the bias fixing.
                 if( apply ) then
                    call sbit(bf_table(5,i),1,1)
                    call sbit(bf_table(5,i),2,1)
                    call sbit(bf_table(5,i),3,1)
                    call sbit(bf_table(5,i),4,1)
                 endif
             endif
         endif
      enddo

****  OK, all of the bias flags that we can set arbitarily as fixed are
*     now set.  Now loop over the bias flags seeing which ones can be
*     resolved
      new = .true.
      call sort_blen(wl_ref_site, bf_sort)
      iter = 2
      do while ( new )
         num_new = 0
         num_resd = 0
         do is = 1, tot_bf
*           See if this bias is already resolved or of it has been
*           marked as un-resolvable.  Get the blen sorted bf_table
*           entry 
            i = bf_sort(is)
            if( i.eq.0 ) then
                write(*,150) i
 150            format('WARNING: BF_SORT entry ',i5,' zero')
                i = is 
            end if
C           i = is
            if( .not. kbit(bf_table(5,i),1)  ) then
*              OK, this bias flag not resolved.  See if we can find another
*              satellite in this interval at this site which is resolved.  If
*              this found find another site with same satellite pair.  bf_table
*              entries are returned. bfs -- Second satellite, bjf second site
*              same reference SV, bfv -- second site second satellite.  
*              EPS and EPE are start and stop epochs of the DD overlap
*              Start at first bias flag at this site: bfs is set to zero
*              either when bias fixed or we run out of flags to test (bfs is
*              is returned zero from find_fix_bfs
               bfs = bf_index(bf_table(2,i))
               bfj = 0
               do while ( bfs.gt.0 )
                  call find_fix_bfs(i, bfs, bfj, bfv, eps, epe)
* 
*                 See if we can resolve WL and NL.  Do not attempt
*                 to resolve biases at are systematic
                  if( bfs.gt.0 .and.
     .               .not. kbit(bf_table(5,i),5) ) then
                     call resolve_wlnl(i, bfs, bfj, bfv, eps, epe, 
     .                     L1r_phs_cse, L2r_phs_cse,
     .                     L1_cyc_cse, L2_cyc_cse, 
     .                     L1r_rng_cse, L2r_rng_cse, 
     .                     ctol_cse, data_flag_cse, azel_cse,
     .                     wl_resd, nl_resd, iter)
                     if( wl_resd ) then
                         num_new = num_new + 1
                         num_resd = num_resd + 1
                         bfs = 0
                         call sbit(bf_table(5,i),1,1)
                     endif
                     if( nl_resd ) then
                         call sbit(bf_table(5,i),4,1)
                     endif 

                  end if
               enddo
            else
               num_resd = num_resd + 1
            end if
         enddo
*        See if we have resolved all biases or if no new ones were resolved
         if( num_resd.eq.tot_bf ) new = .false.
*        If no biases were fixed, then we may still have biases that
*        can not be linked to others so scan the biases to see which
*        ones have not been used in combination with other fixable biases
         if( num_new.eq.0 .and. new ) then
*           Implement later.  Ultimately when code is done new equal
*           to false will be result at end.
            new = .false.
*           Check to see if we have bias flags that could not be fixed
*           because they have never appeared with a fixable double
*           differnce combination.  When such a bias is found, set a
*           new reference bias flag and see how many more we can fix
            longest_ref = 0
            new_ref = 0
            some_unconnect = .false.
            type = 'CON'
            do i = 1, tot_bf
*              Skip any bias marked as unresolvable
               if( .not.kbit(bf_table(5,i),3) .and.
     .             .not.kbit(bf_table(5,i),5)     ) then
*                 Make sure bias is not fixed (it should not be)
                  if( kbit(bf_table(5,i),1) ) then
                      write(*,220) i, bf_table(5,i)
 220                  format('**WARNING** Bias flag ',i4,' is fixed ',
     .                       'but does appear with fixable biases: ',
     .                       'Status ',o3)
                  else 
*                    See how this bias works with other biases
*                    Make sure it overlaps with a resolved satellite
*                    bias
                     some_unconnect = .true.
                     do j = 1, tot_bf
                       if( kbit(bf_table(5,j),1) .and.
     .                     bf_table(2,i).ne.bf_table(2,j) .and.
     .                     bf_table(3,i).eq.bf_table(3,j) ) then
                           overlap = chk_overlap(i, j, 0, 0, 
     .                                          eps, epe)
                           if (overlap.gt.longest_ref ) then
                                new_ref = i
                                longest_ref = overlap
                           endif
                       end if
                     end do  
                  endif
               end if
            end do

****        See if we some un-connected biases that do not appear with
*           any unresolved biases
            if( new_ref.eq.0 .and. some_unconnect ) then
*              Rescan again, but this time do not put on restriction
*              that the bias must overlap with a resolved bias
               type = 'UNC'
               do i = 1, tot_bf
*                 Skip any bias marked as resolvable (and are not systematic)
                  if( .not.kbit(bf_table(5,i),3) .and.
     .                .not.kbit(bf_table(5,i),5)     ) then
*                     Select the longest bias flag that is available
                      if( bf_table(4,i).gt.longest_ref ) then
                          new_ref = i
                          longest_ref = bf_table(4,i)
                      end if
                  endif
               end do
            endif

****        Report the new fix bias.
            if( new_ref.gt.0 ) then
               ns = bf_table(2,new_ref)
               lv = bf_table(3,new_ref)
*              If there is debate about whether a bias can be set
*              as a new reference, the routine below explicitly 
*              checks all possible combinations and output combinations
*              that have resolvable biases (and thus new ref should not
*              be set).
               call verify_newref(new_ref, new)
               write(stat,240) cf_codes(ns), prn_list(lv),
     .                   new_ref, bf_table(4,new_ref), longest_ref, type
 240           format('WL DD New Ref bias ',a4,' PRN_',i2.2,' BF ',i4,
     .          ' Duration ',I4,' Overlap ',i4,' Eps ',a)
               call report_stat('status','autcln','Zero_dd_wl',' ',
     .                           stat,0)
               if( new ) then
                  call sbit(bf_table(5,new_ref),1,1)
               end if
               call sbit(bf_table(5,new_ref),2,1)
               call sbit(bf_table(5,new_ref),3,1)
               new = .true.
            else

****           See if we have more iterations to complete
               iter = iter - 1
               if( iter.gt.0 ) then
                   new = .true.
                   write(*,260) iter
 260               format('Starting Iteration ',i3)
               else
                   new = .false.
               endif
            endif
         endif
      end do
*
*     OK: We are done.  Now output the status of the bias fixing for use
*     solve.  We do this be baseline length to be consistent with what solve
*     expects. (May change later as solve gets more sophisticated).
***  Tell user what we found
      write(*,310) tot_bf, (cf_codes(i), bf_index(i),i=1, num_cfiles)
 310  format('Total of ',i4,' bias flags, By site starts are: ',/,
     .       50(8(a4,1x,I4,1x),:/))

      num_wl_resd = 0
      num_nl_resd = 0     
      do i = 1, tot_bf
         ns = bf_table(2,i)
         lv = bf_table(3,i)
         ep = bf_table(1,i)
         if( 1.eq.2 )
     .   write(*,320) i,cf_codes(ns), prn_list(lv), bf_table(1,i),
     .                bf_table(4,i), bf_table(5,i), 
     .                wl_conf(i), nl_conf(i)
 320     format('BF_ZERO ',i4,1x,a4,1x,'PRN ',i2.2,' EP ',i4,1x,i4,
     .          1x,'S',1x,o4,1x,2F10.2)
         if( kbit(bf_table(5,i),1) ) num_wl_resd = num_wl_resd + 1
         if( kbit(bf_table(5,i),4) ) num_nl_resd = num_nl_resd + 1

      end do
      write(*,340)  tot_bf, num_wl_resd,  num_nl_resd
 340  format('BF_ZERO: Of ',i4,' Biases; ',i4,' WL resolved ',
     .       I4,' NL resolved')

*     Create and output the last bias flag status to the acbias.dat
*     file
      if( acbias_out ) then
         open(205,file='acbias.dat',status='unknown',iostat=ierr)
         call report_error('IOSTAT',ierr,'open','acbias.dat',0,
     .                  'ZERO_DD_WL')
     
         if( ierr.eq.0 ) then
*           Write the reference list and PRN number as first line
            if( single_ref ) then
               write(205,410) wl_ref_svs, prn_list(wl_ref_svs)
            else
               write(205,410) 0,0
            endif
 410        format('REFERENCE ',i2,' PRN ',i2.2)
            do i = 1,tot_bf
               call get_last_bias(i,last_bias)
               if( i.eq.last_bias ) then
                   if( kbit(bf_table(5,i),1) ) then
                       type = 'X'
                   else
                       type = 'R'
                   endif
                   write(205,420,iostat=ierr) i,cf_codes(bf_table(2,i)),
     .                    prn_list(bf_table(3,i)), type
 420               format(i5,1x,a4,' PRN ',i2.2,1x,a1)
               endif
            end do
            close(205)
         end if
      end if

****  Now form and write the double difference biases to the summary file
      call report_stat('status','autcln','Zero_dd_wl',' ',
     .                           'Start generate GAMIT DD',0)
      call gen_gamit_dd(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse)

      call report_stat('status','autcln','Zero_dd_wl',' ',
     .                           'Finised generate GAMIT DD',0)

****  Now check all the double differences
*     code below is still experimental.  Return at this point
      RETURN

c      call report_stat('status','autcln','Zero_dd_wl',' ',
c     .                           'Start checking all DD',0)

***   Allocate memory for normal equations
C     call cma_alloc(cma_data, ineq, ibvec, tot_bf )

C     call clear_bfneq(cma_data(ineq),cma_data(ibvec), tot_bf)

C     tot_dd = 0
C     do i = 1, tot_bf
C        call scan_alldd(i, L1r_phs_cse, L2r_phs_cse,
C    .       L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
C    .       ctol_cse, data_flag_cse, azel_cse, 
C    .       cma_data(ineq),cma_data(ibvec), tot_dd)
C     end do
C     write(stat,510) tot_dd
C510  format('Done checking ',i8,' DD')
C     call report_stat('status','autcln','Zero_dd_wl',' ',stat,0)

****  Now solve the normal equations.  Use the fixed BF_entries to
*     set which bias flags can be set to zero.
C     call solve_bfneq(cma_data(ineq),cma_data(ibvec), 
C    .       L1_cyc_cse, L2_cyc_cse, ctol_cse)
C

C     return
      end

         
CTITLE GET_BLENS

      subroutine get_blens( shortest_pair)

      implicit none

*     Rouitine compute the baseline lengths and return the 
*     pair of stations with the shortest length

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 shortest_pair(2)  ! Two stations with the shortest len

* LOCAL VARIABLES
      integer*4 i,j  ! Loop counters over stations
     .,    ent       ! Entry in lower diagonal matrix of baselens
     .,    bl_ent    ! Function to return baseline entry

      real*8 lat, lng, rad  ! Lat, long (rads) and radius (m) extracted
                     ! from cf_apr_save
     .,    xi(3), xj(3)  ! XYZ of two stations
     .,    bl        ! Length between stations
     .,    min_bl    ! Shortest baseline length

****  Loop over the stations computing and saving lengths.  The 
*     cf_apr_save array is used (compiled from cf_preval).  For 
*     each site the values are lat, long and radius.

      min_bl = 1.d10
      do i = 1, num_cfiles-1
*        Get the XYZ coordinates of this site
         lat = cf_apr_save(1,i)
         lng = cf_apr_save(2,i)
         rad = cf_apr_save(3,i)*1000.d0
         xi(1) = rad*cos(lat)*cos(lng)
         xi(2) = rad*cos(lat)*sin(lng)
         xi(3) = rad*sin(lat)

         do j = i+1, num_cfiles
            lat = cf_apr_save(1,j)
            lng = cf_apr_save(2,j)
            rad = cf_apr_save(3,j)*1000.d0
            xj(1) = rad*cos(lat)*cos(lng)
            xj(2) = rad*cos(lat)*sin(lng)
            xj(3) = rad*sin(lat)
            ent = bl_ent(i,j)
            bl = sqrt((xi(1)-xj(1))**2+(xi(2)-xj(2))**2+
     .                (xi(3)-xj(3))**2)
            baselens(ent) = bl
            orderlen(ent) = bl

*           Now compute the orderlen that depends of the 
*           quality of the MW-widelanes
            if( lc_num(i).eq.0 .or. lc_num(j).eq.0 ) then
                orderlen(ent) = baselens(ent) + 200.d6
            endif
            if( WL_RMS(i).gt.1.d0 .or. WL_RMS(j).gt.1.d0 ) then
                orderlen(ent) = baselens(ent) + 15.d6*
     .             max(WL_RMS(i),WL_RMS(j))
            endif

****        Sort on orderlen
            bl = orderlen(ent)

*           Check to see if this is the smallest
            if( bl.lt.min_bl ) then
                min_bl = bl
                shortest_pair(1) = i
                shortest_pair(2) = j
            end if
         end do
      end do

***** Thats all
      return
      end

CTITLE ENT_BF
      subroutine ent_bf(ent,num_cfiles,pair)

      implicit none

*     Inverse of bl_ent, given baseline entry returns the pair of sites
*     Brute force at the moment could be spead up if needed
      integer*4 ent       ! Baseline entry in lower triangular form
     .,         pair(2)   ! Pair of sites in baseline
* LOCAL
      integer*4 i,j, bl_ent, num_cfiles

*     Brute force search
      do i = 1, num_cfiles-1
         do j = i+1, num_cfiles
             if( bl_ent(i,j).eq.ent ) then
                 pair(1) = i
                 pair(2) = j
             endif
         end do
      end do

***** Thats all
      return
      end

CTITLE FIND_LONG_LAST_BF

      subroutine find_long_last_bf( sites, llbf_sv, llbf_dur,llbf_bf, 
     .                              single_ref)

      implicit none

*     Routine to look at the bais flags for a pair of sites and return
*     the satellite, duration (epochs) and bias flag number of
*     longest last bias flag at the stations that over lap.
*     The single_ref check to see if a single reference satellite (one
*     which is seen by all stations is available.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 sites(2)  ! Site numbers of baseline to check
     .,         llbf_sv  ! Satellite for long last
     .,         llbf_dur(2) ! Duration of long last
     .,         llbf_bf(2)   ! BF_TABLE entry of long last bias

      logical single_ref  ! Set true if single reference satellite

* LOCAL VARIABLES
      integer*4 i,j   ! Loop over bias flags
     .,    max_dur    ! Longest duration bias flag found
     .,    s2_ent     ! Bias table entry for second site
     .,    chk_overlap  ! Function to return overlap between bf_table
                      ! entries
     .,    overlap    ! Overlap between entries
     .,    eps, epe   ! Epoch numbers of overlap.


      logical check  ! Set true if we should duration: Station
                     ! match and either last bias flag or
                     ! statellite changes


****  Loop over the bf_tables finding entries for this station
*     and getting longest, last bias flag
      llbf_sv = 0
      llbf_dur(1) = 0
      llbf_dur(2) = 0
      llbf_bf(1) = 0
      llbf_bf(2) = 0
      max_dur = 1
      do i = 1, tot_bf
*        See if site matches
         if( bf_table(2,i).eq.sites(1) ) then
*            Correct site.  See if this a last bias flag
*            for a satellites.  Check to see of the next
*            entry is last satellite. (Use check construct
*            in case we get to end of table and would check
*            an out-of-bounds entry)
             check = .false.
             if( i.eq.tot_bf ) then
                check = .true.
             elseif( bf_table(3,i+1).ne. bf_table(3,i) .or.
     .               bf_table(2,i+1).ne. bf_table(2,i) ) then
*               Either satellite or station changes
                check = .true.
             endif

*            Make sure that this satellite is visible at all sites and
*            has at least an 30-minutes of data at 30-sec sampling
             if( single_ref ) then
                do j = 1, num_cfiles
                   if( LC_svs_num(bf_table(3,i),j).lt.60 ) then
                      check = .false.
                   end if
                end do
             end if


             if( check ) then

****             Now check the second site and see what the overlap
*                if.
                 s2_ent = 0
                 do j = 1,tot_bf
                     if( bf_table(2,j).eq.sites(2) .and. 
     .                   bf_table(3,j).eq.bf_table(3,i) ) then
*                        Site and satellite match.  See if this
*                        is the last one at the second site for
*                        this satellite
                         if( j.eq.tot_bf ) then
                             check = .true.
                             s2_ent = j
                        elseif ( bf_table(3,j+1).ne. bf_table(3,j) .or.
     .                         bf_table(2,j+1).ne. bf_table(2,j) ) then
*                           Either satellite or station changes
                            check = .true.
                            s2_ent = j
                         endif
                     end if
                 end do  
              endif

*             See if have a mctching entry and that we should check
              if( check .and. s2_ent.gt.0 ) then

*                This is last entry for this satellite at both sites
                 overlap = chk_overlap(i, s2_ent, 0, 0, eps, epe)
 
                 if( overlap.gt.max_dur ) then
                     llbf_sv = bf_table(3,i)
                     llbf_dur(1) = bf_table(4,i)
                     llbf_bf(1) = i
                     llbf_dur(2) = bf_table(4,s2_ent)
                     llbf_bf(2) = s2_ent
                     max_dur = overlap
                 end if
             endif
         endif
      enddo

****  Thats all
      return
      end 

 
 
CTITLE FIND_FIX_BFS

      subroutine find_fix_bfs(bft, bfs, bfj, bfv, eps, epe)

      implicit none

*     Routine to find another bias that is already fixed that
*     we can form double differences with. BFT is the table
*     entry that we are testing.  If a double difference can 
*     be found then the other bf_table entries that form the
*     double difference are returned along with the start
*     and stop epoch numbers of the overlap.

*     If overlap exists but not fixed baises can be found
*     the this is marked (so that we can tell which biases can
*     never be determined)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 bft  ! Test bf_table entry
     .,         bfs  ! Matching bf_table entry at the same
                     ! site but with satelite that is resolved.  
                     ! Current starting point is passed in and
                     ! updated
     .,         bfj  ! Bf_table entry at different site but
                     ! with overlap with initial satellite
     .,         bfv  ! bf_table entry with second site and satellite
     .,         eps, epe  ! Start and stop epoch numbers of overlap.

* LOCAL VARIABLES
      integer*4 s1, s2, v1, v2  ! Sites and satellites in current
                      ! double difference
     .,         i,j, k  ! Loop counters over bias flags
     .,         overlap ! Length of overlap in double difference
     .,         chk_overlap  ! Function to check over lap between bias flags
     .,         start_bf, end_bf  ! Start and end entries in BF_table to check

      logical kbit    ! Tests if bit is set
     .,       more    ! Logical, remains true while more satellites to
                      ! check at station being tested
     .,       cont    ! Remains true more other stations to test


****  Get the site and satellite from test bf_table entry
      s1 = bf_table(2,bft)
      v1 = bf_table(3,bft)

***** Loop over other satellites at this station to see if
*     can find another satellite with biases fixed
      i = bfs - 1
      more = .true.
      do while ( more )
         i = i + 1
*        Save the new entry in case for when we re-enter this routine
         bfs = i
         if( bf_table(2,i).eq.s1 .and. i.ne. bft ) then
*           OK, found an entry at this station.  See if
*           there is over lap in time
            overlap = chk_overlap(i, bft,0,0, eps,epe)

            if( overlap.gt.0 ) then
*               There is overlap, so this one is possible
*               see if can find another staton with this 
*               same pair of satellites
                v2 = bf_table(3,i)
*               Continue the loop from where we left off
                j = bfj
                cont = .true.
                do while ( cont .and. j.lt.tot_bf )
                   j = j + 1
*                  Check different station and if sats match
                   if( bf_table(2,j).ne.s1 .and. 
     .                 bf_table(3,j).eq.v1 ) then
*                      OK, different station, SV matches test entry
                       overlap = chk_overlap(i, j, bft,0, eps, epe)
*                      See if we have overlap.  If so check the
*                      other entries at this second site to see
*                      if satellites match
 
                       if( overlap.gt.0 ) then
                           s2 = bf_table(2,j)
*                          Check BF_table for this site
                           start_bf = bf_index(s2)
                           end_bf = tot_bf
                           if( s2.lt.num_cfiles ) then
                               end_bf = bf_index(s2+1)-1
                           endif
*                          Loop over all
                           k = start_bf - 1
                           do while ( k.lt.end_bf )
                              k = k + 1
*                             We are now scanning for second
*                             satellite at the second station
                              if( bf_table(3,k).eq.v2 ) then
*                                 See if overlap in time
                                  overlap = chk_overlap(k,j,i,bft,
     .                                      eps,epe)
                                  if( overlap.gt.0 ) then
*                                     OK: We have found a double 
*                                     combination which works
                                      bfj = j
                                      bfv = k
*                                     Mark the test bf_table entry
*                                     as having double differences.
*                                     If all the bias flags in this
*                                     double difference are fixable
*                                     then set this one as fixable
                                      if( kbit(bf_table(5,bfj),3)
     .                                    .and.
     .                                    kbit(bf_table(5,bfv),3)
     .                                   .and.
     .                                    kbit(bf_table(5,bfs),3) ) then
                                          call sbit(bf_table(5,bft),3,1)
                                      endif
*                                     See if current combination 
*                                     have resolved biases
                                      if( kbit(bf_table(5,bfj),1) 
     .                                   .and.
     .                                    kbit(bf_table(5,bfv),1)
     .                                   .and.
     .                                    kbit(bf_table(5,bfs),1) ) then
                                          more = .false.
                                          cont = .false.
*                                         Mark bf_table to show that there
*                                         fixed biases which this one could
*                                         resolved against
                                          call sbit(bf_table(5,bft),2,1)
                                          RETURN
                                      end if
                                  endif
                              endif
                           end do
*                          We have looped over all the Bias flags
*                          at this second site and not found anything
                       endif   ! Overlap with second satellite
                   end if      ! This is not a different station or not the same
                               ! satellite
                end do
            end if
         end if
*        We have finished looping over all stations with this second
*        satellite at the reference station.  We now need to move on
*        another satellite at the reference station.
*        Move to the next entry at this station
*        See if we have run out of entries at this station
         if( bf_table(2,i+1).ne.s1 ) then
*            We have run out of entries.  Set bfs = 0 to denote 
*            this
             bfs = 0
             bfj = 0
             more = .false.
         else
*           Reset the scanning station back to zero
            bfj = 0
         end if
      end do

***** Thats all.  If we get to here then we have run out of entries
      bfs = 0
      return
      end

CTITLE RESOLVE_WLNL
      subroutine  resolve_wlnl(bft, bfs, bfj, bfv, eps, epe,
     .                     L1r_phs_cse, L2r_phs_cse,
     .                     L1_cyc_cse, L2_cyc_cse, 
     .                     L1r_rng_cse, L2r_rng_cse, 
     .                     ctol_cse, data_flag_cse, azel_cse,
     .                     wl_resd, nl_resd, iter)

      implicit none

*     Routine to resolve the L2-L1 bias by combination of 
*     MW-WL, L2-L1 phase and L1 and L2 average phase by itself.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 bft ! Current bf_table being tested
     .,         bfs ! bf_table entry at same site but different satellite
     .,         bfj ! bf_table entry at different site but same first satellite
     .,         bfv ! bf_table entry at 2nd site and 2nd satellite
     .,         eps, epe   ! Overlap of epochs in double difference
     .,         iter ! Iteration countdown.  With higher values, criteria are
                     ! more strict (to ensure that base bias fixing is correct)

      logical wl_resd, nl_resd  ! Set true oif WL or NL can be resolved to 
                    ! integer values

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

* LOCAL VARIABLES

      integer*4 i,j
     .,    ep          ! Epoch counter
     .,    ns,lv       ! Site and satellite list number
     .,    ltoc        ! Returns channel number for specific satellite
     .,    ch          ! Channel number being fixed

      real*8 obs(5,4)  ! Obs from each oneway (columns are site
                        ! satellite combinations, rows are data
                        ! types
     .,      wgh(5,4)   ! Weight for each station/satellite
     .,      stats(4,5) ! Statistics: Columns for each data type
                         ! rows for res*wgh res**2*wgh wgh and num
     .,      dL1_cyc, dL2_cyc  ! Change to number of cycles

       logical OK(4)    ! set true for each good station/satellite one-way
       logical kbit

****  Mark the resolved status as false.
      wl_resd = .false.
      nl_resd = .false.

****  Now starting accumulating data that allows us to fix the
*     L2-L1 and maybe the L1 bias.

*     Clear the accumulation arrays
      do i = 1,4  ! Loop over res*wgh, res**2*wgh wgh num 
         do j = 1,5  ! Loop over MW-WL, L2-L1, L1 and L2, LC
             stats(i,j) = 0.d0
         end do
      end do

*    Loop over the epoch range
      do ep = eps,epe
*       Get the residuals for each of one-ways that go into
*       the double difference
         call get_obs(ep, bf_table(2,bft),bf_table(3,bft),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,1), wgh(1,1), OK(1))
         call get_obs(ep, bf_table(2,bfs),bf_table(3,bfs),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,2), wgh(1,2), OK(2))
         call get_obs(ep, bf_table(2,bfj),bf_table(3,bfj),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,3), wgh(1,3), OK(3))
         call get_obs(ep, bf_table(2,bfv),bf_table(3,bfv),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,4), wgh(1,4), OK(4))

****     If all the oneways are good, increment the statistics
         if( OK(1) .and. OK(2) .and. OK(3) .and. OK(4) ) then
              call inc_zero_stats(obs, wgh, stats)
         endif
      end do

****  Finished accumulating statistics.  Get the means and see if
*     can resolve the integer values
      call fin_zero_stats(bft, bfv,stats, dL1_cyc, dL2_cyc,iter)
      if( wl_conf(bft).gt.dchi_wl_tol*iter ) wl_resd = .true.
      if( nl_conf(bft).gt.dchi_wl_tol*iter ) nl_resd = .true.

****  If the wl_resd is reliable update the L2_cyc_cse values
      if ( wl_resd .and. (dL2_cyc.ne.0.d0 .or. dL1_cyc.ne.0.d0) ) then
          ns = bf_table(2,bft)
          lv = bf_table(3,bft)
          if( kbit(status_rep,15) )
     .    write(*,220) cf_codes(ns),prn_list(lv), 
     .          bf_table(1,bft),bf_table(1,bft)+bf_table(4,bft),
     .          dL1_cyc,dL2_cyc
 220      format('ZUPD ',a4,' PRN ',i2.2,' Epochs ',i4,1x,i4,
     .           ' dL12 ',2F6.1)
          do ep = bf_table(1,bft),bf_table(1,bft)+bf_table(4,bft)
             ch = ltoc(ctol_cse(1,ns,ep),lv,actual_max_chan)
             if( ch.gt.0 ) then
                L1_cyc_cse(ch,ns,ep) = L1_cyc_cse(ch,ns,ep)-dL1_cyc
                L2_cyc_cse(ch,ns,ep) = L2_cyc_cse(ch,ns,ep)-dL2_cyc
             endif
          end do
      end if 

***** Thats all        
      return
      end
 

CTITLE CHK_OVERLAP

      integer*4 function chk_overlap(bfi, bfj, bfk, bfl, eps, epe)

      implicit none

*     Function to return the length of the overlap between bf_tables
*     entries bfi and bfj, bfk (if not zero), and bfl (again is not
*     zero)
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED  VARIABLES
      integer*4 bfi, bfj, bfk, bfl  ! Entries in Bf_table to be checked
                       ! bjk and bfl may be zero in which case they are
                       ! not checked
      integer*4 eps, epe  ! Start and stop times of overlap
                          ! between biases (if not overlap then
                          ! eps(tart) will be greater than epe(nd)

****  OK; get the start and ends of the overlap
      eps = max(bf_table(1,bfi), bf_table(1,bfj))
      epe = min(bf_table(1,bfi)+bf_table(4,bfi),
     .          bf_table(1,bfj)+bf_table(4,bfj))
*     OK Now see if there are other entries to check
      if( bfk.gt.0 ) then
          eps = max(eps,bf_table(1,bfk))
          epe = min(epe, bf_table(1,bfk)+bf_table(4,bfk))
*         See if next one should be checked also
          if( bfl.gt.0 ) then
             eps = max(eps,bf_table(1,bfl))
             epe = min(epe, bf_table(1,bfl)+bf_table(4,bfl))
          endif
      endif

*     Finally compute the overlap range  
      chk_overlap = epe - eps + 1

****  Thats all
      return
      end

CTILE GET_OBS
      subroutine get_obs(ep, ns,lv,  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs, wgh, OK)

      implicit none

*     Routine to return the MW-WL resiudal, L2-L1 residual
*     and the L1 and L2 residuals along with a weight and
*     status as to OK or not
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 ep    ! Epoch being extracted
     .,         ns    ! Site number
     .,         lv    ! Satellite list number being extracted


*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

*   obs(5)  - MW-WL, L2-L1, L1, L2, LC residuals in cycles
*   wgh(5)   - Weight of measurement (inversely proporational to
*             1/sin(el)**2 and phase noise on satellite
      real*8 obs(5), wgh(5)

*   OK      - Logical: Set OK if data is OK
      logical OK


* LOCAL VARIABLES
      integer*4 ch    ! Channel number at this particular epoch
     .,      ltoc     ! Function to return channel number of specific
                      ! satellite list number
      real*8 ow(4)     ! Four one way residuals, L1, L2, P1, P2
     .,      rmsr, rmsp    ! RMS values to use for MW-WL and Phase values
     .,      sine      ! Sine Elevation angle

       logical data_OK ! Function which returns true to data passes
                        ! the phs_mask on the data flag

**** OK: Set status bad until measurement found
      OK = .false.
      ch =  ltoc( ctol_cse(1,ns,ep),lv, actual_max_chan)
      if ( ch .gt. 0 ) then
* MOD TAH 050712: Skip codeless data completely
* 
*        Get each the measaurements
         ow(1) = L1r_phs_cse(ch,ns,ep) + L1_cyc_cse(ch,ns,ep)
         ow(2) = L2r_phs_cse(ch,ns,ep) + L2_cyc_cse(ch,ns,ep)
         ow(3) = L1r_rng_cse(ch,ns,ep)
         ow(4) = L2r_rng_cse(ch,ns,ep)

*        Now form the combinations we need: MW-WL, L2-L1,
*        L1 and L2. (MW-WL formed as L2-L1 Phase)
         if( lambda(ch,4,ns).eq.0 ) then
            obs(1) = 0.d0   ! Of no P2 set set, will be detected when
                            ! statistics are finished
         else
            obs(1) = ow(2)-ow(1)+dfsf(lv)*(ow(3)+ow(4))
         endif
         obs(2) = ow(2)-ow(1)
         obs(3) = ow(1)
         obs(4) = ow(2)
         obs(5) = lcf1(lv)*ow(1) + lcf2(lv)*ow(2)
*
*        Get the weight: 1/(RMS**2*sin(el)**2)
         sine = sin(azel_cse(2,ch,ns,ep))
* MOD TAH 080515: Base RMS on WL RMS not phase RMS
C         rms = LC_svs_RMS(lv,ns)
         rmsr = WL_RMS(ns)
         rmsp = max(LC_svs_RMS(lv,ns),0.025d0) 
* MOD TAH 080509: Weight should decrease with smaller elevation
*        angles not increase as with old code
C        wgh = 1.d0/(rms*sine)**2
         wgh(1) = (sine/rmsr)**2  ! MW Range noise
         wgh(2) = (sine/rmsp)**2/2 ! EX
         wgh(3) = (sine/rmsp)**2 ! L1
         wgh(4) = (sine/rmsp)**2 ! L2
         wgh(5) = (sine/rmsp*(lcf1(lv)+lcf2(lv)))**2 ! LC


*        Check the data flags
         OK = data_OK(data_flag_cse(ch,ns,ep),0,phs_mask)
      endif
***   Thats all
      return
      end

CTITLE INC_ZERO_STATS

      subroutine inc_zero_stats(obs,wgh,stats)

      implicit none

*     Routine to from double differences and increment the
*     statistics for the MW-WL, L2-L1, L1 and L2 only
*

* PASSED VARIABLES
      real*8 obs(5,4)  ! Obs from each oneway (columns are site
                         ! satellite combinations, rows are data
                         ! types
     .,       wgh(5,4)    ! Weight for each station/satellite
     .,       stats(4,5) ! Statistics: Columns for each data type
                          ! rows for res*wgh res**2*wgh wgh and num

* LOCAL VARIABLES
      real*8 dd(5)  ! Double differences for the 4 data types
     .,      tot_wgh(5)   ! Total weight for dd (sum of weights)

      integer*4 i   ! Loop counter

      logical noP2   ! Set true if obs(1,x) is exactly zero

* MOD TAH 050713: See if P2 was there
      noP2 = .false.
      do i = 1,4
          if( obs(1,i).eq.0.d0 ) noP2 = .true.
      enddo

****  Form the double differences from each measurement
*     type and the total wght
 
      do i = 1, 5
          dd(i) = (obs(i,1)-obs(i,2))-(obs(i,3)-obs(i,4))
* MOD TAH 080517: Sum inverse weights so low weight will make total low
          tot_wgh(i) = 1.d0/(1.d0/wgh(i,1)+1.d0/wgh(i,2)+
     .                       1.d0/wgh(i,3)+1.d0/wgh(i,4))

      end do

****  Now increment the statistic for each data type
      do i = 1,5
*         If there is no P2 set number do not imcrement the 
*         MW-WL combination (i=1)
          if( noP2 .and. i.eq. 1 ) then
*            Do not doing anything.
          else
             stats(1,i) = stats(1,i) + dd(i)*tot_wgh(i)
             stats(2,i) = stats(2,i) + dd(i)**2*tot_wgh(i)
             stats(3,i) = stats(3,i) + tot_wgh(i)
             stats(4,i) = stats(4,i) + 1
          endif
      end do
         

****  Thats all
       return
       end

CTITLE FIN_ZERO_STATS

      subroutine fin_zero_stats(bft, bfv, stats, dL1_cyc, dL2_cyc, iter)

      implicit none

*     Routine to compute the mean values of DD combinations
*     and to assess the reliability of fixing them to integer
*     values.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 bft,bfv  ! BF_table entry being assessed and the 
                         ! the entry for the second site/satellite
     .,         iter     ! Iteration counter which set criteria

      real*8 stats(4,5)  ! Summation statistics for the 4 DD
                        ! observable types
     .,      dL1_cyc, dL2_cyc  ! Integer values of the offsets
                        ! in L1 and L2 derived from the 4 DD

      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.
                        ! observables
      integer*4 lv      ! Satellite number (bf_table(3,x))

* LOCAL VARIABLES
      integer*4 i   ! Loop counter
      real*8 mean(5), rms(5), sigm(5)  ! Mean, WRMS and sigma of mean
                     ! for the five observables
      real*8 L21_wl, L21    ! Difference between L2 and L1 number of cycles from
                      ! MW widelane and L2-L1
     .,      mean_lc  ! Mean value of LC after correcting for L21 bias
     .,      L1     ! Nearest L1 cycles consistent with mean_lc 

      character*2 status ! String set to XX or RR depending on which
                     ! which biases are fixed.

***** Finish up the statistics.   Loop over each of the
*     observables.  Check that number of measurements is
*     is not zero.
      nl_conf(bft) = 0.d0
      wl_conf(bft) = 0.d0
      dL1_cyc = 0.d0
      dL2_cyc = 0.d0
      L1 = 0.d0
* MOD TAH 050713: Check L1-L2 count since MW-WL may be zero if
*     no P2 data
      if( stats(4,2).eq.0 ) then
         if( 1.eq.2 )
     .   write(*,120) bft, cf_codes(bf_table(2,bft)),
     .        prn_list(bf_table(3,bft)), cf_codes(bf_table(2,bfv)), 
     .      prn_list(bf_table(3,bfv))
 120     format('***WARNING*** No data found for bias ',i4,
     .           ' Zero DD: Site ',a4,' PRN ',i2.2,'(with ',a4,
     .           ' PRN ',i2.2,')')
         RETURN
      endif

*****  For each type get the mean, RMS, and sigma of mean
      do i = 1,5
          if( stats(4,i).gt.0 ) then
             mean(i) = stats(1,i)/stats(3,i)
             rms(i) = sqrt(abs(stats(2,i)/stats(3,i)-mean(i)**2))
             sigm(i) = sqrt(rms(i)**2/stats(4,i))
          else
             mean(i) = 0
             rms(i) = 100.d0
             sigm(i) = 100.d0
          end if
      end do

*****  Now see if we can fix these to integers
* MOD TAH 180320: Remap the frequencies for Glonass.
      lv = bf_table(3,bft)
      fr1 = fL1u(lv)/fL1(lv)
      fr2 = fL2u(lv)/fL2(lv)
      L21_wl = nint(mean(1)*fr1)/fr1
      L21 = L21_wl               
* RWK 150203: This may need to be reformulated to account for differences
*             in frequency for Glonass SVs
      mean_lc = mean(5) - L21_wl*lcf2(1)
      L1 = nint(mean_lc/(lcf1(1)+lcf2(1))*fr1)/fr1
      if( abs(mean(1)-L21_wl).lt.mdev_wl_tol/iter .and. 
     .    sigm(1).lt.msig_wl_tol/iter .and. 
     .    stats(4,1).gt.min_wl_tol*iter ) then
          
          wl_conf(bft) = (1-abs(mean(1)-L21_wl))/sigm(1)
          if( wl_conf(bft).lt.0 ) then
              write(*,998) bft,wl_conf(bft), mean(1), L21_wl, 
     .              sigm(1)
 998          format('ZBFNegWL ',i5,f8.3, 3F10.6)
          end if
****      See if can resolve the narrow lane            
* RWK 150203: This may need to be reformulated to account for differences
*             in frequencies for Glonass SVs
          if( abs(mean_lc-L1*(lcf1(1)+lcf2(1))).lt.mdev_wl_tol .and.
     .        sigm(5).lt.msig_wl_tol/2 ) then
              nl_conf(bft) = (1-abs(mean_lc-L1*(lcf1(1)+lcf2(1))))/
     .             sigm(5)
          endif
      else
*        See if L2-L1 will do it
         L21 = nint(mean(2)*fr1)/fr1
         if( abs(mean(2)-L21).lt.mdev_wl_tol/2/iter .and. 
     .       sigm(2).lt.msig_wl_tol/10/iter  .and. 
     .       stats(4,2).gt.min_wl_tol*iter ) then
* MOD TAH 040430: Changed mean(1) and sigm(1) to the LG values of
*            mean(2) and sigm(2)
             wl_conf(bft) = (1-abs(mean(2)-L21))/sigm(2)
             if( wl_conf(bft).lt.0 ) then
                 write(*,999) bft,wl_conf(bft), mean(2), L21, sigm(2)
 999             format('ZBFNegIO ',i5,f8.3, 3F10.6)
             end if
****         See if can resolve the narrow lane
             mean_lc = mean(5) - L21*lcf2(1)
             L1 = nint(mean_lc/(lcf1(1)+lcf2(1))*fr1)/fr1
             if( abs(mean_lc-L1*(lcf1(1)+lcf2(1))).lt.mdev_wl_tol .and.
     .           sigm(5).lt.msig_wl_tol/2 ) then
                 nl_conf(bft) = (1-abs(mean_lc-L1*(lcf1(1)+lcf2(1))))/
     .                     sigm(5)
             endif
         else
c            L1  = 0.d0
c            L21 = 0.d0
            L21 = L21_wl
            wl_conf(bft) = 0.0d0
            nl_conf(bft) = 0.0d0
         endif
      endif

****  Save the offsets
      dl1_cyc = L1
      dl2_cyc = L1+L21
      if( wl_conf(bft).lt.0 ) then
          write(*,997) bft,wl_conf(bft), mean(1), L21_wl, 
     .           sigm(1), mean(2), L21_wl, sigm(2)
 997      format('ZBFNegAL ',i5,f8.3, 6F10.6)
      endif
 
****  Write out results
      status = 'RR'
      if( wl_conf(bft).gt.dchi_wl_tol*iter ) status(1:1) = 'X'
      if( nl_conf(bft).gt.dchi_wl_tol*iter ) status(2:2) = 'X'

      if( status(1:1).eq.'X' )
     .write(*,220) iter, bft, status, wl_conf(bft), 
     .      cf_codes(bf_table(2,bft)), 
     .      prn_list(bf_table(3,bft)),
     .      cf_codes(bf_table(2,bfv)), 
     .      prn_list(bf_table(3,bfv)), nint(stats(4,1)),
     .     (mean(i),sigm(i),rms(i),i=1,2),
     .     (mean(i),sigm(i),rms(i),i=5,5), dL1_cyc, dL2_cyc

 220  format('ZBF ',i2,1x, i4,1x,a2,1x,F8.1,1x,a4,1x,'PRN ',i2.2,1x,
     .       '(',a4,1x,I2.2,') ',i4,
     .     3(F7.2,1x,f6.2,1x,f6.2),' dL12 ',2F6.1)

****  Thats all
      return
      end

CTITLE GET_LAST_BIAS

      subroutine get_last_bias( in, last)

      implicit none

*     Routine to look at the site and satellite in bf-table entry
*     in and return the last bias table entry associated with this
*     site and satellite
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
      integer*4 in  ! Input bias table entry
     .,         last ! Bf_table entry for the last entry
                     ! associated with this satellite 

* LOCAL VARIABLES
      integer*4 ns, lv     ! Site and satellite being processed


****  Get the site and satellite for this bf_table entry
      ns = bf_table(2,in)
      lv = bf_table(3,in)

****  Set this to be last and then scan to see change
      last = in
      if( last.lt.tot_bf ) then 
         do while ( bf_table(2,last+1).eq.ns .and. 
     .              bf_table(3,last+1).eq.lv .and.
     .              last.lt.tot_bf ) 
            last = last + 1
         end do
      end if

****  Thats all
      return
      end


CTITLE SCAN_WL

      subroutine scan_wl( L1r_phs_cse, L2r_phs_cse,
     .             L1_cyc_cse, L2_cyc_cse, 
     .             L1r_rng_cse, L2r_rng_cse, 
     .             ctol_cse, data_flag_cse, azel_cse)

      implicit none

*     Routine to scan the widelines getting RMS scatter and
*     RMS of differences.  The latter is used to judge the
*     systematic level of the widelane

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

* LOCAL VARIABLES

      integer*4 i,j  ! Loop over bf_table entry
     .,     ns, lv  ! Site and satellite being processed
     .,     ep      ! Epoch loop variables
     .,     nwl, ndwl  ! Number of widelanes and differences
     .,     ltoc    ! Function to go from satellite to channel
     .,     ch      ! Channel number

      real*8 stats(3), dstat(3)  ! Accumulation arrays for status
                    ! and difference status
     .,      ow(4)  ! Four one-way masurements
     .,      mwwl, dmwwl, prev_mwwl ! MW-WL, difference and previous value
     .,      mean, rms, drms   ! Mean, rms and difference rms
     .,      sine, wgh  ! Sin(elevation) and weight

      logical OK, prev_OK     ! Set true for good measurement and
                     ! previous being good
     .,      data_OK ! Function returns true for good data

      character*4 status  ! Set to SYS or RAN depending on nature
                     ! of widelane (ratio to drms to rms).



***** Loop over the bf_table
      write(*,110)
 110  format('BF_SYS  Num Site  PRN      #        RMS      dRMS Status')
      do i = 1, tot_bf
         nwl = 0
         ndwl = 0
         do j = 1,3
            stats(j) = 0.d0
            dstat(j) = 0.d0
         end do
         ns = bf_table(2,i)
         lv = bf_table(3,i)
         prev_OK = .false.
         do ep = bf_table(1,i),bf_table(1,i)+bf_table(4,i)

*           Get the one-way measurements
            ch =  ltoc( ctol_cse(1,ns,ep),lv, actual_max_chan)
            if( ch.gt.0 ) then
              ow(1) = L1r_phs_cse(ch,ns,ep) + L1_cyc_cse(ch,ns,ep)
              ow(2) = L2r_phs_cse(ch,ns,ep) + L2_cyc_cse(ch,ns,ep)
              ow(3) = L1r_rng_cse(ch,ns,ep)
              ow(4) = L2r_rng_cse(ch,ns,ep)
              sine = sin(azel_cse(2,ch,ns,ep))
* MOD TAH 080509: Weight should decrease with lower elevation angle not
*             increase (previous version)
              wgh = sine**2
              OK = data_OK(data_flag_cse(ch,ns,ep),0,phs_mask)
              if( OK ) then
                  mwwl = ow(2)-ow(1)+dfsf(lv)*(ow(3)+ow(4))
                  stats(1) = stats(1) + mwwl*wgh
                  stats(2) = stats(2) + mwwl**2*wgh
                  stats(3) = stats(3) + wgh
                  nwl = nwl + 1
                  if( prev_OK ) then
                      dmwwl = mwwl-prev_mwwl
                      prev_mwwl = mwwl
                      dstat(1) = dstat(1) + dmwwl*wgh
                      dstat(2) = dstat(2) + dmwwl**2*wgh
                      dstat(3) = dstat(3) + wgh
                      ndwl = ndwl + 1
                  endif
                  prev_OK = .true.
              endif
            end if
         end do
*        Finish up statistics
         if( nwl.gt.1 ) then
             mean = stats(1)/stats(3)
             rms = sqrt(stats(2)/stats(3)-mean**2)
             drms = sqrt(dstat(2)/dstat(3))
*            Report and flag if systematic (drms should sqrt(2) larger than
*            rms).  If it is too small we flag
             if( rms.gt.1.0d0 .and. drms/rms.lt.0.40d0 ) then
                 call sbit(bf_table(5,i),5,1)
                 status = 'SYS'
             else
                 status = 'RAN'
             endif
             write(*,120) i, cf_codes(ns), prn_list(lv), nwl, rms, 
     .                    drms, status
 120         format('BF_SYS ',i4,1x,a4,' PRN ',I2.2,1x,i5,
     .              2f10.2,1x,a)
         else
             write(*,140) i, cf_codes(ns), prn_list(lv)
 140         format('**WARNING** No data BF ',i4,1x,a4,' PRN ',i2.2)
         endif
      end do

****  Thats all
      return
      end

CTITLE SORT_BLEN

      subroutine sort_blen(rs, bf_sort)

      implicit none

*     Routine to sort the bf_table by increasing length from the ref_site
*
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 rs  ! Reference site
     .,         bf_sort(max_bf_table)  ! Sorted reference list to bf_table

* LOCAL VARIABLES
      integer*4 i,j,k  ! Loop counters
     .,     st,en      ! Range to copy from bf_table
     .,     nx,ent     ! Site with next length, and baseline length entry
     .,     bl_ent     ! Function to return baseline entry

      logical site_used(max_cfiles)  ! Set to true when a site included in 
                       ! in sorted list

      real*8 min_len, max_len  ! Min and max_length found

****  Loop over the sites
      do i = 1, num_cfiles
         site_used(i) = .false.
      end do
      site_used(rs) = .true.

      min_len = 0.d0
*     Make the first group by the reference site
      st = bf_index(rs)
      if( rs.lt.num_cfiles ) then
         en = bf_index(rs+1)-1
      else
         en = tot_bf
      endif
      do i = st,en
         k = i - st + 1
         bf_sort(k) = i
      end do

*     Now do the rest
      do i = 2, num_cfiles
         max_len = 1.d12
*        Find the next longest line
         do j = 1, num_cfiles
             if( .not.site_used(j) ) then
                ent = bl_ent(j,rs)
                if( baselens(ent).ge.min_len .and.
     .              baselens(ent).lt.max_len ) then
                    max_len = baselens(ent)
                    nx = j
                end if
             end if 
         enddo
****     OK, save the next longest line
         min_len = max_len
         st = bf_index(nx)
         site_used(nx) = .true.
         if( nx.lt.num_cfiles ) then
             en = bf_index(nx+1) - 1
         else
             en = tot_bf
         endif
         do j = st,en
            k = k + 1
            bf_sort(k) = j
         end do
      end do

      write(*,150) tot_bf, (i, bf_sort(i),i=1,tot_bf)
 150  format('BF_SORT Table ',i6,' Table ',/,
     .        1000(10('I ',2I5,1x),/))
      
               
****  Thats all
      return
      end

CTITLE VERFIY_NEWREF

      subroutine verify_newref(rs, new)

      implicit none

*     Routine to verify that the new reference has no overlap
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 rs   ! New reference bias flags

* LOCAL VARIABLES

      integer*4 ns   ! Site number
     .,         start, last  ! Start and end of bias flags for this
                     ! site
     .,         i,j,k  ! Loop variables over other bias flags
     .,         overlap, chk_overlap  ! Overlap and function to check
                       ! overlap
     .,         eps, epe  ! Start and end epochs of overlap
     .,         ks, kv    ! Overlapping site and satellite

      logical kbit     ! check bit
     .,       new      ! Set true unless overlap is found than then set 
                       ! false

***** Get the range of biases for this site
      new = .true.
      ns = bf_table(2,rs) 
      start = bf_index(ns)
      if( ns.lt.num_cfiles) then
          last = bf_index(ns+1)-1
      else
          last = tot_bf
      endif

***** Loop over all bias flags at this site
      do i = start, last
         if( i.ne.rs ) then
*           Get overlap with this satellite and rs satellite
            overlap = chk_overlap(i, rs, 0, 0,  eps, epe)
            if( overlap.gt.0) then

*              Loop for another station with this pair of satellites
*              overlapping
               do j = 1,tot_bf
                  if( bf_table(2,j).ne.ns .and. 
     .               bf_table(3,j).eq.bf_table(3,rs) ) then
*                    See if overlap
                     overlap = chk_overlap(rs,j,i,0,eps,epe)
                     if( overlap.gt.0 ) then
*                       OK see if we can match the 2nd site and 2md 
*                      satellite
                        do k = 1,tot_bf
                           if( bf_table(2,k).eq.bf_table(2,j) .and.
     .                         bf_table(3,k).eq.bf_table(3,i) ) then
*                             See if overlap
                              overlap = chk_overlap(rs,j,i,k,eps,epe)
*                             If there is overlap, see if any of these 
*                             other bias flags is marked as resolvable
                              if( overlap.gt.0 ) then
                                 if( kbit(bf_table(5,i),3) .and.
     .                               kbit(bf_table(5,j),3) .and.
     .                               kbit(bf_table(5,k),3) ) then
                                     ks = bf_table(2,k)
                                     kv = bf_table(3,k)
                                     new = .false.
                                     write(*,120) rs, cf_codes(ns), 
     .                                    prn_list(bf_table(3,rs)),
     .                                    cf_codes(ks),prn_list(kv),
     .                                    eps, epe
 120                                 format('NewRef ',i4,1x,a4,1x,
     .                                   'PRN ',i2.2,' Overlaps ',
     .                                    a4,1x,'PRN ',i2.2,1x, 
     .                                   'EPS',i5,1x,I5)
*                                 Code below is check on non-resolved
*                                 entries.  This list can be long.
C                                 else
C                                    ks = bf_table(2,k)
C                                    kv = bf_table(3,k)
C                                    write(*,140) rs, cf_codes(ns), 
C    .                                    prn_list(bf_table(3,rs)),
C    .                                    cf_codes(ks),prn_list(kv),
C    .                                    eps, epe
C140                                 format('NewRef ',i4,1x,a4,1x,
C    .                                   'PRN ',i2.2,' Non-resolved ',
C    .                                    a4,1x,'PRN ',i2.2,1x, 
C    .                                   'EPOCH ',i5,'-',I5)

                                 end if
                              end if
                           end if
                        end do
                     end if
                  end if
               end do
            end if
         end if
      end do

***** Thats all
      return
      end

CTITLE SCAN_ALLDD

      subroutine scan_alldd(rs, L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse, 
     .    bf_neq, bf_bvec, tot_dd)

      implicit none

*     Routine to scane all double differences with bias rs.  (For a
*     only those greater than the one being considered are scanned.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 rs   ! New reference bias flags
     .,         tot_dd   ! Total number of DDs
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
       real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

      real*8 bf_neq(tot_bf, tot_bf ) ! Normal equations
     .,      bf_bvec(tot_bf,2)  ! Bvector for ML-WL and NL

*
* LOCAL VARIABLES

      integer*4 ns   ! Site number
     .,         start, last  ! Start and end of bias flags for this
                     ! site
     .,         i,j,k  ! Loop variables over other bias flags
     .,         overlap, chk_overlap  ! Overlap and function to check
                       ! overlap
     .,         eps, epe  ! Start and end epochs of overlap


***** Get the range of biases for this site
      ns = bf_table(2,rs) 
      start = rs + 1
      if( ns.lt.num_cfiles) then
          last = bf_index(ns+1)-1
      else
          last = tot_bf
      endif

***** Loop over all bias flags at this site
      do i = start, last
*       Get overlap with this satellite and rs satellite
        overlap = chk_overlap(i, rs, 0, 0,  eps, epe)
        if( overlap.gt.0) then

*          Loop for another station with this pair of satellites
*          overlapping
           do j = 1,tot_bf
              if( bf_table(2,j).ne.ns .and. 
     .            bf_table(3,j).eq.bf_table(3,rs) ) then
*                See if overlap
                 overlap = chk_overlap(rs,j,i,0,eps,epe)
                 if( overlap.gt.0 ) then
*                   OK see if we can match the 2nd site and 2md 
*                   satellite
                    do k = 1,tot_bf
                       if( bf_table(2,k).eq.bf_table(2,j) .and.
     .                     bf_table(3,k).eq.bf_table(3,i) ) then
*                         See if overlap
                          overlap = chk_overlap(rs,j,i,k,eps,epe)
*                         If there is overlap, see if any of these 
*                         other bias flags is marked as resolvable
                          if( overlap.gt.0 ) then
*                             OK: Overlap.  Now get all the doubles
*                             difference statistcs
                              tot_dd = tot_dd + 1
                              call acc_alldd(rs, i, j, k, eps, epe, 
     .                             L1r_phs_cse, L2r_phs_cse,
     .                             L1_cyc_cse, L2_cyc_cse, 
     .                             L1r_rng_cse, L2r_rng_cse, 
     .                             ctol_cse, data_flag_cse, azel_cse,
     .                             bf_neq, bf_bvec, tot_dd)
                          end if
                       end if
                    end do
                 end if
              end if
           end do
        end if
      end do

***** Thats all
      return
      end


      subroutine acc_alldd(bft, bfs, bfj, bfv, eps, epe, 
     .    L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse, bf_neq, bf_bvec, tot_dd)

      implicit none

*     Routine to scane all double differences with bias rs.  (For a
*     only those greater than the one being considered are scanned.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      integer*4 tot_dd   ! Total number of Double Diffs
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

      real*8 bf_neq(tot_bf, tot_bf ) ! Normal equations
     .,      bf_bvec(tot_bf,2)  ! Bvector for ML-WL and NL
  
      integer*4 bft, bfs, bfj, bfv  ! Bias table entries
     .,         eps, epe   ! Overlap of epochs in double difference
      

* LOCAL VARIABLES

      integer*4 i,j
     .,    ep           ! Epoch counter
     .,    bf(4)        ! Four bias table entries

      real*8 obs(5,4)   ! Obs from each oneway (columns are site
                        ! satellite combinations, rows are data
                        ! types
     .,      wgh(5,4)   ! Weight for each station/satellite
     .,      mwgh       ! Weight used for incrementing normal equations
     .,      stats(4,5) ! Statistics: Columns for each data type
                        ! rows for res*wgh res**2*wgh wgh and num
     .,      ap(4)      ! Partials (1 -1 -1 1) 

      real*8 mean(5), rms(5), sigm(5)  ! Mean, WRMS and sigma of mean
                     ! for the five observables

      logical OK(4)    ! set true for each good station/satellite one-way
     .,       kbit     ! Function to check bit status

      character*4 status  ! Four characters as X or R depending on fixed
                          ! or not

      data ap / 1.d0, -1.d0, -1.d0, 1.d0 /
                   
*     Clear the accumulation arrays
      do i = 1,4  ! Loop over res*wgh, res**2*wgh wgh num 
         do j = 1,5  ! Loop over MW-WL, L2-L1, L1 and L2, LC
             stats(i,j) = 0.d0
         end do
      end do

*     Loop over the epoch range
      do ep = eps,epe
*       Get the residuals for each of one-ways that go into
*       the double difference
         call get_obs(ep, bf_table(2,bft),bf_table(3,bft),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,1), wgh(1,1), OK(1))
         call get_obs(ep, bf_table(2,bfs),bf_table(3,bfs),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,2), wgh(1,2), OK(2))
         call get_obs(ep, bf_table(2,bfj),bf_table(3,bfj),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,3), wgh(1,3), OK(3))
         call get_obs(ep, bf_table(2,bfv),bf_table(3,bfv),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,4), wgh(1,4), OK(4))

****     If all the oneways are good, increment the statistics
         if( OK(1) .and. OK(2) .and. OK(3) .and. OK(4) ) then
              call inc_zero_stats(obs, wgh, stats)
         endif
      end do

****  Finish the mean and RMS caculations
      if( stats(4,2).eq.0 ) then
         if( 1.eq.2 )
     .   write(*,120) bft, cf_codes(bf_table(2,bft)),
     .        prn_list(bf_table(3,bft)), cf_codes(bf_table(2,bfv)), 
     .        prn_list(bf_table(3,bfv)), eps, epe
 120     format(' No data found for bias ',i4,
     .           ' Zero DD: Site ',a4,' PRN ',i2.2,'(with ',a4,
     .           ' PRN ',i2.2,') Epochs ',i5,' to ',i5)
         RETURN
      endif

***** For each type get the mean, RMS, and sigma of mean
      do i = 1,5
          if( stats(4,i).gt.0 ) then
             mean(i) = stats(1,i)/stats(3,i)
             rms(i) = sqrt(abs(stats(2,i)/stats(3,i)-mean(i)**2))
             sigm(i) = sqrt(rms(i)**2/stats(4,i))
          else
             mean(i) = 0.d0
             rms(i) = 100.d0
             sigm(i) = 100.d0
          end if

      end do
***** Check the status
      bf(1) = bft
      bf(2) = bfs
      bf(3) = bfj
      bf(4) = bfv
      do i = 1,4
        if( kbit(bf_table(5,bf(i)),1) ) then
            status(i:i) = 'X'
        else
            status(i:i) = 'R'
        end if
      end do

***** Accumulate the normal equations and bvector for MW-WL and NL.  Base
*     weight on sigma of mean unless there are only a couple of observations
*     Weight based on MW-WL (same weight used for NL so that Normal equations
*     are the same).
      if( stats(4,1).gt.2 ) then
          mwgh = 1/(sigm(1)**2)
      else
          mwgh = 1.d0
      end if

      do i = 1,4 
         do j = 1, 4
*           Increment into normal eqautions
            bf_neq(bf(i),bf(j)) = bf_neq(bf(i),bf(j)) + mwgh*ap(i)*ap(j)
         enddo
      end do
*     Increment the bvector
      do i = 1, 4
         bf_bvec(bf(i),1) = bf_bvec(bf(i),1) + mwgh*mean(1)*ap(i)
         bf_bvec(bf(i),2) = bf_bvec(bf(i),2) + mwgh*mean(5)*ap(i)
      end do


****  Output the values
      write(*,220) tot_dd, bf, status, wl_conf(bft), 
     .      cf_codes(bf_table(2,bft)), 
     .      prn_list(bf_table(3,bft)),
     .      cf_codes(bf_table(2,bfv)), 
     .      prn_list(bf_table(3,bfv)),
     .      nint(stats(4,1)), eps,epe,
     .     (mean(i),sigm(i),rms(i),i=1,2),
     .     (mean(i),sigm(i),rms(i),i=5,5)

 220  format('ADD ',i6,1x,4i5,1x,a4,1x,F8.1,1x,a4,1x,'PRN ',i2.2,1x,
     .       1x,a4,1x,'PRN ',I2.2,1x,i4,' Eps ',2i5,
     .     3(F7.2,1x,f6.2,1x,f6.2))

****  Thats all 
      return
      end

CTITLE CMA_ALLOC

      subroutine cma_alloc(cma_data, ineq, ibvec, tot_bf )

      implicit none

*     Routine to allocate memory for bias normal equations

* PASSED VARIABLES
      integer*4 cma_data, tot_bf
      integer*8 ineq, ibvec

* LOCAL VARIABLES
      integer*4 tot_memory  ! Total number of i*4 words needed
      integer*8 offset      ! offset between available start and cma_data
     .,         memassign   ! Integer*8 function to assign memory

      character*128 status

****  Coompute number of i*4 needed
      tot_memory = (tot_bf**2 + 2*tot_bf)*2
      write(status, 110) tot_memory/(1024.d0**2)*4
 110  format('Allocating ',F8.2,' Mbytes for Bias NEQ')
      call report_stat('status','autcln','CMA_ALLOC',' ',
     .                   status,0)

****  Try to allocate
      offset = memassign(tot_memory,1,loc(cma_data))
      if ( offset.eq.0 ) then
          write(*,120) tot_memory/(1024.d0**2)*4
 120      format(' *** DISASTER *** Not enough memory to solve',
     .           ' biases. ',F8.2,' Mbyes needed',/,
     .            ' Either run on a larger machine or reduce',
     .            ' number of cfiles.')
          call report_stat('fatal','autcln','cma_alloc',' ',
     .         'Not enough memory to run program',0)
      end if

****  Get location where we have memory
      ineq = offset
      ibvec = ineq + tot_bf**2*2

****  Thats all
      return
      end

CTTILE CLEAR_BFNEQ

      subroutine clear_bfneq(bf_neq, bf_bvec, tot_bf)

      implicit none

*     Routine to clear the normal equations

* PASSED VARIABLES
      integer*4 tot_bf   ! Number of bias flags
      real*8 bf_neq(tot_bf, tot_bf) ! Normal equations
     .,      bf_bvec(tot_bf,2)  ! Bvector for ML-WL and NL


* LOCAL VARIABLES
      integer*4 i,j      ! Loop counters

****  CLear the memory space
      do i = 1, tot_bf
         do j = 1,2
            bf_bvec(i,j) = 0.d0
         end do
         do j = 1, tot_bf
            bf_neq(j,i) = 0.d0
         end do
      end do

****  Thats all
      return
      end

CTTLE SOLVE_BFNEQ
 
      subroutine solve_bfneq(bf_neq, bf_bvec, 
     .       L1_cyc_cse, L2_cyc_cse, ctol_cse)


      implicit none

*     Routine to solve the normal equations for the bias parameters

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'


* PASSED VARIABLES
      real*8 bf_neq(tot_bf, tot_bf ) ! Normal equations
     .,      bf_bvec(tot_bf,2)  ! Bvector for ML-WL and NL

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
      integer*4 ctol_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)

* LOCAL VARIABLES
      real*8 scale(max_bf_table) ! Scale vector for invert_vis
     .,      avat, dchi
     .,      min_sig
     .,      dL21, dL1, dL2, dLC, dLCm, dLC_obs
     .,      dL21t, dL1t, dL2t, dLCt, dLCs
     .,      dsol, dcov

      integer*4 pivot(max_bf_table) ! Pivot elements
     .,         min_i
     .,         ep, ns, lv, ch
     .,         ltoc        ! 


      real*8 bf_neq_col(max_bf_table), equ_gn(max_bf_table)

      integer*4 i,j
      logical kbit   ! Function to check bit status
     .,       fixing

****  Constrain the bias parameters that have been fixed
      do i = 1, tot_bf
         if( kbit(bf_table(5,i),1) .and. wl_conf(i).eq.0.d0 ) then
            write(*,120) i, bf_neq(i,i), bf_bvec(i,1),bf_bvec(i,2) 
 120        format('Force: BF_TABLE ',i4,' NEQ ',3e12.3)
            bf_neq(i,i) = bf_neq(i,i)+1.d8
            call sbit(bf_table(5,i),6,1)
         else
            bf_neq(i,i) = bf_neq(i,i)+1.d0
         endif
      end do

****  Now solve the system
      do i = 1, tot_bf
         bf_neq_col(i) = bf_bvec(i,1)
      end do

      call invert_vis(bf_neq, bf_bvec, scale, pivot, tot_bf, tot_bf,2)

*     Get chi**2 change
      call dwdot(dchi, bf_neq_col,1, bf_bvec(1,1),1, tot_bf)
      write(*,210) dchi
 210  format('Change in chi**2 ',E16.6)

*     Write out the solution
      do i = 1, tot_bf
         write(*,220) i, cf_codes(bf_table(2,i)), 
     .                prn_list(bf_table(3,i)),bf_table(1,i),
     .                bf_table(1,i)+bf_table(4,i), bf_table(5,i),
     .                bf_bvec(i,1), sqrt(bf_neq(i,i)), bf_bvec(i,2)
 220     format('BF_EST ',I5,1x,a4,' PRN ',I2.2,' Eps ',2I6,1x,o3,
     .          4F12.4)
      end do

****  Next fix to integers.
      fixing = .true.
      do while ( fixing )
         min_sig = 1.d10
         min_i = 0
         fixing = .false.
         do i = 1, tot_bf
            dL21 = nint(bf_bvec(i,1))
            if( .not.kbit(bf_table(5,i),6) .and.
     .          (bf_bvec(i,1)-dL21)**2/bf_neq(i,i)+
     .           bf_neq(i,i)*100.lt.min_sig  ) then
                min_sig = (bf_bvec(i,1)-dL21)**2/bf_neq(i,i)+
     .                     bf_neq(i,i)*100
                min_i = i
            end if
         end do
*        See if any left
         if( min_i.gt.0 ) then
             call sbit(bf_table(5,min_i),6,1)
*            Force this value to an integer
             dLC_obs = bf_bvec(min_i,2)
             dL21 = nint(bf_bvec(min_i,1)) 
* RWK 150203: This may need to be reformulated to account for differences
*             in frequencies for the Glonass SVS
             dLCm = bf_bvec(min_i,2)-dL21*lcf2(1)
             dL1 = nint(dLCm/(lcf1(1)+lcf2(1)))
             dL2 = dL1 + dL21
             dLC = dL1*lcf1(1)+dL2*lcf2(1)

*            Now check out the values to see if L2-L1 could be off
*            by 1 cycle
             if( abs(dLC-dLC_obs).gt.0.1d0 ) then
                 dL21t = dL21+sign(1.d0,bf_bvec(min_i,1)-dL21)
                 dLCt = bf_bvec(min_i,2)-dL21t*lcf2(1)
                 dL1t = nint(dLCt/(lcf1(1)+lcf2(1)))
                 dL2t = dL1t + dL21t
                 dLCs = dL1t*lcf1(1)+dL2t*lcf2(1)
                 if( abs(dLCs-dLC_obs).lt.abs(dLC-dLC_obs) ) then
                    write(*,240) min_i, cf_codes(bf_table(2,min_i)),
     .                      prn_list(bf_table(3,min_i)),
     .                      dL21, dL21t, 
     .                     (dLC_obs-dLCs), (dLC_obs-dLC) 
 240                format('MOD: Changing ',i5,1x,a4,' PRN ',i2.2,
     .                    ' dL21 ',2F6.1, ' DLC ',2F10.3)
                    dL21 = dL21t
                    dLCm = dLCt
                    dL1  = dL1t
                    dL2  = dL2t
                    dLC  = dLCs
                 end if
             end if

                 
                    
*            Now force values
             do i = 1,tot_bf
                bf_neq_col(i) = bf_neq(i,min_i)
             end do

             avat = bf_neq_col(min_i)

             dchi = (dL21-bf_bvec(min_i,1))**2/avat

             write(*,320) min_i, cf_codes(bf_table(2,min_i)),
     .                       prn_list(bf_table(3,min_i)),
     .                       bf_bvec(min_i,1), 
     .                       sqrt(bf_neq(min_i,min_i)),
     .                       bf_bvec(min_i,2), dchi, dL1, dL2,dLC,
     .                       bf_bvec(min_i,2)-dLC
 320         format('FORCE ',i5,1x,a4,' PRN ',i2.2,1x,3F10.3, E10.3, 
     .              4F8.2)

             do i = 1, tot_bf
                equ_gn(i) = bf_neq_col(i)/avat
             end do
*            If the bias seems well fixed, then update the solution
*            Otherwize mark as not fixed and do not propagate
             if( abs(bf_bvec(min_i,2)-dLC).lt.0.1d0 ) then
                 do i = 1,tot_bf
                    dsol = equ_gn(i)*(dL21 - bf_bvec(min_i,1))
                    bf_bvec(i,1) = bf_bvec(i,1) + dsol
                    dsol = equ_gn(i)*(dLC - bf_bvec(min_i,2))
                    bf_bvec(i,2) = bf_bvec(i,2) + dsol
                 end do
                 do i = 1, tot_bf
                    do j = 1, tot_bf
                       dcov = equ_gn(i)*bf_neq_col(j)
                       bf_neq(i,j) = abs(bf_neq(i,j) - dcov)
                    end do
                 end do
                 call sbit(bf_table(5,min_i),1,1)
             else
                 call sbit(bf_table(5,min_i),1,0)
             endif

             ns = bf_table(2,min_i)
             lv = bf_table(3,min_i)
             do ep = bf_table(1,min_i),bf_table(1,min_i)+
     .                                 bf_table(4,min_i)
                ch = ltoc(ctol_cse(1,ns,ep),lv,actual_max_chan)
                if( ch.gt.0 ) then
                   L1_cyc_cse(ch,ns,ep) = L1_cyc_cse(ch,ns,ep)-dL1
                   L2_cyc_cse(ch,ns,ep) = L2_cyc_cse(ch,ns,ep)-dL2
               endif
             end do
             fixing = .true.
         end if
      end do 
*
*     Thats all
      do i = 1, tot_bf
         write(*,420) i, cf_codes(bf_table(2,i)), 
     .                prn_list(bf_table(3,i)),bf_table(1,i),
     .                bf_table(1,i)+bf_table(4,i), bf_table(5,i),
     .                bf_bvec(i,1), sqrt(bf_neq(i,i)), bf_bvec(i,2)
 420     format('BF_UPD ',I5,1x,a4,' PRN ',I2.2,' Eps ',2I6,1x,o3,
     .          4F12.4)
      end do

      return
      end

CTITLE gen_gamit_dd
 
      subroutine gen_gamit_dd(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse)

      implicit none
 
*     Routine to write out the double difference ambiquities to 
*     be used by solve.  Scheme here is to generate all the baselines
*     form a single station and then to replace longer baselines with
*     shorter ones.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)
 

* LOCAL VARIABLES

      integer*4 nblen   ! number of baselines
     .,   bfi((max_cfiles*(max_cfiles+1))/2) ! Index of baselines, 
                        ! set 0 when baseline included in sort list
     .,   bfj((max_cfiles*(max_cfiles+1))/2) ! Sort list of baseline lengths
     .,   i,j,k         ! Loop countres
     .,   pair(2)       ! Sites in current baseline
     .,   last          ! BF table entries of the last bias on each satellite
                        ! in entries being tested
     .,   numout, numfxd  ! Number of DD biases writen and of these the number
                        ! fixed.
     .,   last_bfs(max_gprn,max_cfiles)  ! List of last bias flags in the bf_table each
                        ! site and satellite
     .,   lbf_status(max_gprn,max_cfiles)  ! Status of last bias flags, Set -1 if no 
                        ! data, 0 initially and then 1 once it hs been used
     .,   nnew, nuse    ! Number of satellites new for baseline (not used at either
                        ! site) and used (used at both sites)
     .,   svnew(max_gprn),svuse(max_gprn) ! List of satellites new for baseline and 
                        ! used for baselined.  Only nnew-1 ambiguities are added at
                        ! each baseline, but differences can be formed with nnew+nuse-1
                        ! satellites
     .,   js, ks, jt    ! Satellite numbers
     .,   num_last_bfs  ! Number of last bias flags
     .,   bl_ent        ! Function to return entry in baseline table for sites i,j
     .,   longcom_dur, longcom_bf ! Longest duration and bf_table number for satellite
                        ! visible at all sites
     .,   longsng_dur, longsng_bf ! Longest duration and bf_table number for satellite
                        ! that is visible only at some stations
     .,   longone_dur, longone_bf ! Longest duration and bf_table number for satellite
                        ! that is visible only at some stations
     .,   bf            ! Generic bf_table entry number
*    .,   namb          ! Entry in ambiquity tables (computed from site and satellite 
*                       ! number)
*    .,   amb_tab(5,max_cfiles*max_gprn)  ! Ambiquity table.  Entries in table are:
*                       ! 1 -- Site 1 (starts as ref_site)
*                       ! 2 -- Site 2 (Second site in baseline)
*                       ! 3 -- SVS  1 (starts as Ref_svs)
*                       ! 4 -- SVS  2 (second satellite in DDiff)
*                       ! 5 -- Status: 1 if amb OK, -1 if reference site/svs
*                       !      -2 No data on satellite
     .,   kbf(4)        ! The four bf_table entries that make up double differnce
                        ! ambiquity being output
     .,   cbl           ! Baseline number
     .,   chk_overlap   ! Function to check overlap of double diffs
     .,   overlap, max_ovl, max_sv  ! Overlap, max value and satellite associated
                        ! maximum value
     .,   eps, epe      ! Epoch start and stop of overlap
     .,   numovl        ! Number of values in overlap for DDs
     .,   DD(2,max_cfiles-1)  ! Baselines that go into ambiquities for
                        ! each satellite.  Site entries set negative if not needed
     .,   DDBL(max_cfiles-1,max_gprn)  ! Baseline numbers associated with DD and
                        ! saved by satellite
     .,   usedDD(max_cfiles-1)  ! Flag to show baseline has been used.
     .,   svs_list(max_gprn)    ! Sorted list of amount of satellite data for each
                        ! baseline.
     .,   num_list(max_gprn)    ! Number of data per satellite
     .,   expected_namb ! Expected number of ambiguities
     .,   nzero         ! Number of zeroes
     .,   usedSite(max_cfiles)  ! Status of site in forming baselines
                        ! 0 - Site not used yet, -1 no data on PRN, used
     .,   minblindx     ! Minimum baseline index in sorted list (proxy for length)
     .,   ib, kb, nb    ! Sites with shortest baselines
 

      real*8  WLAmb         ! Value for widelane ambiquitiy


      real*8 maxlen     ! Maximum length of currently shortest baseline
     .,   confid        ! Sum of confidence for dd wl.

      real*8 mean(5), rms(5), sigm(5)  ! Mean, WRMS and sigma of mean
                     ! for the five observables (MW-WL, L1-L2, L1, L2 and LC)
  
      logical kbit      ! check bit status
     .,   fixd          ! Set true is all oneways in dd are fixed
     .,   OK_svs        ! Set true to indicate satellite present at all sites.
     .,   OK_site       ! Set true is site has last bias on all satellites or
                        ! all stations are missing that same satellite
     .,   bldone        ! Set true when baselines all filled

      character*1 type  ! Set X for fixed, R for Free bias

      write(uns,120)
 120  format('ACBIAS',4x,'double_difference biases:')

*     First sort the baselines into increasing length
      nblen = (num_cfiles)*(num_cfiles-1)/2
*     Form an index array of all baselines, when one is included
*     in the sorted list, this entry is set zero so that it 
*     won't be used again.  This covers the cases when some baselines
*     have exactly the same length.
      do i = 1, nblen
          bfi(i) = i
      end do

****  Check for sites that have been deleted and move these to the
*     end if of the length list
      do i = 1, num_cfiles-1
         do j = i+1, num_cfiles
            k = bl_ent(i,j)
            if( lc_num(i).eq.0 .or. lc_num(j).eq.0 ) then
                baselens(k) = baselens(k) + 20.d6
            endif
         end do
      end do

      do i = 1, nblen
          k = 0
          maxlen = 1.d10
          do j = 1,nblen
             if( baselens(j).le.maxlen .and. bfi(j).ne.0 ) then
                 maxlen = baselens(j)
                 k = j
             endif
          end do
          bfj(i) = k
          bfi(k) = 0
      end do

***** Generate the list of last bias flags for each site and satellites
      num_last_bfs = 0
      do i = 1, num_cfiles
         do j = 1, num_sat
             last_bfs(j,i) = 0  ! Clear the point for each site/sv combination
                                ! (If value remains zero then no data on sv)
          end do
      end do
      do i = 1, tot_bf
         call get_last_bias(i,last)
         if( last.eq.i ) then
             last_bfs(bf_table(3,i),bf_table(2,i)) = i
             num_last_bfs = num_last_bfs + 1
         end if
      end do
      write(*,140) tot_bf, num_last_bfs
 140  format('Of ',i5,' Total biases, there are ',i5,' last biases')

***** Find the satellite that is present at all stations and is the longest 
*     length.  We will use this station and satellite as the base
      longcom_dur = 0
      longsng_dur = 0
      longone_dur = 0
      do i  = 1, num_cfiles
         write(*,150) cf_codes(i),(last_bfs(j,i),j=1,num_sat)
* MOD MAF 210701: Updated 32(I5) to 50(I5) to allow for 45 Beidou satellites
c150     format('LBF ',a4,32(I5))
 150     format('LBF ',a4,50(I5))

* MOD TAH 050410: Only consider sites that have data on all biases
        OK_site = .true.
        do j = 1, num_sat
*           Station is missing a satellite, exclude from consideration
            if( last_bfs(j,i).eq.0 ) then
               OK_svs = .false.
*              Check to see if any other stations have the same satellite
               do k = 1, num_cfiles
                  if( last_bfs(j,k).gt.0 ) OK_svs = .true.
               end do
*              If some other sites have this satellite then exclude this
*              station.
               if( OK_svs ) OK_site = .false.
            end if
        end do
*       If site excluded, report now, otherwise use it for check
        if( .not. OK_site ) then
            write(*,160) cf_codes(i)
 160        format('Removing ',a4,' as possible reference site due ',
     .             'to missing satellite')
* MOD TAH 050618: Keep track of longest here in case no station can 
*       see all satellites
            do j = 1, num_sat
               bf = last_bfs(j,i)
               if( bf.gt.0 ) then
                   if( bf_table(4,bf).gt.longone_dur ) then
                       longone_dur = bf_table(4,bf)
                       longone_bf  = bf
                   end if
                endif
             end do

        ELSE
*          Check for longest sequence
           do j = 1, num_sat
               bf = last_bfs(j,i)

               OK_svs = .true.
               if( bf.gt.0 ) then
                  do k = 1, num_cfiles
*                    Only check sites that have at least some data
                     if( last_bfs(j,k).eq.0 .and. LC_num(k).gt. 0 ) 
     .                                              OK_svs = .false.
                  end do
                  
*                 If this satellite is everywhere, check the length
                  if( OK_svs ) then
                      if( bf_table(4,bf).gt.longcom_dur ) then
                          longcom_dur = bf_table(4,bf)
                          longcom_bf  = bf
                      end if
                  else     ! Save the non-all present longest one
                      if( bf_table(4,bf).gt.longsng_dur ) then
                          longsng_dur = bf_table(4,bf)
                          longsng_bf  = bf
                      end if
                  endif
               endif
           end do
        endif
      end do
*
*     Report the longest last bias flag and save the referernces 
      if( longcom_dur.gt.0 ) then
         write(*,200) cf_codes(bf_table(2,longcom_bf)), 
     .                prn_list(bf_table(3,longcom_bf)),
     .                bf_table(1,longcom_bf), bf_table(4,longcom_bf)
 200     format('Longest com BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longcom_bf)
         ref_svs  = bf_table(3,longcom_bf)
      else if( longsng_dur.gt.0 ) then
         write(*,210) cf_codes(bf_table(2,longsng_bf)), 
     .                prn_list(bf_table(3,longsng_bf)),
     .                bf_table(1,longsng_bf), bf_table(4,longsng_bf)
 210     format('Longest sng BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longsng_bf)
         ref_svs  = bf_table(3,longsng_bf)
      else if( longone_dur.gt.0 ) then
         write(*,220) cf_codes(bf_table(2,longone_bf)), 
     .                prn_list(bf_table(3,longone_bf)),
     .                bf_table(1,longone_bf), bf_table(4,longone_bf)
 220     format('Longest one BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longone_bf)
         ref_svs  = bf_table(3,longone_bf)
      else
         write(*,230)
 230     format('NO REF SITE and SVS FOUND: Stop')
         stop 'AUTCLN: No ref site and satellite found'
      endif

****  Initialize the last bias flag status (lbf_status)
      do i = 1, num_cfiles
         do j = 1, num_sat
             if( last_bfs(j,i).eq.0 ) then
*                There is not last bias flag for this 
*                combination
                 lbf_status(j,i) = -1
             else
                 lbf_status(j,i) = 0
                 if( LC_SVS_num(j,i).eq.0 ) then
                    write(*,320) cf_codes(i), prn_list(j),
     .                 lbf_status(j,i) 
 320                format('GAMITDD Warning no data for ',a4,
     .                     ' PRN',i2.2,' Last BF ',i5)
                 endif
             endif
         end do
      end do

      namb = 0
      do j = 1, num_cfiles
         write(*,360) cf_codes(j),(lbf_status(k,j),k=1,num_sat)
         do k = 1, num_sat
            if( lbf_status(k,j).eq.-1 ) then
                namb = namb + 1
            endif
         end do
      end do

*     Now account for missing stations and satellites
      ib = 0
      do j = 1, num_cfiles
*         See if all satellites missing
          ib = 0
          do k = 1, num_sat
             if( lbf_status(k,j).ne.-1 ) ib = ib + 1
          end do
*         Decrement namb if no last bias flags at all
          if( ib.eq.0 ) namb = namb - 1
      enddo
*     Repeat for satellites
      do k = 1, num_sat
         ib = 0
         do j = 1,num_cfiles
            if( lbf_status(k,j).ne.-1 ) ib = ib + 1
          end do
*         Decrement namb if no last bias flags at all
          if( ib.eq.0 ) namb = namb - 1
      enddo

      expected_namb = (num_cfiles-1)*(num_sat-1)-namb
      write(*,325) (num_cfiles-1)*(num_sat-1), namb, 
     .             expected_namb
 325  format('GAMITDD: All biases ',i4,' Missing SVS ',I4,
     .       ' Expected ambiguities ',i4)

C     do i = 1, nblen
C        call ent_bf(bfj(i),num_cfiles, pair)
C        write(*,365) i, bfj(i), pair, cf_codes(pair(1)), 
C    .                cf_codes(pair(2))
C 365    format('BL ',i4,i5, 2i4, 1x, a4,1x,a4)
C     end do

****  Set up baselines to be used for each satellite
      nzero = 0
      do j = 1, num_sat
*        Initilize the DD array with simple baselines
C        k = 0
C        OK_site = .false.
C        do while ( k.lt. num_cfiles-1 .and. .not.OK_site )
C           k = k + 1
C           if( lbf_status(j,k).eq.0 ) then
C              OK_site = .true.
* MOD TAH 080414: Looks like index should start at k+1 not k.
C              do i = k+1, num_cfiles
C                 if( lbf_status(j,i).eq. 0 ) then
C                    DD(1,i-1) = k
C                    DD(2,i-1) = i
C                    usedDD(i-1) = 0
C                 else
*                    No data for this baseline so leave out
C                    DD(1,i-1) = -k
C                    DD(2,i-1) = -i
C                    usedDD(i-1) = -1
C                 endif
C              end do
C           else
C              DD(1,k) = -1
C              DD(2,k) = -k
C              usedDD(k) = -1
C           end if
C        end do
C        write(*,327) prn_list(j),(DD(1,k),DD(2,k),k=1,num_cfiles-1)
C327     format('DD PRN',i2.2,' DDPairs ',50(2I3,2x))
C              
*        OK Now fill out the baseline in ascending order.
*        Thus algorithm will skip shorter baselines but
*        should always generate independent baselines.
C        do i = 1,nblen
C           call ent_bf(bfj(i),num_cfiles, pair) 
C           if( usedDD(pair(2)-1).eq.0  ) then
*               Make sure that site 1 has data on this satellite
C               if( lbf_status(j,pair(1)).eq.0 ) then
C                  DD(1,pair(2)-1) = pair(1)
C                  usedDD(pair(2)-1) = 1
C               endif
*           Try reversing the baselines
C            elseif( usedDD(pair(1)-1).eq.0  ) then
*               Make sure that site 1 has data on this satellite
C                if( lbf_status(j,pair(2)).eq.0 ) then
C                   DD(1,pair(1)-1) = pair(2)
C                   usedDD(pair(1)-1) = 1
C                endif
C           endif
C        end do
*        Now switch any reversed baselines
C         do k = 1, num_cfiles -1
C            if( DD(1,k).gt.DD(2,k) ) then
C                cbl = DD(1,k)
C                DD(1,k) = DD(2,k)
C                DD(2,k) = cbl
C            endif
C         enddo

*     Now generate the list of independent baselines in ascending order
*     Start with the shortest baseline and then get closest stations
*     to end and work out from there finding next closest station.
*
         do i = 1, num_cfiles-1
            DD(1,i) = 0
            DD(2,i) = 0
            usedDD(i) = 0
         end do
         do i = 1, num_cfiles
            usedSite(i) = 0
            if( lbf_status(j,i).eq.-1 ) usedSite(i) = -1
         end do

*        Now fill the DD array with the shortest baselines that can see
*        this satellite
         bldone = .false.
*        Find the first baseline with data on this satellite
         k = 0
         do while ( k.lt. nblen .and. .not. bldone )
             k = k + 1 
             call ent_bf(bfj(k),num_cfiles, pair)
             if( lbf_status(j,pair(1)).eq.0 .and.
     .           lbf_status(j,pair(2)).eq.0 ) then
*                Baseline has good data at both sites so start
*                here.
                 DD(1,1) = pair(1)
                 DD(2,1) = pair(2)
                 usedDD(1) = 1
                 bldone = .true.
                 usedSite(pair(1)) = 1
                 usedSite(pair(2)) = 1
             endif
         end do
         if( .not.bldone) then
             bldone = .true.  ! No use searching for more baselines
         else
             bldone = .false.
         endif

*****    OK Now fill out all the baseline
         nb = 1
         do while ( .not. bldone )
             minblindx = nblen+1
             ib = 0
             do i = 1, num_cfiles
                if( usedSite(i).gt.0 .and.
     .              lbf_status(j,i).eq. 0 ) then
*                  Site has been used and has data on satellite.
*                  Check lengths to unused used site that have 
*                  good data
                   do k = 1, num_cfiles                      
                      if ( usedSite(k).eq.0 .and.
     .                     lbf_status(j,k).eq. 0 ) then
                         cbl = bl_ent(i,k)
*                        Find where this baseline ranks
                         do ks = 1, nblen 
                            if( bfj(ks).eq.cbl ) js = ks
                         end do
                         if( js.lt.minblindx ) then
*                            This is a shorter baseline so save
                             minblindx = js
                             kb = k
                             ib = i
                         end if
                      end if
                   end do
                end if
             end do
*            OK: If ib and kb have been set then we have
*            found a good baseline so add
             if( ib.gt. 0 ) then
                 nb = nb + 1
*                See if we need to reverse the baseline.
                 if( kb.lt.ib ) then
                    k = kb 
                    kb = ib
                    ib = k
                 endif
                 DD(1,nb) = ib
                 DD(2,nb) = kb
                 usedDD(nb) = 1
                 usedSite(ib) = usedSite(ib)+1
                 usedSite(kb) = usedSite(kb)+1
             else
*                Nothing left so we are done
                 bldone = .true.
             endif
         end do

         if( 1.eq.2 )
     .   write(*,328) prn_list(j),(DD(1,k),DD(2,k),usedDD(k), 
     .                             k=1,num_cfiles-1)
 328     format('DD PRN',i2.2,' DDTrip ',50(2I3,2x,i2,1x))

*        Now assign baseline numbers to the DD entries and save
*        by satellite
         do i = 1,num_cfiles-1
            if( DD(1,i).gt.0 ) then
                cbl = bl_ent(DD(1,i),DD(2,i))
                DDBL(i,j) = cbl
            else
                DDBL(i,j) = 0   ! No baseline with this site
            end if
         end do

         if( 1.eq.2 )
     .   write(*,330) prn_list(j),(DDBL(i,j),i=1,num_cfiles-1)
 330     format('DDBL PRN',i2.2,' BL ',50I5)
         do i = 1,num_cfiles-1
            if ( DDBL(i,j).eq.0 ) then
              nzero = nzero + 1
            endif
         end do

      end do
C     write(*,335) nzero
C335  format('GAMITDD: Nzeros in DDBL ',i4)                     

            
****  Now scan down the baselines in increasing length to
*     set up double difference biases
      namb = 0          
      do i = 1, nblen
*        Get the sites in this baseline
         call ent_bf(bfj(i),num_cfiles, pair)

*        OK: For this baseline, generate sorted list of satellites
         do j = 1,num_sat
             kbf(1) = last_bfs(j,pair(1))
             kbf(2) = last_bfs(j,pair(2))
             if( kbf(1).gt.0 .and. kbf(2).gt.0 ) then
                overlap = chk_overlap(kbf(1), kbf(2),0,0,eps, epe)
             else
                overlap = -9999
             end if
*            See where this fits in list
             OK_svs = .false.
             do js = 1, j-1
                 if( overlap.lt.num_list(js) .and. .not.OK_svs ) then
*                    Move list up so we can insert
                     do ks = j-1,js, -1
                        num_list(ks+1) = num_list(ks)
                        svs_list(ks+1) = svs_list(ks)
                     end do
                     num_list(js) = overlap
                     svs_list(js) = j
                     OK_svs = .true.
                  endif
              end do
              if( .not. OK_svs ) then   ! Add to end
                  num_list(j) = overlap
                  svs_list(j) = j
                  OK_svs = .true.
              endif
         end do
         if( 1.eq.2 ) 
     .   write(*,337) pair,cf_codes(pair(1)),
     .                cf_codes(pair(2)),
     .                (svs_list(j),num_list(j),j=1,num_sat)
 337     format('GAMITDD: BSRT ',2i3,1x,a4,1x,a4,1x,50(I2.2,1x,i5,2x)) 
        
*        Now find all new and used satellite pairs for this
*        baseline
         nnew = 0
         nuse = 0
         do jt = 1, num_sat
*           Get the sorted satellite number.
            j = svs_list(jt)
*           Set if not used and not -1
c            if( (lbf_status(j,pair(1)).eq.0 .or.
c     .          lbf_status(j,pair(2)).eq.0 ).and.
c     .          (lbf_status(j,pair(1)).ge.0 .and.
c    .           lbf_status(j,pair(2)).ge.0)   ) then
*               Satellite at one of these stations has
*               not been used yet
c                nnew = nnew + 1
c                svnew(nnew) = j
c            endif
*           See if we use this baseline for this satellite
            do k = 1, num_cfiles - 1
               if( DDBL(k,j).eq.bfj(i) ) then
*                  This is baseline for this site so add
                   nnew = nnew + 1
                   svnew(nnew) = j
               end if
            end do 

*           See if satellites used before but not one that we just
*           added.
            if( lbf_status(j,pair(1)).gt.0 .and.
     .          lbf_status(j,pair(2)).gt.0  .and.
     .          svnew(nnew).ne. j   ) then
*               Satellite has been used and so can be
*               for double differences
                nuse = nuse + 1
                svuse(nuse) = j
            endif
         end do


****     OK: Now see if we have any new satellites
         if( nnew.gt.0 ) then
*            For each new satellite, find longest overlap
*            others.
* NOTE: Could replace here with code that generates smallest sigma
*       for the double differences!.
             max_ovl = 0
             max_sv = 0
             ks = 0
             do j = 1, nnew
                max_ovl = 0
                js = svnew(j)
                kbf(1) = last_bfs(js,pair(1))
                kbf(2) = last_bfs(js,pair(2))
                do k = j+1,nnew
                    ks = svnew(k)
                    kbf(3) = last_bfs(ks,pair(1))
                    kbf(4) = last_bfs(ks,pair(2))
                    if( kbf(1).gt.0 .and. kbf(2).gt.0 .and. 
     .                  kbf(3).gt.0 .and. kbf(4).gt.0 ) then
                        overlap = chk_overlap(kbf(1), kbf(2), 
     .                        kbf(3), kbf(4),eps, epe)
                        if( overlap.gt. max_ovl ) then
                            max_ovl = overlap
                            max_sv  = ks
                        end if
                    end if
                end do 

*               Now test against reference list to see if longer
                do k = 1,nuse
                    ks = svuse(k)
                    kbf(3) = last_bfs(ks,pair(1))
                    kbf(4) = last_bfs(ks,pair(2))
                    if( kbf(1).gt.0 .and. kbf(2).gt.0 .and. 
     .                  kbf(3).gt.0 .and. kbf(4).gt.0 ) then
                        overlap = chk_overlap(kbf(1), kbf(2), 
     .                        kbf(3), kbf(4),eps, epe)
                        if( overlap.gt. max_ovl ) then
                            max_ovl = overlap
                            max_sv  = ks
                        end if
                    end if
                end do

*****           It is possible that there is no overlap on
*               this baseline for this new satellite.  This can
*               happen on early segements of data
                if( max_ovl.eq. 0 .and. js.ne.ks) then
                    max_ovl = -1
                    if( ks.gt.0 ) then
                       max_sv  = ks  ! Use last satellite?
                    else
                       max_sv = svs_list(num_sat)
* MOD TAH 070410: Make sure not same SVS
                       if( max_sv.eq.js ) then
                           max_sv = svs_list(num_sat-1)
                       endif
                    end if
                endif

                if( max_ovl.ne.0 ) then
                    namb = namb + 1
                    amb_tab(1,namb) = pair(1)
                    amb_tab(2,namb) = pair(2)
                    amb_tab(3,namb) = js
                    amb_tab(4,namb) = max_sv
                    amb_tab(5,namb) = max_ovl
                    lbf_status(js,pair(1))     = 
     .                                   lbf_status(js,pair(1))+1
                    lbf_status(max_sv,pair(1)) = 
     .                                   lbf_status(max_sv,pair(1))+1
                    lbf_status(js,pair(2))     = 
     .                                   lbf_status(js,pair(2))+1
                    lbf_status(max_sv,pair(2)) = 
     .                                   lbf_status(max_sv,pair(2))+1
                    write(*,340) namb, cf_codes(pair(1)),
     .                           cf_codes(pair(2)),
     .                           prn_list(js),prn_list(max_sv), max_ovl
 340                format('GAMITDD Amb ',i4,1x,a4,1x,a4,1x,I2.2,
     .                     1x,I2.2,1x,' MAX_OVL ',i5)
                endif
             enddo 
         endif      
      enddo 

      do i = 1, num_cfiles
         write(*,360) cf_codes(i),(lbf_status(j,i),j=1,num_sat)
* MOD MAF 210701: Updated 32(I5) to 50(I5) to allow for 45 Beidou satellites
c360     format('GAMITDD ',a4,1x,32(I5))
 360     format('GAMITDD ',a4,1x,50(I5))
      end do          

*
***** OK: All done, write out table
      numout = 0
      numfxd = 0
      do i = 1,namb

         if( amb_tab(5,i).ne.0 ) then
             fixd = .false.
             type = 'R'
             kbf(1) = last_bfs(amb_tab(3,i),amb_tab(1,i))
             kbf(2) = last_bfs(amb_tab(3,i),amb_tab(2,i))
             kbf(3) = last_bfs(amb_tab(4,i),amb_tab(1,i))
             kbf(4) = last_bfs(amb_tab(4,i),amb_tab(2,i))
 
             if( kbf(1).gt.0 .and. kbf(2).gt.0 .and. 
     .           kbf(3).gt.0 .and. kbf(4).gt.0 ) then 
                if( kbit(bf_table(5,kbf(1)),1) .and.
     .              kbit(bf_table(5,kbf(2)),1) .and.
     .              kbit(bf_table(5,kbf(3)),1) .and.
     .              kbit(bf_table(5,kbf(4)),1) ) then
                    fixd = .true.
                    type = 'X'
                    numfxd = numfxd + 1
                endif
                confid = 1.d0/(1.d0/wl_conf(kbf(1)) + 
     .                         1.d0/wl_conf(kbf(2)) +
     .                         1.d0/wl_conf(kbf(3)) +
     .                         1.d0/wl_conf(kbf(4)) )
                overlap = chk_overlap(kbf(1), kbf(2), 
     .                        kbf(3), kbf(4),eps, epe)
             else
                overlap = 0
                write(*,250) cf_codes(i), prn_list(j),
     .             (amb_tab(k,i),k=1,5), kbf
 250            format('NOKBF: ',a4,' PRN ',i2.2,' AMB_TAB ',
     .              5i4,' KBF ',4i6)
             end if
   
****         Do final calc on wide lane
             WLAmb = 0
             numovl = 0
             if( overlap.gt.0 ) then
                 call check_gamitdd(kbf(1),kbf(2), kbf(3), 
     .                kbf(4), eps, epe, 
     .                L1r_phs_cse, L2r_phs_cse, L1_cyc_cse,
     .                L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .                ctol_cse, data_flag_cse, azel_cse, 
     .                mean, rms, sigm, numovl )

****             See if we should adjust value
                 cbl = bl_ent(amb_tab(1,i),amb_tab(2,i))
                 if( numovl.gt.min_wl_tol ) then
                     call eval_WLAmb(mean,
     .                    sigm, baselens(cbl), WLAmb, type )
                 endif
                     
             else
                 mean(1) = 0.d0
                 sigm(1) = 0.d0
                 mean(2) = 0.d0
                 sigm(2) = 0.d0
             endif

             if ( overlap.lt.0 ) then 
*               check if any over lap
                max_ovl = 0
                do k = 1, num_cfiles
                   if( k.ne.amb_tab(2,i) ) then
                      kbf(3) = last_bfs(amb_tab(4,i),k)
          
                      overlap = chk_overlap(kbf(3), kbf(4), 
     .                       0, 0,eps, epe)
                      if( overlap.gt.max_ovl ) max_ovl = overlap
                   end if
                end do
                if( max_ovl.eq.0 ) then
                   overlap = -1
                else
                   overlap =  0
                endif
             endif
             k = bl_ent(amb_tab(1,i),amb_tab(2,i))
* MOD TAH 070405: Flip the sign of the WLAmb value written out
* Seems to be long standing bug. (WLAmb replaced with -WLAmb).
* Flip the signs of the means as well to avoid confusion.
             write(uns,420) cf_codes(amb_tab(1,i)),
     .                      cf_codes(amb_tab(2,i)),
     .                      prn_list(amb_tab(3,i)),
     .                      prn_list(amb_tab(4,i)), -WLAmb,
     .                      confid, type, overlap, numovl,
     .                      -mean(1), sigm(1), 
     .                      -mean(2), sigm(2),
     .                      baselens(k)/1000.d0

 420         format('ACBIAS',4x,'B1L2 ',a4,1x,a4,1x,
     .                i2.2,1x,i2.2,1x,F6.1,1x,
     .                F10.2,1x,a1,4x,' ! ',2i5,
     .              ' WL Mean ',F7.3,' +- ',F7.3,
     .              ' L1-L2 ',F7.3,' +-  ',F7.3,
     .              ' LEN ',F8.2,' km')

             numout = numout + 1
         else 
             write(*,480) i,cf_codes(amb_tab(1,i)),
     .                      cf_codes(amb_tab(2,i)),
     .                      prn_list(amb_tab(3,i)),
     .                      prn_list(amb_tab(4,i)),
     .                      amb_tab(5,i)
 480         format('GAMITDD Warning: No data on ambiguity ',
     .               i4,1x,a4,1x,a4,1x,' PRN ',2I3.2,1x,i5)
         endif
      end do
*
***** Tell user how many fixed
      write(uns,'(a)') 'ACBIAS    '
      write(uns,520) numout, num_last_bfs, numfxd,
     .               expected_namb
 520  format('NBIAS out ',2i6,' Number fixed ',i6,' Expected ',i6)

****  Thats all
      return
      end

CTITLE BL_ENT

      integer*4 function bl_ent(i,j)

      implicit none

*     Function to return the entry between stations i and j

* PASSED VARIABLES
      integer*4 i,j

*     Standard lower diagonal calculation
      if( i.lt.j ) then
         bl_ent = ((j-2)*(j-1))/2 + i
      elseif ( j.lt.i ) then
         bl_ent = ((i-2)*(i-1))/2 + j
      else
         write(*,120) i,j
 120     format('**WARNING** BL_ENT called with same station'
     .         ,' pair ',2i4)
         bl_ent = 0
      endif

****  Thats all
      return
      end


CTITLE Check_used

C     logical function check_used(refbf, last, bf_baseline, 
C    .                            bfj, nbl,lv, bf_dep)

C     implicit none

*     Routine to check to see if the one ways have been used with
*     a common station previously.  If they have then they are 
*     marked as used, other wise used is set false.

* INCLUDES
C
C     include '../includes/kalman_param.h'
C     include '../includes/const_param.h'
C     include 'ctogobs_com.h'

* PARAETERS 

C     integer*4 max_stack  ! Maximum stack size while searching
C    .,         max_bls    ! Maximum number of baselines

C     parameter ( max_stack = max_cfiles )
C     parameter ( max_bls = (max_cfiles*(max_cfiles+1))/2)
C
* PASSED VARIABLES
C     integer*4 refbf(2) ! Entries in the bf_table for the reference
C                        ! satellite
C    .,  last(2)  ! Entries in bf_table for the new oneways
C              ! being considered
C     integer*4 bf_baseline(max_bls) ! Baselines that
C                 !  have been used.  Bit set for each satellite as it is used
C    .,  bf_dep(max_bls) ! Bit set if baseline is 
C                 ! a dependent one.  These baselines are not check when
C                 ! when other dependent/independent baselines are searched.

C    .,   bfj(max_bls) ! Sort list of baseline lengths
C    .,   nbl     ! Baseline number from the bfj list being checked.
C    .,   lv      ! Satellite being processed

* LOCAL VARIABLES
C     integer*4 i,j      ! looop counter
C    .,  pair(2)         ! Pairs of sites in baselines
C    .,  s1, s2          ! The two sites being considered
C    .,  bls(max_bls) ! Stack of baselines that need to be further
C                        ! checked
C    .,  rfs(max_bls) ! Stack of reference sites for the baselines.
C                        ! When one reference site is found, it is pushed on the stack and
C                        ! the second site in the baseline becomes the reference.
C    .,  nb              ! Number of current baseline being processed.
C    .,  ns              ! Number of entries on current stack
C    .,  rs              ! Current reference site for baseline
C    .,  it              ! Iterations though code

C     logical done       ! Logical set true when all baselines have been scanned or a 
C                        ! linkage of baselines between sites has been found.
C    .,   s1_used, s2_used  ! Set true if the sites have been used in an eariler
C                        ! baseline.
C    .,   onstack        ! Set true if the current baseline we are condsidering is already
C                        ! on the stack.
C    .,   kbit           ! Checks bit status
C
***** First check the baselines we have processed already and
*     see if one or both of these sites have never been used 
*     before.  If one of them has not been used, then clearly
*     independent.
C     s1 = bf_table(2,refbf(1))   ! S1 and S2 are the two sites in the baseline
C     s2 = bf_table(2,refbf(2))
C     s1_used = .false.
C     s2_used = .false.

C     do i = 1, nbl - 1  ! Test up to the baseline before this one
C         call ent_bf(bfj(i), num_cfiles, pair)
C         if( pair(1).eq.s1 .or. pair(2).eq.s1 ) s1_used = .true.
C         if( pair(1).eq.s2 .or. pair(2).eq.s2 ) s2_used = .true.
C     end do

*     If either site has not been used, then this is indepednent baseline
C     if( .not.s1_used .or. .not.s2_used ) then
C         check_used = .false.
C         RETURN
C     end if

****  OK: Now the tricky part.  Both sites have been used but there may not
*     be squence of baselines that connect them.  So now we need to search
*     down each baseline tree starting from S1 and see if we can get to S2
*     through the baselines.  The algoritm is basically recursive but is implemented
*     through a stack approach
C     rs = s1      ! Initial reference site (first in baseline)
C     done = .false.
C     nb = 0
C     ns = 0    ! Number of entries on stack
C     check_used = .false.
C     it = 0

C     do while ( .not. done )
*        Move to the next baseline
C        nb = nb + 1
C        it = it + 1
C        if( it.gt.64000 ) then
C           print *,' Iteration overflow S1, S2 ',s1,s2,nb, ns, 
C    .               (i,bls(i),i=1,ns)
C           check_used = .true.
C           return
C        endif
C        if( nb.lt.nbl ) then    ! OK, if this is before the baseline 
C                                ! we are considering
*           Check to see if this baseline is already on the stack, if it is 
*           not then continue checking 
C           onstack = .false.
C           do j = 1, ns
C              if( nb.eq.bls(j) ) onstack = .true.
C           end do
*           See if this is dependent baseline.  If it is then do not
*           search (indicate that is on the stack)
C           if( kbit(bf_dep(bfj(nb)),lv) ) onstack = .true.

C           if ( .not. onstack ) then
*               See if one of the sites matchs the current reference sites
C               call ent_bf(bfj(nb),num_cfiles, pair)
*               Switch the baseline order is necessary to make pair(1) be the reference site
C               if( pair(2).eq.rs ) then
C                   pair(2) = pair(1)
C                   pair(1) = rs
C               end if
C               if( pair(1).eq.rs .and. 
C    .              kbit(bf_baseline(bfj(nb)),lv) ) then    ! First site 
*                                                     matches and SV used
*                   Check to see if second site is actually the end of 
*                   the baseline chain
C                   if( pair(2).eq.s2 ) then
*                       Connection has been made to second site so this is not
*                       an independent baseline.  Get out of test
C                       check_used = .true.
C                       done = .true.
C                   else     ! Site does not match end yet, so push the 
C                            ! entry onto thestack and start
C                            !  searching baselines from this point
C                       ns = ns + 1
C                       if( ns.gt.nbl ) then
C                          write(*,999) ns, nbl, onstack, 
C    .                                  (i,bls(i),i=1,ns)
C999                       format('Stack overflow ',2i5,L3,/,
C    .                            (100(2i5,1x),/))
C                          stop 
C                       endif
C                       bls(ns) = nb
C                       rfs(ns) = rs
C                       nb = 0     ! Start searching baseline again
C                       rs = pair(2)
C                   end if
C               end if
C           end if       ! Not on the stack already
C        else    ! We have exceeded the baseline count, so pop the stack and follow
C                ! another path to see where it leads
C           if( ns.ge.1 ) then
C               nb = bls(ns)
C               rs = rfs(ns)
C               ns = ns - 1
C               if( nb.ge.nbl-1 ) done = .true.
C           else
C              done = .true.
C           end if
C        end if
C     end do
*
*     Scanned all used baselines or we done a path
C     return
C     end

CTITLE check_gamitdd

      subroutine check_gamitdd(bft, bfs, bfj, bfv, eps, epe, 
     .    L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse, mean, rms, sigm, numovl )

      implicit none

*     Routine to check the WL average value for an amb_tab entry
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
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
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

      real*8 mean(5), rms(5), sigm(5)  ! Mean, WRMS and sigma of mean
                     ! for the five observables
  
      integer*4 bft, bfs, bfj, bfv  ! Bias table entries
     .,         eps, epe   ! Overlap of epochs in double difference
     .,         numovl     ! Number of valeus in overlap
      

* LOCAL VARIABLES

      integer*4 i,j
     .,    ep           ! Epoch counter

      real*8 obs(5,4)   ! Obs from each oneway (columns are site
                        ! satellite combinations, rows are data
                        ! types
     .,      wgh(5,4)   ! Weight for each station/satellite
     .,      stats(4,5) ! Statistics: Columns for each data type
                        ! rows for res*wgh res**2*wgh wgh and num


      logical OK(4)    ! set true for each good station/satellite one-way

                   
*     Clear the accumulation arrays
      do i = 1,4  ! Loop over res*wgh, res**2*wgh wgh num 
         do j = 1,5  ! Loop over MW-WL, L2-L1, L1 and L2, LC
             stats(i,j) = 0.d0
         end do
      end do

*     Loop over the epoch range
      do ep = eps,epe
*       Get the residuals for each of one-ways that go into
*       the double difference
         call get_obs(ep, bf_table(2,bft),bf_table(3,bft),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,1), wgh(1,1), OK(1))
         call get_obs(ep, bf_table(2,bfs),bf_table(3,bfs),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,2), wgh(1,2), OK(2))
         call get_obs(ep, bf_table(2,bfj),bf_table(3,bfj),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,3), wgh(1,3), OK(3))
         call get_obs(ep, bf_table(2,bfv),bf_table(3,bfv),  
     .       L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .       L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .       azel_cse, obs(1,4), wgh(1,4), OK(4))

****     If all the oneways are good, increment the statistics
         if( OK(1) .and. OK(2) .and. OK(3) .and. OK(4) ) then
              call inc_zero_stats(obs, wgh, stats)
         endif
      end do

****  Finish the mean and RMS caculations
      numovl = stats(4,1)
      if( stats(4,2).eq.0 ) then
         if( 1.eq.2 )
     .   write(*,120) bft, cf_codes(bf_table(2,bft)),
     .        prn_list(bf_table(3,bft)), cf_codes(bf_table(2,bfv)), 
     .        prn_list(bf_table(3,bfv)), eps, epe
 120     format(' No data found for bias ',i4,
     .           ' Zero DD: Site ',a4,' PRN ',i2.2,'(with ',a4,
     .           ' PRN ',i2.2,') Epochs ',i5,' to ',i5)
          do i = 1,5
             mean(i) = 0
             rms(i) = 100
             sigm(i) = 100
          end do

         RETURN
      endif

***** For each type get the mean, RMS, and sigma of mean
      do i = 1,5
          if( stats(4,i).gt.0 ) then
             mean(i) = stats(1,i)/stats(3,i)
             rms(i) = sqrt(abs(stats(2,i)/stats(3,i)-mean(i)**2))
             sigm(i) = sqrt(rms(i)**2/stats(4,i))
          else
             mean(i) = 0.d0
             rms(i) = 100.d0
             sigm(i) = 100.d0
          endif
      end do

****  OK Output results
      write(*,220) cf_codes(bf_table(2,bft)),cf_codes(bf_table(2,bfv)), 
     .        prn_list(bf_table(3,bft)), prn_list(bf_table(3,bfv)),
     .        eps, epe, mean(1),sigm(1), mean(2), sigm(2),
     .        mean(5), sigm(5), numovl
 220  format('GAMITDD ',a4,1x,a4,1x,i2.2,1x,i2.2,' Ep ',i5,1x,i5,1x,
     .       'MWWL ',F7.3,' +- ',F7.3,' L1-L2 ',F7.3,' +- ',F7.3,
     .       ' LC  ',F7.3,' +- ',F7.3,' Num ',i5)

      return
      end

CTITLE eval_WLAmb

      subroutine eval_WLAmb(mean,sigm,
     .           blen, WLAmb, type )

      implicit none

*     Routine to make final adjustment to the WL ambiguity
*     this time using the ddiff passed to solve.
* MOD TAH 070312: Added some additional conditions to check quality of
*     fix and to weight L1-L2 values.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
      real*8 WLAmb     ! Returned value of the WL ambiquity

      real*8 mean(5)   ! Mean values of MW WL, L1-L2 plus
     .,      sigm(5)   ! Sigmas of mean values
     .,      blen      ! Length of baseline (m)

      character*(*) type ! First character X/R if fixed or not.

* LOCAL VARIABLES
      real*8    L21_wl   ! estimate of the DD widelane value
     .,         L21      ! estimate from L2-L1
     .,         err_mwwl, err_exwl  ! Error in MWWL and EXWL
     .,         sig_mwwl, sig_exwl  ! Assinged sigmas for MWWL and EXWL

      real*8 L21_scale   ! Scaling for L2-L1 estimate based 
                         ! on length
      logical fix_OK      ! Set true if ambiquity can be fixed

***** See how much we trust L1-L2 (should be good on short baselines)
      fix_OK = .false.
      WLAmb = 0.0d0

      L21_wl = nint(mean(1))
      L21    = nint(mean(2))
      err_mwwl = mean(1)-L21_wl 
      sig_mwwl = sigm(1)
*     See if 1/2 cycle allowed
      if( mean(1).eq. 0.0d0 ) then
         L21 = nint(mean(2)*2)/2.d0
      endif

* MOD TAH 070312: For length of 100 km, EXWL sigma is four times
*     larger and increases quadratically
      L21_scale = (200.d0+blen/1000.d0)/200.d0
      err_exwl = mean(2)-L21
      sig_exwl = sigm(2)*L21_scale

* MOD TAH 070312: See if L21 sigma is small enough that we should
*     use value
      if( sig_exwl.lt.0.10d0 ) then
*         We will use the EXWL. If is it inconsisent with MWWL, then
*         force sigma to be smaller
          if( L21_wl.eq.L21 ) then
             if( abs(err_exwl).lt.mdev_wl_tol .or. 
     .          (abs(err_mwwl).lt.mdev_wl_tol .and. 
     .           sig_mwwl.le. msig_wl_tol)         ) then
                WLAmb = L21
                fix_OK = .true.
             end if
          else
*            Only fix using EXWL if very close to integer and
*            small sigma
             if( abs(err_exwl).lt.mdev_wl_tol/2 .and.
     .           sig_exwl.lt.0.05d0 ) then
                 WLAmb = L21
                 fix_OK = .true.
             else if( abs(err_mwwl).lt.mdev_wl_tol .and.
     .           sig_mwwl.lt.msig_wl_tol/2.0 ) then
*                See if MWWL is "very good"
                 WLAmb = L21_wl
                 fix_OK = .true.
             endif
          endif
      else
*         We need to rely on the MWWL because the EXWL not well enough
*         deterimined.
          if( abs(err_mwwl).lt.mdev_wl_tol .and. 
     .        sig_mwwl.lt.msig_wl_tol ) then
             WLAmb = L21_wl
             Fix_OK = .true.
          end if
      endif

*     Set final values
      if( Fix_OK ) then
         type(1:1) = 'X'
      else
         type(1:1) = 'R'
      endif




c OLD Algorithm
C     if( abs(mean(1)-L21_wl).lt.mdev_wl_tol .and. 
C    .    sigm(1).lt.msig_wl_tol ) then
C         type(1:1) = 'X'
C         WLAmb = L21_wl
C     else if( abs(mean(2)-L21).lt.mdev_wl_tol/2 .and. 
C    .       sigm(2).lt.msig_wl_tol/10 .and. 
C    .       L21_scale.lt.1000.d0 ) then
C           type(1:1) = 'X'
C           WLAmb = L21
C     else
*        Neither test was passed
C        if( L21_wl .ne. L21 .and. L21_scale.lt.1000.d0 ) then
C           WLAmb = 0
C           type(1:1) = 'R'
C        end if

C     endif

****  Thats all
      return
      end

CTITLE scan_nodd
 
      subroutine scan_nodd(ctol_cse, data_flag_cse, bf_type_cse)

      implicit none
 
*     This routine scans all data and flags any data with no-double differences
*     Bias flags need to be pushed when this is done.
 
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

 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep)

* LOCAL VARIABLES
      integer*4 i,j         ! Site and satellite counter
     .,         ep          ! Epoch counter
     .,         ch          ! Channel for SV at epoch 
     .,         num_dded    ! Number of data edited for Site/SV
                            ! due to no dd's
     .,         c2          ! second channel at site to form SD
     .,         j2          ! Satellite associated with c2
     .,         k           ! Second station loop
     .,         l           ! Loop over second channel number
     .,         p2          ! SV numbers at second site
     .,         ddfnd       ! Counter for SV's at second site 
                            ! (needs to be 2 if double diff found)
     .,         ndd         ! Number off DD founds (ends do while
                            ! when first one found).
     .,         ltoc        ! Function that returns channel number
                            ! for specific SV (zero if not present)
     .,         ntot        ! Counts for number of data at site
     .,         nsvs        ! Cpunts for number of data per satellite

      logical data_OK       ! Function returns true for good data
     .,       good_bf       ! Function returns true for bias flag
                            ! on good data point.

      character*80 oneddwarn ! Warning when just 1 dd


****  Set the editing maskes
      call set_phs_mask( phs_mask, phs_bias_mask )

****  Loop over sites and satellites and epochs
      call report_stat('status','autcln','scan_nodd',' ',
     .                     'Starting scan',0)
      do i = 1, num_cfiles
         do j = 1, num_sat
            num_dded = 0
            do ep = 1, num_ep

*              OK: See if have data on this satellite at the epoch
               ch = ltoc(ctol_cse(1,i,ep),j, actual_max_chan)
               if( ch.gt.0 ) then   ! We have data on this satellite
                  if( data_OK(data_flag_cse(ch,i,ep),0, phs_mask) ) then
*                     We have good data on this site/satellite
*                     Scan all other ocmbinations
                      ndd = 0
                      do c2 = 1,actual_max_chan
                         if( c2.ne.ch .and.
     .                       data_OK(data_flag_cse(c2,i,ep),0, 
     .                       phs_mask) ) then
*                            Get the list number for this second SV
                             j2 = ctol_cse(c2,i,ep)
****                         see if any other station has this combination
*                            of two satellites.  If it does count DD
                             k = 0
                             do while( k.lt.num_cfiles )
                                k = k + 1
                                if( k.ne.i ) then  ! Scan
                                  ddfnd = 0
                                  do l = 1,actual_max_chan
*                                    Get SV number for this channel at this site
                                     p2 = ctol_cse(l,k,ep)
                                     if( (p2.eq.j .or.p2.eq.j2) .and.
     .                                  data_OK(data_flag_cse(l,k,ep),0, 
     .                                  phs_mask) ) ddfnd = ddfnd + 1
                                  end do
                                  if( ddfnd.eq.2 ) then
                                      ndd = ndd + 1
                                      k = num_cfiles  ! Found a DD so quit
                                  end if
                                end if
                             end do
                         endif
                      end do
*                     Scanned over all data.  See if we found any double
*                     differences
                      if( ndd.eq.0 ) then   ! No-double differences
*                        Flag data and push any biases
                         if( good_bf(data_flag_cse(ch,i,ep),0, 
     .                               phs_mask) ) then
*                            Push bias flag to next good epoch 
                              call push_ddbf( ep, i, j, data_flag_cse, 
     .                                        bf_type_cse,  ctol_cse)
                         endif
                         call sbit(data_flag_cse(ch,i,ep),23,1)
                         num_dded = num_dded + 1
                      end if
                  end if
               end if
            end do
            if( num_dded.gt.0 .and. 1.eq.2 ) then
               write(*,240) num_dded, cf_codes(i),prn_list(j)
  240          format('NODD Removed ',i5,' epoch on ',a4,' PRN',I2.2)
            end if

         end do
      end do
*     Now re-count number of data that we have (if not-double difference 
*     then data not counted

      do i = 1, num_cfiles
         ntot = 0
         do j = 1, num_sat
            nsvs = 0
            do ep = 1, num_ep
*              OK: See if have data on this satellite at the epoch
               ch = ltoc(ctol_cse(1,i,ep),j, actual_max_chan)
               if( ch.gt.0 ) then   ! We have data on this satellite
                  if( data_OK(data_flag_cse(ch,i,ep),0, phs_mask) ) then
*                     We have good data on this site/satellite
                      nsvs = nsvs + 1
                  end if
               endif
            end do

*           See if we have a case with just 1 double difference.  If so remove the 
*           data point.
            if( nsvs.eq. 1 ) then
                write(oneddwarn,320) cf_codes(i),prn_list(j)
 320            format('Only 1 double diff on ',a4,' PRN_',i2.2)
                call report_stat('warning','autcln','scan_nodd',' ',
     .                            oneddwarn,0)
                do k = 1, num_ep
                   ch = ltoc(ctol_cse(1,i,k),j, actual_max_chan)
                   if( ch.gt.0 ) then   ! We have data on this satellite
                      if( data_OK(data_flag_cse(ch,i,ep), 
     .                            0, phs_mask) ) then
                         call sbit(data_flag_cse(ch,i,k),23,1)
                         call sbit(data_flag_cse(ch,i,k),31,0)
                      endif 
                   end if
                enddo
                nsvs = 0
             endif

             LC_SVS_num(j,i) = nsvs
             ntot = ntot + nsvs
         enddo
         LC_num(i) = ntot
      end do

*     Write out number of values in estimates.
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
      write(*,340) (prn_list(k),k=1,num_sat)
 340  format('NODD: Number of data by site and satellite',/,
     .       'NODD  Site   All',50I4.2)
 
      do j = 1, num_cfiles
         write(*,350) cf_codes(j), LC_num(j), 
     .               (LC_svs_num(k,j), k=1,num_sat)
 350     format('NODD ',1x,a4,i6,50I4)
      end do


      call report_stat('status','autcln','scan_nodd',' ',
     .                     'Finishing scan',0)

      return
      end

      
 
                  
                     
