CTITLE FLAT_DD
 
      subroutine flat_dd(opass, L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse,L1r_rng_cse, L2r_rng_cse,
     .    ctol_cse, data_flag_cse, bf_type_cse, 
     .    params_cse, par_flag_cse  )

      implicit none
 
*     This routine will use a sequential technique to try and flatten
*     the double differences.  The basic procedure is a list
*     of all bias flags and their duration is compiled.  The longest
*     duration one-way bias flag is choosen as reference.  All biases
*     are fixed relative to reference set.  (A second or higher reference
*     may be needed if their is no connection to some bias flags from
*     the reference set).
 
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
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES

*  i -- Loop counter
*  ns, lv, js, kv  -- Site and satellite numbers used in double diffs.
*  epn, dur  -- Epoch and duration of bias being resolved
*  entk  -- Entry for the second satellite at this station (in bf_table)
*  num_resolved  -- Number of biases resolved in current pass
*  iter  -- Iteration counter through loop 

      integer*4 i, ns, lv, js, kv, epn, dur, entk, num_resolved,
     .          iter 

                          
*  all_resolved -- Logical to indicate that all baises have been resolved
*      or that there are too many bias flags.
*  svs_OK, dd_OK  -- Satellite and double diff OK logicals
*  svs_resolved   -- Set true when we find a satellite at the current
*     station that can be used to determine double diffs.

      logical all_resolved, svs_OK, dd_OK, svs_resolved

****  Check to see the type of data we have (ie. dual frequency or L1 only)
      do i = 1, num_cfiles
         if( lambda(1,2,i).eq.0 .and. .not.nol1only ) fdd_L2_fact = 0
      end do
      if( fdd_L2_fact.eq.0 ) then
          call report_stat('status','autcln','flat_dd',' ',
     .                 'L1 only data detected, L1 clocks generated',0) 
      end if

* MOD TAH 031221: Reset all of the cycle offsets to integer values
      call reset_cyc( L1_cyc_cse,L2_cyc_cse )       


***** Start: scan the complete data set and get the list of
*     all biases and their duration.  Once this is done, check to
*     see which is the longest.
      call get_biases( data_flag_cse, ctol_cse, all_resolved, 'REP' )

*     Now scan to get longest sequence.  All resolved returns true
*     if there are no unresolved bias flags.
      call scan_longest( all_resolved )

***** This is first case, so we can mark all bias flags at this
*     site as fixed.  We do this for the bias flags which overlap
*     the reference site bias flag.  Also mark this satellite at
*     all other stations as fixed.
      if( .not. all_resolved ) call set_ref_bfs

****  Now start the loop, until we have resolved all bias flags
      iter = 0

      do while ( .not. all_resolved )

*        Work down the list of bias flags.  If the duration is negative
*        then bias is resolved 
         num_resolved = 0
         iter = iter + 1

         do i = 1, tot_bf

*           See if resolved
            if( bf_table(4,i).gt.0 ) then      

*               This bias flag is not resolved.   See if we can resolve
*               it by forming double differences
                epn = bf_table(1,i)
                ns  = bf_table(2,i)
                lv  = bf_table(3,i)
                dur = bf_table(4,i)

*               See if we can find a satellite at this station which is resolved
                svs_resolved = .false.
                kv = 0
                do while ( .not.svs_resolved .and. kv.le. num_sat)

                    call fdd_check_svs( i, ns, lv, kv, entk,
     .                   data_flag_cse, ctol_cse, svs_OK)

*                   If we found a satellite, see if we can find another
*                   station with this pair of satellites.
                    if( svs_OK ) then

*                       Check the other stations to see if match can
*                       be found and the mean bias set to zero.
                        call fdd_check_dd( i, ns, lv, js, kv, 
     .                        entk,  L1r_phs_cse, L2r_phs_cse,
     .                        L1_cyc_cse, L2_cyc_cse,
     .                        L1r_rng_cse, L2r_rng_cse, 
     .                        ctol_cse, data_flag_cse,  
     .                        params_cse, par_flag_cse, dd_OK )


*                       OK, Double difference is good.  Now find the
*                       mean bias and remove
                        if( dd_OK ) then
                            call fdd_rm_bias( i, ns, lv, 
     .                           L1r_phs_cse, L2r_phs_cse,
     .                           L1_cyc_cse, L2_cyc_cse,
     .                           ctol_cse, data_flag_cse )

*                           Mark this bias as being resolved
                            bf_table(4,i) = -bf_table(4,i)
                            svs_resolved  = .true.
                            num_resolved  = num_resolved + 1
                         end if
                    end if
                 end do
*            end test if bias needs resolving
             end if
         end do

*        See if we still baises that need resolving
         call fdd_check_resolved( all_resolved )

*        See if all are not yet resolved and yet we did not fix any
*        new ones
         if( .not. all_resolved .and. num_resolved.eq.0 ) then
*            We need to set another reference site/sv combination
             call scan_longest( all_resolved )

*            Mark this series as set and continue loop
             bf_table(4,ref_ent) = -bf_table(4,ref_ent)
         end if

*        Report current status
         write(*,310) iter, num_resolved, cf_codes(ref_site),
     .                prn_list(ref_svs), ref_ent
 310     format('FDD: Iter ',i5,': ',i4,' Biases resolved ',
     .          a4,' PRN',i2.2,' Entry ',i4,' Set Ref')

      end do

****  Thats all; all double differences should now be flat
      return
      end


CTITLE GET_BIASES

      subroutine get_biases( data_flag_cse, ctol_cse, all_resolved, rep)

      implicit none

*     Routine to scan all bias flags and make a list of the entries.
*     The all_resolved logical is passed in case there are too many bias
*     flags.  By setting this true, the rest of sequence will exit.

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
     .          ctol_cse(num_chan, num_cfiles, num_ep)

*   all_resolved -- Logical to indicate that all baises have been resolved
*      or that there are too many bias flags.

      logical all_resolved 

*   rep  -- String passed as 'REP' to report the flags
      character*(*) rep

* LOCAL VARIABLES:

*   i,j,k       - Loop counters
*   last_good_bf  -- Epoch number of last good bias flag
*   last_good_ep  -- Epoch number of last good epoch of data
*   ch          -- Channel number of current PRN
*   ltoc        - function to return the channel number
*               - for a particular observation.  If the
*               - satellite is not being oberved then -1
*               - is returned 

      integer*4 i,j,k, last_good_bf, last_good_ep, ch, ltoc

*  good_bf  -- Logical function returns is valid bias flag on this
*     observation
*  data_OK  -- Logical returns true if data point OK

      logical good_bf, data_OK 

***** Loop over all the data, getting all the bias flags
      tot_bf = 0
      all_resolved = .false.

      call set_phs_mask( phs_mask, phs_bias_mask)

      do i = 1, num_cfiles
         bf_index(i) = tot_bf + 1
         do j = 1, num_sat
            last_good_bf = 0
            last_good_ep = 0
            do k = 1, num_ep

*              Check if this is bias flag.  Get channel number for
*              this satellite if it is observed at this time.
               ch = ltoc( ctol_cse(1,i,k),j,actual_max_chan) 
               if( ch.gt.0 ) then

*                  See if this is a good bias flag
                   if( good_bf(data_flag_cse(ch,i,k),0,phs_mask) ) then

*                      Save the duration of the last bias (if there
*                      was one already)
                       if( last_good_bf.gt.0 ) then
                           bf_table(4,tot_bf) = last_good_ep -
     .                                          last_good_bf
                       endif
*                      Save information about this bias flag
                       tot_bf = tot_bf + 1
                       if( tot_bf.gt.max_bf_table ) then
                           call report_stat('warning','autcln',
     .                          'get_biases',' ',
     .                          'Too many bias flags',0)
                           all_resolved = .true.
                           RETURN
                       endif
                       bf_table(1,tot_bf) = k
                       bf_table(2,tot_bf) = i
                       bf_table(3,tot_bf) = j
                       bf_table(5,tot_bf) = 0
                       last_good_bf = k
                       last_good_ep = k
                   end if

*                  If data point is OK, then save last_good_ep
                   if( data_OK(data_flag_cse(ch,i,k),0,phs_mask) ) then
                       last_good_ep = k
                   end if
               end if
            end do

*           Save the information about the last bias flag
            if( last_good_bf.gt.0 ) then
                bf_table(1,tot_bf) = last_good_bf
                bf_table(2,tot_bf) = i
                bf_table(3,tot_bf) = j
                bf_table(4,tot_bf) = last_good_ep-last_good_bf
                bf_table(5,tot_bf) = 0
            end if
         end do
      end do

****  Tell user what we found
      write(*,310) tot_bf, (cf_codes(i), bf_index(i),i=1, num_cfiles)
 310  format('Total of ',i4,' bias flags, By site starts are: ',/,
     .       50(8(a4,1x,I4,1x),:/))

      if( rep.eq.'REP' ) then     
         do i = 1, tot_bf
c           write(*,320) i,(bf_table(j,i), j=1,4)
c320        format('BF_TABLE ',i4,' ENTS ',4i5)
            write(*,320) i, bf_table(1,i), cf_codes(bf_table(2,i)),
     .            prn_list(bf_table(3,i)),bf_table(4,i),
     .            bf_table(1,i)+bf_table(4,i)
 320        format('BF_TABLE ',i4,' Start EP ',i5,1x,a4,' PRN',i2.2,
     .             ' NumEPs ',i4,' EndEp ',I5) 
         end do
      endif


****  Thats all
      return
      end

CTITLE SCAN_LONGEST

      subroutine scan_longest( all_resolved )

      implicit none

*     Routine to check the bf_table for unresolved bias flags
*     and select the reference site and SVS which has the largest
*     amount of data.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
* all_resolved  -- Set true if all the biases have been resolved

      logical all_resolved

* LOCAL VARIABLES
* i  -- Loop counter
* max_dur  -- Longest duration bias flag
* num_unresolved -- Number of unresolved biases left
* max_overlap    -- Longest over lap between a non-resolved
*     and resolved bias.
* ent_max_dur    -- Entry with maxium duration
* ent_max_ovr    -- Entry with maximum overlap (if exists is used in
*                   preference to max_dur)
* ent_max_loc    -- Entry for the max overlap partner.
* loc_max_ovr, loc_max_j -- Local max overlap and entry while scanning
*     an individual station.
* over_lap       -- Amount of overlap between resolved and unresolved 
*                   bias flags.

      integer*4 i,j, max_dur, num_unresolved,  max_overlap, ent_max_dur,
     .          ent_max_ovr, loc_max_ovr, loc_max_j, over_lap,
     .          ent_max_loc

* check          -- Logical set true if we have overlapping segments of
*                   data

      logical check

* type   -- Type of max selected (OVER for overlap, DUR for duration)
      character*4 type


****  Scan over the bf_table and see what we have
      all_resolved = .true.
      max_dur = 0
      max_overlap = 0
      ent_max_dur = 0
      ent_max_ovr = 0
      ent_max_loc = 0
      num_unresolved = 0

      do i = 1, tot_bf
         if( bf_table(4,i).gt.0 ) then
*            OK, unresolved bias flag
             all_resolved = .false.
             num_unresolved = num_unresolved + 1
*  MOD TAH 031220: Normalize by the range noise (in cycles)
             if( bf_table(4,i)/rng_noise(bf_table(2,i)).gt. 
     .           max_dur ) then
                 max_dur   = bf_table(4,i)/rng_noise(bf_table(2,i))
                 ent_max_dur = i
             end if

*            Now check the overlap between this entry and
*            all the resolved entries.  Make sure that either
*            satellite or station overlaps.
             loc_max_ovr = 0
             do j = 1, tot_bf
                if( bf_table(4,j).lt.0 .and. 
     .             (bf_table(2,j).eq.bf_table(2,i) .or. 
     .              bf_table(3,j).eq.bf_table(3,i) )   ) then

*                   Make sure we overlap:
                    check = .false.
                    if( bf_table(1,j).le.bf_table(1,i) .and. 
     .                  bf_table(1,j)-bf_table(4,j).gt.
     .                      bf_table(1,i) )  check = .true.
                    if( bf_table(1,j).gt.bf_table(1,i) .and.
     .                  bf_table(1,i)+bf_table(4,i).gt.
     .                      bf_table(1,j) ) check = .true.
                    if( check ) then
                        over_lap = min(bf_table(1,i)+bf_table(4,i),
     .                                 bf_table(1,j)-bf_table(4,j)) -
     .                             max(bf_table(1,i),bf_table(1,j))
                        if( over_lap/rng_noise(bf_table(2,i)).gt.
     .                      loc_max_ovr ) then
                            loc_max_ovr = over_lap/
     .                                    rng_noise(bf_table(2,i))
                            loc_max_j   = j
                        end if
                    end if
                end if
             end do

****         See if this is the largest
             if( loc_max_ovr.gt.max_overlap ) then
                 max_overlap = loc_max_ovr
                 ent_max_loc = loc_max_j
                 ent_max_ovr = i         
             end if
          end if
      end do

****  Now see what we have:
      if( max_overlap.gt.0 ) then
          ref_ep    = bf_table(1,ent_max_ovr)
          ref_site  = bf_table(2,ent_max_ovr) 
          ref_svs   = bf_table(3,ent_max_ovr)
          ref_ent   = ent_max_ovr
          type      = 'OVER'
      else
          ref_ep    = bf_table(1,ent_max_dur)
          ref_site  = bf_table(2,ent_max_dur) 
          ref_svs   = bf_table(3,ent_max_dur)
          ref_ent   = ent_max_dur
          type      = 'DUR'
      end if

***** Thats all
      write(*,310) tot_bf, num_unresolved, ref_ent,
     .             cf_codes(ref_site),  prn_list(ref_svs), ref_ep,
     .             bf_table(4,ref_ent), type, ent_max_loc
 310  format('SCAN_LONGEST: ',i4,' total bf, ',i4,' unresolved ',
     .       i4,' set as reference ',
     .       a4,' PRN ',i2.2,' Epoch ',i4,' Dur ',i4,' is longest ',a,
     .          ' with ',i5)

      return
      end 
 
CTITLE SET_REF_BFS

      subroutine set_ref_bfs

      implicit none

*     This rouitine will mark the bias flags on the reference site
*     and satellite at the reference epoch as being resolved.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* LOCAL VARIABLES

* i,j,k   - Loop variables
* last_sv, last_site, last_ent -- Last satellite, site and entry in
*        bf_table marked.  Used to make sure we don't mark more
*        than one satellite or site relative to the reference site
* bgn, enn  -- Beginning and ending epochs for the overlap region of
*        a bias flag for the new bias being checked
* bgo, eno  -- Same as above but for the old bias flags already marked.

      integer*4 i,j, last_sv, last_site, last_ent, bgn, enn,
     .          bgo, eno 

* mark  -- Logical set true if we should mark a bias flag as resolved
      logical mark 



****  Loop over all the bias flags and see which can be marked as
*     resolved
      last_sv = 0
      last_site = 0
      
      do i = 1, tot_bf

*        See if this is reference site or satellite
         if( bf_table(2,i).eq. ref_site ) then

*            See if reference satellite and epoch.  If it is mark as
*            resolved
             if( i.eq.ref_ent ) then
                 bf_table(4,i) = -abs(bf_table(4,i))
             else

*                See of epoch and duration match for overlap.
*                The 6 below ensure that we have at least
*                6 overlapping epochs in double difference.
*                (Without testing data_flag we can't be sure
*                that there is data.  If this is a problem then
*                later we may need to test to be sure).
                 mark = .false.
                 if( bf_table(1,i).le.ref_ep .and. 
     .               bf_table(1,i)+abs(bf_table(4,i)).gt.
     .                   ref_ep+6 )  mark = .true.
                 if( bf_table(1,i).gt.ref_ep .and.
     .               ref_ep+abs(bf_table(4,ref_ent)).gt.
     .                   bf_table(1,i)+6 ) mark = .true.

                 if( mark ) then
*                   OK, now we need to check if we already
*                   marked as resolved this satellite.   If
*                   have then, use the one with longest
*                   overlap
                    if( bf_table(3,i).ne.last_sv ) then
*                       We are OK, this is a different satellite
                        bf_table(4,i) = -abs(bf_table(4,i))
                        last_sv = bf_table(3,i)
                        last_ent = i
C                       write(*,300) 'Mark', i, (bf_table(k,i),k=1,4)
C300                    format(a,' Set of Ent, SITE, SVS, DUR',6i5)
                    else
*                       We have already marked this satellite
*                       as resolved.  See if the current segment
*                       is a better choice.
                        bgn = max(bf_table(1,i),ref_ep)
                        enn = min(bf_table(1,i)+abs(bf_table(4,i)),
     .                            ref_ep+abs(bf_table(4,ref_ent)) )
                        j = last_ent
                        bgo = max(bf_table(1,j),ref_ep)
                        eno = min(bf_table(1,j)+abs(bf_table(4,j)),
     .                            ref_ep+abs(bf_table(4,ref_ent)))

*                       See if this segment is better
                        if( enn-bgn.gt.eno-bgo ) then
*                           Yes, it is better.  Un-mark the old one
*                           mark the new one
C                       write(*,300) 'UnMark', j, (bf_table(k,j),k=1,4)
C                       write(*,300) 'Mark  ', i, (bf_table(k,i),k=1,4)
                            bf_table(4,j) = abs(bf_table(4,j))
                            bf_table(4,i) = -abs(bf_table(4,i))
                            last_sv = bf_table(3,i)
                            last_ent = i
                        end if
                    end if
                 end if
             end if
         else if( bf_table(3,i).eq.ref_svs ) then

*            This is another site but for the correct satellite.
*            See if epochs overlap and if we have alreay done this
*            site
             mark = .false.
             if( bf_table(1,i).le.ref_ep .and. 
     .           bf_table(1,i)+abs(bf_table(4,i)).gt.
     .               ref_ep+6 )  mark = .true.
             if( bf_table(1,i).gt.ref_ep .and.
     .           ref_ep+abs(bf_table(4,ref_ent)).gt.
     .               bf_table(1,i)+6 ) mark = .true.
             if ( mark ) then

*                Check to see if already have done this site
                 if( last_site.ne. bf_table(2,i)) then
                     bf_table(4,i) = -abs(bf_table(4,i))
                     last_ent = i
                     last_site = bf_table(2,i)

C                    write(*,300) 'Mark', i,(bf_table(k,i),k=1,4) 
                 else
*                    See if this segment is a better choice for
*                    the reference. 
                     bgn = max(bf_table(1,i),ref_ep)
                     enn = min(bf_table(1,i)+abs(bf_table(4,i)),
     .                        ref_ep+abs(bf_table(4,ref_ent)) )
                     j = last_ent
                     bgo = max(bf_table(1,j),ref_ep)
                     eno = min(bf_table(1,j)+abs(bf_table(4,j)),
     .                        ref_ep+abs(bf_table(4,ref_ent)))

*                    See if this segment is better
                     if( enn-bgn.gt.eno-bgo ) then
*                        Yes, it is better.  Un-mark the old one
*                        mark the new one
                         bf_table(4,j) = abs(bf_table(4,j))
                         bf_table(4,j) =  abs(bf_table(4,j))
                         bf_table(4,i) = -abs(bf_table(4,i))
C                       write(*,300) 'UnMark', j, (bf_table(k,j),k=1,4)
C                       write(*,300) 'Mark  ', i, (bf_table(k,i),k=1,4)
                         last_site = bf_table(2,i)
                         last_ent = i
                     end if
                 end if
             end if
         end if
      end do

****  Thats all
      return
      end

CTITLE FDD_CHECK_SVS

      subroutine fdd_check_svs( entn, ns, lv, kv, entk,
     .           data_flag_cse, ctol_cse, svs_OK)

      implicit none

*     This routine will check to see if we can find another
*     satellites whose bias is resolved to single difference
*     with the one we have.

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
     .          ctol_cse(num_chan, num_cfiles, num_ep)

*   entn   -- Entry number for the bias flag we are trying to resolve
*   ns, lv -- Site number and satellite number trying to be resolved
*   kv     -- Candidate satellite to difference with
*   entk   -- Entry number of the bias flag for matching satellite

      integer*4 entn, ns, lv, kv, entk

*   svs_OK  -- set true if we find an appropriate satellite

      logical svs_OK

* LOCAL VARIABLES

*  ltoc    -- Function which returns channel number for a statellite
*  c1, c2  -- Test channel numbers for satellites being checked
*  bge, ene -- Epoch numbers for over lap of bias flags
*  num_sd   -- Number of single differences entries. 
*  k        -- Loop counter

      integer*4 ltoc, c1, c2, bge, ene, num_sd, k


*  done    -- Indicates that we have done looking for next entry
*             (RETURN if the station number changes) 
*  test    -- Logical set true if the time ranges overlap
*  data_OK  -- Function returns true is data is good

      logical done, test, data_OK

****  Check up the next group of satellites to see if can find a
*     good differencing satellite.  If this is first try set the
*     entry at the start of the station.
      if ( kv.eq.0 ) entk = bf_index(ns) - 1
     
      svs_OK = .false.
      done   = .false.
      do while ( .not. done )

*        Move to next entry in bf_tables for this site
         entk = entk + 1
*        If this entry is for another station then we have run out
*        entries
         if( bf_table(2,entk).ne.ns ) then
*            We are out of entry
             kv = num_sat+1
             RETURN
         end if

*        See if we can find over lap.  Make sure we don't check our
*        own entry that we are trying to find a match for.  Check
*        only the resolved biases.
         if( entk.ne.entn .and. bf_table(4,entk).lt.0 ) then 

C            write(*,210) ns, lv, kv, entn,  entk
C210         format('FDD_TEST S1/C1 C2 ',3i4,' Entries ',2i6)
*            See if epochs seem to overlap
             test = .false.
             if( bf_table(1,entk).le.bf_table(1,entn) .and.
     .           bf_table(1,entk)+abs(bf_table(4,entk)).gt.
     .               bf_table(1,entn)+6 ) test = .true.
             if( bf_table(1,entk).gt.bf_table(1,entn) .and.
     .           bf_table(1,entn)+abs(bf_table(4,entn)).gt.
     .               bf_table(1,entk)+6 ) test = .true.

*            If we should test, then check to see we really have
*            overlapping data.
             if( test ) then
                 bge = max(bf_table(1,entk),bf_table(1,entn))
                 ene = min(bf_table(1,entk)+abs(bf_table(4,entk)),
     .                     bf_table(1,entn)+abs(bf_table(4,entn)))
                 num_sd = 0
                 kv = bf_table(3,entk) 

*                Now loop over the overlap region getting the
*                count of the number of good data.
                 do k = bge, ene
                    c1 = ltoc(ctol_cse(1,ns,k),lv,actual_max_chan)
                    c2 = ltoc(ctol_cse(1,ns,k),kv,actual_max_chan)
*                   See if both measurements available at this epoch
                    if( c1.gt.0 .and. c2.gt. 0 ) then
                        if( data_OK(data_flag_cse(c1,ns,k),0,phs_mask)
     .                                      .and.
     .                      data_OK(data_flag_cse(c2,ns,k),0,phs_mask))
     .                                       then
                            num_sd = num_sd + 1 
                         end if
                    end if
                 end do
*                Finally check that we have enough over lapping data
                 if( num_sd.gt.6 ) then
                    svs_OK = .true.
                    done   = .true.
C                   write(*,310) entk, kv, ns, lv, bf_table(1,entn),
C    .                           bf_table(1,entk), num_sd
C310                format('CHECK_SVS: Found ',i4,' SVS ',i2.2,
C    .                     ' For ',i3,' SVS ',i3,' Ep ',2i6,
C    .                     ' Total overlap ',i4)
                 end if
             end if
         end if
      end do

****  Thats all
      return 
      end

CTITLE FDD_CHECK_DD

      subroutine fdd_check_dd( entn, ns, lv, js, kv, entk,
     .       L1r_phs_cse, L2r_phs_cse,
     .       L1_cyc_cse, L2_cyc_cse,L1r_rng_cse, L2r_rng_cse,
     .       ctol_cse, data_flag_cse,  
     .       params_cse, par_flag_cse, dd_OK  )     

      implicit none
 
*     Routine to find a station that can be double differenced
*     with the current single difference that has been form.

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
     .          ctol_cse(num_chan, num_cfiles, num_ep), 
     .          par_flag_cse(num_param, num_ep)
*   entn   -- Entry number for the bias flag we are trying to resolve
*   ns, lv -- Site number and satellite number trying to be resolved
*   js, kv -- Candidate site and satellite to difference with
*   entk   -- Entry number of the bias flag for matching site

      integer*4 entn, ns, lv, js, kv, entk

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
     .    params_cse(num_param, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   dd_OK  -- set true if we find an appropriate site

      logical dd_OK

* LOCAL VARIABLES

*  j, k  -- Loop counters over the bf_table entries
*  ends  -- Last bf_table entry for the second site being considered  

      integer*4 j,k, ends

*  freq_OK  -- Logical to make sure we only match dual frequency
*     with each other
*  test     -- Logical set true if we should continue 
 
      logical freq_OK, test


***** Start scanning up the bf_table to see what we have available
      fdd_max_dd = 0
      dd_OK = .false.
  
      do j = 1, tot_bf

****     Make sure that site does not overlap and that this is a
*        resolved bias flag.  Also check to see that if this is
*        dual frequency then we only compare with dual frequnecy 
         freq_OK = .true.
*        Removed test: With mixed L1 only, LC data we will use only
*        L1 data.
C        if( fdd_L2_fact.ne.0 .and. 
C    .       lambda(bf_table(3,j),2,bf_table(2,j)).eq.0 ) 
C    .                                         freq_OK = .false.
         if( bf_table(2,j).ne.ns .and. bf_table(4,j).lt.0 .and.
     .       freq_OK ) then

*            See if it matches one of the satellites that we need and
*            that the time range overlaps
             test = .false.
             if( bf_table(1,j).le.bf_table(1,entn) .and.
     .           bf_table(1,j)+abs(bf_table(4,j)).gt.
     .               bf_table(1,entn)+6 ) test = .true.
             if( bf_table(1,j).gt.bf_table(1,entn) .and.
     .           bf_table(1,entn)+abs(bf_table(4,entn)).gt.
     .               bf_table(1,j)+6 ) test = .true.

             if( (bf_table(3,j).eq.lv .or. bf_table(3,j).eq.kv) .and.
     .           test        ) then
*                Matches one of the satellites.  Now loop up through
*                the rest of the table for this station finding the 
*                other satellite
                 if( bf_table(2,j).lt.num_cfiles ) then
                      ends = bf_index(bf_table(2,j)+1) - 1
                 else
                      ends = tot_bf
                 end if
                 do k = j+1, ends

*                   Make sure that bias is resolved
                    if( bf_table(4,k).lt.0 ) then
*                       See if it is the other satellite and then check
*                       time range
                        js = bf_table(2,j)
                        if( bf_table(3,j).ne.bf_table(3,k) .and.
     .                      (bf_table(3,k).eq.lv .or. 
     .                       bf_table(3,k).eq.kv     )     ) then 
                            call fdd_form_dd( entn, entk, j, k, 
     .                           ns, lv, js, kv, 
     .                           L1r_phs_cse, L2r_phs_cse,
     .                           L1_cyc_cse, L2_cyc_cse,
     .                           L1r_rng_cse, L2r_rng_cse,
     .                           ctol_cse, data_flag_cse,  
     .                           params_cse, par_flag_cse  )
                        end if
                    end if
                 end do
             end if
         end if
      end do

****  See if managed to find double differences
      if( fdd_max_dd.gt.0 ) then
          dd_OK = .true.
      end if

****  Thats all
      return
      end
         
CTITLE FDD_FORM_DD

      subroutine fdd_form_dd( entn, entk, j, k, ns, lv, js, kv, 
     .      L1r_phs_cse, L2r_phs_cse,
     .      L1_cyc_cse, L2_cyc_cse,L1r_rng_cse, L2r_rng_cse,
     .      ctol_cse, data_flag_cse,  
     .      params_cse, par_flag_cse  )

      implicit none

*     Routine to form the double differences for the selected
*     station satellite pair and to compute the mean and rms
*     If this combination is the largest number of double differences
*     then this offset is saved.


* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*   entn  -- bf_table entry for one-way being resolved
*   entk  -- bf_table entry for the second satellite at the main site.
*   j,k   -- bf_table entries for second station and pair of
*            satellites being check.
*   ns, lv, js, kv -- Site and satellites being used for one-way and
*            the matching double difference

      integer*4 entn, entk, j,k, ns, lv, js, kv  

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
     .    params_cse(num_param, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)


* LOCAL VARIABLES
* bge, ene  -- Epoch range to form double differences.  These
*     are the latest start time for a bias flag, and the
*     eariliest end time.
* num_dd    -- Number of double differences found
* l         -- Loop counter
* c11, c12, c21, c22 -- Channel numbers for DD entries
* ltoc      -- Functio to return channel number

      integer*4 bge, ene, num_dd, l,  c11, c12, c21, c22, ltoc


* sum_dd(2)    -- Sum of double differences for LC and MW-WL
* ssm_dd(2)    -- Sum of double differences squared.

* r11, r12, r21, r22 -- Phase residuals for one-ways.
* mw11, mw12, mw21, mw22 -- MW-WL for each widelane
* dd_res    -- Double difference residual
* rms_dd, mean_dd -- Current calcs for rms and mean of double 
*     differences for LC and MW-WL
* phs_omc   -- Phase residual funnction

      real*8 sum_dd(2), ssm_dd(2),  r11, r12, r21, r22,
     .        mw11, mw12, mw21, mw22, 
     .       dd_res(2), rms_dd(2), mean_dd(2), phs_omc

*  test  -- Used to see if should continue
*  d11, d12, d21, d22 -- Set true if data good for each one-way
*  data_OK  -- Function returns true if data is OK

      logical test, d11, d12, d21, d22, data_OK, kbit, ft

****  Start, the final formation of the double differences.  We
*     have second entry with correct satellite, now finally check
*     that have overlap in time.  If this is the case compute
*     mean double difference and make sure we have data.

*     OK, we have found the correct satellite
*     pair.  Check to see if times overlap.
      test = .false.
      if( bf_table(1,k).le.bf_table(1,entn) .and.
     .    bf_table(1,k)+abs(bf_table(4,k)).gt.
     .                bf_table(1,entn)+6 ) test = .true.
      if( bf_table(1,k).gt.bf_table(1,entn) .and.
     .    bf_table(1,entn)+bf_table(4,entn) .gt.
     .                bf_table(1,k)+6 ) test = .true.

*     See if we should test
      if( test ) then
*         Yes, no see how many data point actually
*         over lap in the double differences
          bge = max(bf_table(1,entn), bf_table(1,entk),
     .              bf_table(1,j),  bf_table(1,k) )
          ene = min(bf_table(1,entn) +  bf_table(4,entn),
     .              bf_table(1,entk) +  abs(bf_table(4,entk)),
     .              bf_table(1,j)    +  abs(bf_table(4,j)),
     .              bf_table(1,k)    +  abs(bf_table(4,k)) )

****      Initialize
          num_dd = 0
          sum_dd(1) = 0
          ssm_dd(1) = 0
          sum_dd(2) = 0
          ssm_dd(2) = 0

*         Loop over the epoch range given
          ft = .false.
          do l = bge, ene

*            Check to see if all the data is present
*            at this epoch
             c11 = ltoc(ctol_cse(1,ns,l),lv,actual_max_chan)
             c12 = ltoc(ctol_cse(1,ns,l),kv,actual_max_chan)
             c21 = ltoc(ctol_cse(1,js,l),lv,actual_max_chan) 
             c22 = ltoc(ctol_cse(1,js,l),kv,actual_max_chan)
*            See if all present
             if( c11.gt.0 .and. c12.gt.0 .and. c21.gt.0 .and.
     .           c22.gt.0 ) then

*                Make sure all data is good
                 d11 = data_OK(data_flag_cse(c11,ns,l),0,phs_mask) 
                 d12 = data_OK(data_flag_cse(c12,ns,l),0,phs_mask)                                                 
                 d21 = data_OK(data_flag_cse(c21,js,l),0,phs_mask)
                 d22 = data_OK(data_flag_cse(c22,js,l),0,phs_mask)
* MOD TAH 050712: Make sure we L2 range data 
                 if( lambda(c11,4,ns).eq.0 ) d11 = .false.
                 if( lambda(c12,4,ns).eq.0 ) d12 = .false.
                 if( lambda(c21,4,js).eq.0 ) d21 = .false.
                 if( lambda(c22,4,js).eq.0 ) d22 = .false.

                 if( d11 .and. d12 .and. d21 .and. d22 ) then

*                    Now compute the phase residuals so that we can 
*                    form the double difference.  Check if dual or
*                    single frequnecy needed.
*                    phs_omc now detects the fdd_L2_fact and so single
*                    call is needed.
                     r11 = phs_omc(L1r_phs_cse(c11,ns,l),
     .                     L2r_phs_cse(c11,ns,l),     
     .                     L1_cyc_cse(c11,ns,l), L2_cyc_cse(c11,ns,l),     
     .                     apr_clk_val(ns),apr_clk_val(num_cfiles+lv),
     .                     fL1(lv), fL2(lv) )
                     r12 = phs_omc(L1r_phs_cse(c12,ns,l),
     .                     L2r_phs_cse(c12,ns,l),     
     .                     L1_cyc_cse(c12,ns,l), L2_cyc_cse(c12,ns,l),     
     .                     apr_clk_val(ns),apr_clk_val(num_cfiles+kv),
     .                     fL1(kv), fL2(kv) ) 
                     r21 = phs_omc(L1r_phs_cse(c21,js,l),
     .                     L2r_phs_cse(c21,js,l),     
     .                     L1_cyc_cse(c21,js,l), L2_cyc_cse(c21,js,l),     
     .                     apr_clk_val(js),apr_clk_val(num_cfiles+lv),
     .                     fL1(lv), fL2(lv) ) 
                     r22 = phs_omc(L1r_phs_cse(c22,js,l),
     .                     L2r_phs_cse(c22,js,l),     
     .                     L1_cyc_cse(c22,js,l), L2_cyc_cse(c22,js,l),     
     .                     apr_clk_val(js),apr_clk_val(num_cfiles+kv),
     .                     fL1(kv), fL2(lv) )
*                    Now compute the MWWL for each one
C                    mwwl = ow(2)-ow(1)+dfsf(lv)*(ow(3)+ow(4))
                     mw11 = (L2r_phs_cse(c11,ns,l)+L2_cyc_cse(c11,ns,l))
     .                     -(L1r_phs_cse(c11,ns,l)+L1_cyc_cse(c11,ns,l))
     .                     +dfsf(lv)*( L1r_rng_cse(c11,ns,l)+
     .                             L2r_rng_cse(c11,ns,l) )
                     mw12 = (L2r_phs_cse(c12,ns,l)+L2_cyc_cse(c12,ns,l))
     .                     -(L1r_phs_cse(c12,ns,l)+L1_cyc_cse(c12,ns,l))
     .                     +dfsf(kv)*( L1r_rng_cse(c12,ns,l)+
     .                             L2r_rng_cse(c12,ns,l) )
                     mw21 = (L2r_phs_cse(c21,js,l)+L2_cyc_cse(c21,js,l))
     .                     -(L1r_phs_cse(c21,js,l)+L1_cyc_cse(c21,js,l))
     .                     +dfsf(lv)*( L1r_rng_cse(c21,js,l)+
     .                             L2r_rng_cse(c21,js,l) )
                     mw22 = (L2r_phs_cse(c22,js,l)+L2_cyc_cse(c22,js,l))
     .                     -(L1r_phs_cse(c22,js,l)+L1_cyc_cse(c22,js,l))
     .                     +dfsf(kv)*( L1r_rng_cse(c22,js,l)+
     .                             L2r_rng_cse(c22,js,l) )


*                    Now form dd_res
                     dd_res(1) = (r11-r12)-(r21-r22)
                     dd_res(2) = (mw11-mw12) - (mw21-mw22)
                     num_dd = num_dd + 1
                     sum_dd(1) = sum_dd(1) + dd_res(1)
                     ssm_dd(1) = ssm_dd(1) + dd_res(1)*dd_res(1)
                     sum_dd(2) = sum_dd(2) + dd_res(2)
                     ssm_dd(2) = ssm_dd(2) + dd_res(2)*dd_res(2)
                 end if
             end if
          end do

****      See what we have now that we have loop over everything
*         Increased overlap from 6 to 50 for robustness
          if( num_dd.gt.60 ) then
              mean_dd(1) = sum_dd(1)/num_dd
              rms_dd(1) = sqrt((ssm_dd(1)-num_dd*mean_dd(1)**2)/num_dd)
              mean_dd(2) = sum_dd(2)/num_dd
              rms_dd(2) = sqrt((ssm_dd(2)-num_dd*mean_dd(2)**2)/num_dd)
	      if( kbit(status_rep,14) )
C    .        write(*,310) ns, lv, js, kv, bge, ene, num_dd, mean_dd,
C    .                     rms_dd, fdd_max_dd
C310          format('FDD_DD S1/C1 S2/C2 ',4i3,' EPS ',2i5,' Num ',
C    .               I4,' Mean, RMS ',2F10.2,' Max ',I4)
     .        write(*,315) cf_codes(ns), prn_list(lv), cf_codes(js),
     .                prn_list(kv),bge, ene, num_dd, mean_dd(1),
     .                rms_dd(1),mean_dd(2),rms_dd(2) , fdd_max_dd 
 315          format('FDD_DD ',a4,' PRN ',i2.2,' (',a4,' PRN ',i2.2,
     .               ')  EPS ',2i5,' Num ',
     .               I4,' LC Mean, RMS ',2F10.2,
     .                  ' WL Mean, RMS ',2f10.2,' Max ',I4)
*             See if we have more data than the current max
              if( num_dd.gt. fdd_max_dd ) then
                  fdd_max_dd = num_dd
                  fdd_max_mean(1) = mean_dd(1)
                  fdd_max_rms(1)  = rms_dd(1)
                  fdd_max_mean(2) = mean_dd(2)
                  fdd_max_rms(2)  = rms_dd(2)
              end if
          else
C             write(*,310) ns, lv, js, kv, bge, ene, num_dd
          end if
      end if
C     write(*,320) ns, lv, js, kv, entn, entk, j, k
C320  format('FDD_No Match S1/C1 S2/C2 ',4i4,' Entries ',4i6)
      

****  Thats all
      return
      end

CTITLE FDD_RM_BIAS

      subroutine fdd_rm_bias( entn, ns, lv, 
     .                        L1r_phs_cse, L2r_phs_cse,
     .                        L1_cyc_cse, L2_cyc_cse,
     .                        ctol_cse, data_flag_cse ) 

      implicit none

****  This routine will apply the cycle offset to make the 
*     double differences zero mean.


* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES 

*   entn  -- Entry in bf_table being resolved.
*   ns, lv  -- Station and satellite

      integer*4 entn, ns, lv
 
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
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES

* i   -- Loop counter over epochs
* bgn, enn  -- Beginning and ending epochs
* ltoc      -- function that return channel number for satellite number
* ch        -- Channel for current observation

      integer*4 i, bgn, enn, ltoc, ch

* dL1c, dL2c  -- Change in the number of cycles
* dL21   -- Difference in L2-L1 cycles (from MWWL)
* dL1, dL2 -- Changes in dL1 and dL2 (correction is opposite sign)

      real*8 dL1c, dL2c, dL21, dL2, dL1

      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.
      logical kbit
                            
*     Get the channel number
      ch = ltoc(ctol_cse(1,1,i),lv, actual_max_chan)


****  OK, compute the number of cycles that need to be applied
      if( fdd_L2_fact.ne.0 ) then
C          dL1c = -fdd_max_mean/(lcf1(lv)+lcf2(lv))
C          dL2c = dL1c
* MOD TAH 031221: Compute MWWL and LC offsets needed
           fr1 = fL1u(lv)/fL1(lv)
           fr2 = fL2u(lv)/fL2(lv)
           dL21 = nint(fdd_max_mean(2)*fr1)
           dL1 = nint((fdd_max_mean(1)*fr1-
     .                dL21*lcf2(lv))/(lcf1(lv)+lcf2(lv)))
           dL2  = dL1 + dL21
           dL1c = -dL1/fr1
           dL2c = -dL2/fr2
      else
*         Single frequency results
          dL1c = -fdd_max_mean(1)
          if( lambda(lv,2,ns).ne.0 ) then
             dL2c = dL1c
          else
             dL2c = 0 
          end if
      end if

*     Now get the epoch range for cycles to be removed
      bgn = bf_table(1,entn) 
      if( entn.lt.tot_bf ) then
*         See if next entry is same sation and satellite
          if( bf_table(2,entn+1).eq.bf_table(2,entn).and.
     .        bf_table(3,entn+1).eq.bf_table(3,entn) ) then
*             Go to epoch before next bias flag
              enn = bf_table(1,entn+1)-1
          else
*             Either sation or satellite has changed.  Go to
*             end of data
              enn = num_ep
          end if
      else
*         Last bias flag, so go to end of data
          enn = num_ep
      end if

****  Tell user what is happening and adjust the number of cycles
      if( kbit(status_rep,14) )
     .write(*,210) cf_codes(ns), prn_list(lv), entn, bgn, enn,
     .             dL1c, dL2c
 210  format('FDD_RM_BIAS: ',a4,' PRN',i2.2,' Entry ',i4,' EPS ',
     .       2i5,' dL1/dL2 ',2f10.2)
      do i = bgn, enn

*        Get the channel number
         ch = ltoc(ctol_cse(1,ns,i),lv, actual_max_chan)
         if( ch.gt.0 ) then
             L1_cyc_cse(ch,ns,i) = L1_cyc_cse(ch,ns,i) + dL1C
             L2_cyc_cse(ch,ns,i) = L2_cyc_cse(ch,ns,i) + dL2C
         end if
      end do

****  Thats all
      return
      end

CTITLE FDD_CHECK_RESOLVED

      subroutine fdd_check_resolved( all_resolved )

      implicit none

*     Routine to see if all the bias flags are resolved yet

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
* all_resolved  -- Set false if any bias remained unresolved.

      logical all_resolved

* LOCAL VARIABLES
* i  -- Loop variable

      integer*4 i

****  Check the complete list
      all_resolved = .true.
      do i = 1, tot_bf
          if( bf_table(4,i).gt.0 ) all_resolved = .false.
      end do

****  Thats all
      return
      end
















