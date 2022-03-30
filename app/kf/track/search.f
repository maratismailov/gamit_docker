CTITLE SELECT_BEST

      subroutine select_best(ib)

      implicit none

*     Routine to look at the best RMS choices for the ambiquities
*     and select the best based the one that performs best overall.
                               
      include 'track_com.h'

* PASSED VARIABLES
* ib  -- Index to the position of the best choice of ambiquities
      integer*4 ib

* LOCAL VARIABLES
* i,j,k -- Counters
* nr -- Number of entry being checked to see if new amgiquity combination
*       can be found that should be included in the ranking
* cr -- Current number of ranked ambiquity combinations
* na -- Ambiquity number for marking values as fixed
* nm -- Number of matching entries in the ranking.
* saved_epoch(max_amb_save)  -- Epoch number used in saving the ranked
*       list.  Each epoch can only be used once in forming the ranking.

      integer*4 i,j,k,l, nr, cr, na, tot, nm, 
     .          saved_epoch(max_amb_save)

* used(max_amb_save) -- Logical set true once a set of ambiquity has
*       been added to the ranked list
* done   -- Set true once we have all ambiquity sets ranked
* new_cr -- Set true when a new amquitity sequence is found
* match  -- Set true is ambiquity set matches the combination we are
*           trying to rank.

      logical used(max_amb_save), done, new_cr, match, test

      character*1024 line
      integer*4 trimlen

****  Choose the ambiquity combination which has the lowest RMS
*     and rank those with the same combination.  The ranking is
*     done by the sum of the inverse of RMS squared with a lower
*     RMS of rank_min_rms
      rank_min_rms = 0.001

*     Get the first combination
      rank_nu = save_nu(1)
      test = .false.
      if( test ) then 
         do i = 1, min(num_amb_save,3)
            write(line,50) i, save_nu(i), save_ep(i), 
     .                     save_rms(i)*1000.d0,
     .                    (save_na(j,i),j=1,save_nu(i))
 50         format(' Saved ',i4,i3,i6,F8.1,20i5)
            write(line(trimlen(line)+2:),60) (save_amb(1,j,i),
     .          save_amb(2,j,i),j=1,save_nu(i)) 
 60         format(' AMB ', 50I3)
            write(*,'(a)') line(1:trimlen(line))
         end do
      end if

****  Save the ambiquity numbers being ranked              
      do i = 1, save_nu(1)
         rank_na(i) = save_na(i,1)
      end do
     
      nr = 0
      cr = 0
      done = .false.
      do i = 1, num_amb_save
         used(i) = .false.
      end do

*     Loop through the saved ambiquities, finding each combination
*     and ranking them
      do while ( .not. done )

*        Goto the next ambiquity and see if it has already been used
*        in the ranking.
         nr = nr + 1
         tot = 0
*        If we are on the last entry, set done true
         if ( nr.ge. num_amb_save-1 ) done = .true.  

*        See if have ranked this combination before
         new_cr = .false.
         if( .not.used(nr) ) then
*            OK, we have not used the combination before so create
*            new combination to rank.  Make sure it matches the
*            set of ambiquity numbers we are trying to rank
*            See if we have all the desired ambiquities in this saved
*            value
             nm = 0
             do k = 1, save_nu(nr)
                do l = 1, rank_nu
                   if( save_na(k,nr).eq.rank_na(l) ) nm = nm + 1
                end do
             end do
             match = .false.
             if( nm.eq.rank_nu ) match = .true.
             
****         If we have a match then we can creat a new entry to
*            rank
             if( match ) then
                cr = cr + 1
                new_cr = .true.
                comb_ivar(cr) = 0.d0
                comb_ion(cr)  = 0.d0
                do j = 1, save_nu(nr)
*                  Loop over the ambiquities we are trying to resolve
*                  can copy the values to the combination array.
                   do k = 1, rank_nu
                      if( save_na(j,nr).eq.rank_na(k) ) then
                         comb_amb(1,k,cr) = save_amb(1,j,nr)
                         comb_amb(2,k,cr) = save_amb(2,j,nr)
                      end if
                   end do
                end do

*               Clear the list of epochs used in this save
                do j = 1, num_amb_save
                   saved_epoch(j) = 0
                end do

             end if
         end if

****     If we have added to new combination to check, scan the list
*        now and get its rank
         if( new_cr ) then
            do j = 1, num_amb_save
*              See if this entry matches the combination we are ranking
*              Sum up all the matching entries and make sure it matches
*              the number we are trying to rank.
               nm = 0
               do l = 1, save_nu(j)

*                 Has the correct number, see if they match
                  do k = 1, rank_nu
                     if( save_na(l,j).eq.rank_na(k) .and.   ! Ambiquity number matches
     .                   save_amb(1,l,j).eq.comb_amb(1,k,cr) .and.
     .                   save_amb(2,l,j).eq.comb_amb(2,k,cr) ) then
*                        Ok, all the elements match.  Do nothing in this
*                        case.  If not all matched then match set to false
                         nm = nm + 1
                     end if
                  end do
               end do
               match = .false.
               if( nm.eq.rank_nu ) match = .true.
               if( match ) used(j) = .true.

*              Make sure we have not used this epoch before
               do k = 1, tot
                  if( save_ep(j).eq.saved_epoch(k) ) match = .false.
               end do

*              If we have a match, then accumulate 1/rms**2
               if( match ) then
                   comb_ivar(cr) = comb_ivar(cr) + 
     .                      1.d0/(save_rms(j)**2+rank_min_rms**2)
                   comb_ion(cr)  = comb_ion(cr) + 1.d0/save_rmsi(j)
                   tot = tot + 1
                   used(j) = .true.
                   saved_epoch(tot) = save_ep(j)
               end if
            end do

*           Now add to ranking list.  List sorted from largest to smallest
            comb_tot(cr) = tot
            call add_to_rank(cr)
         end if
      end do

****  Write out the ranking 
      num_ranked = cr 
      write(*,520) num_ranked,rank_nu,(rank_na(j),j=1,rank_nu)
 520  format('Ranked ',i4,' entries for ',i3,' ambiquities',/,
     .       '  Occ Sum 1/Var  RelRank   Sum 1/ion  RelRank ',
     .       50(' #',I4,' L1 L2'))    
  
      do i = 1, min(num_ranked,10)
         write(*,540) comb_tot(i), comb_ivar(i), 
     .                comb_ivar(1)/comb_ivar(i),
     .                comb_ion(i), comb_ion(1)/comb_ion(i),
     .                ((comb_amb(k,j,i),k=1,2),j=1,rank_nu)
 540     format(i4, 2(F11.2,1x,F8.2,1x),40I6)
      end do

*     Set the choice as the top ranked entry
      ib = 1

*     Now if constrast good enough set ambiquity as resolved
      if( num_ranked.eq.1 ) comb_ivar(2) = 1
      
      if( comb_ivar(1)/comb_ivar(2).ge. relrank_limit .and. 
     .    num_ranked.gt.0 .and.  comb_tot(1).gt. num_amb_samp/5 .and.
     .    comb_ivar(1).gt.comb_tot(1)*1000 ) then
          do i = 1, rank_nu
             na = rank_na(i)
             call sbit(bf_ents(5,na),2,1)
             num_tot_resolved = num_tot_resolved + 1
          end do
      else
*         See if this is better than any previous ones
          do k = 1, rank_nu
             if( saved_rank(rank_na(k)).lt. 
     .           comb_ivar(1)/comb_ivar(2)  ) then
*                Saved rank is less than current value.  Update
*                saved rank.
                 if( num_ranked.gt.1 ) then
                    saved_rank(rank_na(k)) = comb_ivar(1)/comb_ivar(2)
                 else
                    saved_rank(rank_na(k)) = 1.0
                 end if
*                Update all the trials to reflect these values in 
*                case to do single ambiguity search
                 do i = 1, num_amb_save
                    do j = 1, save_nu(i)
                       if( save_na(j,i).eq.rank_na(k) ) then
                           save_amb(1,j,i) = save_amb(1,j,i) -
     .                           comb_amb(1,k,1)
                           save_amb(2,j,i) = save_amb(2,j,i) -
     .                           comb_amb(2,k,1)
                       end if
                    end do
                 end do

             else
*                This is not as good as previous ones, so reset the
*                ambiquities to zero.
                 comb_amb(1,k,1) = 0
                 comb_amb(2,k,1) = 0
             end if
          end do
      end if

****  Thats all
      return
      end
 

CTITLE ADD_TO_RANK

      subroutine add_to_rank( cr )

      implicit none

*     Routine to merge the current combination of ambiguities into
*     the ranked list.  This involves finding its rank and merging
*     it into the ranking

      include 'track_com.h'

* PASSED VARIABLES
* cr -- The current number of entries in the ranked list.  The
*       last entry is the one we just computed and it needs to
*       be moved to the correct position in the list
      integer*4 cr

* LOCAL VARIABLES
* i,j  -- Loop counters
* ent  -- Entry position of the next highest ranked one above  the 
*         current combination (i.e, current entry is added as ent+1
* copy_amb(2,max_chan)  -- Copy of the current entry which is moved
*         from the end of the list to its ranked position
* copy_tot  -- Copy of total value to allow for swapping entries
      integer*4 i,j, ent, copy_amb(2,max_chan), copy_tot

* copy_ivar -- Copy of the rank of the current entry
* copy_ion  -- Copy of the ion rank for the current entry

      real*4 copy_ivar, copy_ion

****  OK, find where the new entry to should be added into the
*     list
      ent = 0
      do i = 1, cr-1
         if( comb_ivar(cr).lt.comb_ivar(i) ) ent = i
      end do
             
****  See we need to do anything.  If ent = cr-1 then the
*     current entry is at the end of the list so we don't need to
*     anything
      if( ent.eq.cr-1 ) RETURN

****  Now copy the cr entry, and move all the ones of lower rank
*     down in the list and add the new entry
      copy_ivar = comb_ivar(cr)
      copy_ion  = comb_ion(cr)
      copy_tot  = comb_tot(cr)
      do j = 1, rank_nu
         copy_amb(1,j) = comb_amb(1,j,cr)
         copy_amb(2,j) = comb_amb(2,j,cr)
      end do

*     Now move the list down from ent+1 to cr-1
      do i = cr-1,ent+1,-1
         comb_ivar(i+1) = comb_ivar(i)
         comb_ion(i+1)  = comb_ion(i)
         comb_tot(i+1)  = comb_tot(i)
         do j = 1, rank_nu
            comb_amb(1,j,i+1) = comb_amb(1,j,i)
            comb_amb(2,j,i+1) = comb_amb(2,j,i)
         end do
      end do

*     Now add the copy into its ranked position
      comb_ivar(ent+1) = copy_ivar
      comb_ion(ent+1)  = copy_ion
      comb_tot(ent+1)  = copy_tot
      do j = 1, rank_nu
         comb_amb(1,j,ent+1) = copy_amb(1,j)
         comb_amb(2,j,ent+1) = copy_amb(2,j)
      end do   

****  Thats all
      return
      end

CITILE CHECK_TRIAL
 
      subroutine check_trial( ep, sent, sol_upd, data_type, sn, rms )

      implicit none

*     Routine to check the trial values of the ambiquities to see
*     if we get better residuals than those in sol_upd from the
*     nominal set of ambiquities.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep  -- Epoch number
* sent(max_chan) -- Array containing pointer to which
*     entries in the search arrays that we actually doing at the moment
* sn  -- Stattion number 

      integer*4 ep, sent(max_chan), sn

* sol_upd(max_obs)  -- Current set of residuals from the nominal set
*     of ambiquities
* rms  -- RMS scatter of residuals for this trial (999.d0 if bad residuals
*     are found first)

      real*8 sol_upd(max_obs), rms

* data_type  -- String with type of data used

      character*(*) data_type

* LOCAL VARIABLES
* i,j -- Loop counters
* trial_amb(max_obs) -- Trial ambiquity values
* ent  -- Entry in the saved ambiquitity tables for the point before
*     where we should insert the current entry.
* nd   -- Number of dd actually used

      integer*4 i,j, trial_amb(max_obs), trimlen, ent, nd

* dres -- Change in residual
* sum_sqr  -- Sum of residuals squared
* sol_trial(max_obs) -- Final residuals for combination which looks
*      good
* rmsi   -- ION RMS

      real*8 dres, sum_sqr, sol_trial(max_obs), rmsi

* outline  -- String to contain output line

      character*512 outline

***** For the set of ambiquities we are trying, find out the changes in
*     the double difference data used in the nominal solution.
*     Generate the set of L1,L2 ambiquity values
      do i = 1, num_sent
         do j = 1, 2
            trial_amb((i-1)*2+j) = bf_search(j,i,sent(i))
         end do
      end do

****  OK, now start multipling to see how the residuals look
      sum_sqr = 0.d0
      do i = 1, num_dd
         dres = 0.d0
         nd = 0
         if( ow_des(1,dd_des(1,i)).eq.sn ) then
            do j = 1, num_sent*2
               dres = dres + amb_to_res(i,j)*trial_amb(j)
            end do
            if( abs(sol_upd(i)-dres).gt.max_dd_elev(i) ) then
                rms = 999.d0
                RETURN
            else
                nd = nd + 1
                sol_trial(i) = sol_upd(i)-dres
                sum_sqr = sum_sqr + (sol_upd(i)-dres)**2
            end if
         end if
      end do

****  Ok, Each of the residuals was OK.  See what RMS we have
      if( nd.eq.0 ) then
          rms = 999.d0
          RETURN
      end if
      rms = sqrt(sum_sqr/nd)
      write(outline,320) ep, rms, (trial_amb(j),j=1, 2*num_sent)
 320  format('SOK RMS ',i8,1x,F6.3,' Trial ',40i5)
      write(outline(trimlen(outline)+1:),340) 
     .                            (sol_trial(j), j=1,num_dd)
 340  format(' RES ',50F8.3)
      
      if( 1.eq.2 ) 
     .write(*,'(a)') outline(1:trimlen(outline))

****  OK, Now search the ionospheric constraint
      sum_sqr = 0.d0
      do i = 1, num_ddion
         dres = 0.d0
         nd = 0
         if( ow_des(1,dd_des(1,i)).eq.sn ) then
            nd = nd + 1 
            do j = 1, num_sent*2
               dres = dres + amb_to_ion(i,j)*trial_amb(j)
            end do
            sol_trial(i) = (ion_obs(i)-dres)*vel_light/fR1
            sum_sqr = sum_sqr + (sol_trial(i))**2
         end if
      end do

****  Now add in the ion RMS

      rmsi = sqrt(sum_sqr/nd)
      rms = sqrt(rms**2+ion_wght*sum_sqr/nd)

      write(outline,420) ep, rms, rmsi, (trial_amb(j),j=1, 2*num_sent)
 420  format('SOK RMS ',i8,1x,2F6.3,' Trial ',40i5)
      write(outline(trimlen(outline)+1:),340) 
     .                            (sol_trial(j), j=1,num_ddion)
      
      if( 1.eq.2 ) 
     .write(*,'(a)') outline(1:trimlen(outline))

****  Now add this one to list
      ent = 0
      do i = 1, num_amb_save
         if( rms.gt.save_rms(i) ) ent = i
      end do

*     Now push down the values to make room for current entry
       
      do i = num_amb_save, ent+1, -1
         save_rms(i+1) = save_rms(i)
         save_rmsi(i+1) = save_rmsi(i)
         save_ep(i+1)  = save_ep(i)
         save_nu(i+1)  = save_nu(i)
         do j = 1, save_nu(i)
            save_na(j,i+1) = save_na(j,i)
            save_amb(1,j,i+1) = save_amb(1,j,i)
            save_amb(2,j,i+1) = save_amb(2,j,i)
         end do
      end do

****  OK, now add our new entry at ent+1 position
      save_rms(ent+1) = rms
      save_rmsi(ent+1) = rmsi
      save_ep(ent+1)  = ep
      save_nu(ent+1)  = num_sent
      do j = 1, num_sent
         save_na(j,ent+1) = amb_search(j)
         save_amb(1,j,ent+1) = bf_search(1,j,sent(j))
         save_amb(2,j,ent+1) = bf_search(2,j,sent(j))
      end do

****  Now update the number of saved values
      num_amb_save = num_amb_save + 1
      if( num_amb_save.ge.max_amb_save ) num_amb_save = max_amb_save-1
                           
      return
      end

CTITLE FIND_BFS

      subroutine find_bfs( ep, sn )

      implicit none

*     Routine to find the bias flags which are used at this epoch
*     (all bias flags entries are returned.  Only those biases to
*     site ns are returned (ns is the site with the original bias
*     we are trying to find)
 
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* sn   -- Site number whose bias is being searched

      integer*4 ep, sn

* LOCAL VARIABLES
* nb  -- Short name for num_bf (number of bias flags)
* i   -- Loop counter
* ns, pn -- Site and PRN in double difference

      integer*4 nb, i, j,k, ns, pn

* in_dd -- set true if bias flag we are testing is the double differences

      logical in_dd

****  Search accross the table of ambiquities and see which ones
*     apply
      nb = 0
      do i = 1, num_ambs
         if( ep.ge.bf_ents(3,i) .and. ep.le. bf_ents(4,i) .and.
     .       sn.eq.bf_ents(1,i) ) then
*            This epoch in inside time range.  Now check that this
*            bias flag was used in the double differences that we
*            are testing.
             in_dd = .false.
             do j = 1, num_dd
                do k = 1, 4

*                  Get PRN and site of this one-way
                   ns = ow_des(1,dd_des(k,j))
                   pn = ow_des(2,dd_des(k,j))
                   if( bf_ents(1,i).eq.ns .and. bf_ents(2,i).eq.pn )
     .                                               in_dd = .true.
                end do
             end do      
             if( in_dd ) then
                 nb = nb + 1
                 bf_flags(nb) = i 
             end if
         end if
      end do

      num_bf = nb

****  Thats all
      return
      end

CTITLE SELECT_TRIALS

      subroutine select_trials( ep, data_type )

      implicit none

*     Routine to the select the ranges of trial values needed for
*     each of bias flags.  Algorithm based on the mean uncertainties
*     for the MW- and extra-wide lanes. 

      include '../includes/const_param.h' 
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number

      integer*4 ep

* data_type -- character string with type of data being used (L1L2LC)

      character*(*) data_type

* LOCAL VARIABLES
* i, j, k
* na    -- Ambiquity number
* nt    -- Counter of number of L1 and L2 pairs of searchs to run
* ne    -- counter for number of entries for each search
* dL1, dL2 -- Range of values to be tried for L1 and L2 (+- these values
*          checked
* tot   -- Total number of combinations to be tried. Make real*8 to aviod
*          overflow problems
* dd_sng(4) -- Signs used in the double differences (1 -1 -1 1)

      integer*4 i,j,k, l, na, nt, ne, dL1, dL2, dd_sgn(4)
      real*8 tot 

* wls_rng -- Range of MW widelanes that could be needed 
* range_L1 -- Estimate of range of L1 values needed.
* amb_to_dd(2,max_obs) -- Pair of columns which map a particular bias
*     flag ambiquity to double difference observations. 

      real*8 wls_rng, range_L1, amb_to_dd(2,max_obs)

* kbit  -- Bit check logical function
* ion_OK -- Set true is selected biases satisify the ionospheric 
*     contraint

      logical kbit, ion_OK

      data dd_sgn / 1, -1, -1, 1 /

****  Loop over each of bias flags to get the range of values
*     to test.  

      nt = 0
      do i = 1, num_bf

*        See if this bias is already known
         na = bf_flags(i)
         if( wls_ref(1,na).le.0 ) then
*            This bias flag has been set as a reference choice
*            so we don't need to search.  Set the fully resolved
*            bit.
             call sbit(bf_ents(5,na),2,1)
          end if

*         Now check to see if this is a fully resolved bias
          if( .not.kbit(bf_ents(5,na),2) ) then
*             OK, not fully resolved add this to trial values.
*             Find out range of L1 and L2 values needed
              if( 2*wls_sig(1,na).gt.max_wlrng ) then
                  wls_rng = max_wlrng/2
              else
                  wls_rng = wls_sig(1,na)
                  if( abs(wls_all(1,na)).gt.0.30 ) wls_rng = 2*wls_rng
              end if

*             Select the sigma based on MW-WL range and the variance of
*             EX-wide lane
              range_L1 = (77.d0/17.d0)*wls_rng + 
     .                   (60.d0/17.d0)*wls_sig(2,na)
              dL1 = int(range_L1) + 1
              dL2 = int(wls_rng)

*             Now generate the range of values we need.  This is
*             done for L1 and L2
              ne = 0
              nt = nt+1

****          Precalculate the allowable size of the ionospheric
*             delay. 
              call precalc_ion(ep, na, nt )

              do j = -dL1, dL1
                 do k = -dL2, dL2

*                   Now see how this combination will work with the
*                   ionospheric delay constraint.  Here we use the
*                   baseline length, ionospheric activity and the
*                   elevation angle to see how far we can go.  OK is
*                   returned true is ionospheric constraint satisfied.
                    call check_ion(ep, nt, j, k+j, ion_OK )
                    if( ion_OK ) then
                       ne = ne + 1
                       bf_search(1,nt,ne) = j
                       bf_search(2,nt,ne) = k+j
                    end if
                 end do
              end do
*             Save the number of search values for L1 and L2
              num_search(nt) = ne
*             Save the ambiquity number that this search points to
              amb_search(nt)   = na
          end if
      end do

****  Save the total number of search values
      num_sent = nt 

*     Compute the allowable limits on the size of residials allowed
*     due to the effects of the atmosphere on the residuals.
      do i = 1, num_dd
         call precalc_atm( ep, i)
      end do

      if( num_sent.eq.0 ) RETURN
      
      if( 1.eq.2 ) then
         do i = 1, num_sent
            write(*,230) i,'L1',num_search(i),(bf_search(1,i,j),
     .                     j=1,num_search(i))
            write(*,230) i,'L2',num_search(i),(bf_search(2,i,j),
     .                     j=1,num_search(i))

 230        format('SENT: ',i3,1x,a,1x,i4,1x,40(':',2i4))
         end do
      end if

      tot = 1
      do i = 1, num_sent
         tot = tot*num_search(i)
      end do

      if( tot.lt.1.d10 ) then
          write(*,220) ep, num_sent, tot, (amb_search(i), num_search(i),
     .             i=1, num_sent)
 220      format('Epoch ',i6,i4,' Sets, Total Search ',f12.0,
     .          ' Ambiquities: ',100(I4,' # ',I3,2x))
      else
          write(*,225) ep, num_sent, tot, (amb_search(i), num_search(i),
     .             i=1, num_sent)
 225      format('Epoch ',i6,i4,' Sets, Total Search ',d12.6,
     .          ' Ambiquities: ',100(I4,' # ',I3,2x))
      end if      

      if( tot.gt.max_tot_search ) then
         write(*,*) 'Total too large; Limit ',max_tot_search 
         num_sent = 0
         RETURN
      end if

****  Now form the operator which will map the trial ambiquity values
*     to the double difference data used.
      do i = 1, num_sent

*        See which double differences this pair of ambiquities are used
*        in
         do j = 1, num_dd
            amb_to_dd(1,j) = 0
            amb_to_dd(2,j) = 0
         end do

*        Generate one column each for L1 and L2 of the effects of this 
*        ambiquity set on the double difference.  These columns we 
*        multiply by the obs_to_res matrix to get a pair of columns of the 
*        amb_to_res matrix.  Loop with number of data types since these
*        are paired possibly if L1L2 is used. (P1P2PC should not
*        by used in part of the code).
         do j = 1, num_dd, num_dtype
*           See if this ambiquity is used and get the correct sign
            do k = 1, 4
*              Check if the site and satellite match between the
*              ow_des and bf_ents for the ambiquity
               if( ow_des(1,dd_des(k,j)).eq.bf_ents(1,amb_search(i)) 
     .                               .and.
     .             ow_des(2,dd_des(k,j)).eq.bf_ents(2,amb_search(i)) )
     .                                                       then
*                  OK, Site and satellite match and do this ambiquity
*                  is used in the double difference.  Generate the
*                  entries for each this pair.  
* NOTE: There are complications here.  For the moment assume LC only,
*        Solved by using counter l to put entries in correct slot.
                   l = 0 
                   if( index(data_type,'L1').gt.0 ) then
                       amb_to_dd(1,j) = dd_sgn(k)*vel_light/fR1
                       l = l+1
                   end if
                   if( index(data_type,'L2').gt.0 ) then
                       amb_to_dd(2,j+l) = dd_sgn(k)*vel_light/fR2
                       l = l+1
                   end if
                   if( index(data_type,'LC').gt.0 ) then
                       amb_to_dd(1,j+l) = dd_sgn(k)*lcf1*vel_light/fR1
                       amb_to_dd(2,j+l) = dd_sgn(k)*lcf2*vel_light/fR1
                       l = l+1
                   end if
               end if
            end do
         end do

****     We have generate the column pair for ambiquities to double
*        differences.  Now, generate the column pair for ambiquities
*        to residuals remembering to subtract the direct contribution
*        from changing the ambiquity.
         do j = 1, num_dd
            amb_to_res(j,(i-1)*2+1) = -amb_to_dd(1,j)
            amb_to_res(j,(i-1)*2+2) = -amb_to_dd(2,j)
            do k = 1, num_dd
               amb_to_res(j,(i-1)*2+1) = amb_to_res(j,(i-1)*2+1) +
     .              obs_to_res(j,k)*amb_to_dd(1,k) 
               amb_to_res(j,(i-1)*2+2) = amb_to_res(j,(i-1)*2+2) +
     .              obs_to_res(j,k)*amb_to_dd(2,k) 
            end do
         end do

      end do 

****  Now do the ionosphere amb to res matrix:
      do i = 1, num_sent

*        See which double differences this pair of ambiquities are used
*        in
         do j = 1, num_ddion
            amb_to_dd(1,j) = 0
            amb_to_dd(2,j) = 0
         end do

*        Generate one column each for L1 and L2 of the effects of this 
*        ambiquity set on the double difference.  These columns we 
*        multiply by the obs_to_res matrix to get a pair of columns of the 
*        amb_to_res matrix.  Loop with number of data types since these
*        are paired possibly if L1L2 is used. (P1P2PC should not
*        by used in part of the code).
         j = 0
         do l = 1, num_dd, num_dtype
            j = j + 1
*           See if this ambiquity is used and get the correct sign
            do k = 1, 4
*              Check if the site and satellite match between the
*              ow_des and bf_ents for the ambiquity
               if( ow_des(1,dd_des(k,l)).eq.bf_ents(1,amb_search(i)) 
     .                               .and.
     .             ow_des(2,dd_des(k,l)).eq.bf_ents(2,amb_search(i)) )
     .                                                       then
*                  OK, Site and satellite match and do this ambiquity
*                  is used in the double difference.  Generate the
*                  entries for each this pair.  
                   amb_to_dd(1,j) = dd_sgn(k)
                   amb_to_dd(2,j) = -(fR1/fR2)*dd_sgn(k)
               end if
            end do
         end do

****     We have generate the column pair for ambiquities to double
*        differences.  Now, generate the column pair for ambiquities
*        to residuals remembering to subtract the direct contribution
*        from changing the ambiquity.
         do j = 1, num_ddion
            amb_to_ion(j,(i-1)*2+1) = -amb_to_dd(1,j)
            amb_to_ion(j,(i-1)*2+2) = -amb_to_dd(2,j)
            do k = 1, num_ddion
               amb_to_ion(j,(i-1)*2+1) = amb_to_ion(j,(i-1)*2+1) +
     .              ion_to_res(j,k)*amb_to_dd(1,k) 
               amb_to_ion(j,(i-1)*2+2) = amb_to_ion(j,(i-1)*2+2) +
     .              ion_to_res(j,k)*amb_to_dd(2,k) 
            end do
         end do

*        Test write:
         if( 1.eq.2 ) then
             write(*,310) i,'L1',(amb_to_dd(1,j), j =1,num_ddion)
             write(*,310) i,'L2',(amb_to_dd(2,j), j =1,num_ddion)
 310         format('IA_DD  ',i4,1x,a,1x,20(F10.4))
             write(*,330) i,'L1',(amb_to_ion(j,(i-1)*2+1),j=1,num_ddion)
             write(*,330) i,'L2',(amb_to_ion(j,(i-1)*2+2),j=1,num_ddion)
 330         format('I_RES ',i4,1x,a,1x,20(F10.4))
         end if
      end do
      
*     Thats all
      return
      end
        
CTITLE INC_SENT

      subroutine inc_sent( sent, num, finished )

      implicit none

*     Routine to increment the search entries.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* sent(max_chan) -- The search entry array.  Note: Arrays are
*     one dimensional with L1/L2 searchs linked together
* num(max_chan) --  Number of search entries per combination (num_search
*     is the actual array).

      integer*4 sent(max_chan), num(max_chan)

* finished   -- Set true when all combinations are searched 

      logical finished

* LOCAL VARIABLES
* ns  -- Counter for position in search entry array

      integer*4 ns

* done  -- Logical to indicate that we have incremented all the the
*     array values needed (basically once one search ends, we need to
*     increment the previous one and so on)

      logical done

****  Increment the last entry and check for "rollover"
      done = .false.
      ns = num_sent
      do while ( .not.done )
         sent(ns) = sent(ns) + 1
         if( sent(ns).le.num(ns) ) then
             done = .true.
         else
*            Reset this value and set ns to increment the
*            previous entry
             sent(ns) = 1
             ns = ns - 1

*            See if we have run out of entries completely
             if( ns.le.0 ) then
                 done = .true.
                 finished = .true.
             end if
         end if
      end do

****  Thats all
      return
      end

CTITLE CHECK_ION

      subroutine check_ion(ep, nt, dL1, dL2, ion_OK ) 

      implicit none

*     Routine to check the impact of the changes in the L1 and L2
*     ambiquities on the extra-wide lane.

      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* nt   -- Search entry bias parameter number
* dL1, dL2 -- Trial changes in the L1 and L2 ambiquities

      integer*4 ep, nt, dL1, dL2 

* ion_OK   -- Set true if combination passes the Extra-wide lane
*             test.

      logical ion_OK

* LOCAL VARIABLES
* dexwl -- Change in the extra-wide implied by the choice of ambiquities
      
      real*8 dexwl


***** Compute the change in the extra-wide lane implied by this
*     combination
      dexwl = dL1 - (77.d0/60.d0)*dL2
      ion_OK = .false.
      if( dexwl.ge.ion_lim(1,nt) .and. dexwl.le.ion_lim(2,nt) ) then
          ion_OK = .true.
      end if

****  Thats all
      return 
      end


CTITLE PRECALC_ION

      subroutine precalc_ion(ep, na, nt ) 

      implicit none

*     Routine to pre-calculate the allowable magnitude of the changes
*     in the extra wide lane. The ionospheric algorithm based separation 
*     of sites, the size of the ionosphere and the elevation angles (of the
*     two satellites in the difference).  

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* na   -- Ambiquity number
* nt   -- Search entry number for this bias flag to be checked 

      integer*4 ep, na, nt

* LOCAL VARIABLES
* pn, pn_ref  -- PRNS for main station and reference station
* ns, ns_ref  -- Site numbers for main and reference station
* i           -- Loop counter

      integer*4 pn, pn_ref, ns, ns_ref, i

* basel  -- Distance between sites (m)
* elev(2) -- The elevation angles to the two satellites (deg)
* ion_zd(2) -- Zenith angle to satellite at the ion layer height (rad)
* ion_mag   -- Computed possible magnitude of ionospheric delay (cycles)

      real*4 basel, elev(2), ion_zd(2), ion_mag


****  OK, Start with the ionospheric delay limits. 
*     Get the reference satellite and other site for this  ambiquity
      pn     = bf_ents(2,na)
      ns     = bf_ents(1,na)
     
      pn_ref = bf_ents(2,wls_ref(1,na))
      ns_ref = bf_ents(1,wls_ref(2,na))

*     Compute the baseline length between the two sites
      basel = sqrt((curr_site_xyz(1,ns)-curr_site_xyz(1,ns_ref))**2 +
     .             (curr_site_xyz(2,ns)-curr_site_xyz(2,ns_ref))**2 +
     .             (curr_site_xyz(3,ns)-curr_site_xyz(3,ns_ref))**2 )

*     Get elevation angles at our site to the two satellites
      elev(1) = 0.0
      elev(2) = 0.0
      do i = 1, num_chan_se(ns,ep)
         if( ctop_cse(i,ns,ep).eq.pn ) then
             elev(1) = elev_cse(i,ns,ep)
         else if( ctop_cse(i,ns,ep).eq.pn_ref ) then
             elev(2) =  elev_cse(i,ns,ep)
         end if
      end do

*     The elevation angles may not be present if a satellite has
*     not been observed at this epoch.  Just set value to 10 deg
*     to allow the tolerances to be computed.
      if( elev(1).le.0 .or. elev(2).le.0 ) then 
          if( elev(1).le.0 ) elev(1) = 10.0
          if( elev(2).le.0 ) elev(2) = 10.0
      end if

****  Now compute the expected magnitude of the ionospheric noise
      do i = 1,2
         ion_zd(i) = asin(earth_rad/(earth_rad+ion_ht)*
     .                 sin((90-elev(i))*pi/180.d0))
      end do

      ion_mag = ion_ppm*(basel/(vel_light/fR1))*(1.d0/cos(ion_zd(1))+
     .                                           1.d0/cos(ion_zd(2)))

****  Now compute minumum and maximum change allowed to the extra-
*     widelane so that the values will be +- ion_mag.  Put lower limit
*     of 0.1 cycles on ionosphere.
      ion_lim(1,nt) = -(ion_mag+0.3+wls_all(2,na)) - wls_sig(2,na)
      ion_lim(2,nt) =  (ion_mag+0.3-wls_all(2,na)) + wls_sig(2,na)

*     OK, tell user what is happening
      if( 1.eq.2 ) 
     .write(*,220) ep, na, nt, elev, basel/1000.d0, 
     .             (ion_lim(i,nt),i=1,2), ion_mag
 220  format('EX-WL constraint at ep ',i6,' Amb ',2i3,' Elevs ',
     .       2F7.2,' Length ',F9.3,' km; Bounds ',3F7.3,' cycles')

****  Thats all
      return 
      end 
         
CTITLE PRECALC_ATM

      subroutine precalc_atm(ep, nd ) 

      implicit none

*     Routine to pre-calculate the allowable magnitude of the changes
*     in the atmospheric delay for each double difference.  

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* nd   -- Double difference number 

      integer*4 ep,  nd

* pn, pn_ref  -- PRNS for main station and reference station
* ns, ns_ref  -- Site numbers for main and reference station
* i           -- Loop counter

      integer*4 pn, pn_ref, ns, ns_ref, i

* basel  -- Distance between sites (m)
* elev(2) -- The elevation angles to the two satellites (deg)
* atm_lim -- Limit on size of atmospheric delay allowed in search (m)
* amap(2) -- Approximate mapping function to map atm to elevation

      real*4 basel, elev(2), atm_lim, amap(2), tolbas

****  OK, Start with the atmospheric delay size calculation
*     Get the reference satellite and other site for this double difference
      ns     = ow_des(1,dd_des(1,nd))
      pn     = ow_des(2,dd_des(1,nd))
     
      ns_ref = ow_des(1,dd_des(4,nd))
      pn_ref = ow_des(2,dd_des(4,nd))

*     Compute the baseline length between the two sites
      basel = sqrt((curr_site_xyz(1,ns)-curr_site_xyz(1,ns_ref))**2 +
     .             (curr_site_xyz(2,ns)-curr_site_xyz(2,ns_ref))**2 +
     .             (curr_site_xyz(3,ns)-curr_site_xyz(3,ns_ref))**2 )

*     Get elevation angles at our site to the two satellites
      elev(1) = 0.0
      elev(2) = 0.0
      do i = 1, num_chan_se(ns,ep)
         if( ctop_cse(i,ns,ep).eq.pn ) then
             elev(1) = elev_cse(i,ns,ep)
         else if( ctop_cse(i,ns,ep).eq.pn_ref ) then
             elev(2) =  elev_cse(i,ns,ep)
         end if
      end do

      if( elev(1).le.0 .or. elev(2).le.0 ) then 
          write(*,120) nd, ns, ns_ref, pn, pn_ref, elev
 120      format('**ERROR** Problem with elevation angles for',
     .           ' Double Difference ',i3,' Sites ',2i3,' PRNS ',
     .           2(i3.2,1x),' Elevations ',2F6.2,' deg')
          if( elev(1).le.0 ) elev(1) = 10.0
          if( elev(2).le.0 ) elev(2) = 10.0
      end if


****  Now handle the atmospheric delay size.  Here we use the distance
*     and elevation angles.  (Maybe in the future the height would also
*     be a useful parameterization.)
      stratm_dcorr = 0.01   ! Units: m/sqrt(D km)
      stratm_base  = 0.01
      tolbas = stratm_dcorr*sqrt(basel/1000.d0) 
      atm_lim = max(stratm_base,tolbas)
      amap(1) = 1./(sin(elev(1)*pi/180.)+0.030)
      amap(2) = 1./(sin(elev(2)*pi/180.)+0.030)
      max_dd_elev(nd) = atm_lim*(amap(1)+amap(2))
      if( 1.eq.2 )
     .write(*,220) ep, nd, pn, pn_ref, elev, basel/1000.d0, 
     .             max_dd_elev(nd)
 220  format('ATM   constraint at ep ',i6,' DD  ',3i3,' Elevs ',
     .       2F7.2,' Length ',F9.3,' km; Limit ',F7.3,' m')

****  Thats all
      return
      end 
       
CTITLE VERIFY_SEARCH

      subroutine verify_search( an, ist, iend)

      implicit none

*     Routine to verify that we can fix this bias with a search by
*     checking the times it is above the elevation cutoff and that
*     a fixed bias is also above the elevation cutoff.

      include 'track_com.h'

* PASSED VARIABLES
* an -- Ambiquity number being tested 
* ist, iend -- Start and stop epoch numbers for the search.  Returned
*       as zero if search not possible

      integer*4 an, ist, iend

* LOCAL VARIABLES
* i,j  -- Loop counter
* id_ist, id_iend -- Start and stop of good data with ambiquity
* ir_ist, ir_iend -- Start and stop of fixed bias parameters
* ns, na  -- site and ambiquity numbers

      integer*4 i, j, id_ist, id_iend, ir_ist, ir_iend, ns, na

* data_OK  -- Function to test is data is good.
* kbit     -- Check status of bit
* ddbf_OK  -- Checks if all oneways in the definition of ambiquity are
*             good.

      logical data_OK, kbit, ddbf_OK

****  OK, initialize every to zero
      id_ist  = 0
      id_iend = 0
      ir_ist  = 0
      ir_iend = 0 

****  Now loop over duration of bias parameter (not all data may be
*     present due to elevation cut-off and bad data)
      ns = bf_ents(1,an)
      do i = bf_ents(3,an), bf_ents(4,an)

*        See if we have data start
         do j = 1, num_chan_se(ns,i)
            if( data_OK(data_flag_cse(j,ns,i),data_mask) ) then
*               Good data point found.  See which satellite and
*               see if it is either the one we are looking for or 
*               if its amabiquity is fixed.
                if( ctop_cse(j,ns,i).eq.bf_ents(2,an) ) then
*                   This is our satellite.  If start time not
*                   set, set it now
                    if( ddbf_OK(i,an) ) then
                       if( id_ist.eq.0 ) id_ist = i
                       id_iend = i 
                    end if
                else
*                   See if ambiquity has been fixed. Get ambiquity number
                    na = amb_point_cse(j,ns,i)
                    if( kbit(bf_ents(5,na),2) ) then
*                       OK: Fixed ambiquity mark start and end times
                        if( ir_ist.eq.0 ) ir_ist = i
                        ir_iend = i
                    end if
                end if
            end if
         end do
      end do

****  OK, Now see what we have.  The start time will be max of start times
*     and end will be the min of the times
      ist  = max(id_ist,  ir_ist)
      iend = min(id_iend, ir_iend)

      if( ist.eq.0 .or. iend.eq.0 ) then
          ist = 0
          iend = -1
      else
C         write(*,220) an, ist, iend, id_ist, id_iend, ir_ist, ir_iend
C220      format('AMB ',i3,' Start/End ',2i7,' Data ',2i7,' Ref ',2i7)
      end if

*     Thats all
      return
      end

CTITLE DDBF_OK

      logical function ddbf_OK(ep, an)

      implicit none

*     Logical function that returns true if the data for the sites
*     and satellites at this epoch for this ambiquity number is all
*     OK. 


      include 'track_com.h'

* PASSED VARIABLES
* ep -- Epoch number
* an -- Ambiquity number being tested
 
      integer*4 ep, an

* LOCAL VARIABLES
* i1, i2  -- Two stations in the definition of the bias flag
* s1, s2  -- Two satellites in the definition of the bias flag
* na      -- Ambiquity number of the other ambiquities
* j       -- Loop counter

      integer*4 i1, i2, s1, s2, na, j 

* data_OK  -- Logical function that returns true if data OK

      logical data_OK 

****  Get each of the oneways involved and see if they are OK
      i1 = bf_ents(1,an)
      s1 = bf_ents(2,an)

*     Now for the other pair of bias flags in the double difference
*     get the satellite and site.      
      na = wls_ref(1,an)
      s2 = bf_ents(2,na)
      na = wls_ref(2,an)
      i2 = bf_ents(1,na)

*     OK, now test if all these one ways are good.
      ddbf_OK = .true.
      do j = 1, num_chan_se(i1,ep)
         if( ctop_cse(j,i1,ep).eq.s1 .and. 
     .       .not.data_OK(data_flag_cse(j,i1,ep),data_mask) ) 
     .                                              ddbf_OK = .false.
         if( ctop_cse(j,i1,ep).eq.s2 .and. 
     .       .not.data_OK(data_flag_cse(j,i1,ep),data_mask) ) 
     .                                              ddbf_OK = .false.
      end do

*     Check the second site
      do j = 1, num_chan_se(i2,ep)
         if( ctop_cse(j,i2,ep).eq.s1 .and. 
     .       .not.data_OK(data_flag_cse(j,i2,ep),data_mask) ) 
     .                                              ddbf_OK = .false.
         if( ctop_cse(j,i2,ep).eq.s2 .and. 
     .       .not.data_OK(data_flag_cse(j,i2,ep),data_mask) ) 
     .                                              ddbf_OK = .false.
      end do

****  Thats all
      return
      end

CITLE REORDER

      subroutine reorder

      implicit none

*     This routine will loop over the data at each station and
*     reorder the channels so that they are in decreasing order
*     of elevation angle.  When double differences are formed
*     later, this will ensure that the highest elevation angle
*     one-way will be used as the reference.

      include 'track_com.h'

* PASSED VARIABLES
* None

* LOCAL VARIABLES
* ns, ep  -- Station number and epoch counter
* j       -- Channel loop counter
* in      -- Index to sorted elevations
* nc      -- number of channels (short name)
* elev_ind(max_chan) -- Pointer to sorted elevation list

      integer*4 ns, ep, i, j, in, nc, elev_ind(max_chan)

* elev_org(max_chan) -- Original elevation angles which need to
*     be sorted. 

      real*8  elev_org(max_chan)

* L1o_srt(max_chan) - Sorted L1 data
* L2o_srt(max_chan) - Sorted L2 data
* P1o_srt(max_chan) - Sorted P1 data
* P2o_srt(max_chan) - Sorted P2 data

      real*8 L1o_srt(max_chan), L2o_srt(max_chan), P1o_srt(max_chan),
     .       P2o_srt(max_chan)
                  
* elev_srt(max_chan) - Sorted elevation angles

      real*8 elev_srt(max_chan)
                  
* data_flag_srt(max_chan) -- Sorted data flag
* ctop_srt(max_chan) -- Sort channel to PRN array

      integer*4 data_flag_srt(max_chan), ctop_srt(max_chan)

      real*4 tec_srt(max_chan)  ! Copy of tec_los value for re-ordering

****  OK loop over each station.  (est_pos must have been called 
*     before this)
      do ns = 1, num_site

*        Loop over all the epochs at this site
         do ep = 1, num_epochs

*           Now sort the elevation angles by size
            do j = 1, num_chan_se(ns,ep)
               elev_org(j) = elev_cse(j,ns,ep)
               elev_ind(j) = j 
            end do
            nc = num_chan_se(ns,ep)

*           Now sort the elevation angles returning the index
*           array with the sorted values
            call sort_elev(nc, elev_org, elev_ind, max_chan)   

*           OK, now re-arrange the data so that they are in
*           decreasing elevation order.
            do i = 1, nc
*              Copy all of the data into the correctly ordered
*              positions and then copy back
               in = elev_ind(i)
               L1o_srt(i) = L1o_all_cse(in,ns,ep)
               L2o_srt(i) = L2o_all_cse(in,ns,ep)               
               P1o_srt(i) = P1o_all_cse(in,ns,ep)
               P2o_srt(i) = P2o_all_cse(in,ns,ep) 

               elev_srt(i) = elev_cse(in,ns,ep)
    
               data_flag_srt(i) = data_flag_cse(in,ns,ep)
               ctop_srt(i) = ctop_cse(in,ns,ep)
               tec_srt(i) = tec_los_cse(in,ns,ep)
            end do

*           The _srt arrays are now sorted.  Move them back into
*           data space.
            do i = 1, nc
               L1o_all_cse(i,ns,ep) = L1o_srt(i)      
               L2o_all_cse(i,ns,ep) = L2o_srt(i)                 
               P1o_all_cse(i,ns,ep) = P1o_srt(i)      
               P2o_all_cse(i,ns,ep) = P2o_srt(i)      
                                                        
               elev_cse(i,ns,ep)    = elev_srt(i)     
                                                        
               data_flag_cse(i,ns,ep) = data_flag_srt(i)
               ctop_cse(i,ns,ep)      = ctop_srt(i)
               tec_los_cse(i,ns,ep)   = tec_srt(i)
            end do

         end do
      end do

****  Thats all
      return
      end
     
CTITLE SORT_ELEV

      subroutine sort_elev(nd, data, s, na)

      implicit none

*     Routine to sort an array and return an index array
*     from maximum to minimun returning the array s with the
*     indices for the sorted order/
*     Gang Chen MIT 970701
*          simple sort method. will update to more effecient form when
*          i have time later

* nd, na  -- Number of data points and the dimension of the array
      integer*4 nd, na

*      s(na) --- the index order of data array
      real*8 data(na)
      integer*4 s(na)
           
* LOCAL VARIABLES 

      integer*4 i, j

****  Do a brute force sort since there are not many values.
      s(1)= 1
      do i = 2, nd
         j = i
         do while ((j.gt.1))
           if (data(i).gt.data(s(j-1))) then
            s(j)=s(j-1)
            j=j-1
           else
             GOTO 27
           end if
         enddo
27       s(j) = i
      enddo

****  Thats all
      return
      end       

CTITLE BASE_LEN

      real*8 function base_len(ep, s1,s2)

      implicit none

*     Routine to return baseline length in meters bewteen sites
*     i and j

      include 'track_com.h'

* PASSED VARIABLES
* s1, s2  -- The two site numbers
* ep    -- Epoch number

      integer*4 s1, s2, ep

* LOCAL VARIABLES
* i,j -- Loop counters

      integer*4 i,j

* xyz_s1(3), xyz_s2(3) -- Coordinates of the two sites

      real*8 xyz_s1(3), xyz_s2(3)

* s1_known, s2_known -- Logical to indicate if coordinates are known
      logical s1_known, s2_known


****  Get the coordinates of the two sites.  If site is not-kinematic
*     use the apriori, other wise get the coordinates from the kine_pos
*     array.
*     See if same site numbers (length = 0)
      if( s1.eq.s2 ) then
          base_len = 0.d0
          RETURN
      end if

*     See if site is kinematic 
      s1_known = .false.
      s2_known = .false.

      do i = 1, num_kine
         if( kine_to_site(i).eq.s1 ) then
             do j = 1, 3
                xyz_s1(j) = kine_xyz(j,i,ep)
             end do
             s1_known = .true.
         else if ( kine_to_site(i).eq.s2 ) then
             do j = 1, 3
                xyz_s2(j) = kine_xyz(j,i,ep)
             end do 
             s2_known = .true.
         end if
      end do

*     If either coordinate is not known then get value from apriori
      if( .not. s1_known ) then
          do j = 1, 3
             xyz_s1(j) = site_apr(j,s1)
          end do
      end if
      if( .not. s2_known ) then
          do j = 1, 3
             xyz_s2(j) = site_apr(j,s2)
          end do
      end if

****  OK, now compute the length
      base_len = sqrt((xyz_s1(1)-xyz_s2(1))**2 +
     .                (xyz_s1(2)-xyz_s2(2))**2 +
     .                (xyz_s1(3)-xyz_s2(3))**2 ) 


****  Thats all
      return
      end
           
CTITLE GET_OW_ION

      subroutine get_ow_ion(ep, i, j, ni, elev ) 

      implicit none

*     Routine to get the one-way L1 and L2 values at epoch ep, for
*     site i, channel j.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* i, j -- Site and channel numbers
* ni   -- Running sum of number of ion oneways

      integer*4 ep, i,j, ni

* elev  -- Elevation angle (deg)

      real*8  elev


* LOCAL VARIABLES
* na   -- Number of current ambiquity
* k    -- Loop counter

      integer*4 na, k

* ion_zd  -- Zenith distance to satellite at ionosheric layer
*     height

      real*8 ion_zd

****  OK, loop over each data type
*     OK, based on data types, compute the o-minus-c values
      k = 0 

      ni = ni + 1
      na = amb_point_cse(j,i,ep)
      ow_ion(ni) = (L1o_all_cse(j,i,ep)+ambiq_all(1,na)) -
     .     (fR1/fR2)*(L2o_all_cse(j,i,ep)+ambiq_all(2,na))     
      ow_iov(ni) = data_var(1,ctop_cse(j,i,ep))

      do k = 1, num_site
         ow_ionp(k,ni) = 0.d0
      end do

****  Now compute the partial
      ion_zd =  asin(earth_rad/(earth_rad+ion_ht)*
     .                 sin((90-elev)*pi/180.d0))
      ow_ionp(i,ni) = 1.d0/cos(ion_zd) 


****  Thats all for the moment, add more later

      return
      end

CTITLE RESOLVE_SLIP

      subroutine resolve_slip(ep, pe, i,j, lv, kv, dexwl)

      implicit none

*     Routine to resolve whether a slip on sites i-j and satellites
*     lv-kv ocurrs on satellite lv or kv.  Once determinated, the
*     slip is marked on correct one-way.  (The j-site will be
*     marked).

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep  -- Epoch number
* pe  -- Previous epoch to which comparisons are being made
* i,j -- Two site numbers
* lv, kv -- Two satellite PRN numbers

      integer*4 ep, pe, i,j, lv, kv

* dexwl  -- Change in the extra-widelane

      real*8 dexwl

* LOCAL VARIABLES
* m, mv -- Third satellite index and satellite number
* sv    -- Final selection of bad satellite

      integer*4 m, mv, sv

* exwl  -- Estimate of the extra wide-lane (also known in gamit as LG)
* exwlp -- Estimate from previous epoch
* data(4,4) -- Estimates of the phase and range for each of the one-ways
*              in the double difference.  The first index is over phase
*              and range (L1,L2,P1,P2); the second index is by site
* dexnew     -- Change in widelane from third satellite choice
* elev(4)   -- Four elevation angles

      real*8 exwl, exwlp, data(4,4), dexnew, elev(4)

* dd_OK     -- Set true is double difference should be OK
* OK(4)     -- The status of each of the oneways.

      logical dd_OK, OK(4), kbit

* df(4) -- Data flags from get_ow
      integer*4 df(4)

****  OK, form the combinations with the other satellites that are
*     visible at this and the previous epoch.

      do m = 1, num_prn
         if( m.ne.lv .and. m.ne.kv ) then
*           OK, this is a different satellite.  Form the change
*           in widelane with lv and see if matches the change
*           with kv (if it does then lv is the problem).
            mv = prn_used(m)
            call get_ow(pe, i, mv, data(1,1), elev(1), OK(1),df(1))
            call get_ow(pe, i, lv, data(1,2), elev(2), OK(2),df(2))
            call get_ow(pe, j, mv, data(1,3), elev(3), OK(3),df(3))
            call get_ow(pe, j, lv, data(1,4), elev(4), OK(4),df(4))           
  
            dd_OK = OK(1) .and. OK(2) .and. OK(3) .and. OK(4)
            if( dd_OK ) then
                exwlp =   (data(1,1)-data(1,2)) -
     .                     (data(1,3)-data(1,4))   - (fR1/fR2)*
     .                    ((data(2,1)-data(2,2)) -
     .                     (data(2,3)-data(2,4)) )

*               Get the next epoch where slip occurrs
                call get_ow(pe, i, mv, data(1,1), elev(1),OK(1),df(1))
                call get_ow(pe, i, lv, data(1,2), elev(2),OK(2),df(2))
                call get_ow(pe, j, mv, data(1,3), elev(3),OK(3),df(3))
                call get_ow(pe, j, lv, data(1,4), elev(4),OK(4),df(4))
                dd_OK = OK(1) .and. OK(2) .and. OK(3) .and. OK(4)
                if( dd_OK ) then
                    exwl =   (data(1,1)-data(1,2)) -
     .                       (data(1,3)-data(1,4))   - (fR1/fR2)*
     .                       ((data(2,1)-data(2,2)) -
     .                       (data(2,3)-data(2,4)) ) 
                    dexnew = exwl - exwlp
                    if( abs(dexnew).gt.max_ion_jmp ) then

*                       OK, we see a slip here as well so the slip 
*                       must be on lv, so mark.  Otherwize it must
*                       on kv so mark that.
                        sv = lv
                    else
                        sv = kv
                    end if

*                   OK, Mark the slip.  See if already marked.
                    do mv = 1, num_chan_se(j,ep)
                       if( ctop_cse(mv,j,ep).eq. sv .and.
     .                     .not.kbit(data_flag_cse(mv,j,ep),4) ) then
                          write(*,200) ep, j, sv, dexwl, data_mask
 200                      format('Marking EX-WL slip at epoch ',i6,
     .                           ' Site ', I2,' PRN ',i3.2,' Size ',
     .                           F20.2,' cycles Mask ',i4)
                          call mark_slip( ep, j, sv, 4, dexwl)
                       end if
                    end do
                end if
            end if
         end if
      end do

****  Thats all
      return 
      end

CTITLE FIND_GAP

      subroutine find_gap(ep,  i,j, lv, kv)

      implicit none

*     Routine to find which one-way sequence (j/lv or j/kv) has the
*     gap in it.  We do this by checing the previous epoch

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep  -- Epoch number
* i,j -- Two site numbers
* lv, kv -- Two satellite PRN numbers

      integer*4 ep, i,j, lv, kv

* LOCAL VARIABLES
* data(4) -- Estimates of the phase and range for each of the one-ways
*              in the double difference.  The first index is over phase
*              and range (L1,L2,P1,P2); the second index is by site
* elev    -- Elevation angle

      real*8 data(4), elev

* OK     -- The status of each of the oneways.

      logical OK

* df     -- Data flag
      integer*4 df

****  OK, Check the previous epoch to see where gap is

      call get_ow(ep-1, i, lv, data, elev, OK, df)
      if( .not.OK ) then
c         write(*,200) ep, j, lv
c200      format('Marking Gap bias flag at epoch ',i6,' Site ',
c    .           I2,' PRN ',i2.2,' first')
          call mark_slip( ep, j, lv, 4, -1.d0)
      end if

*     Check the other satellite
      call get_ow(ep-1, i, kv, data, elev, OK, df)
      if( .not.OK ) then
c         write(*,220) ep, j, kv
c220      format('Marking Gap bias flag at epoch ',i6,' Site ',
c    .           I2,' PRN ',i2.2,' second')
c         write(*,200) ep, j, kv
          call mark_slip( ep, j, kv, 4, -1.d0)
      end if

****  Thats all
      return
      end

CTITLE UPDATE_WLS_ALL

      subroutine update_wls_all(na, dL1, dL2, solupdate)

      implicit none

*     Routine to updates the average values of the MW and EX wide
*     lanes based on changing ambiquity number na by dL1 and dL2.
*     Since the ambiquities are formed as double differences, check
*     all double differences that depend on na and update appropriately,

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* na  -- Ambiquity number being updated
      integer*4 na

* dL1, dL2 -- changes to the number of cycles
      integer*4 dL1, dL2

* solupdate -- Set true to update solution vector as well

      logical solupdate, kbit

* LOCAL VARIABLES
* i  -- Loop counter
      integer*4 i

      integer*4 PtoL   ! PRN to list
     .,         sv     ! List number for primary satellite in ambiguity

****  Update the one which is being changed
      sv = PtoL(bf_ents(2,na))

C     ambiq_all(1,na) = ambiq_all(1,na) + dL1*fR1/fL1(sv)
C     ambiq_all(2,na) = ambiq_all(2,na) + dL2*fR2/fL2(sv)
      ambiq_all(1,na) = ambiq_all(1,na) + dL1
      ambiq_all(2,na) = ambiq_all(2,na) + dL2

C     wls_all(1,na) = wls_all(1,na) - (dL1-dL2)*fR1/fL1(sv)
C     wls_all(2,na) = wls_all(2,na) - (dL1-(fR1/fR2)*dL2)*fR1/fL1(sv)
      wls_all(1,na) = wls_all(1,na) - (dL1-dL2)
      wls_all(2,na) = wls_all(2,na) - (dL1-(fR1/fR2)*dL2)
      if( solupdate ) call update_sol(na, dL1, dL2)

c      print *,'update_wls_all', na, dL1, dL2

*     Now loop over all other ambiquities to see which one depend
*     on this one (Make sure bias is not fixed already.  Can happen
*     when user inputs ambiqiuties)

      do i = 1, num_ambs
         if( wls_ref(1,i).eq.na .and. i.ne.na .and. 
     .       .not. kbit(bf_ents(5,i),2)  ) then
*            This depends on na, adjust the mean value for the
*            wide lanes
!            ambiq_all(1,i) = ambiq_all(1,i) - dL1
!            ambiq_all(2,i) = ambiq_all(2,i) - dL2
             wls_all(1,i) = wls_all(1,i) + (dL1-dL2)
             wls_all(2,i) = wls_all(2,i) + (dL1-(fR1/fR2)*dL2)
!             print *,'NA2 ', i, wls_all(2,i),dL1,dL2
c             print *,'update_sol1 ', i, dL1, dL2
!             if( solupdate) call update_sol(i, -dL1,-dL2)

         end if 

*        Check the second entry in the double differences:
*        Normally the next two are reference station entries
*        and these are unlikely to change except for multi-station
*        processng.
         if( wls_ref(2,i).eq.na  .and. i.ne.na .and. 
     .       .not. kbit(bf_ents(5,i),2) ) then
*            This depends on na, adjust the mean value for the
*            wide lanes
!            ambiq_all(1,i) = ambiq_all(1,i) + dL1
!            ambiq_all(2,i) = ambiq_all(2,i) + dL2
             wls_all(1,i) = wls_all(1,i) - (dL1-dL2)
             wls_all(2,i) = wls_all(2,i) - (dL1-(fR1/fR2)*dL2)
!             print *,'NA3 ', i, wls_all(2,i),dL1,dL2
!             if( solupdate) call update_sol(i, +dL1,+dL2)

         end if

*        Check the second entry in the double differences
         if( wls_ref(3,i).eq.na  .and. i.ne.na.and. 
     .        .not. kbit(bf_ents(5,i),2) ) then
*            This depends on na, adjust the mean value for the
*            wide lanes
!            ambiq_all(1,i) = ambiq_all(1,i) - dL1
!            ambiq_all(2,i) = ambiq_all(2,i) - dL2
             wls_all(1,i) = wls_all(1,i) + (dL1-dL2)
             wls_all(2,i) = wls_all(2,i) + (dL1-(fR1/fR2)*dL2)
!             print *,'NA4 ', i, wls_all(2,i),dL1,dL2
!             if( solupdate) call update_sol(i, -dL1,-dL2)

         end if
        
      end do

****  Thats all
      return
      end

CTITLE UPDATE_SOL

      subroutine update_sol( na, dL1, dL2 )

      implicit none

*     Routine to update sol_vec and float estimate of LC or L1+L2 for
*     a change in ambiguity of one-ways in double difference.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* na        -- Ambiguity to be updated
* dL1, dL2  -- Updates changes to the number of cycles
      integer*4 na, dL1, dL2

* LOCAL VARIABLES

      integer*4 ina  ! Index in the estimated ambiguties
     .,         ne   ! Number of data types

****  Update the values for the ambig_float (m) and sol_vec (m)
      ne = 0
      if( index(float_type,'LC').gt.0 )  ne = 1
      if( index(float_type,'L1').gt.0 .and.
     .    index(float_type,'L2').gt.0 )  ne = 2
      ina = amb_parn(1,na)
      if( ina.gt.0 ) then

*         See if L1+L2 or LC
          if( ne.eq. 1 )  then  ! LC solution
              sol_vec(ina) = sol_vec(ina) + (lcf1*dL1 + lcf2*dL2)*
     .                                       lambda(1)
              ambiq_float(1,na) = sol_vec(ina)
c              print *,'Update ',na, dL1, dL2, ina, sol_vec(ina), 
c     .                lcf1*dL1 + lcf2*dL2
          else
              sol_vec(ina) = sol_vec(ina) + dL1*lambda(1)
              ambiq_float(1,na) = sol_vec(ina)
              ina = amb_parn(2,na)
              if( ina.gt.0 ) then
                 sol_vec(ina) = sol_vec(ina) + dL2*lambda(2)
                 ambiq_float(2,na) = sol_vec(ina)
              endif
          endif
      endif

*     Thats all
      return
      end


CTITLE OUT_RESID

      subroutine out_resid(ep, ks)

      implicit none

*     Routine to write out the residuals from the current epoch of
*     data.  The files to which the residuals are written are 
*     assumed to have been opened in the main program

      include 'track_com.h'

* PASSED VARIABLES
* ep -- Epoch number
* ks -- Kinematic site being processed

      integer*4 ep, ks

* LOCAL VARIABLES
* i,j,k  -- Loop counters
* luo    -- Unit number to write to
* num(max_dtype) -- Number of measurements of each data type
* pf_prn(max_obs) -- Satellite numbers for each residual
* pf_site(max_obs) -- Site number for residuals
* pf_amb(max_obs)   -- Ambiuity number pointed to be measurements
* pf_bf5(max_obs)   -- Bias flag entry 5 (fixed status) for each measurement
* date(5)  -- Date of measurement
* doy      -- Day of year 
* num_res  -- Number of single-difference residuals to output
* pn -- Short version of pf_prn(k)
* eb -- Elevation bin number (1-18 in 5 degress bins)

      integer*4 i,j,k, luo, num(max_dtype), pf_prn(max_obs), 
     .          pf_site, pf_amb(max_obs), pf_bf5(max_obs),
     .          date(5), doy, num_res, pn, eb, is

      integer*4 ref_prn  ! Reference PRN for Double Diffs

* pf_elev(max_obs)  -- Elevation angles  

      real*4 pf_elev(max_obs), pf_az(max_obs)

* sum(max_dtype) -- Sum for computing mean
* pf_res(max_dtype, max_obs) -- Postfit residuals (mm)
* sectag -- Seconds tag
* sec_in_day -- Seconds into day
* curr_mjd -- MJD of measurement
* fday     -- Fractional doy of year
* mean     -- Mean residual

      real*8 sum(max_dtype), pf_res(max_dtype, max_obs), sectag,
     .       curr_mjd, fday, mean, sec_in_day

* format  -- Format for writing out residuals (varies because of number of
*            data types possible). 
      character*80 format

* pf_ref(max_obs) -- Character to show which is the reference satellite
      character*1 pf_ref(max_obs)

      integer*4 UtoL   ! Function to return the number in list of used
                       ! satellites (num_prn and prn_used).

*
****  Loop over the number of double differences we have and convert
*     back to single difference
      pf_site = kine_to_site(ks)
      ref_prn = 0
      do j = 1, num_dtype
         sum(j) = 0
         num(j) = 0
         do i = 1, num_dd, num_dtype
*            if( ep.ge.debug_start .and. ep.le. debug_end ) then
*               print *,'Site Test ',i, num_dd, pf_site, dd_des(:,i),
*    .             'OW_DES ', ow_des(:,dd_des(1,i))
*            endif
             if( pf_site.eq.ow_des(1,dd_des(1,i))) then
                num(j) = num(j) + 1
                sum(j) = sum(j) + sol_obs(i+j-1)
                if( ref_prn.eq.0 ) then
                    ref_prn = ow_des(2,dd_des(4,i))
                else
                    if( ref_prn.ne.ow_des(2,dd_des(4,i)) ) then
                        write(*,220) ep, ref_prn, i, dd_des(:,i),
     .                               ow_des(:,dd_des(:,i))
 220                    format('REFPRN Change Ep ',I6,' DD ',i4,
     .                         'DD_DES ',4i4,4(' Site/PRN ',2I3))
                    endif
                endif
             end if
*            if( ep.ge.debug_start .and. ep.le. debug_end ) then
*               print *,'Site Sum  ',i, j, num_dd, pf_site, num(j)
*            endif

         end do
      end do

****  OK, Get the mean and remove from all data and assign to
*     reference satellite
      num_res = 0
*     if( ep.ge.debug_start .and. ep.le. debug_end ) then
*        print *,'NumberRes ',num_res, num(:), 'NDR ',num_dtype
*     endif
      do j = 1, num_dtype
          if( num(j).gt.0 ) then
C            mean = -sum(j)/(num(j)+1)
             mean = -sum(j)/num(j)
             k = 0
             do i = 1, num_dd, num_dtype
*               Routine is called for each kinematic site.  Make
*               sure this double difference if for this site
                if( pf_site.eq.ow_des(1,dd_des(1,i))) then
                   k = k + 1
                   pf_res(j,k) = (sol_obs(i+j-1) + mean)*1000.d0
                   pf_prn(k) = ow_des(2,dd_des(1,i))
                   pf_ref(k) = ' '
                end if
             end do
*            Finally assign the mean the reference satellite
             num_res = k+1
             pf_res(j,k+1) = mean*1000.d0
             pf_prn(k+1) = ref_prn    !   ow_des(2,dd_des(4,1))
             pf_ref(k+1) = 'R'
          end if
      end do
      if( ep.ge.debug_start .and. ep.le. debug_end ) then
          write(*,810) num_dtype, num_dd, num_res,ks, pf_site,
     .                (k, pf_res(1,k),pf_prn(k), k= 1,num_res)
 810      format('POSTFIT RES NumT ',i4,' NumDD ',i4,' NumRes ',3i4,
     .              50(i3,F10.3,' P',i3.3))
      end if

****  Now get the information on the site and satellite
c      do k = 1, num_res-1
c         i = (k-1)*num_dtype+1
c         pf_prn(k) = ow_des(2,dd_des(1,i))
c         pf_site(k) = ow_des(1,dd_des(1,i))
c         pf_ref(k) = ' '
c      end do

*     Get the information on the reference site (assumed to be the
*     same reference satellite for all data), and the ambiquity and
*     satellite number
c      pf_prn(num_res) = ow_des(2,dd_des(4,1))
c      pf_site(num_res) = ow_des(1,dd_des(1,1))
c      pf_ref(num_res) = 'R'
    
      do k = 1, num_res
         do j = 1, num_chan_se(pf_site,ep)
            if( ctop_cse(j, pf_site,ep).eq.pf_prn(k) ) then
               pf_elev(k) = elev_cse(j,pf_site,ep)
               pf_az(k) = az_cse(j,pf_site,ep)
               pf_amb(k) = amb_point_cse(j,pf_site,ep)
               pf_bf5(k) = bf_ents(5,pf_amb(k))
            end if
         end do
      end do
 

****  Now generate the output format that we need (based on number of
*     data types)
      call gen_res_format(num_dtype,'data', format)

*     Now output
      curr_mjd = ref_start + ((ep-1)*usr_interval +
     .  	 sec_offset(pf_site,ep))/86400.d0
      call mjd_to_ymdhms(curr_mjd, date, sectag)
      call ymd_to_doy( date, doy)
      sec_in_day = date(4)*3600.d0 + 
     .  	   date(5)*60.d0 + sectag

      fday = doy*1.d0 + date(4)/24.d0 + date(5)/1440.d0 +
     .  	sectag/86400.d0 
      if( write_res ) then
         is = kine_to_site(ks)
         do k = 1, num_res
             luo = 120+(is-1)*max_sat+pf_prn(k)
*            OK write the values out
             write(luo,format) date,sectag, 
     .  	  (pf_res(j,k),j=1,num_dtype),
     .  	   pf_elev(k), pf_az(k), pf_prn(k), 
     .  	   pf_amb(k), pf_ref(k),
     .  	   pf_bf5(k),  fday, ep
              if( ep.ge.debug_start .and. ep.le. debug_end) 
     .        write(*,format) date,sectag, 
     .  	  (pf_res(j,k),j=1,num_dtype),
     .  	   pf_elev(k), pf_az(k), pf_prn(k), 
     .  	   pf_amb(k), pf_ref(k),
     .  	   luo,  fday, ep
         end do
      end if

****  Now accumulate the statistics 
      do k = 1,num_res
****     Get elevation dependence
         eb = pf_elev(k)/5+1
         if( eb.lt.1 ) eb = 1
         if( eb.gt.18 ) eb = 18
         pn = UtoL(pf_prn(k))
         do j = 1, num_dtype
            sum_res(j,ks,pn) = sum_res(j,ks,pn)+pf_res(j,k)
            sum_var(j,ks,pn) = sum_var(j,ks,pn)+pf_res(j,k)**2
            sum_num(j,ks,pn) = sum_num(j,ks,pn)+1

*           Summ elevation dependence
            elv_res(j,ks,eb) = elv_res(j,ks,eb)+pf_res(j,k)
            elv_var(j,ks,eb) = elv_var(j,ks,eb)+pf_res(j,k)**2
            elv_num(j,ks,eb) = elv_num(j,ks,eb)+1
         end do
      end do

            

****  Thats all
      return
      end

CTITLE GEN_RES_FORMAT
     
      subroutine gen_res_format(num_dtype,type, format)

      implicit none

*     Generates the output format either for the data (type data) or
*     for the header (type label)

* PASSED VARIABLES
* num_dtype  -- Number of data types

      integer*4 num_dtype

* type  -- Type of format to generate
* format  -- The format string

      character*(*) type, format

* LOCAL VARIABLES
* j    -- Loop counter

      integer*4 j

****  See the type we have
      if( type.eq.'data' ) then
          write(format,120) '(i5,4i3,1x,f7.3,', num_dtype,
     .        'F12.2,1x,2f9.3,2x,i3.3,1x, i4,a1,1x,i4,1x,F15.9,1x,i6)'
 120      format(a,i2,a)
      end if     
      
      if( type.eq.'label' ) then
          write(format,140) '*YYYY MM DD HR MN   Sec   ',
     .      ('    Res (mm)',j=1,num_dtype),
     .      '    Elev    Azim    PRN   BF  Fixd     Fract DOY    Epoch' 
 140      format(40a)
      end if


****  Thats all
      return
      end
      
CTITLE GET_NUM_DTYPE

      subroutine get_num_dtype(type, num_dtype, dtypes)

      implicit none

*     Routine to return the number of data types in an analyis

* PASSED VARIABLES
* num_dtype -- Number of data types
      integer*4 num_dtype

* type  -- String with type of analysis
      character*(*) type
      character*(*) dtypes(*)

****  OK, count the types
      num_dtype = 0
      if( index(type,'L1').gt.0 ) then
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'L1'
      end if
      if( index(type,'L2').gt.0 ) then
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'L2'
      end if
      if( index(type,'LC').gt.0 ) then 
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'LC'
      end if 
      if( index(type,'P1').gt.0 ) then 
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'P1'
      end if
      if( index(type,'P2').gt.0 ) then
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'P2'
      end if
      if( index(type,'PC').gt.0 ) then
          num_dtype = num_dtype + 1
          dtypes(num_dtype) = 'PC'
      end if

      return
      end

CTITLE clear_stats

      subroutine clear_stats

      implicit none

*     Routine to clear statistics for each analysis type

      include 'track_com.h'

      integer*4 i,j,k

      do i = 1, max_dtype
         do j = 1, max_kine
            do k = 1, max_sat
               sum_res(i,j,k) = 0.d0
               sum_var(i,j,k) = 0.d0
               sum_num(i,j,k) = 0
            end do
            do k = 1, 18
               elv_res(i,j,k) = 0.d0
               elv_var(i,j,k) = 0.d0
               elv_num(i,j,k) = 0
            end do
         end do
      end do

      return
      end


CTITLE phase_stats

      subroutine phase_stats(un, type)

      implicit none

*     Routine to finish up and write the phase statistics from
*     the run

      include 'track_com.h'

* PASSED VARIABKES
      integer*4 un   ! Unit numebr for output
      character*(*) type  ! Analysis type

* LOCAL VARIABLES
      integer*4 i,j,k
      integer*4 pn(max_sat)    ! PRN's with data
      integer*4 pnn            ! number of PRNs with adta

      real*8 res_mean(max_dtype,max_sat), 
     .       res_rms(max_dtype,max_sat),
     .       sum_all_var(max_dtype),
     .       all_rms(max_dtype)

      integer*4 res_num(max_dtype,max_sat),
     .          sum_all_num(max_dtype)

      character*2 dtypes(max_dtype)
    

***** Write heaser
      write(un,120) type
 120  format(/,'TRACK Residual statistics for Analysis ',a)
      
      call get_num_dtype(type,num_dtype,dtypes)

****  Finishup statistics
      do j = 1, num_kine
         do i = 1, num_dtype
            pnn = 0
            sum_all_var(i) = 0
            sum_all_num(i) = 0
            do k = 1, num_prn
               if( sum_num(i,j,k).gt.1 ) then
                  pnn = pnn + 1
                  res_mean(i,pnn) = sum_res(i,j,k)/sum_num(i,j,k)
                  res_rms(i,pnn) = sqrt((sum_var(i,j,k) -
     .                             res_mean(i,pnn)**2*sum_num(i,j,k))/
     .                             (sum_num(i,j,k)-1))
                  res_num(i,pnn) = sum_num(i,j,k)
C                 pn(pnn) = k
                  pn(pnn) = prn_used(k)
               
                  sum_all_var(i) = sum_all_var(i) + sum_var(i,j,k)
                  sum_all_num(i) = sum_all_num(i) + sum_num(i,j,k) 
               end if
            end do
            all_rms(i) = sqrt(sum_all_var(i)/sum_all_num(i))
         end do

****     Write out values
         write(un, 210) site_names(kine_to_site(j)),type,
     .                  (pn(i),i=1,pnn)
 210     format('STATISTICS BY PRN for site ',a,' DataTypes ',a,/,
     .          'TYPE     Site  DT    ALL',100(1x,I5.2))
         do i = 1, num_dtype
            write(un,220) site_names(kine_to_site(j)),dtypes(i),
     .                    0.d0,(res_mean(i,k),k=1,pnn)
 220        format('PMEAN ',A8,1x,a2,F7.1,100(F6.1))
         end do
         do i = 1, num_dtype
            write(un,240) site_names(kine_to_site(j)),dtypes(i),
     .                    all_rms(i),(res_rms(i,k),k=1,pnn)
 240        format('PRMS  ',A8,1x,a2,F7.1,100(F6.1))
         end do
         do i = 1, num_dtype
            write(un,260) site_names(kine_to_site(j)),dtypes(i), 
     .                    sum_all_num(i),(res_num(i,k),k=1,pnn)
 260        format('PNUM  ',A8,1x,a2,1x,I6,100(1x,i5))
         end do
      end do

****  Now do statistics by elevation angles
      do j = 1, num_kine
         do i = 1, num_dtype
            pnn = 0
            sum_all_var(i) = 0
            sum_all_num(i) = 0
            do k = 1, 18
               if( elv_num(i,j,k).gt.1 ) then
                  res_mean(i,k) = elv_res(i,j,k)/elv_num(i,j,k)
                  res_rms(i,k) = sqrt((elv_var(i,j,k) -
     .                             res_mean(i,k)**2*elv_num(i,j,k))/
     .                             (elv_num(i,j,k)-1))
                  res_num(i,k) = elv_num(i,j,k)
               else
                  res_mean(i,k) = 0.0
		  res_rms(i,k) = 0.0
                  res_num(i,k) = 0
               end if
            end do
         end do

****     Write out values
         write(un, 310) site_names(kine_to_site(j)),type,
     .                  ((i-1)*5,i*5,i=1,18)
 310     format('STATISTICS BY Elevation for site ',a,' DataTypes ',a,/,
     .          'TYPE     Site  DT',5x,18(1x,i2.2,'-',i2.2))
         do i = 1, num_dtype
            write(un,320) site_names(kine_to_site(j)),dtypes(i),
     .                   (res_mean(i,k),k=1,18)
 320        format('EMEAN ',A8,1x,a2,5x,100(F6.1))
         end do
         do i = 1, num_dtype
            write(un,340) site_names(kine_to_site(j)),dtypes(i),
     .                   (res_rms(i,k),k=1,18)
 340        format('ERMS  ',A8,1x,a2,5x,100(F6.1))
         end do
         do i = 1, num_dtype
            write(un,360) site_names(kine_to_site(j)),dtypes(i),
     .                    (res_num(i,k),k=1,18)
 360        format('ENUM  ',A8,1x,a2,5x,100(1x,i5))
         end do
      end do

****  Thata all
      return
      end

CTITLE UtoL 

      integer*4 function UtoL( prn )

      implicit none  

*     Routine to returm the list number from the PRN_USED array.
*     (similar to PtoL which uses the list of PRNs in the SP3 file)

      include 'track_com.h'

* PASSED 
      integer*4 prn   ! PRN number offset for the GNSS constellation.

* LOCAL
      integer*4 i     ! Loop counter

      UtoL = -1
      do i = 1, num_prn
         if( prn.eq.prn_used(i) ) then
           UtoL = i
           exit
         end if
      end do

****  Thats all
      return
      end


