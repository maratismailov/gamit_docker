CTITLE Check_used

      logical function check_used(refbf, last, bf_baseline, bfj, nbl)

      implicit none

*     Routine to check to see if the one ways have been used with
*     a common station previously.  If they have then they are 
*     marked as used, other wise used is set false.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
      integer*4 refbf(2) ! Entries in the bf_table for the reference
                         ! satellite
     .,  last(2)  ! Entries in bf_table for the new oneways
         	  ! being considered
     .,  bf_baseline((max_cfiles*(max_cfiles+1))/2) ! Baselines that
                  !  have been used.  Bit set for each satellite as it is used

     .,   bfj((max_cfiles*(max_cfiles+1))/2) ! Sort list of baseline lengths
     .,   nbl     ! Baseline number from the bfj list being checked.

* LOCAL VARIABLES
      integer*4 i,j      ! looop counter
     .,  pair(2)         ! Pairs of sites in baselines
     .,  s1, s2          ! The two sites being considered
     .,  bls((max_cfiles*(max_cfiles+1))/2) ! Stack of baselines that need to be further
                         ! checked
     .,  rfs((max_cfiles*(max_cfiles+1))/2) ! Stack of reference sites for the baselines.
                         ! When one reference site is found, it is pushed on the stack and
                         ! the second site in the baseline becomes the reference.
     .,  nb              ! Number of current baseline being processed.
     .,  ns              ! Number of entries on current stack
     .,  rs              ! Current reference site for baseline

      logical done       ! Logical set true when all baselines have been scanned or a 
                         ! linkage of baselines between sites has been found.
     .,   s1_used, s2_used  ! Set true if the sites have been used in an eariler
                         ! baseline.
     .,   onstack        ! Set true if the current baseline we are condsidering is already
                         ! on the stack.


***** First check the baselines we have processed already and
*     see if one or both of these sites have never been used 
*     before.  If one of them has not been used, then clearly
*     independent.
      s1 = bf_table(2,refbf(1))   ! S1 and S2 are the two sites in the baseline
      s2 = bf_table(2,refbf(2))
      s1_used = .false.
      s2_used = .false.

      do i = 1, nbl - 1  ! Test up to the baseline before this one
          call enf_bf(bfj(i), num_cfiles, pair)
          if( pair(1).eq.s1 .or. pair(2).eq.s1 ) s1_used = .true.
          if( pair(1).eq.s2 .or. pair(2).eq.s2 ) s2_used = .true.
      end do

*     If either site has not been used, then this is indepednent baseline
      if( .not.s1_used .or. .not.s2_used ) then
          check_used = .false.
          RETURN
      end if

****  OK: Now the tricky part.  Both sites have been used but there may not
*     be squence of baselines that connect them.  So now we need to search
*     down each baseline tree starting from S1 and see if we can get to S2
*     through the baselines.  The algoritm is basically recursive but is implemented
*     through a stack approach
      rs = s1      ! Initial reference site (first in baseline)
      done = .false.
      nb = 0
      ns = 0    ! Number of entries on stack
      check_used = .false.

      do while ( .not. done )
*        Move to the next baseline
         nb = nb + 1
         if( nb.lt.nbl ) then    ! OK, if this is before the baseline we are considering
*           Check to see if this baseline is already on the stack, if it is not then
*           continue checking 
            onstack = .false.
            do j = 1, ns
               if( nb.eq.bls(i) ) onstack = .true.
            end do
            if ( .not. onstack ) then
*               See if one of the sites matchs the current reference sites
                call enf_bf(bfj(nb),num_cfiles, pair)
*               Switch the baseline order is necessary to make pair(1) be the reference site
                if( pair(2).eq.rs ) then
                    pair(2) = pair(1)
                    pair(1) = rs
                end if
                if( pair(1).eq.rs ) then    ! First site matches
*                   Check to see if second site is actually the end of the baseline chain
                    if( pair(2).eq.s2 ) then
*                       Connection has been made to second site so this is not
*                       an independent baseline.  Get out of test
                        check_used = .true.
                        done = .true.
                    else     ! Site does not match end yet, so push the entry onto the
                             ! stack and start searching baselines from this point
                        ns = ns + 1
                        bls(ns) = nb
                        rfs(ns) = rs
                        nb = 0     ! Start searching baseline again
                        rs = pair(2)
                    end if
                end if
            end if       ! Not on the stack already
         else    ! We have exceeded the baseline count, so pop the stack and follow
                 ! another path to see where it leads
            nb = bls(ns)
            rs = rfs(ns)
            ns = ns - 1
*           Have we run out of stack?
            if( ns.eq.0 ) done = .true.
         end if
      end do
*
*     Scanned all used baselines or we done a path
      return
      end

    


