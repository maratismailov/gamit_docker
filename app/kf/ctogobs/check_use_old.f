CTITLE Check_used

      logical function check_used(refbf, last, bf_baseline)

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


* LOCAL VARIABLES
      integer*4 i,j      ! looop counter
     .,  lv              ! Satellite number
     .,  pair(2), next(2)   ! Pairs of sites in baselines
     .,  nblen           ! Number of baselines
     .,  s1, s2          ! The two sites being considered

      logical kbit       ! Checks bit status
     .,   used           ! Logical which set status


****  First see if oneways have been used all previously
      used = .false.
      if( kbit(bf_table(5,refbf(1)),5) .and.
     .    kbit(bf_table(5,refbf(2)),5) .and.
     .    kbit(bf_table(5,last(1)), 5) .and.
     .    kbit(bf_table(5,last(2)), 5) ) then
          used = .true.
      end if

****  See if it is clear that it is not used
      if( .not.used ) then
          check_used = used
          RETURN
      endif
      return
      s1 = bf_table(2,refbf(1))
      s2 = bf_table(2,refbf(2))

***** If some have not been used before, then OK. If they have
*     all been used, check to see if there is common site
*     between them.  If there isn't then we set this combination
*     as not used. Eg. baseline 1-2 and 3-4, 2-3 is still OK even
*     though all oneways have been used before
      nblen = (num_cfiles)*(num_cfiles-1)/2
      if( used ) then
          used = .false.
*         OK, see if common site
          lv = bf_table(3,last(1))
*         Loop over the baselines getting those that
*         use the reference site.  When one is found
*         loop again over the second site to see if
*         there is a common site (in which case this 
*         is dependent baseline and should not be used.
          do i = 1, nblen
              if( kbit(bf_baseline(i),lv) ) then
                 call ent_bf(i,num_cfiles, pair)
*                See if match on one of our sites
                 if( pair(1).eq.bf_table(2,refbf(1)) .or.
     .               pair(2).eq.bf_table(2,refbf(1))) then

*                    see if the second site in this baseline
*                    has been used with the other site
                     do j = 1,nblen
                        if( kbit(bf_baseline(j),lv) .and. 
     .                      i.ne.j) then
                            call ent_bf(j,num_cfiles, next)
*                           Do this baseline match our second
*                           site?
                            if( next(2).eq.bf_table(2,refbf(2)) .or.
     .                          next(1).eq.bf_table(2,refbf(2)) ) then
*                               OK: We have another baseline that
*                               used the second site.  See if the
*                               other site is these baselines is
*                               if the same
                                if( next(1).eq.pair(1) .or. 
     .                              next(1).eq.pair(2) ) used = .true.
                                if( next(2).eq.pair(1) .or. 
     .                              next(2).eq.pair(2) ) used = .true.
*                               If used is true, I could return at
*                               this point
                                if( used ) then
                                    check_used = used
                                    RETURN
                                end if
                            end if
                        end if
                     end do
                 end if
              end if
          end do
      end if

****  Set the status and return
      check_used = used
      return
      end   


 
