CTITLE MERGE_NONSEC

      subroutine merge_nonsec(save_num_nonsec, save_param_nonsec,
     .                  save_val_nonsec)

      implicit none 

*     Routine to merge the non-secular terms from the orginal values
*     used in the globk solution and the values to be used as
*     apriori's in the glorg solution.  Really only valid for 
*     single epoch globk combinations

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'

* PASSED VARIABLES
* save_num_nonsec -- Original number of non-secular terms
* save_param_nonsec(2,max_nonsec) -- Parameter pointers
      integer*4 save_num_nonsec, save_param_nonsec(2,max_nonsec)
* save_val_nonsec(8,max_nonsec) -- Actual parameters of non-secular 
*     terms
      real*8 save_val_nonsec(8,max_nonsec) 

* LOCAL VARIABLES
* i,j,k -- Loop counters
* is    -- Site number
* norg  -- Original number of nonsecular terms (used so that do
*          loop end variable is not changed)
* jel   -- Parameter pointer

      integer*4 i,j,k, is, norg, jel

* dsol_new(3) -- Non-secular change for new values
* dsol_old(3) -- Non-secular change for old vales (The old value is
*          removed to be replaced with new value)
* old_log(3), new_log(3) -- old and new values of the log terms
*     needed to adjust the parameter estimates when apriori changed.

      real*8 dsol_new(3), dsol_old(3), old_log(3), new_log(3)

* duplicate -- Set true if duplicate between new and old list
* test      -- Set true if one value duplicates in merge
* remove    -- Set true if non-secular term should be removed 
*              becuase it can not updated.
* same      -- Set true if duplicate non-secular terms are same

      logical duplicate, test, remove , kbit, same

* fatal  -- Fatal error message if too many terms
* mess   -- Warning message 
      character*128 fatal, mess

      character*16 non_sec_type

****  OK: First scan over the list of non-secular terms just read
*     and see which ones can actually be changed.  The changable 
*     ones are log terms when log terms are estimated (and not too
*     tightly constrained).
      i = 0
      do while ( i.lt. num_nonsec )
         i = i + 1
*        See if this is type 4 (log term) and that we have an 
*        estimate to go with it.
         remove = .true. 
         if( param_nonsec(2,i).eq.4 ) then
*            OK log term, see if estimate go with it (assumes 
*            that if one component is estimated, all components
*            are estimated).
             is = param_nonsec(1,i)
             jel = parn_log(1,is)
             if( jel.gt.0 ) remove = .false.
         end if

*        See if we should remove
         if( remove ) then

* NOD TAH 070924: Add two checks before reporting:
*            (1) Is this site actually used and
*            (2) Are we trying to change the value (Do later)

             same = .false. 
             test = .false.
             do j = 1,i-1
                if( param_nonsec(1,i).eq. param_nonsec(1,j) .and.
     .              param_nonsec(2,i).eq. param_nonsec(2,j) ) then
*                   See if values match
                    test = .true.
                    do k = 1,8
                        if( apr_val_nonsec(k,i).ne. 
     .                      apr_val_nonsec(k,j) ) then 
                            test = .false.
                        endif
                    enddo
                 endif 
*                We have found a matching entry so no need to report.
                 if ( test ) then
                    same = .true.
                    exit
                 endif
              end do
* MOD TAH 190531: Now see if we can find this in saved non-secular
*             terms 
              if (.not.same ) then 
                 do j = 1,save_num_nonsec
                   if( param_nonsec(1,i)  .eq. 
     .                 save_param_nonsec(1,j) .and.
     .                 param_nonsec(2,i).eq. 
     .                 save_param_nonsec(2,j) ) then
*                      See if values match
                       test = .true.
                       do k = 1,8
                           if( apr_val_nonsec(k,i).ne. 
     .                         save_val_nonsec(k,j) ) then 
                               test = .false.
                           endif
                       enddo
                    endif 

*                   We have found a matching entry so no need to report.
                    if ( test ) then
                       same = .true.
                       exit
                    endif
                 end do
              endif
 
* MOD TAH 190513: Only report if no duplicate entry found
             if( kbit(guse_site,is) .and. .not. same ) then
                 non_sec_type = 'EQ_LOG'
                 if( param_nonsec(2,i).eq.2 ) non_sec_type = 'PERIODIC'
                 if( param_nonsec(2,i).eq.4 ) non_sec_type = 'EQ_EXP'

                 write(mess,120) trim(non_sec_type), gsite_names(is)
 120             format('Extended ',a,' for ',a8,' not estimated, ',
     .                  'no change to apriori allowed in glorg')
                  call report_stat('WARNING','GLORG','merge_nonsec',
     .                     ' ',mess,0)
             endif

             do j = i+1, num_nonsec
        	 do k = 1,2
        	    param_nonsec(k,j-1) = param_nonsec(k,j)
        	 end do
        	 do k = 1,8
        	    apr_val_nonsec(k,j-1) = apr_val_nonsec(k,j)
        	 end do
             end do
             i = i - 1
             num_nonsec = num_nonsec - 1
          end if
      end do  

      do i = 1, num_nonsec

*        OK, Loop over the ones already used.
         j = 1
         duplicate = .false.
         do while ( j.lt.save_num_nonsec .and. 
     .              .not.duplicate )
            j = j + 1 
            duplicate = .true.
            do k = 1,2
               if( param_nonsec(k,i).ne.save_param_nonsec(k,j))
     .             duplicate = .false.
            end do

*           Check if a duplicate
            if( duplicate .and. param_nonsec(1,j).eq.4 ) then
*              OK: This term was used in globk solution.  Compute the
*              change needed for the new value and save in the parameter
*              update array
               is = param_nonsec(1,i)
 
               call eval_nonsec(is, gepoch_out,1,param_nonsec(1,i),
     .              apr_val_nonsec(1,i), dsol_new,0)
               call eval_nonsec(is, gepoch_out,1,
     .              save_param_nonsec(1,j),save_val_nonsec(1,j), 
     .              dsol_old,0)

*              Compute the change in the parameter estimates
               do k = 1,3
                  if( parn_site(k,1,is).ne.0 ) then
                      jel = parn_site(k,1,is)
                      parm_change(jel) = parm_change(jel) +
     .                   dsol_new(k) - dsol_old(k)
                  end if
               end do

*              Now see if we have an explicit estimate of log term 
*              associated with this non-secular term
               if( save_param_nonsec(2,j).eq.4 ) then
                   call nonsec_convert('TONEU',1,save_val_nonsec(3,j),
     .                   old_log, apr_val_site(1,1,is))
                   call nonsec_convert('TONEU',1,apr_val_nonsec(3,i),
     .                   new_log, apr_val_site(1,1,is))
                   do k =1,3
                      if( parn_log(k,is).ne.0 ) then
                          jel = parn_log(k,is)
                          parm_change(jel) = parm_change(jel) +
     .                        new_log(k) - old_log(k)
                      endif
                   enddo
               endif  

            end if
         end do
****     If no duplicate was found then this is a new non-secular term that
*        was not used in the globk run.  The only time this is valid to do
*        is for log terms and only then if logs are estimated. (The only
*        updates that get through are for log terms with estimates so we
*        can assume this is the case.
         if( .not.duplicate ) then
            is = param_nonsec(1,i)
            call get_nonlog(is, new_log)
            do k =1,3
               if( parn_log(k,is).ne.0 ) then
            	   jel = parn_log(k,is)
            	   parm_change(jel) = parm_change(jel) + new_log(k) 
               endif
            enddo
         end if 
      end do

****  OK, now copy over the non-secular terms that were used in
*     globk but have not been updated.
      norg = num_nonsec
      do j = 1, save_num_nonsec
         duplicate = .false.
         do i = 1, norg
            test = .true.
            do k = 1,2
               if( param_nonsec(k,i).ne.save_param_nonsec(k,j))
     .             test = .false.
               if( apr_val_nonsec(k,i).ne.save_val_nonsec(k,j))
     .             test = .false.
            end do
            if( test ) duplicate = .true.
         end do

*        Check if a duplicate
         if( .not.duplicate ) then
            num_nonsec = num_nonsec + 1
            if( num_nonsec.gt. max_nonsec) then
                write(fatal,220) max_nonsec
 220            format('Too many nonsecular terms during merge. ',
     .                  'Max allowed is ',i6)
                call report_stat('FATAL','GLOBK',
     .                  'merge_nonsec',' ',fatal,0)
            end if
            do k = 1,2
               param_nonsec(k,num_nonsec) = save_param_nonsec(k,j)
            end do
            do k = 1,8
               apr_val_nonsec(k,num_nonsec) =  save_val_nonsec(k,j)
            end do
         end if
      end do

****  Thats all
      return
      end
    

