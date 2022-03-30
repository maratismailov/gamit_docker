CTITLE MUL_PMU_PRED

      subroutine mul_pmu_pred( cov_parm, sol_parm )

      implicit none 

*     Routine to update the process noise to step forward or backward to
*     next block of multi-pmu estimates.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*   sol_parm(num_glb_parn)              - Solution vector
 
      real*8 cov_parm(num_glb_parn,num_glb_parn),
     .    sol_parm(num_glb_parn)

* LOCAL VARIABLES

*  i,j,k,l  -- Loop variables
*  ni(2)  -- Parameter pointers to offset and rate
*  of     -- Offset in index to account for the time interval between
*            current day and last multi-day PMU value (approximate for
*            continuous data)

      integer*4 i,j,k, ni(2), of

*  pa(2)  -- Parameter numbers for j offset and rate EOP
*  pb(2)  -- Parameter numbers for k offser and rate EOP

      integer*4 pa(2), pb(2)
 
*  cov_init(2,2) -- Initial variance covariance of offset and rate
*  cov_irw(2,3)  -- Variance of offset and rate and covariance between
*                   them
*  mar_rw, mar_irw -- Markov process noise to use.

      real*8 cov_init(2,2), cov_irw(2,2)

      real*4 mar_rw, mar_irw
 
      logical kbit 

****  See if we need to update the process noise (ie., do we have a
*     new block of mul_pmu entries), or 
*     if we are treating the values as independent, or
*     there is no need to do because it was done in init_cov.
      if( .not. kbit(mul_pmu_opt,32) .or.
     .          kbit(mul_pmu_opt, 1) .or.
     .    .not. kbit(mul_pmu_opt,31)     ) RETURN

*     If we are not doing multiday polar motion return also
      if( num_mul_pmu.eq.0 ) RETURN
 
****  Now see if should propagate the system forward or backward
*     Depending on the direction, use either the first or last 
*     offset and rate sigma to set the apriori  

      of = nint(abs((deltat*365.25d0)/spacing_mul_pmu)) 

      do i = 1,3   !  Loop over X Y pole and UT1
*        Based on type of EOP get the statistics that we need.
         if( i.eq.1 ) then
             mar_rw        = mar_wob(1)
             mar_irw       = mar_wob(3)
         else if( i.eq.2 ) then
             mar_rw        = mar_wob(2)
             mar_irw       = mar_wob(4)
         else if ( i.eq.3 ) then
             mar_rw        = mar_ut1(1)
             mar_irw       = mar_ut1(2)
         end if

         if ( deltat.gt.0 ) then
             ni(1) = parn_mul_pmu(1,i, num_mul_pmu)
             ni(2) = parn_mul_pmu(2,i, num_mul_pmu)
*            Propagate the offset and rate offsets into the state vector
             do j = 1, num_mul_pmu
                pa(1) = parn_mul_pmu(1,i,j)
                pa(2) = parn_mul_pmu(2,i,j)
                if( ni(2).gt.0 ) then
                    sol_parm(pa(1)) = sol_parm(ni(1)) + 
     .                     (deltat*365.25d0+(j-1)*spacing_mul_pmu)*
     .                      sol_parm(ni(2))
                    sol_parm(pa(2)) = sol_parm(ni(2))
                else
                    sol_parm(pa(1)) = sol_parm(ni(1))
                end if
             end do 
         else 
             ni(1) = parn_mul_pmu(1,i, 1)
             ni(2) = parn_mul_pmu(2,i, 1)
*            Propagate the offset and rate offsets into the state vector.
*            (This time in reverse order)
             do j =  num_mul_pmu, 1, -1 
                pa(1) = parn_mul_pmu(1,i,j)
                pa(2) = parn_mul_pmu(2,i,j)
                if( ni(2).gt.0 ) then
                    sol_parm(pa(1)) = sol_parm(ni(1)) + 
     .                     (deltat*365.25d0-(j-1)*spacing_mul_pmu)*
     .                      sol_parm(ni(2))
                else
                    sol_parm(pa(1)) = sol_parm(ni(1))
                end if
             end do
                
         end if
         do j = 1, 2
            do k = 1, 2
               if( ni(j).gt.0 .and. ni(k).gt.0 ) then
                   cov_init(j,k) = cov_parm(ni(j),ni(k))
               else
                   cov_init(j,k) = 0.d0
               end if
            end do
         end do
                
*        Now compute and add the process noise.
         do j = 1, num_mul_pmu
            do k = j, num_mul_pmu
*              Get the parameter numbers that we need for the
*              offset and rate terms
               pa(1) = parn_mul_pmu(1,i,j)
               pa(2) = parn_mul_pmu(2,i,j)  
               pb(1) = parn_mul_pmu(1,i,k)
               pb(2) = parn_mul_pmu(2,i,k) 

*              Now get the covariance matrix contribution for these
*              two values separated in time.
               if( deltat.gt.0 ) then
                   call irw_cov(cov_irw, cov_init,  mar_rw, 
     .                          mar_irw, j+of, k+of,  
     .                          spacing_mul_pmu) 
               else
                   call irw_cov(cov_irw, cov_init,  mar_rw, 
     .                          mar_irw, num_mul_pmu-j+of, 
     .                          num_mul_pmu-k+of,  
     .                          spacing_mul_pmu)
               end if 

*              Now add to covariance matrix
               if( pa(1).gt.0 .and. pb(1).gt.0 )
     .            cov_parm(pa(1),pb(1)) = cov_parm(pa(1),pb(1))
     .                                    + cov_irw(1,1)
               if( pa(1).gt.0 .and. pb(2).gt.0 )
     .            cov_parm(pa(1),pb(2)) = cov_parm(pa(1),pb(2))
     .                                    + cov_irw(1,2)
               if( pa(2).gt.0 .and. pb(1).gt.0 )
     .            cov_parm(pa(2),pb(1)) = cov_parm(pa(2),pb(1))
     .                                    + cov_irw(2,1)
               if( pa(2).gt.0 .and. pb(2).gt.0 )
     .            cov_parm(pa(2),pb(2)) = cov_parm(pa(2),pb(2))
     .                                    + cov_irw(2,2)

*              Save into the lower block as well
               if( j.ne.k ) then
                   if( pa(1).gt.0 .and. pb(1).gt.0 )
     .                cov_parm(pb(1),pa(1)) = cov_parm(pa(1),pb(1))
                   if( pa(1).gt.0 .and. pb(2).gt.0 )
     .                cov_parm(pb(2),pa(1)) = cov_parm(pa(1),pb(2))
                   if( pa(2).gt.0 .and. pb(1).gt.0 )
     .                cov_parm(pb(1),pa(2)) = cov_parm(pa(2),pb(1))
                   if( pa(2).gt.0 .and. pb(2).gt.0 )
     .                cov_parm(pb(2),pa(2)) = cov_parm(pa(2),pb(2))
               end if
            end do 
         end do 

      end do 

      call sbit(mul_pmu_opt,31,0 )

****  Thats all
      return 
      end 





      

