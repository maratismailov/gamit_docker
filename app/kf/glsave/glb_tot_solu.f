CTITLE GLB_TOT_SOLUTION

      subroutine glb_tot_solution ( sol_parm )

      implicit none  
 
*     Routine to add the apriori values into the estimated parameter
*     adjustments so that we save the total value in the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k   - Loop counters
*   np      - Parameter number
 
      integer*4 i,j,k, np
 
*   sol_parm(1) - The solution vector from globk
*   dt          - Time difference between experiment and
*               - reference epoch of solution
*   dnon_sec(3) - Contributions from the non-secular terms
 
      real*8 sol_parm(1), dt, dnon_sec(3)

***** Loop over site positions
 
      do i = 1, gnum_sites
          dt = (gepoch_out - site_epoch(i))/365.25
* MOD TAH 0509927: Compute the non-secular terms
          call eval_nonsec(i, gepoch_out, num_nonsec, param_nonsec,
     .                    apr_val_nonsec, dnon_sec,0)
          do j = 1, 3
              np = parn_site(j,1,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_site(j,1,i)
     .                         + apr_val_site(j,2,i)*dt
     .                         + dnon_sec(j)
              end if
          end do
      end do
 
****  Do site velocities
 
      do i = 1, gnum_sites
          do j = 1, 3
              np = parn_site(j,2,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_site(j,2,i)
              end if
          end do
      end do

* MOD TAH 0301615: Generate full LOG solution
      do i = 1, gnum_sites
         call get_nonlog(i,apr_val_log(1,i))
         do j = 1, 3
            np = parn_log(j,i)
            if( np.gt.0 ) then
                sol_parm(np) = sol_parm(np) + apr_val_log(j,i)
            end if
         end do
      end do 

* MOD TAH 040703: Atmospheric delays
      do i = 1, gnum_sites
         np = parn_atm(i)
         if( np.gt.0 ) then
             sol_parm(np) = sol_parm(np) + apr_val_atm(i)
         end if
      end do
 
*     Axis offset
      do i = 1, gnum_sites
          dt = (gepoch_out - axo_epoch(i))/365.25
*                             ! Loop over value and rate
          do j = 1,2
              np = parn_axo(2,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_axo(j,i)
                  if( j.eq.2 ) then
                        sol_parm(np) = sol_parm(np)
     .                               + apr_val_axo(2,i)*dt
                  end if
              end if
          end do
      end do
 
***** Loop over source positions
      do i = 1, gnum_sources
          dt = (gepoch_out - source_epoch(i))/365.25
          do j = 1, 2
              np = parn_source(j,1,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_source(j,1,i)
     .                         + apr_val_source(j,2,i)*dt
              end if
          end do
      end do
 
*     Loop over source rates
      do i = 1, gnum_sources
          do j = 1, 2
              np = parn_source(j,2,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_source(j,2,i)
              end if
          end do
      end do

*     Loop over the satellite ephermides parameters
      do i = 1, gnum_svs
         do j = 1, max_svs_elem
            np = parn_svs(j,i)
            if( np.gt.0 ) then
                sol_parm(np) = sol_parm(np) + apr_val_svs(j,i)
            end if
         end do
      end do
 
***** Polar motion/UT1 values
*                     ! X/Y pole
      do i = 1,4
          np = parn_wob(i)     
          if( np.gt.0 ) then
              sol_parm(np) = sol_parm(np) + apr_val_wob(i)
* MOD TAH 030130: Increase index test (covers rates as well)
              if( i.le.4 ) then
                  sol_parm(np) = sol_parm(np) + gwob_apr(i)
              end if
          end if
      end do
 
*                     ! UT1 model
      do i = 1,2
          np = parn_UT1(i)
          if( np.gt.0 ) then
              sol_parm(np) = sol_parm(np) + apr_val_ut1(i)
              if( i.le.2 ) then
                  sol_parm(np) = sol_parm(np) + gut1_apr(i)
              end if
          end if
      end do

* MOD TAH 981104: Add back the apriori for multi-day polar motion
      do i = 1,3
         do j = 1, 2
            do k = 1, num_mul_pmu
               if( parn_mul_pmu(j,i,k).ne.0 ) then
                   np = parn_mul_pmu(j,i,k)
                   sol_parm(np) = sol_parm(np) + apr_val_mul_pmu(j,i,k)
               end if
            end do
         end do
      end do

 
***** Nutation angle
*                     ! DPSI, Deps and seasonal model
      do i = 1,8
          np = parn_nut_ang(i)
          if( np.gt.0 ) then
              sol_parm(np) = sol_parm(np) + apr_val_nut_ang(i)
              if( i.le.2 ) then
                  sol_parm(np) = sol_parm(np) + gnut_ang_apr(i)
              end if
          end if
      end do
 
***** Nutation series coefficients
      do i = 1, max_nut_coeff
*                             ! In and out of phase
          do j = 1,2
              np = parn_nut_coeff(j,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_nut_coeff(j,i)
              end if
          end do
      end do
 
***** Tides (standard model)
      if( glb_glb_tides ) then
          k = 1
      else
          k = gnum_sites
      end if
 
*     Save values (Check order of the apriori values)
      do i = 1, k
*                                         ! h and l values
          do j = 1, 2
              np = parn_tid(j,i)
              if( np.gt.0 ) then
                  sol_parm(np) = sol_parm(np) + apr_val_tid(j+1,i)
              end if
          end do
 
*                                         ! Lag angle
          if( parn_tid(3,i).gt.0 ) then
              sol_parm(parn_tid(3,i)) = sol_parm(parn_tid(3,i)) +
     .                                  apr_val_tid(1,i)
          end if
      end do
 
***** Gamma
      if( parn_gamma.gt.0 ) then
          sol_parm(parn_gamma) = sol_parm(parn_gamma) + apr_val_gamma
      end if
 
***** Thats all
      return
      end
 
