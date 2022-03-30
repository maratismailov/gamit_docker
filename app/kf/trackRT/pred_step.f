      subroutine pred_step(ep)

      implicit none

*     Routine to step the KF state vector and covarance matrix
*     forward one step and add process noise

      include '../includes/const_param.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep    ! Epoch counter (process noise scaled to
                      ! to be per-epoch).

* LOCAL 
      integer*4 i,j   ! counters
     .,         np    ! Parameter number

****  Step the state vector forward (all process noises are
*     random walks so easy)
      do i = 1, num_parm
         sol_vecm(i) = sol_vecp(i)
         do j = 1, num_parm
            cov_parmm(i,j) = cov_parmp(i,j)
         end do
      end do

****  Now add process noise
      do i = 1, num_site
         do j = 1, 3
            if( site_parn(j,i).gt. 0 .and. mar_site(j,i).gt.0 ) then
               np = site_parn(j,i)
               cov_parmm(np,np) = cov_parmm(np,np)+mar_site(j,i)*kf_step
            end if
         end do
         if( atm_parn(i).gt.0 .and. mar_atm(i).gt.0 ) then
            np = atm_parn(i)
            cov_parmm(np,np) = cov_parmm(np,np)+mar_atm(i)*kf_step
         end if
      end do

****  Thats all
      return
      end

  
