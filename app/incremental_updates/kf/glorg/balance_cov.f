CTITLE BALANCE_COV

      subroutine balance_cov ( cov_parm )

      implicit none

*     Routine to used the number of times at site
*     have been used to re-scale covariances matrix.  Only applied
*     for single days.  Mostly used for IGS combination.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED 
      real*8 cov_parm(num_glb_parn,num_glb_parn)   ! Solution covarince 
                ! matrix 

* LOCAL
      real*8 scale(max_glb_parn)  ! scale factor; sqrt(number of times
                                  ! a site is used) and 1 for other
                                  ! parameters
      real*8 pcoscale ! Scale factor for Satellite PCOs.
      integer*4 i, j  ! Loop counters
      integer*4 np    ! Parameter number
     

***** Compute the scaling factors
      scale(1:num_glb_parn) = 1.d0
      do i = 1,gnum_sites
         do j = 1,3
            np = parn_site(j,1,i)
* MOD TAH 210510: Due to book keeping count is 3 times actual
*           usage.
            if( np.gt.0 ) scale(np) = sqrt(times_used(i)/3.d0)
         end do
      end do

***** Satellite PCO's
      pcoscale = sqrt(dble(num_glb_sol))
      do i = 1, gnum_svs
         do j = max_svs_elem-2, max_svs_elem
            np = parn_svs(j,i)
            if( np.gt.0 ) scale(np) = pcoscale
         enddo
      enddo

***** Now scale the matrix; scale rows
      do i = 1, num_glb_parn
         call dwsmy( scale(i), cov_parm(i,1), num_glb_parn, 
     .                         cov_parm(i,1), num_glb_parn,
     .                         num_glb_parn)
      end do
c
***** Scale the columns
      do i = 1, num_glb_parn 
         call dwsmy( scale(i), cov_parm(1,i),1, 
     .                         cov_parm(1,i),1, num_glb_parn)
      end do

***** Thats all
      return
      end





