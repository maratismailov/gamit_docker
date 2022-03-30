CTITLE GLB_AVERAGE
 
      subroutine glb_average(cov_sav_parm, sol_sav_parm, jmat,
     .      cov_parm, sol_parm )

      implicit none 
 
c
c     This routine will average the forward and back solutions
c     and return the averaged values in cov_sav_parm and
c     sol_sav_parm.  The algorithm used is based on the standard
c     kalman filter estimation algorithm.
c
c Include files
c -------------
      include '../includes/kalman_param.h'
c
      include '../includes/globk_common.h'
c
c Variables
c ---------
c cov_sav_parm -- the covariance matrix from the forward running
c     solution.  This matrix is read from disk
c sol_sav_parm -- the solution vector from the forward running
c     solution.  This vector is read from disk.
c jmat -- the kalman gain matrix for the averaging process
c cov_parm  -- the covariance matrix from the backward running solution
c sol_parm  -- the solution vector from the backward running solution.
c sav_fin_parm -- the parameter estimates from the final solution
c
      real*8 cov_sav_parm(num_glb_parn,num_glb_parn),
     .    sol_sav_parm(num_glb_parn), jmat(num_glb_parn,num_glb_parn),
     .    cov_parm(num_glb_parn,num_glb_parn), sol_parm(num_glb_parn)
 
c
c
c Local variables
c ---------------
c ident -- idenification for use of common
c scale -- a scaling vector used by invert_vis.  This variable is also
c     used as a scratch vector in various calculations as well.
c ipivot -- the pivot vector used by invert vis.
c
      real*8 scale(max_glb_parn)
 
c
 
      integer*4 ident, ipivot(max_glb_parn)
 
*   i,j,k   - Loop counters
*   im      - Index to row of matrix
*   end,start   - End and start row numbers when we are computing
*               - covariance matrix.  If postfit residuals are to
*               - to be calculated we compute whole matrix
*   kbit    - Bit checking function.

      integer*4 i,j, im, end,start

      logical kbit

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters
 
c
c
      common ipivot, scale
 
      data I8 / 1 /

c.... Set the idenification for use of common
c
      ident = 9
c
c
c.... First add the covariance matrices (dwadd: v +v =v )
c                                                1  2  3
      call dwadd8(cov_sav_parm,1, cov_parm,1, cov_sav_parm,1,
     .   (I8*num_glb_parn)*num_glb_parn)
c
c.... Now invert this sum in place.
      call invert_vis( cov_sav_parm, sol_sav_parm, scale, ipivot,
     .   num_glb_parn, num_glb_parn, 0)
c
c.... Now form the kalman gain matrix. Note: we form the transpose
c     of this matrix, so that we can always access it by column
c     rather than by row.  (For ema manipulation the column access
c     is much more efficient than row access)
c     Only form columns for the markov elements
*                                 ! row
      do i = 1, num_glb_parn
*                                 ! column
         do j = 1, num_glb_parn
c
            call dwdot(scale(j), cov_parm(1,i),1, cov_sav_parm(1,j),1,
     .         num_glb_parn)
c
         end do
c
c....    Move the column into jmat
         call dwmov(scale,1,jmat(1,i),1, num_glb_parn)
c
      end do
c
c.... Now get the averaged estimates of the markov elements
c     Form difference sol_sav_parm-sol_parm and save in sol_sav_parm
      call dwsub(sol_sav_parm,1, sol_parm,1,
     .           sol_sav_parm,1,num_glb_parn)
c
c.... Multiply difference by jmat
c     Only do the markov elements
      do im = 1, num_glb_parn
c
         call dwdot(scale(im),jmat(1,im),1, sol_sav_parm,1,
     .      num_glb_parn)
c
      end do
c
c.... Now add correction to sol_parm, for the markov or copy earlier
c     results for deterministic
      do im = 1, num_glb_parn
            sol_sav_parm(im) = sol_parm(im) + scale(im)
      end do
c
c
c.... Now complete the covariance matrix. (Only do Marvkov)
      do i = 1,num_glb_parn
 
*        If we are going to do postfit residuals, compute whole
*        matrix, otherwize just do the diagonal+3 (allow for NEU calculation)
* MOD TAH 910403: Compute full matrix if baseline length components
*        are to be shown as well
* MOD TAH 950926: The difference below makes very little difference so
*        compute the out to 6 to cover orbital elements 
* MOD TAH 161104: Expand futher to allow for 15 orbit parameters or
*        upto 3*num_num_pmu (used to be 6)
* MOD TAH 190528: Increased off diagonal for ECOMC orbit model from 
*        15 to 19 (additional cos/sin 2U and 4U direct)
         if( compute_glb_res .or. kbit(bak_opts,2) ) then
            start = i
            end   = num_glb_parn
         else
            start = i
            end   = i+max(19,3*num_mul_pmu)
            if( end.gt.num_glb_parn ) end = num_glb_parn
         end if
            start = i
            end   = num_glb_parn
 
 
****     Now do calculations for upper diagonal
         do j = start,end
 
            call dwdot(scale(j),jmat(1,j),1, cov_parm(1,i),1,
     .         num_glb_parn)
 
            cov_sav_parm(j,i) = cov_parm(j,i) - scale(j)
            cov_sav_parm(i,j) = cov_sav_parm(j,i)
 
         end do
 
      end do
c
c.... Thats all
      return
      end
 
