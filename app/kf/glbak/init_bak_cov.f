CTITLE INIT_BAK_COV
 
      subroutine init_bak_cov ( cov_parm, sol_parm, cov_sav_parm,
     .                          sol_sav_parm )

      implicit none  
 
*     Routine to initialize the variances for the bak kalman filter.
*     The values used are 1000 times larger than the final sigma
*     for the forward solution.  This should leave the solution
*     "fairly" unconstrainted .i.e. about 0.1% effect on varaiances.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
*   i       - Loop counter
 
      integer*4 i
 
*   cov_parm(num_glb_parn,num_glb_parn) - Matrix to be initialized
*   cov_sav_parm(num_glb_parn,num_glb_parn) - Final solution from
*                                       - forward solution
 
*   multiplier              - Multiplier for get apriori variances
*                           - from end variances
 
*   sol_parm(num_glb_parn)  - Solution vector to be cleared
*   sol_sav_parm(num_glb_parn)  - Final solution vector
 
 
      real*8 cov_parm(num_glb_parn,num_glb_parn),
     .    cov_sav_parm(num_glb_parn,num_glb_parn), multiplier,
     .    sol_parm(num_glb_parn), sol_sav_parm(num_glb_parn)
 
 
 
***** Loop over columns in cov_parm and set sol_parm
 
      multiplier = 3000.d0
 
C     call dwsmy(0.d0, sol_parm,1, sol_parm,1, num_glb_parn)
      call dwint(0.d0, sol_parm,1, num_glb_parn)
 
      do i = 1, num_glb_parn
 
C         call dwsmy(0.d0, cov_parm(1,i),1, cov_parm(1,i),1,
C    .               num_glb_parn)
          call dwint(0.d0, cov_parm(1,i),1, num_glb_parn)
* MOD TAH 190528: Changed dmax1 to max (generic routine) and
*         put upper limit of 1.d4 
          cov_parm(i,i) = min(max(multiplier*cov_sav_parm(i,i),
     .                        dble(cov_apr(i))),1.d4)
      end do
 
***** Thats all
      return
      end
 
