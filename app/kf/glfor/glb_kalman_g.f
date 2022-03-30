CTITLE GLB_KALMAN_GAIN
 
      subroutine glb_kalman_gain( cov_parm, a_part, part_pnt, cov_obs,
     .                            temp_gain, kgain )

      implicit none  
 
*     Routine to compute the Kalman filter gain matrix.
*     The gain is computed using
*
*              m    T       m   T -1
*     KGAIN = C    A (Q + AC   A )
*              m+1          m+1
*
*            m
*     where C    is the predicted covariance matrix at this step
*            m+1                                       (cov_parm)
*
*           A    is the partials matrix (a_part)
*           Q    is the covariance matrix of the data (cov_obs)
*
*      m   T
*     C   A   is saved in temp_gain for speed of execution.
*      m+1
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'
      include '../includes/glb_hdr_def.h'
 
*   ir,ic       - Row and column counters
*   ipivot(max_glb_parn)    - Pivot values for invert VIS
*   part_pnt(2,max_glb_deriv,cnum_parn) - Pointer to which
*                   - parameters as used in the partials matrix
 
 
      integer*4 ir,ic, ipivot(max_glb_parn),
     .    part_pnt(2,max_glb_deriv,cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn)     - Compressed partials
*                   - matrix
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix of the
*                   - parameters being estimated
*   cov_obs(cnum_parn,cnum_parn)        - Covariance matrix of the
*                   - parameters from the SOLVK solution
*   kgain(cnum_parn,num_glb_parn)       - Kalman filter gain
*                   - matrix
*   scr_real(max_glb_parn)              - Temporary storage for
*                   - rows of matrix during calculations
*   temp_gain(cnum_parn,num_glb_parn)   - Temporay matrix as
*                   - decribed above
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn,cnum_parn), kgain(cnum_parn,num_glb_parn),
     .    scr_real(max_glb_parn), temp_gain(cnum_parn,num_glb_parn)
 
      real*8 min_cor, cor
      integer*4 irmc, icmc
      logical kbit 

*     Find the max negative correlation
      min_cor = 0      
      do ir = 1, cnum_used
         do ic = 1,ir-1
            cor = cov_obs(ir,ic)/sqrt(cov_obs(ir,ir)*cov_obs(ic,ic))
            if( cor.lt.min_cor ) then
                min_cor = cor
                irmc = ir
                icmc = ic
            end if
         end do
      end do

      if( kbit(crt_opts,16) )
     .write(log_unit,100) ' Cov_obs before update', irmc,icmc, min_cor
 100  format(a,' Row, Col ',i5,i5,' Min cor ',f10.5)
 
*           m   T
***** Form C   A  by rows
*           m+1
*                                 ! Loop over rows
      do ir = 1, num_glb_parn
*                                 ! Loop over columns
          do ic = 1, cnum_used
   
              call glb_dot_part( scr_real(ic), cov_parm(1,ir),1,
     .             a_part(1,ic), part_pnt(1,1,ic), indx_pnt(ic) )
 
          end do
 
*         Now move the column into temp_gain
          call dwmov( scr_real,1, temp_gain(1,ir),1, cnum_used)
      end do
      
*                          m   T
*     Now finish forming AC   A   and add to Q
*                          m+1
 
      do ic = 1, cnum_used
*                                 ! Only do lower diagonal and copy
          do ir = ic, cnum_used
 
              call glb_dot_part( scr_real(ir),
     .             temp_gain(ic,1),cnum_parn, a_part(1,ir),
     .             part_pnt(1,1,ir), indx_pnt(ir) )
          end do
 
*         Add this contribution to cov_obs and then move values to
*         upper diagonal
 
          do ir = ic, cnum_used
              cov_obs(ir,ic) = cov_obs(ir,ic) + scr_real(ir)
          end do
 
*         Move
          call DWMOV(cov_obs(ic+1,ic),1, cov_obs(ic,ic+1),cnum_parn,
     .               cnum_used-ic)
      end do
 
*     Find the max negative correlation
      min_cor = 0      
      do ir = 1, cnum_used
         do ic = 1,ir-1
            cor = cov_obs(ir,ic)/sqrt(cov_obs(ir,ir)*cov_obs(ic,ic))
            if( cor.lt.min_cor ) then
                min_cor = cor
                irmc = ir
                icmc = ic
            end if
         end do
      end do

      if( kbit(crt_opts,16) )
     .write(log_unit,100) ' Cov_obs after  update', irmc,icmc, min_cor

*     Now invert the covariance matrix
 
      call invert_vis( cov_obs, kgain, scr_real, ipivot, cnum_used,
     .                 cnum_parn, 0)
 
*     Now complete the calcualtion of the Kalman Gain
 
      do ir = 1, num_glb_parn
          do ic = 1, cnum_used
 
              call DWDOT(scr_real(ic), temp_gain(1,ir),1,
     .                   cov_obs(1,ic),1, cnum_used)
 
          end do
 
*         Move scr_real into the Kalman gain matrix
          call dwmov( scr_real(1), 1, kgain(1,ir),1, cnum_used)
 
*         DEBUG
C         write(*,100) ir,(kgain(ic,ir),ic=1,cnum_used)
C 100     format(' Param ',i4,10( 5(f12.6,1x),:/,11x))
 
      end do
 
***** Thats all
      return
      end
 
