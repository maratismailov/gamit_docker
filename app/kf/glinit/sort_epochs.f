CTITLE SORT_EPOCHS
 
      subroutine sort_epochs( expt_names, sepoch_expts, expts_var)
 
 
*     This routine will take the epochs of the global solutions
*     and sort them into time order (either ascending or descending)
*     depending on sort direction. The file names are also moved.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
 
*   i,j,k                           - Loop counters
*   expt_names(max_glb_sol )        - Names of the global
*                                   - solutions
*   ctemp                           - Storage for name switching
 
      integer*4 i,j
      character*(sort_recl) expt_names(max_glb_sol), ctemp
 
*   sepoch_expts( max_glb_sol )     - Epochs of each solution.
*                                   - List will be sorted based
*                                   - on these values
*   temp_epoch                      - Temporary storage for
*                                   - switching epochs.
*   expts_var(max_glb_sol)          - Variances to be given to
*                                     each experiment (to make chi**2
*                                     unity).
*   temp_var                        - Temporary storage for shifting
*                                     variance
*   apr_values( max_glbapr_sz )   - Storage for covariance
*                                   - matrix and solution vector
 
      real*8 sepoch_expts( max_glb_sol ), expts_var(max_glb_sol)

      real*8 temp_epoch, temp_var
 
***** Loop doing a bubble sort
 
      do i = 1, num_glb_sol - 1
 
          do j = 1, num_glb_sol - i
 
              if( sort_direction*sepoch_expts(j).gt.
*                                                         ! switch values
     .            sort_direction*sepoch_expts(j+1) ) then
 
                  call switch_8(sepoch_expts(j),sepoch_expts(j+1),
     .                          temp_epoch)

                  call switch_8(expts_var(j),expts_var(j+1),temp_var)
                  call switch_ch(expt_names(j),expt_names(j+1),ctemp )
              end if
          end do
      end do
 
*     Save epochs
      if( sort_direction.eq.1 ) then
* MOD TAH 020323: start and end now done in glinit before call to save_epoch
C         gepoch_start = sepoch_expts(1)
          gepoch_expt  = gepoch_start
          gepoch_prev  = gepoch_start
C         gepoch_end   = sepoch_expts(num_glb_sol)
      else
C         gepoch_start = sepoch_expts(num_glb_sol)
C         gepoch_end   = sepoch_expts(1)
          gepoch_expt  = gepoch_end
          gepoch_prev  = gepoch_expt
      end if
 
***** Thats all
      return
      end
 
