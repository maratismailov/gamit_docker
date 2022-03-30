CTITLE CLOSE_SOLVK
 
      subroutine close_solvk(ierr)

      implicit none
 
c     Close Solvk file
c MOD TAH 870128 Modified to use FmpClose
 
 
      include '../includes/kalman_param.h'
      include '../includes/solvk_cntl.h'
 
*   FmpClose            - HP function to close file, return error code
 
      integer*4 ierr, FmpClose
 
      ierr = FmpClose(isoldc)
 
      call report_error('FmpClose',ierr,'clos','SOLVK common',0,
     .                  'CLOSE_SOLVK')
 
      control_read = .false.
 
      end
 
