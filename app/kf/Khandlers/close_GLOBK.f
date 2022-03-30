CTITLE CLOSE_GLOBK
 
      subroutine close_GLOBK(ierr)

      implicit none
 
c     Close Globk file
c MOD TAH 870128 Modified to use FmpClose
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
*   FmpClose            - HP function to close file, return error code
 
      integer*4 ierr, FmpClose
 
      ierr = FmpClose(glb_com_dcb)
 
      call report_error('FmpClose',ierr,'clos','GLOBK common',0,
     .                  'CLOSE_GLOBK')
 
      glb_con_read = .false.
 
      end
 
