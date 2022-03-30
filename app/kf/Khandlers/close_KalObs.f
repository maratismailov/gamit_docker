CTITLE CLOSE_KALOBS
 
      subroutine close_KalObs(ierr)

      implicit none
 
c     Close KalObs
c MOD TAH 870128 Modified to use FmpClose
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*   FmpClose            - HP function to close file, return error code
 
      integer*4 ierr, FmpClose
 
      ierr = FmpClose(ko_dcb)
 
      call report_error('FmpClose',ierr,'clos','KalObs file',0,
     .                  'CLOSE_KALOBS')
 
      values_read = .false.
 
      end
 
