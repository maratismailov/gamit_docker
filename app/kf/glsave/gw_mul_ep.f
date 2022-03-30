CTITLE GW_MUL_EP
 
      subroutine gw_mul_ep

      implicit none  
 
*     Routine to write the epochs of the multi-epoch parameter
*     and apriori values.
*
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
 
 
      integer*4 ierr
 
***** Now write out record(s)
      if( cnum_par_ep.gt.0 ) then
         call writd(cglb_dcb,ierr,mul_par_ep,128*cnum_par_ep,
     .               crec_par_ep)
          call report_error('FmpWrite',ierr,'writ',
     .                      'Parameter multi-epochs',1,
     .                      'qw_CODES')
      end if
      if( cnum_apr_ep.gt.0 ) then
          call writd(cglb_dcb,ierr,mul_par_ep,128*cnum_apr_ep,
     .               crec_apr_ep)
          call report_error('FmpWrite',ierr,'writ',
     .                      'Apriori multi-epochs',1,
     .                      'qw_CODES')
      end if
 
 
***** Thats all
      return
      end

