CTITLE GW_CODES
 
      subroutine gw_codes

      implicit none  
 
*     Routine to write the list of parameter and aproiri codes
*     to the global file.
*
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
 
 
      integer*4 ierr
 
***** Now write out record(s)
 
      call writd(cglb_dcb,ierr,glb_codes,128*cnum_par_types,
     .           crec_par_types)
      call report_error('FmpWrite',ierr,'writ','Global codes',1,
     .                  'GW_CODES')
 
      call writd(cglb_dcb,ierr,apr_codes,128*cnum_apr_types,0)
      call report_error('FmpWrite',ierr,'writ','Aproiri codes',1,
     .                  'GW_CODES')
 
***** Thats all
      return
      end
 
