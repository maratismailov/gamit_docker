CTITLE GW_APRIORIS
 
      subroutine gw_aprioris 

      implicit none 
 
*     Routine to make the list of apriori values and write these
*     to the global file.
*
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   type        - Type of apriori value
*   indx        - Index to the array pointed to the code
 
      integer*4 i, ierr, type, indx
 
***** Make up the list.  Here we decode the codes we have to find the
*     apriori to be saved
 
      do i = 1, num_apr_codes

*         Use codes with original site numbers in them 
          call decode_code( org_apr_codes(i), type, indx )
          call glb_save_apr( type, indx, apr_list(i) )
 
      end do
 
***** Now write out record(s)
 
      call writd(cglb_dcb,ierr,apr_list,128*cnum_apr_vals,crec_apr_vals)
      call report_error('WWRIT',ierr,'writ','aproiri vals',1,
     .                  'GW_aproiris')
 
***** Thats all
      return
      end
 
