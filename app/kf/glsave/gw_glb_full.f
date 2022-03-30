
CTITLE gw_glb_full
 
      subroutine gw_glb_full

      implicit none  
 
*     Routine to create the list of site and source names and save
*     these in the global file.
*
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   max_names   - Maximum number of names we might get
 
      integer*4 max_names
 
      parameter ( max_names = max_glb_sites + max_glb_sources )
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   full_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
 
      integer*4 i, ierr, full_list(max_glb_sites*8)
 
*   cnames_list(max_names)  - List of site and sources names
 
      character*32 cnames_list(max_glb_sites)
 
      equivalence ( cnames_list, full_list )
 
***** Make the list of site and source names
 
      do i = 1, cnum_sites
          cnames_list(i) = gsite_full(i)
      end do
 
***** Now write out record(s)
 
      call writd(cglb_dcb,ierr, full_list,128*cnum_full ,crec_full)
      call report_error('FmpWrite',ierr,'writ','Full  block',1,
     .                  'gw_glb_full')
 
***** Thats all
      return
      end
 
