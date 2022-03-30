CTITLE GW_GLB_NAMES
 
      subroutine gw_glb_names

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
*   names_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
 
      integer*4 i, ierr, names_list(1), j
      logical kbit
 
*   cnames_list(max_names)  - List of site and sources names
 
      character*8 cnames_list(max_names)
 
      equivalence ( cnames_list, names_list )
 
***** Make the list of site and source names
      j = 0 
      do i = 1, gnum_sites
          if( kbit(guse_site,i) ) then
              j = j + 1
              cnames_list(j) = gsite_names(i)
          end if
      end do
 
      do i = 1, cnum_sources
          cnames_list(i+cnum_sites) = gsource_names(i)
      end do

      do i = 1, cnum_svs
          cnames_list(i+cnum_sites+cnum_sources) = gsvs_names(i)
      end do
 
***** Now write out record(s)
 
      call writd(cglb_dcb,ierr, names_list,128*cnum_names,crec_names)
      call report_error('FmpWrite',ierr,'writ','Names block',1,
     .                  'GW_glb_Names')
 
***** Thats all
      return
      end
 
