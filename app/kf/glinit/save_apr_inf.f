CTITLE save_apr_INF
 
      subroutine save_apr_inf ( apr_values , crunjd, grunjd)

      implicit none
 
 
*     Routine to get and save the names of sites, sources, and their
*     positions and axis offsets (for sites).
*     We also save any apriori values for the tidal parameters (both
*     standard model and the extended model values).
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   ierr            - File reading error flag
 
      integer*4 ierr
 
*   apr_values(65642j)  - EMA storage for values read from global
*                   - file (Use a latge number to convince compiler
*                   - to use I*4 index)
 
      real*8 apr_values(*)
 
      real*8 crunjd, grunjd  ! Run time of current hfile being
              ! read (crunjd), run time of latest solution in this
              ! epoch of data (grunjd) and sectag for conversion.
              ! Values are used to see if SV antenna offsets should
              ! by updated. 
 
* MOD TAH 150825: Now check the satellite antenna offset values.  Here we 
*     read the records and update if this hfile has a run-time greater than
*     any previously processed ones. 
*     Read over the satellite information records
      call read_svinf_rec 

***** Firstly get the site and source names for this experiment.
*     Read the names one block at a time.  After this call we have
*     all the site and source names added to our list.  There are also
*     arrays ltog_sites and ltog_sources which convert from local to
*     global site and source numbers (for this experiment)
 
      call rw_names_block('R')
 
***** Now read the apriori values, and get as many of the site and
*     source values as we can.  We also use the solution itself
*     to provide apriori values for those parameters which were
*     estimated

      call get_aprioris( ierr, apr_values, crunjd, grunjd )

****  Now see if we need to change site names due to Earthquakes
*     or due to renaming sites.

      call eq_name_change('UPD')
 
***** Thats all
      return
      end
 
