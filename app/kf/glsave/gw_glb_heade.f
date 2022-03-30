CTITLE GW_GLB_HEADER

      subroutine gw_glb_header

      implicit none  
 
*     Routine to create and write the global header record
*     for the global solution file.  This routine will also purge
*     and create the file to which the solution will be written.
*
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
      include '../includes/globk_common.h'
      include '../includes/obs_header.h'
 
*   AddressOf   - Gets the address of a variable (used to find
*               - out how large the global header records are)
*   i           - Loop counter
*   ierr        - File error flag
*   sln_len     - Length of the solution records

* MOD TAH 190511: Changed to I*8 for 64-bit memory 
      integer*8 AddressOf
      integer*4 i, ierr, sln_len, trimlen
* MOD TAH 190511: Introduced when number of parameters > 32767.
      integer*8 i8  !  i8 value of 1 to force I8 calculations
 
*   jsize       - Size of the file in blocks
 
      integer*4 jsize
      logical kbit
 
*   Full_file_name  - Full name of the global file (with size
*                   - Information)
 
      character*64 Full_file_name
 
***** First compute the record boundaries of all the blocks for the
*     global file
*     Get number of used sites
      i8 = 1
      cnum_sites = 0
      do i = 1, gnum_sites
         if( kbit(guse_site,i) ) cnum_sites = cnum_sites + 1
      end do
 
      cnum_header = ( AddressOf(last_glb_header_wrd) -
     .                AddressOf(glb_header) )/(4*128)  + 1
      sln_len     = ( AddressOf(sln_dummy(1)) -
     .                AddressOf(sdelete_count(1))-4)/(4*128) + 1
      cnum_names  = ( cnum_sites+gnum_sources+gnum_svs )*2/128 + 1
      cnum_full   = ( cnum_sites*8 - 1)/128 + 1
*     Add additional solution record for each of the global
*     combined solutions used.      
*     This total number of solution records in all the input
*     global files.  (These records are added in glinit).
* MOD TAH 190621: Only update if not BACK type.
      if( .not.back_type ) then 
          cnum_soln_recs  = gnum_soln_recs

*          If we have combined more than input file, then add an
*          extra solution record for the combined solution.  This also
*          cause and extra site info record to be written.
           if( num_glb_sol.gt.1 ) cnum_soln_recs = cnum_soln_recs + 1
      endif

*     Compute number of records needed for the site information 
*     for each solution.  The site information is allocated 
*     256 bytes per station.  The solution records themselves
*     contain the record number of the station information.  
*     Initialize number of site info records with the number
*     needed for this combined solution.
      cnum_sinf = ((cnum_sites-1)/2+1)
* MOD TAH 981222: Records per global are not saved so code
*     below commented out:
c     do i = 1, gnum_soln_recs 
c        if( ns_by_sol(i).gt.0 ) then
c            cnum_sinf = cnum_sinf + ((ns_by_sol(i)-1)/2+1)
c        end if
c     end do 
* MOD TAH 050622: Save the number of satellite information records
      cnum_svinf = ((gnum_svs-1)/4+1) 

****  Save the number of combined solution headers there are
      cnum_comb = gnum_comb 
      if( num_glb_sol.gt.1 ) cnum_comb = cnum_comb + 1

*     Now add in the number of stations needed by combined 
*     global solutions.           
              
      cnum_par_types = (num_glb_parn   -1)/128 + 1
      cnum_apr_types = (num_apr_codes  -1)/128 + 1
      cnum_apr_vals  = (2*num_apr_codes-1)/128 + 1
C     cnum_par_vals  = (2*num_glb_parn*(num_glb_parn+1)-1)/128 + 1 
      cnum_par_vals  = (2*num_glb_parn*(num_glb_parn+i8)-1)/128 + 1 

* MOD TAH 981104: Add in the multi-day parameter epochs
      cent_par_ep    = ent_par_ep
      cent_apr_ep    = ent_par_ep
      cnum_par_ep    = (cent_par_ep - 1)/64 + 1
      cnum_apr_ep    = (cent_apr_ep - 1)/64 + 1

****  Save the number of entries on the apriori constraints.  This
*     is num_glb_parn + off diagonals for NEU constraints (3 per one)
      cnum_acvc      = 4*(num_glb_parn+3*gnum_off_diag)
      cnum_acvc_recs = (cnum_acvc - 1)/128 + 1
 
      crec_names     = 1              + cnum_header
      crec_full      = crec_names     + cnum_names
      crec_solutions = crec_full      + cnum_full 
      crec_sinf      = crec_solutions + cnum_soln_recs*sln_len
*     Save the record for the solution record that referrs to the
*     combined global file.
      crec_comb_soln = crec_sinf      - sln_len
      crec_svinf     = crec_sinf      + cnum_sinf
* MOD TAH 050622: Add records for satellite information
      crec_par_types = crec_svinf     + cnum_svinf
      crec_apr_types = crec_par_types + cnum_par_types

* MOD TAH 981020: Add the records for the multi-epoch polar motion
      crec_par_ep    = crec_apr_types + cnum_apr_types
      crec_apr_ep    = crec_par_ep    + cnum_par_ep
      
      crec_apr_vals  = crec_apr_ep    + cnum_apr_ep
      crec_par_vals  = crec_apr_vals  + cnum_apr_vals

*     Compute record number for start of apriori constraint block.
*     The apriori constraint values and the number of records will
*     computed later as the records are written.
      crec_acvc      = crec_par_vals +  cnum_par_vals

*     jsize not really needed san so we don't need to be exact here.
      jsize  = crec_acvc      + cnum_acvc_recs - 1
 
***** Purge any file with the same name as we will create
 
      call clean_file( glb_out_file )
 
      call FullFileName( glb_out_file, 1, jsize, 128, full_file_name)
 
      call FmpOpen( cglb_dcb, ierr, Full_File_Name, 'wc', 1)
      call report_error('FmpOpen',ierr,'creat',Full_file_name,1,
     .                  'GW_glb_header')
 
***** Now set up the rest of the header blocks
* MOD TAH 190916: Only update these values if we are not using
*     BACK type solution creating a GLX file in glbak.
      cfile_type   = 'GLOBAL'
      if( trimlen( gexpt_title).le.0 ) then
          cexpt_title  = 'COMBINED GLOBAL from ' // list_file
      else
          cexpt_title  = gexpt_title
      end if
      cnum_apr_codes = num_apr_codes
      cnum_parn    = num_glb_parn
*     Cnum_sites is computed above accounting for the sites that
*     are actually used.  So dont resave here.  (The suppress sites
*     will have 0 codes and should be ignored later).
C     cnum_sites   = gnum_sites
      cnum_sources = gnum_sources
      cnum_svs     = gnum_svs
      cgamit_mod   = ggamit_mod 
 
      do i = 1, 7
          crun_time(i) = grun_time(i)
      end do
 
      ctai_utc     = gtai_utc
      cgpst_utc    = ggpst_utc
 
      do i = 1, max_edit_types
          cdelete_count(i) = gdelete_count(i)
      end do

* MOD TAH 190628: Only update value if not back_type (BACK glsave)
      if( .not. back_type ) cnum_obs = gnum_obs

* MOD TAH 960812: Save the chi**2 valued for the solution
*     in the combined binary (Useful later for checking 
*     quality and for Sinex Ver 1.0)      
      if( sum_chi_num.gt.0 ) then
          cchisq   = sum_chi/sum_chi_num 
      else
          cchisq   = 0
      end if
      
      if( .not. back_type ) then  
         cepoch_end   = gepoch_end
         cepoch_expt  = gepoch_out
         cepoch_start = gepoch_start
      endif

      csvs_epoch   = svs_epoch(1)
 
      do i = 1,3
          cetd_apr(i) = apr_val_tid(i,1)
      end do
 
      do i = 1,8
          cnut_ang_apr(i) = gnut_ang_apr(i)
      end do
 
      do i = 1,6
          cut1_apr(i) = gut1_apr(i)
      end do
 
      do i = 1,8
          cwob_apr(i) = gwob_apr(i)
      end do

****  Add the new values

      cglb_vers = glbf_version
      ccons_type = gcons_type
      csys_type  = gsys_type
      call get_institute( cowner, default_institute)
      call get_institute( ccreator, default_institute)
      cprog_gen = '=GLK'
      call get_anal_type( canal_type )

*     Copy over new variables MOD TAH 050622:
      cload_mod = gload_mod
* MOD TAH 130419: See if need to update the loads applied
*     based app_modl command options.
      if( kbit(appload_mod,1) ) then   ! Command has been used, see what the
*                                        status is
*         Test atmospheric load status (see if change in model)
          if( kbit(appload_mod,2) ) then
              call sbit(cload_mod,9,1) 
          else   ! correction removed so unset
              call sbit(cload_mod,9,0)
          end if
*         Test Hydrologic load
          if( kbit(appload_mod,3) ) then
              call sbit(cload_mod,25,1) 
          else   ! correction removed so unset
              call sbit(cload_mod,25,0)
          end if
      end if
*
      
      cspeopmod = gspeopmod 
      cetidemod = getidemod  
      cotidemod = gotidemod  
      coatmlmod = goatmlmod  
      catmtdmod = gatmtdmod  
      chydromod = ghydromod 
      cgnut     = ggnut
      cggrav    = gggrav
* MOD TAH 140327: Added Earth radiation and antenna thrust.
      ceradmod  = ggeradmod
      cantradmod = ggantradmod
* MOD TAH 140403: Added atm modeling options
      cdryzen   = ggdryzen 
      cwetzen   = ggwetzen 
      cdrymap   = ggdrymap 
      cwetmap   = ggwetmap 
      cionsrc   = ggionsrc 
      cmagfield = ggmagfield 

***** Now write out the header records to the file
 
      call writd( cglb_dcb, ierr, glb_header, 128*cnum_header, 1)
      call report_error('FmpWrite',ierr,'writ','global header',1,
     .                  'GW_glb_header')
 
***** Thats all
      return
      end

CTITLE GET_ANAL_TYPE

      subroutine get_anal_type ( anal_type )

      implicit none 

*     Generates the string with the analysis type

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* anal_type  - Type of analysis

      character*(*) anal_type

      integer*4 trimlen, i
      logical   estimated

*     Set station positions assuming that we have these
      anal_type = 'S'

***** See if estimated velocities
      estimated = .false.
      do i = 1, gnum_sites
         if( parn_site(1,2,i).gt.0 ) estimated = .true. 
      end do
      if( estimated ) then
         anal_type(trimlen(anal_type)+2:) = 'V'
      end if

*     Check EOP parameters
      estimated = .false.
      do i = 1, 4
         if( parn_wob(i).gt.0 ) estimated = .true. 
      end do
      do i = 1,2
         if( parn_ut1(i).gt.0 ) estimated = .true. 
      end do

*     Check multiday PMU
      if( parn_mul_pmu(1,1,1).gt.0 ) estimated = .true.

      if( estimated ) then
         anal_type(trimlen(anal_type)+2:) = 'E'
      end if

*     Orbit parameters
      estimated = .false.
      do i = 1, gnum_svs
         if( parn_svs(max_svs_elem,i).gt.0 ) estimated = .true.
      end do
      if( estimated ) then
         anal_type(trimlen(anal_type)+2:) = 'A'
      end if
 
****  Thats all
      return
      end
