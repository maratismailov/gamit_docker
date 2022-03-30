      subroutine glout(ms_type, outfile, opts)

      implicit none  
 
*     This is the output program for GLOBK.  It reads the GLOBK
*     common, and the last covariance matrix from the GLFOR
*     solution, and prints out the results
*
*     The specifics of the output are determined by the options
*     variable passed in the runstring.
*
*     The runstring for GLOUT is:
*     CI> GLOUT lu option glb_com_file.
*     lu is output lu,
*     option is the bit mapped options (See globk_control.ftni)
*     glb_com_file is the name of the global solution control file.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/sd_common.h'

* PASSED variables

* ms_type - sets if called from main progrom (MAIN) or
*           subroutine
* outfile - Output file name
* opts    - Options for output

      character*(4) ms_type
      character*(*) outfile
      integer*4 opts
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   ierr        - SEGLD error flag
*   iout        - Output LU
*   options     - options passed through runstring
*   i           - Loop counter
 
      integer*4 ema_data(1), iout, options, ierr, trimlen, i

*   appy_nut    - Logical indicating that we should apply nutation series 
*                 corrections
*   apply_sd    - Indicates we are apply SD corrections

      logical apply_nut, apply_sd, kbit
 
*   scratch(max_glb_parn)   - Scratch common area to be used
*               - for miscellaneous calculations
 
      real*8 scratch(max_glb_parn)

      common / progcon / iout, options
 
      common / globc_ema / ema_data
 
      common scratch
 
***** The first segment decodes the runstring, and opens and reads the
*     files we will need.
      options = opts
      call GLOUO(ms_type, outfile)

* MOD TAH 070823: See if mid-point option is selected
      if ( kbit(options,29) ) then
           gepoch_out = (gepoch_start+gepoch_end)/2.d0
      else
           gepoch_out = gepoch_expt
      endif

*     Get any updates to the satellites orbits and possibly polar motion
* MOD TAH 070926: Added explicit satellite epoch

      call glb_upd_apr( gepoch_out, gepoch_expt, .false., 
     .                  ema_data(1), .false. )

*     Read in the new nutation series.  To apply planetary nutation need
*     also to nutation series as well.
      if( trimlen(nut_inp_file).gt.0 ) then
          apply_nut = .true.
          call read_nut_series(110, nut_inp_file, plan_inp_file,ierr )
          if( ierr.ne.0 ) apply_nut = .false.
      else
          apply_nut = .false.
      end if

*     Read the short period UT1 corrections
      call init_sd( sd_inp_file, apply_sd, .false., gsite_names,
     .              gnum_sites, gepoch_out )

*     Now compute nutation series corrections.
      if( apply_nut ) then
          call compute_std( gepoch_out, gnut_ang_apr )
      end if

      call get_sd_mid( gepoch_out, apply_sd, gnum_sites)

*     Compute the contributions to the postition at the epoch of
*     output for the nonsecular terms
      do i = 1, gnum_sites
         call eval_nonsec(i, gepoch_out, num_nonsec, param_nonsec,
     .                    apr_val_nonsec, cont_nonsec(1,i),0)
      end do

* MOD TAH 050927: Get the long site names from the head.snx file (either
*     local or in $HELP_DIR or gg/tables)

      call get_long_names

 
*     The next segment stablizes the solution (i.e., sets translation
*     origin, and orientation origin
 
      call GLOUS
 
*     Now output the solution to the output device
 
      call GLOUW
 
***** Thats all
      end
 
 
