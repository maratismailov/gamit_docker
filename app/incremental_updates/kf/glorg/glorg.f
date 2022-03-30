      subroutine glorg( ms_type, cmdfile, outfile, opts)
 
      implicit none 
 
*     This is the output program for GLOBK.  It reads the GLOBK
*     common, and the last covariance matrix from the GLFOR
*     solution, and prints out the results
*
*     The specifics of the output are determined by the options
*     variable passed in the runstring.
*
*     The runstring for GLORG is:
*     CI> GLORG <lu> <option> <command file> <glb_com_file>.
*     lu is output lu,
*     option is the bit mapped options (See globk_control.ftni)
*     glb_com_file is the name of the global solution control file.
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/sd_common.h'
      include '../includes/glorg_common.h'
c      external glorg_bd 
     
* PASSED variables

* ms_type - sets if called from main progrom (MAIN) or
*           subroutine
* cmdfile - Command file name
* outfile - Output file name
* opts    - Options for output

      character*(4) ms_type
      character*(*) cmdfile, outfile
      integer*4 opts
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   ident(10)   - Ten words of identification for use of
*               - scratch common (not used in this program)
*   ierr        - SEGLD error flag
*   iout        - Output LU
*   options     - options passed through runstring
*   trimlen     - Length of string
*   nu          - Counter for number of sites used in plate
*                 roatation determination
*   i           - Loop counter
*   jcov_obs    - Start address in ema_data for cov_obs
*   jsol_obs    - Start adresss in ema_data for sol_obs
*   jplate_partials - Start address in ema_data for plate partials
*   jatv        - Start address in ema_data for atv
*   jerr, FmpPurge - IOSTAT error for deleting file with FmpPurge routine
 
      integer*4 ema_data(1), ident(10), iout,
     .    options, trimlen, ierr, i, nu,
     .    dcb_sol(16), jerr, FmpPurge

* MOD TAH 070112: Made the memory address variables integer*8 to be
*     consistent with 64-bit machines.
      integer*8 jcov_obs, jsol_obs, jplate_partials, jatv


*   appy_nut    - Logical indicating that we should apply nutation series
*                 corrections
*   apply_sd    - Indicates we are apply SD corrections

      logical apply_nut, apply_sd
    
      logical kbit ! Function to test if bit is set (count 1-32)

 
*   scratch(max_glb_parn)   - Scratch common area to be used
*               - for miscellaneous calculations
 
      real*8 scratch(max_glb_parn)
 
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
 
      common ident, scratch
 
***** The first segment decodes the runstring, and opens and reads the
*     files we will need.
      options = opts

      last_glb_ema = 1
      call GLORO(ms_type, cmdfile, outfile)

*     Get any updates to the satellites orbits and possibly polar motion
* MOD TAH 070926: Added explicit satellite epoch
      call glb_upd_apr( gepoch_out, gepoch_expt, .false., 
     ,                  ema_data(1), .false. )
 
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
                 
***** Set up parameter names

      call set_param_names

***** See if we are applying equates and forces before
*     stabalization
      if( first_eqf .and. num_stab_iter.le.1 ) then
          call report_stat('status','globk','glorg',' ',
     .                 'Applying equates: First',0)
          call GLOREF
      else 
          if( first_eqf ) then
             write(iout,150)
 150         format('**WARNING** FIRST_EQ can not be used when ',
     .              'stabablization iterated.  Equates will be ',
     .              'done later')
             first_eqf = .false.
          end if
      end if

*     The next segment stablizes the solution (i.e., sets translation
*     origin, and orientation origin
      call report_stat('status','globk','glorg',' ',
     .                 'Frame Realization',0)
  
      call GLORS

***** See if we are applying equates and forces after
*     stabalization
      if( .not.first_eqf ) then
          call report_stat('status','globk','glorg',' ',
     .                 'Applying equates',0)
          call GLOREF
      end if

****  Now see if we are to estimate plate rotation vectors
      if( num_plates.gt.0 ) then
          call report_stat('status','globk','glorg',' ',
     .                 'Fitting plates',0)

****      GEt number of sites in plate determiantions
	  nu = 0
	  do i = 1, gnum_sites
	     if( plate_number(i).gt.0 ) then
		 nu = nu + 1
             end if
          end do
	  num_plate_sites = nu
*         Get memory allocation
	  jcov_obs = isol_parm+2*num_glb_parn
	  jsol_obs = jcov_obs +2*9*num_plate_sites*num_plate_sites
	  jplate_partials = jsol_obs + 2*3*num_plate_sites
* MOD TAH 030513: Added 1 to num_plates to allow for translation
	  jatv     = jplate_partials + 
     .                   2*9*num_plate_sites*(num_plates+1)

	  call est_plate(iout, ema_data(icov_parm), 
     .            ema_data(isol_parm), 
     .            ema_data(jcov_obs), ema_data(jsol_obs),
     .            ema_data(jplate_partials), ema_data(jatv) )
      end if

***** Now check that we covariance matrix is OK
      write(iout, 200)
 200  format(' Checking covariance matrix after equate and force')
      call check_covar(ema_data(icov_parm), num_glb_parn, num_glb_parn)

* MOD TAH 210509: Added BALN org_opt check to re-scale covariance
*     matrix based on the number of times a site has been used.
      if( kbit(org_opts,32) .and. gepoch_end-gepoch_start.lt.2.0 ) then
         call report_stat('status','globk','glorg',' ',
     .                 'Balancing numbers of stations',0)
         call balance_cov( ema_data(icov_parm))
      endif
 
*     Now output the solution to the output device
*     Update the name of apriori file
      call report_stat('status','globk','glorg',' ',
     .                 'Outputting solution',0)

      call GLORW
      
****  See if we are going to write a new solution file
      if( trimlen(glr_sol_file).gt.0 ) then
* MOD TAH 060814: If del_scratch option set, then remove the old
*         sol and comfiles before creating new ones
          if( del_scratch ) then
             if( trimlen(glb_com_file).gt.0 ) then 
                jerr = FmpPurge(glb_com_file)
                call report_error('IOSTAT',jerr,'delet',
     .                          glb_com_file,0,'globk') 
             endif
             if( trimlen(glb_sol_file).gt.0 ) then 
                jerr = FmpPurge(glb_sol_file)
                call report_error('IOSTAT',jerr,'delet',
     .                          glb_sol_file,0,'globk')
             endif
          endif 
          call report_stat('status','globk','glorg',glr_sol_file,
     .                 'Updating sol_file',0)

          glb_sol_file = glr_sol_file
          call rw_glb_covar('C', dcb_sol, ema_data(icov_parm))
          call rw_glb_covar('W', dcb_sol, ema_data(icov_parm))
* MOD TAH 051211: Create a new common and write this out.  It is needed 
*         if the apriori files have changed
          glb_com_file = glr_sol_file(1:trimlen(glr_sol_file)) // '.COM'
          call rw_globk_common('N')
      end if
      call report_stat('status','globk','glorg',' ',
     .                 'Normal finish',0)
 
***** Thats all
      end
 
CTITLE GLORO
 
C                        ,Open files segments
      subroutine gloro(ms_type, cmdfile, outfile )

      implicit none  
 
*     This segment decodes the runstring and opens the files
*     we will need
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glorg_common.h'
      
* PASSED variables
* ms_type  - Sets the type of program call (MAIN for program)
* cmdfile  - Command file

      character*(*) ms_type, cmdfile, outfile
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   ierr        - SEGLD error flag
*   iout        - Output LU
*   options     - options passed through runstring
 
      integer*4 ema_data(1), iout, options, trimlen, ierr
      integer*4 num
      
*   kbit         - Check bit status
      logical kbit       

      character*12 full_ver, hsver
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
****  Decode the runstring
      com_orgopt = ' '
      if( ms_type(1:4).eq.'MAIN' ) then
 
          call decode_GLORG_run( iout, options, glb_org_file,
     .         glorg_command_file, glb_com_file, com_orgopt )
          newoutfile = glb_org_file
      else
         iout = 200
         glorg_command_file = cmdfile

*        See if option give to erase file first.
         if( kbit(options,17) ) then
             call open_lu(iout,outfile,ierr,'unknown') 
         else
              call open_lu(iout,outfile,ierr,'append') 
         end if
         newoutfile = outfile
         
         call report_error('IOSTAT',ierr,'open_lu',outfile,
     .          0,'GLORG/GLORO')
         if( ierr.ne.0 ) iout = 6
      endif 
        
*     Now open the GLOBK common and read
      if( ms_type(1:4).eq.'MAIN' ) then 
          call rw_globk_common('R')     
          istart_vma = 0
      end if

* MOD TAH 070823: See if mid-point option is selected
      if ( kbit(options,29) ) then
           gepoch_out = (gepoch_start+gepoch_end)/2.d0
      else
           gepoch_out = gepoch_expt
      endif


****  See if new command file option has been passed
      write(iout,'(a)') gdescription(1:trimlen(gdescription))

      full_ver = hsver(glorg_version)
      write(iout,100) full_ver(1:trimlen(full_ver))  
 100  format(/,' +++++++++++++++++++++++++++++++++++++',
     .       /,' + GLORG                 Version ',a,' +',
     .       /,' +++++++++++++++++++++++++++++++++++++',/)
 
****  See if new command file option has been passed
      if( trimlen(com_orgopt).eq.0 ) then
          com_orgopt = comopt
      endif
      if( trimlen(com_orgopt).gt.0 ) then
          write(*,110) com_orgopt(1: trimlen(com_orgopt))
          if( iout.ne.6 )
     .    write(iout,110) com_orgopt(1: trimlen(com_orgopt))
 110      format('COMOPT: Line starting with ',a,
     .           ' will be interpretted',/) 
      endif
*     Now open and read the covariance matrix and solution vector
      if( ms_type(1:4).eq.'MAIN' ) then
          if ( most_cparn_num.lt.35 ) then
              print *,'GLORG Memory WARNING: most_cparn_num ',
     .                ' to small ',most_cparn_num,'. Reset to 35'
     .              
              most_cparn_num = 35
          end if
          call glfor_mem( ema_data ) 
      end if

      call glb_out_map

*     This is needed on the HP
      inquire(file=glb_sol_file,iostat=ierr, number=num)
      if( num.gt.0 ) close(num)

      call rw_glb_covar('O', dcb, ema_data(icov_parm) )
      call rw_glb_covar('L', dcb, ema_data(icov_parm) )

*     Now read the glorg commands

      call read_glorg_commands( ema_data(isol_parm), options )
 
***** Thats all
      return
      end
 
CTITLE GLORS
 
C                        ,Stablizing segment
      subroutine glors

      implicit none  
 
*     This segment stablizes the solution.  At the moment it just
*     fixes the translation, and the RA origins
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glorg_common.h'
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   gnum_used   - Number of global sites used in origin
*   i,j,k       - Loop counter
*   iout        - Output LU
*   options     - options passed through runstring
* MOD TAH 980417: Introduced iterated stabalization.
*   it          - Iteration count for stabalization.
 
      integer*4 ema_data(1), i,j, iout,
     .    options, it
 
*   pmu_parts(3,3,max_glb_sites)  - Polar motion UT1 partial derivatives

      real*8 pmu_parts(3,3,max_glb_sites)
 
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
***** Stabilize the translation origin (position and velocity)
 
*     Stabilize the positions
 
      do i = 1, gnum_sites
*                     ! XYZ
          do j = 1,3
              call pmu_part(i, pmu_parts(1,j,i), apr_val_site(1,1,i),
     .            j, gut1_apr(1))
          end do
      end do
 
***** See if new(new) orientation system:
*     Initialize the site variance to be all the same.
      do i = 1, gnum_sites
         stab_site_var(i,1) = 1.d0
         stab_site_var(i,2) = 1.d0
      end do

      do i = 1, max_glb_site_wrds
         use_pos(i) = use_sites(i)
         use_rat(i) = use_sites(i)
      end do

      write(iout, 300) (1-stab_rel)*100, stab_rel*100, stab_nsig,
     .                 cnd_hgt_var, 
     .                 (stab_min_dh(i), stab_min_rms(i), 
     .                  stab_min_dne(i), i=1,2),
     .                  use_ratio(1), use_ratio(2) 
 300  format(/,' Stabilization with ',F5.1,'% constant, ',F5.1,
     .     '% site dependent weighting.',/,
     .     ' Delete sites with ',F5.1,'-sigma condition.',/,
     .     ' Height variance factor ',F10.2,' Position,',F10.2,
     .     ' Velocity',/,
     .     ' For Position: Min dH sigma ',F6.4,' m;    Min RMS ',
     .       F6.4,' m,    Min dNE sigma ',F7.5,' m'/,
     .     ' For Velocity: Min dH sigma ',F6.4,' m/yr; Min RMS ',
     .       F6.4,' m/yr, Min dNE sigma ',F7.5,' m/yr',/,
     .     ' Sigma Ratio to allow use: Position ',F6.2,
     .     ' Velocity ',F6.2)

      do it = 1, num_stab_iter
          if( cnd_pos_bits.ne.0 ) then
              call apply_cond_full( iout, parn_site(1,1,1), 2,
     .               pmu_parts, gnum_sites, ema_data(icov_parm), 
     .               ema_data(isol_parm), apr_val_site(1,1,1),
     .               cond_var(1,1),  num_glb_parn, 
     .               ema_data(isol_parm+2*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+14*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+28*num_glb_parn),
     .               use_pos  , cnd_pos_bits, cnd_hgt_var(1),
     .               use_ratio(1), 'Position', gsite_names, 
     .               stab_site_var(1,1), stab_site_err,
     .               stab_nsig, stab_rel, stab_min_dh(1), 
     .               stab_min_rms(1), stab_min_dne(1), 
     .               it, options, list_file)
          end if

*         Do the rates
          if( cnd_rat_bits.ne.0 ) then

              call apply_cond_full( iout, parn_site(1,2,1), 2, 
     .               pmu_parts, gnum_sites, ema_data(icov_parm), 
     .               ema_data(isol_parm), apr_val_site(1,1,1),
     .               cond_var(1,2),  num_glb_parn, 
     .               ema_data(isol_parm+2*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+14*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+28*num_glb_parn),
     .               use_rat  , cnd_rat_bits, cnd_hgt_var(2),
     .               use_ratio(2), 'Velocity', gsite_names,
     .               stab_site_var(1,2), stab_site_err,
     .               stab_nsig, stab_rel, stab_min_dh(2), 
     .               stab_min_rms(2), stab_min_dne(2), 
     .               it, options, list_file)
          end if

*         If we are not at the last iteration, re-read the loose 
*         covarinace matrix and re-update the parameter changes.
          if ( it.lt.num_stab_iter ) then
              call rw_glb_covar('L', dcb, ema_data(icov_parm) )
              call glorg_upd_apr( ema_data(isol_parm) )
          end if
      end do
 
***** Now check that we covariance matrix is OK
 
      call check_covar(ema_data(icov_parm), num_glb_parn, num_glb_parn)
 
***** Thats all
      return
      end
 
CTITLE GLORW
 
C                        ,Stablizing segment
      subroutine glorw

      implicit none  
 
*     This segment writes out the final solution to the output
*     device with the options passed through the runstring
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   i           - Loop counter
*   iout        - Output LU
*   options     - options passed through runstring
 
      integer*4 ema_data(1), iout, options
 
*   kbit        - Bit checking function
 
      logical kbit
 
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
***** Start,
 
      call write_glb_header( iout, options,
     .       ema_data(icov_parm), ema_data(isol_parm))
 
*     Now the parameter estimates
 
*     Write out summary velocity field
      IF( KBIT(OPTIONS,5) ) THEN
      call write_glb_sum( iout, options, use_rat , ema_data(icov_parm),
     .                                   ema_data(isol_parm),
     .                                   list_file )
      END IF

*     Write out summary position field
      IF( KBIT(OPTIONS,14) ) THEN
      call write_glb_pos( iout, options, use_pos , ema_data(icov_parm),
     .                                   ema_data(isol_parm),
     .                                   list_file )
      END IF

****  Write log estimates
      call write_log_sum(iout, options, use_rat, 
     .                   ema_data(icov_parm), ema_data(isol_parm) )

*     Now write out the full solution

      call write_glb_params( iout, options, ema_data(icov_parm),
     .                                      ema_data(isol_parm) )
 
*     Write out baseline lengths
 
      IF( KBIT(OPTIONS,2) ) THEN
      call write_glb_basel( iout, options, 1, ema_data(icov_parm),
     .                                        ema_data(isol_parm) )
      call write_glb_bcomp( iout, options, 1, ema_data(icov_parm),
     .                                        ema_data(isol_parm) )
      END IF
 
*     Write out baseline rate information
      if( kbit(options,3) ) then
          call write_glb_bcomp( iout, options, 2, ema_data(icov_parm),
     .                                            ema_data(isol_parm) )
      end if
 
*     Now write out correlations
 
      call write_glb_corel( iout, options,  ema_data(icov_parm),
     .                                      ema_data(isol_parm) )

*     Write a blank line at the end of output
      write(iout,'(1x)')
 
***** Thats all
      return
      end
 
CTITLE DECODE_GLORG_RUN
 
      subroutine decode_glorg_run( iout, options, outfile, 
     .           glorg_command_file, glb_com_file, com_orgopt )
 

      implicit none 
 
*     Routine to decode runstring
 
*   ierr            - Error during conversion
*   iout            - Ouput Lu
*   len_run         - Length of runstring
*   LogLu           - HP function for users LU
*   options         - Options for output
*   rcpar           - HP function to read runstring
*   lenfn           - Length of output file name
 
      integer*4 ierr, iout, len_run, LogLu, options,
     .    rcpar, indx, lenfn
     
*   kbit         - Check bit status
      logical kbit       
 
*   glb_com_file    - Name of the common file
*   glorg_command_file  - OPtional name of command file
*   outfile         - Output file name 
*   com_orgopt  - Optional string at start of lines 
 
      character*(*) outfile, glb_com_file, glorg_command_file,
     .              com_orgopt
 
*   runstring       - Runstring for integer values
      character*128 runstring
 
***** Get the Lu for output
 
      lenfn = rcpar(1, outfile)
 
*     Get options
      len_run = rcpar(2, runstring)
      if( len_run.gt.0 ) then
          indx = 1
          call decode_prt_opt(runstring, indx, options )
      else
          options = 0
      end if
      
* MOD TAH 970909: Moved file open to after runstring to see
*     if we need to erase first (ERAS option).      
      if( lenfn.gt.0 ) then
          iout = 200
          if( kbit(options,17) ) then
              call open_lu(iout, outfile, ierr, 'unknown')
          else
              call open_lu(iout, outfile, ierr, 'append')
          endif 
          call report_error('IOSTAT',ierr,'open_lu',outfile,
*                                                 ! Dont kill
     .                      0,'Decode_GLORG_run')
          if( ierr.ne.0 ) then
              iout = LogLu( ierr )
          end if
*                     ! use user's terminal
      else
          iout = LogLu(ierr)
      end if
 
*     Get optional command file name
      len_run = rcpar(3, glorg_command_file)
      if( len_run.eq.0 ) then
*         Stop running.
          call proper_runstring('glorg.hlp','glorg',-1)
      end if
 
*     Get optional comfile name
      len_run = rcpar(4, glb_com_file)
      if( len_run.eq.0 ) then
          glb_com_file = ' '
      end if

****  See if command read option passed
      len_run = rcpar(5, com_orgopt)


 
***** Thats all
      return
      end
 
CTITLE APPLY_ROT
 
      subroutine apply_rot( parn, dim, sol_parm, pmu_parts,
     .           gnum_sites, pmu_changes )
 
      implicit none  
 
*     Routine to apply a rotational change
 
*   i,j,k       - Loop counter
*   iel         - position in sol_parm
*   dim         - Dim of parn
*   gnum_sites  - Number of sites
*   parn(3,dim,1)   - Parmeter numbers
 
      integer*4 i,j, iel, dim, gnum_sites, parn(3,dim,1)
 
*   corr        - Correction to rate
*   sol_parm(1) - Parameter estimates
*   pmu_parts(3,3,1)  - PMU partials
*   pmu_changes(3)  - Changes to pmu values (mas)
 
      real*8 corr, sol_parm(1), pmu_parts(3,3,1), pmu_changes(3)
 
 
*     Loop over sites
      do i = 1,gnum_sites
*                     ! Components
          do j = 1,3
              if( parn(j,1,i).ne.0 ) then
                  iel = parn(j,1,i)
                  corr        = pmu_parts(1,j,i)*pmu_changes(1)
     .                        + pmu_parts(2,j,i)*pmu_changes(2)
     .                        + pmu_parts(3,j,i)*pmu_changes(3)
                  sol_parm(iel) = sol_parm(iel) - corr
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
 
CTITLE APPLY_TRANS
 
      subroutine apply_trans( parn, dim, sol_parm, gnum_sites,
     .                        rate_trans )

      implicit none  
 
*     Routine to apply a translation change
 
*   i,j,k       - Loop counter
*   iel         - position in sol_parm
*   dim         - Dim of parn
*   gnum_sites  - Number of sites
*   parn(3,dim,1)   - Parmeter numbers
 
      integer*4 i,j, iel, dim, gnum_sites, parn(3,dim,1)
 
*   corr        - Correction to rate
*   sol_parm(1) - Parameter estimates
*   rate_trans(3)   - Translation
 
      real*8  sol_parm(1), rate_trans(3)
 
*     Loop over sites
      do i = 1,gnum_sites
*                     ! Components
          do j = 1,3
              if( parn(j,1,i).ne.0 ) then
                  iel = parn(j,1,i)
                  sol_parm(iel) = sol_parm(iel) - rate_trans(j)
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
 
CTITLE GLOREF
 
C                        ,Equate/force segment
      subroutine gloref

      implicit none  
 
*     This segment equates and forces adjustments to parameters
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glorg_common.h'
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   options     - options passed through runstring
 
      integer*4 ema_data(1), iout, options
 
      common / progcon / iout, options
 
      common / globc_ema / ema_data
 
 
****  See if we need to rotate covariance matrix into NEU before
*     applyiing condintions

      if( equate_loc ) then 
          write(iout,'(a)')'Rotating into local coordinates for equates'
          call cov_xyz_neu( 1, ema_data(icov_parm),
     .                       ema_data(isol_parm))
      end if

*     See if have any parameters to be equated.
      if( num_equates.gt.0 ) then
          call equate(iout, ema_data(icov_parm), ema_data(isol_parm))
      end if

*     See if we have parameters we want to force.
      if( num_force.gt.0 ) then
          call force(iout, ema_data(icov_parm), ema_data(isol_parm))
      end if


*     if we equate in local frame, rotate back to XYZ
      if( equate_loc ) then
          call cov_xyz_neu( -1, ema_data(icov_parm),
     .                       ema_data(isol_parm))
      end if

****  Thats all
      return
      end

      subroutine wr_solp ( seg, sol_parm )

      implicit none 

      real*8 sol_parm(3)
      character*(*) seg

      write(*,'(a,1x,3F16.4)') seg, sol_parm(1),sol_parm(2),sol_parm(3)
      end  

*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include 'glorg_bd.f' 
