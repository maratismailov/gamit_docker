      subroutine glbak(ms_type)
 
      implicit none 

*     This the back global Kalman filter program.  It is made
*     up of three segments plus the main:
*     GLBAO -- Opens/reads/closes files
*     GLBAP -- Computes O-C and partial derivatives
*     GLBAF -- Actually does the Kalman filtering.
*     GLBAR -- Computes postfit residials (if option selected)
*     GLBAW -- Writes out the back solution results
*
*                                  2:01 PM  TUE., 18  AUG., 1987

* MOD TAH 190525: Added frame stabilization to back solution when
*     glorg is used in the globk solution.  Run string extended to
*     pass glorg_command file and the COM_org_opt string
*     Runsting:
*     % glbak [com_file] <glorg command file> <option string>
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
      include '../includes/glb_hdr_def.h'

* ms_type   - Character string which indicates that glfor is a program
*          (MAIN) or a subroutine (SUBR)

      character*(4) ms_type

* MAIN PROGRAM VARIABLES 
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
 
*   iscr_real(max_glb_parn) - Scracth common area
 
*   len_read    - Actual length of record read
*   pcontrol    - Segment control used by GLFOO.
*               - pcontrol = 1, open the files we need
*               - pcontrol = 2, open/read new global solution
*               - pcontrol = 3, open/read new global solution, and
*               -               read previous covariance matrix
*               - pcontrol = 4, re-read current global solution
*               - pcontrol = 5, open and read solution and save
*               -               GLOBK common
 
 
      integer*4 cov_dcb(16), ema_data(1), i, j, ierr,
     .    iscr_real(max_glb_parn), pcontrol, trimlen
 
*   timr        - LIBHS function to return run time
*   running_time  - Incremental running time for program
 
      real*4 timr, running_time

*   ephem_prev  - Previous satellite empheris JD
*   secs        - USed in systime call.
*   glb_var     - REscaling factor on covaraince matrix
*   glb_diag    - Diagonal scaling 

      real*8 ephem_prev, secs, glb_var, glb_diag 

*   apply_sd   - Indicates that we should apply SD corrections
*   apply_nut  - Indicates that we should apply nutation corrections

      logical apply_sd, apply_nut

* glb_used - Logical to denote if global file used
* MOD TAH 190528: Addded option
* first_soln -- Set true if first solution (no average formed)
 
      logical glb_used,  first_soln
 
 
*   scr_real(max_glb_parn)  - Scratch common area
*   final_gep_expt - Actual referr to date for the analysis.
*                    This is gepoch_expt at the end of the
*                    forward solution which we save here.
 
      real*8 scr_real(max_glb_parn), final_gep_expt
      
* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file.
*   for_chi, bak_chi -- Forward and backwarsd chi**2/f
      real*4 for_chi, bak_chi
 
*   cname       - Name of next solution equivalanced to iname
 
      character*(sort_recl) cname
 
*   cr          - Contains carriage return
 
      character*1 cr
      
*   outline     - Line to be output to the status report
*   progname    - Name of program

      character*256 outline, progname
      integer*4 lenprog, rcpar, FmpClose

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters


      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
*     scr_real, iscr_real     - Unlabeled common
 
      common scr_real, iscr_real

*                        approximately 1900.
      data  ephem_prev / 2415020.d0 /
      data I8 / 1 /

 
***** Get the runstring, and open the GLOBK common, and the sort file
 
      cr = char( 13 )

      running_time = timr(0)
      lenprog = rcpar(0,progname)
 
      pcontrol = 1
      call GLBAO (ms_type) 
 
*     Open up the loglu
      log_unit = 200
      call open_lu(log_unit, glb_log_file, ierr, 'append')

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
     .              gnum_sites, gepoch_expt )

*     Set to show co-seismic offset not applied yet
      do i = 1, num_eq
          eq_co_applied(i) = .false.
      end do

*     Set bit 32 of the print options to show the baseline
*     routines that baselines should only be written if the
*     station is really used in this experiment.

      call sbit(bak_opts,32,1)
      call sbit(mul_pmu_opt,32,1)
 
***** Now loop over the experiments in the global solution, adding
*     each one to the solution.

      final_gep_expt = gepoch_expt
 
      do i = num_glb_sol, 1, -1

          if( i.eq.num_glb_sol) first_soln = .true.
          if( i.ne.num_glb_sol) first_soln = .false.

           
          read(100, iostat=ierr, rec=i) cname, glb_var,
     .                                for_chi, bak_chi
          call report_error('IOSTAT',ierr,'read',sort_file,0,
     .                      'GLFOR')


* MOD TAH 970706: Check if we are still using this file
          if( ierr.eq.0 ) then
             if( ichar(cname(1:1)).gt.128 ) then 
                 ierr = -2

*                Read the previous global solution records but
*                check first if the file has been opened yet. 
                 if( i.eq.num_glb_sol )  then
                     call rw_glb_covar('O', cov_dcb, 
     .                             ema_data(icov_sav_parm))
                     call rw_glb_covar('L', cov_dcb, 
     .                             ema_data(icov_sav_parm))
                 else   
                     call rw_glb_covar('P', cov_dcb, 
     .                             ema_data(icov_sav_parm))
                 end if 
              end if 
          end if

* MOD TAH 970706: Breakup the scaling factors
          if( glb_var.lt.0 ) then
              glb_diag = (abs(glb_var*1.d3) - int(abs(glb_var*1.d3)))*
     .                   1000.d0
              glb_var = abs(glb_var) - glb_diag/1000.d3
              if( glb_diag.lt.1.d0 ) glb_diag = 1.d0
          else  if ( glb_var.eq.0 ) then
              glb_var = 1.d00
              glb_diag = 1.d0
          else
              glb_diag = 1.d0
          end if
 
*                                     ! Continue, OK
          if( ierr.eq.0 ) then
*                                     ! Save name for next segment
              glb_inp_file = cname
 
              if( i.eq.num_glb_sol ) then
*                                         ! Open the next global solution and
                  pcontrol = 2
*                                         ! read last covariance matrix
              else
*                                         ! Open the next global solution and
                  pcontrol = 3
*                                         ! read previous covariance matrix
              end if
              call GLBAO (ms_type) 
 
*****         Get epoch change between experiments (gepoch_epxt
              gepoch_prev = gepoch_expt
              gepoch_expt = cepoch_expt
              if( i.eq.num_glb_sol ) gepoch_prev = gepoch_expt
 
*                                                          ! Years
              deltat = (gepoch_expt-gepoch_prev)/365.25d0
*             We need a small offset so that the earthquake modules
*             know which direction time is going.
              if( sort_direction.eq. 1.and.deltat.eq.0 ) deltat=-1.d-12

              deltaephem = (csvs_epoch  - ephem_prev ) /365.25d0

*             Make sure we have satelites in this particular global
*             before updating the ephem epoch.  This allows VLBI
*             and GPS solutions to be processed together.
              if( cnum_svs.gt.0 ) then
                  ephem_prev = csvs_epoch
              else
                  deltaephem = 0.d0
              end if
              if( abs(deltaephem).lt.1.d-6 ) then
                  deltaephem = 0.d0
              else
                  if( gnum_svs.gt.0 .and. cnum_svs.gt.0 ) then
                      write(*,90) deltaephem
  90                  format(/' Updating SV ephemeris epoch by ',
     .                         f10.4,' years')
                  end if
              end if
 
              running_time = timr(-1)
 
              write(*,100) i,(isol_obs+2*cnum_parn-icov_sav_parm+128)*
     .                     4.0/1024./1024., running_time, glb_var,
     .                     glb_diag 
              if( log_unit.ne.6 ) then
                  write(log_unit,100) i,
     .                     (isol_obs+2*cnum_parn-icov_sav_parm+128)*
     .                     4.0/1024./1024., running_time, glb_var,
     .                     glb_diag 
              end if
             
  100         format(' Global ',i4,' using ',f7.1,' Mb. Running time ',
     .                f8.2,' Scaling by ',F10.3,F11.8 )
     
              write(outline,110) i,(isol_obs+2*cnum_parn-icov_parm+128)*
     .                    4.0/1024./1024., glb_var, glb_diag 
  110         format('File ',i4,' using ',f8.2,
     .               ' Mb. Scaled by ',F10.3,F11.8)
              call report_stat('status',progname,'glbak',glb_inp_file,
     .                         outline,0)

*****         Now scale the covariance matrix

              if( glb_var.ne.1.0d0 ) then
                  call dwsmy(glb_var, ema_data(icov_obs), 1,
     .                 ema_data(icov_obs), 1, (I8*cnum_parn)*cnum_parn)
              end if

*             Run down the diagonal scaling only the diagonal.
              if( glb_diag.ne.1.d0 ) then
                  call dwsmy(glb_diag, ema_data(icov_obs), cnum_parn+1,
     .                  ema_data(icov_obs), cnum_parn+1, cnum_parn)
              end if 
 
*****         Get the new apriori's for this experiment.  (Currently it
*             is satellite orbits, but later could be polar motion UT1 values
*             as well.
* MOD TAH 070926: Added explicit satellite epoch.
              call glb_upd_apr( cepoch_expt, cepoch_expt,  .true. , 
     .                          ema_data(isol_obs), .true.)

*****         Get any time markov elements
              call get_mar_svs( cepoch_expt )

*             Now compute nutation series corrections.
              if( apply_nut ) then
                  call compute_std( cepoch_expt, vnut_ang_apr )
              else
                  vnut_ang_apr(1) = cnut_ang_apr(1)
                  vnut_ang_apr(2) = cnut_ang_apr(2)
              end if

              call get_sd_mid( cepoch_expt, apply_sd, gnum_sites)

* MOD TAH 190528: Get the non-secular contribution at this tme
*             Compute the contributions to the postition at the epoch of
*             output for the nonsecular terms
              do j = 1, gnum_sites
                  call eval_nonsec(j, cepoch_expt, num_nonsec,
     .                    param_nonsec, apr_val_nonsec, 
     .                    cont_nonsec(1,j),0)
              end do

* MOD TAH 050927: Get the long site names from the head.snx file (either
*     local or in $HELP_DIR or gg/tables)
              call get_long_names
*
*****         Now get the O-C and partials
 
              call GLBAP
 
*****         Now actually do the Kalman filtering
 
              call GLBAF(glb_used, first_soln)

*             Now save the chi**2 values              
*                 If the global file was not used, then update the list file
*                 to denote this
* MOD TAH 980519: Write out the chi**2/f to the srt_file
                  if( .not.glb_used ) then
C -- wrong order of brackets I think? (simon) cname(1:1) = char(ichar(cname(1:1)+128))
                      cname(1:1) = char(ichar(cname(1:1))+128)
                  endif
                  if( glb_diag.ne.1.d0 ) then
                      glb_var = -(glb_var+glb_diag/1000.d3)
                  end if

* MOD TAH 980519: Get the forward chi**2                   
                  bak_chi = dchi_save
* MOD TAH 190621: Save for_chi, bak_chi for use in glsave
*                 writing GLX file during back solution/
                  for_chi_save = for_chi
                  bak_chi_save = bak_chi
                  
                  write(100,iostat=ierr,rec=i) cname, glb_var,
     .                                for_chi, bak_chi
                  
                  call report_error('IOSTAT',ierr,'writ','SRT file'
     .                                 ,0,'glfor') 

*****         See if we need to do frame realization before output.
*             Unit 202: Output bak_file.
              if( num_stab_iter.gt.0 ) then

                  call GLBAS(202)     ! Frame realozaton code from glorg

              endif



 
*****         Now re_read input global solution so that we can compute
*             post-fit residuals
              if( compute_glb_res .and. glb_used ) then
                  pcontrol = 4
                  call GLBAO(ms_type) 
 
*                                            ! Redo prefit residuals
                  call GLBAP 
 
*                                            ! Compute residuals
                  call GLBAR
 
              end if

*****         Save the run time for this solution.
        
              call systime( grun_time, secs)
              grun_time(6) = secs
              grun_time(7) = 0
       
*****         Write out the current solution.  Force the multiday
*             polar motion estimates to be printed on the first solution.
              if( i.eq. num_glb_sol ) call sbit(mul_pmu_opt,31,1)
              call GLBAW 
 
              call GLBAV

*             Out put the PMU estimates.
              call print_pmu( 6 ,
     .         ema_data(icov_parm), ema_data(isol_parm),
     .                 num_glb_parn)
              if( log_unit.ne.6 ) call print_pmu(log_unit,
     .             ema_data(icov_parm), ema_data(isol_parm),
     .             num_glb_parn)

* MOD TAH 190620: If we ran glorg, see if we are are saving GLX file.
              if( num_stab_iter.gt.0 ) then

* MOD TAH 190619: See if we now will save this solution as a GLX file.
                  if( trimlen(glb_out_file).gt.0 ) then
*                     Save the GLX file for this day
                      write(*,220) gepoch_expt
 220                  format('Saving Epoch ',F10.2,' to binary hfile')
                      call glsave('BACK','NO')
                   endif
              endif
 
*                         ! No error getting next name
          end if
*                         ! Looping over the eperiments
      end do
 
***** Now go back and read the last solution in the list so that
*     we can save the wobble,UT1 and nutation angles for the correct
*     experiment
 
      read(100, iostat=ierr, rec=num_glb_sol ) cname
      if( ichar(cname(1:1)).gt.128 ) then
          cname(1:1) = char(ichar(cname(1:1))-128)
      end if 
      glb_inp_file = cname 

*     Re-instate the actual epoch for the solution.
      gepoch_expt = final_gep_expt

*                             ! Re-Read solution and save common
      pcontrol = 5
      call GLBAO (ms_type) 

* MOD TAH 961028: close the cov_dcb file (reopened later)
      ierr = FmpClose( cov_dcb )

***** Thats all
      close(202)
      close(100)
      if( log_unit.ne.6 ) close(log_unit)
      end
 
 
      subroutine GLBAO(ms_type) 

      implicit none 
 
*     This segment of the back global Kalman filter program will:
*     GLFOO -- Opens/reads/closes files
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
      include '../includes/glb_hdr_def.h'
 
* ms_type   - Character string which indicates that glfor is a program
*          (MAIN) or a subroutine (SUBR)

      character*(4) ms_type

* VARIABLES 

*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   len_read    - Length read with VREAD
*   pcontrol    - Segment control used by GLFOO.
*               - pcontrol = 1, open the files we need
*               - pcontrol = 2, open and read new global solution
*               - pcontrol = 3, open and read new global solution and
*               -               read previous covariance matrix
*               - pcontrol = 4, re-read current global solution
*               - pcontrol = 5, open/read solution and save globk
*               -               common
*   vma_bytes   - Number of bytes needed.
* MOD TAH 150524: Introduced integer*8 version for large 
*     solutions >20000 parameters
 
      integer*4 cov_dcb(16), ema_data(1), ierr,
     .    len_read, pcontrol, trimlen, fmpclose
      integer*8 vma_bytes

      logical kbit

      common / progcon / pcontrol, cov_dcb
 
      common / globc_ema / ema_data
 
***** Depending on pcontrol do the desired operation
 
*                               ! Get the runstring, and
      if( pcontrol.eq.1 ) then
*                               ! open the GLOBK common, and the sort file
 
***       Decode the runstring
          com_orgopt = ' '
          glorg_command_file = ' '
          if( ms_type(1:4).eq.'MAIN' ) then
              call decode_glbak_run( glb_com_file, glorg_command_file,
     .                           com_orgopt ) 
 
***           Open the globk common, and read the segments
              call open_globk( glb_com_file, ierr)
              call report_error('FmpOpen',ierr,'open',glb_com_file,1,
     .                      'Segment GLBAO')
 
              call rw_globk_block('R','CONTROL',glb_control,ierr)
              call report_error('FmpRead',ierr,'read','control block',1,
     .                      'Segment GLBAO')
              istart_vma = 0
              call rw_globk_block('R','MARKOV', glb_markov ,ierr)
              call report_error('FmpRead',ierr,'read','markov block',1,
     .                      'Segment GLBAO')
 
              call rw_ema_block('R','EMA',glb_ema,ierr)
              call report_error('VREAD',ierr,'read','ema block',1,
     .                      'Segment GLBAO')

* MOD TAH 190525: See if we need to update comopt after reading the common files.
              if( trimlen(glorg_command_file).eq.0 ) then
                   glorg_command_file = glr_cmd_file ! use name from globk common
              endif
              if( trimlen(com_orgopt).eq.0 ) then  ! Use OPTION string from
*                                                    globk common
                   com_orgopt = comopt          
              endif 
          else    ! Use  glr_cmd_file ! use name from globk common     
                   glorg_command_file = glr_cmd_file ! use name from globk common
                   com_orgopt = comopt          
          end if

* MOD TAH 190525: If the glorg_command file name has been passed, open and
*         read the commands associated with frame realization
          num_stab_iter = 0   ! This value will be used to test if frame 
*                             ! realization should be done.
          if( trimlen(glorg_command_file).gt.0 ) then
              call read_glbak_orgcmd(6) 
          endif
 
***       Now open the sort file
* MOD TAH 980519: Increased the record length by 8 bytes for for and
*         bak chi**2/f Saving
          open(100,file=sort_file, iostat=ierr, status='old',
     .         access='direct', recl=sort_recl+8+8)
          call report_error('IOSTAT',ierr,'open',sort_file,1,
     .                      'Segment GLFOO')
 
***       Initalize the covariance matrix for the solution
 
          call get_max_glb_deriv
          call glfor_mem( ema_data )
 
          call glb_bak_map( cnum_parn, vma_bytes ) 
 
***       Open the output file for the temporary covariance matrices
          call rw_glb_covar('O', cov_dcb, ema_data(icov_sav_parm))
 
****      Read the last solution, and then use 1000. times the variances
*         as the apriori variances for this solution
          call rw_glb_covar('L', cov_dcb, ema_data(icov_sav_parm) )
 
          call init_bak_cov ( ema_data(icov_parm), ema_data(isol_parm),
     .         ema_data(icov_sav_parm), ema_data(isol_sav_parm))
 
****      Open the back solution output file
          if( kbit(bak_opts,17) ) then
              call open_lu(202,glb_bak_file, ierr,'unknown') 
          else
              call open_lu(202,glb_bak_file, ierr,'append')
          endif
          call report_error('IOSTAT',ierr,'open',glb_bak_file, 1,
     .                      'GLBAO')
 
*****     Now write out the description of this solution
          write(202,'(a)', iostat=ierr) gdescription
!         call dump_cov_parm( 'PCONTROL 1', ema_data(icov_parm), 
!    .          num_glb_parn)
 
      end if
 
 
***** See if global solution file operation.
      if( pcontrol.eq.2 .or. pcontrol.eq.3  ) then
 
*****     Open the current global solution file, and read the header, and
*         names
          call FmpOpen(cglb_dcb,ierr,glb_inp_file,'RO',1)
          call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'Segment GLFOO')
 
***       Read the header
          call rw_glb_header('R',ierr)
          call report_error('FmpRead',ierr,'read',glb_inp_file,0,
     .                      'Segment GLFOO')

*****     See if we need to set up a new group of multi-day polar
*         motion/UT1 values

          call setup_mul_pmu

 
****      Save wobble, ut1 and nutation angle
          call save_glb_apr
 
*         MOD TAH 150825: Get the PRN to SVN information
          call read_svinf_rec

****      Get the names and local to global conversions
          call rw_names_block('R')

****      Update names due to earthquakes
          call eq_name_change('NO') 
 
****      Get the parameter codes for this solution
 
          call get_param_codes( gpar_codes )
 
*                         ! See if nutation/PMU estimated
          call scan_parm
 
*****     Map the ema area usage and read the solution and covariance
*         matrix into ema
 
          call get_max_glb_deriv
 
          call glb_bak_map( cnum_parn, vma_bytes )
 
          call readd(cglb_dcb,ierr,ema_data(icov_obs),
     .               128*cnum_par_vals, len_read, crec_par_vals)
 
          call report_error('VREAD',ierr,'read','global solution',0,
     .                      'Segment GLFOO')
 
 
****      Now get the covariance matrix from the forward solution
          if( pcontrol.eq.2 ) then
 
C             Read last was already done above
C             call rw_glb_covar('L', cov_dcb, ema_data(icov_sav_parm))
          else
              call rw_glb_covar('P', cov_dcb, ema_data(icov_sav_parm))
!         call dump_cov_parm( 'PCONTROL 3', ema_data(icov_parm), 
!    .          num_glb_parn)
          end if
      end if
 
****  See if we should re-read the input (observation) covariance matrx
 
      if( pcontrol.eq.4 ) then
 
*         Get the parameter codes back
 
          call get_param_codes( gpar_codes )
 
          call readd(cglb_dcb,ierr,ema_data(icov_obs),
     .               128*cnum_par_vals, len_read, crec_par_vals)
 
          call report_error('VREAD',ierr,'read','global solution',0,
     .                      'Segment GLFOO')
          ierr = fmpclose(cglb_dcb)
 
      end if
 
***** See if open, read solution and save common
 
      if( pcontrol.eq.5 ) then
 
***       Open global solution file
          call FmpOpen(cglb_dcb,ierr, glb_inp_file,'RO', 1)
 
***       Read header
          call rw_glb_header('R',ierr)
 
***       Now save aprioris
          call setup_mul_pmu 
          call save_glb_apr
          call read_pmu_inp

          ierr = fmpclose(cglb_dcb)
 
***       Write out updated parts of the globk common if main
          if( ms_type(1:4).eq.'MAIN' ) then 
              call rw_globk_block('W','CONTROL',glb_control,ierr)
              call report_error('FmpWrite',ierr,'writ','control block',
     .                      1,'Segment GLBAO')
              call rw_globk_block('W','MARKOV', glb_markov ,ierr)
              call report_error('FmpWrite',ierr,'writ','markov block',
     .                      1,'Segment GLBAO')
 
              call rw_ema_block('W','EMA',glb_ema,ierr)
              call report_error('VWRIT',ierr,'writ','ema block',1,
     .                      'Segment GLBAO')
          end if

*         If this not MAIN, see we should the common
          if( ms_type(1:4).ne.'MAIN' .and. 
     .        trimlen(glb_com_file).gt.0 ) then
             call rw_globk_common('W')
          end if
 
      end if
 
***** Thats all
      return
      end
 
CTITLE GLBAP
 
      subroutine glbap
 
      implicit none 

 
*     This segment of the forward global Kalman filter program:
*     GLBAP -- Computes O-C and partial derivatives
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Contols segment GLFOO (not used here)

*   Looks_OK    - Returns True if all OK so far
 
      integer*4 cov_dcb(16), ema_data(1), pcontrol

      logical looks_ok
 
      common / progcon / pcontrol, cov_dcb
 
      common / globc_ema / ema_data
 
***** Compute the O-C values (Note the ema usage was mapped in GLFOO
*     when this global solution was first read)
 
      call glb_o_minus_c( ema_data(isol_obs), looks_ok )
 
***** Now compute the partial derivatives which we need
 
      call glb_partials( ema_data(ipart_pnt), ema_data(ia_part) )
      
***** Now check the sol_obs array and see if we should suppress
*     bad elements

      call check_param( ema_data(isol_obs), cnum_parn, indx_pnt )
  
***** Thats all
      return
      end
 
CTITLE GLBAF
 
      subroutine glbaf(glb_used, first_soln) 
 
      implicit none 

 
*     This segment of the forward global Kalman filter program:
*     GLBAF -- Does the Kalman filtering
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(1) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Controls segment GLFOO (not used here)
 
      integer*4 cov_dcb(16), ema_data(1),  pcontrol

* glb_used - Logical to denote if global file used
* first_soln - Set true for first day so no averaging to get cov_sav_parm 
      logical glb_used,  first_soln 
  
 
      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
***** Call the routines to add in the next solution into the
*     Kalman filer solution
 
*     First compress the parameters from the current SOLVK solution
*     which were not used in the solution
 
      call remove_params ( ema_data(icov_obs),   ema_data(isol_obs),
     .                     ema_data(ipart_pnt),  ema_data(ia_part)  )
 
      call glb_bak_filter( ema_data(icov_sav_parm),
     .                     ema_data(isol_sav_parm), ema_data(ijmat),
     .                     ema_data(icov_parm),     ema_data(isol_parm),
     .                     ema_data(irow_copy),
     .                     ema_data(ipart_pnt),     ema_data(ia_part),
     .                     ema_data(itemp_gain),    ema_data(ikgain),
     .                     ema_data(icov_obs),      ema_data(isol_obs),
     .                     glb_used,  first_soln)
 
 
***** Thats all
      return
      end
 
 
CTITLE GLBAW
 
      subroutine glbaw
 
      implicit none 


*     This segment of the forward global Kalman filter program:
*     GLBAW -- writes out the results
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i,j         - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Controls segment GLFOO (not used here)
 
      integer*4 cov_dcb(16), ema_data(1), pcontrol
 
      common / progcon / pcontrol, cov_dcb
 
      common / globc_ema / ema_data
 
***** write solution
 
      call glb_bak_writ1( ema_data(icov_sav_parm),
     .                    ema_data(isol_sav_parm),
     .                    ema_data(icov_obs), ema_data(isol_obs),
     .                    ema_data(irow_copy) )
 
***** Thats all
      return
      end
 
 
CTITLE GLBAV
 
      subroutine glbav
 
      implicit none 

 
*     This segment of the forward global Kalman filter program:
*     GLBAW -- writes out the results
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i,j         - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Controls segment GLFOO (not used here)
*   kbit        - Checks status of bits
 
      integer*4 cov_dcb(16), ema_data(1), pcontrol

      logical kbit
 
      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
***** write solution
 
      call glb_bak_writ2( ema_data(icov_sav_parm),
     .                    ema_data(isol_sav_parm),
     .                    ema_data(icov_obs), ema_data(isol_obs),
     .                    ema_data(irow_copy) )

*     Now write out baseline lengths and components
      if( kbit( bak_opts,2) ) then
          call write_glb_basel( 202 ,bak_opts, 1, 
     .                          ema_data(icov_sav_parm),
     .                          ema_data(isol_sav_parm))
          call write_glb_bcomp( 202, bak_opts, 1, 
     .                          ema_data(icov_sav_parm),
     .                          ema_data(isol_sav_parm))
      end if
 
***** Thats all
      return
      end
 
CTITLE GLBAR
 
      subroutine glbar
 
      implicit none 

 
*     This segment of the forward global Kalman filter program:
*     GLBAR -- Compute postfit residuals
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i,j         - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   loc_parn(3,max_glb_sites)  - local parameter number array
*               - generated for the site origin.
*   pcontrol    - Controls segment GLFOO (not used here)
*   used(8)     - Indicates to use the site (bit mapped) in the
*               - origin (Currently defaults to using all sites)
 
      integer*4 cov_dcb(16), ema_data(1), i,j, 
     .    loc_parn(3,max_glb_sites), pcontrol, used(8)
 
*   origin_gain(max_glb_parn)   - Kalman gain vector for the
*               - translation origin.
*   translation(3)  - Translation of coordinate system to minimize
*                   - the horizontal displacements
 
      real*8 origin_gain(max_glb_parn), translation(3)

*   all_est  -  True if all site positions estimated (one for each coordinate)
      logical all_est(3)
 
 
      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
***** If needed compute postfit resiuduals
      if( compute_glb_res ) then
 
*         Remove parameters from the solvk solution which are not
*         not being used being used here.
          call remove_params( ema_data(icov_obs), ema_data(isol_obs),
     .                        ema_data(ipart_pnt), ema_data(ia_part) )
 
*         Now compute the residuals
          call glb_residuals( ema_data(icov_sav_parm),
     .                        ema_data(isol_sav_parm),
     .                        ema_data(ia_part), ema_data(ipart_pnt),
     .                        ema_data(irow_copy),
     .                        ema_data(icov_obs),ema_data(isol_obs),
     .                        ema_data(itemp_gain) )
 
*         Get the corrections to the aprioris to allow the "observed" values
*         to be calculated.  The observed local parmameters are those values
*         with translations and rotations removed.
 
          do i = 1, cnum_used
 
              call get_bak_apr( i, gpar_codes(i),
     .             ema_data(isol_sav_parm), ema_data(irow_copy) )
 
          end do
 
*         Now stabilize the results  (Site positions).  Fix the orgin first
*         to minimize the horizintal displacements
 
*                             ! Set to use all sites in orgin. Will be
          do i = 1,8
*                             ! Changed later.
              used(i) = -1
          end do
 
          do i = 1,3
*                                     ! Clear parn array first
              do j = 1, cnum_sites
                  loc_parn(i,j) = 0
              end do
 
              call codes_to_parn( i+6, gpar_codes, cnum_used,
     .             loc_parn(i,1),3)
          end do

*         scan to see if site coordinates estimated
          do i = 1,3
             all_est(i) = .false.
             do j = 1, cnum_sites
                if( loc_parn(i,j).eq.0 ) all_est(i) = .false.
             end do
          end do
 
C         We should not try to fix the orgin here.  It is already known
C         call fix_origin( 0.d0, apr_val_site, 2, 'EMA', translation,
C    .        ema_data(isol_obs), loc_parn, 1, used, cnum_sites,
C    .        ema_data(itemp_gain), ltog_sites)
 
 
          do i = 1,3
 
              translation(i) = 0.d0
              if( all_est(i) ) then 
                 call glb_stabilize( loc_parn(i,1),  cnum_sites, 3,
     .                               ema_data(icov_obs), cnum_parn,
     .                               cnum_used, ema_data(isol_obs),
     .                               origin_gain, used, translation(i))
              end if
          end do
 
*         Stabilize the RA orgin if needed
          do j = 1, cnum_sources
              loc_parn(1,j) = 0
          end do
 
          call codes_to_parn( 11, gpar_codes, cnum_used, loc_parn, 1)
          all_est(1) = .true.
          do j = 1, cnum_sources
             if( loc_parn(1,j).eq.0 ) all_est(1) = .false.
          end do

          if( cnum_sources.eq.0 ) all_est(1) = .false.

          if( all_est(1) ) then
             call glb_stabilize( loc_parn, cnum_sources, 1,
     .                           ema_data(icov_obs), cnum_parn,
     .                           cnum_used, ema_data(isol_obs),
     .                           origin_gain, used, translation)
          end if
 
*****     Now increment the post-fit chi**2 (Rigorous)
 
          call glb_inc_postfit(ema_data(icov_obs), ema_data(isol_obs),
     .                         ema_data(itemp_gain) )
          call inc_bak_stats ( ema_data(icov_obs), ema_data(isol_obs))
 
      end if
 
***** Thats all
      return
      end
 
CTITLE GLBAS
 
      subroutine glbas(iout) 

      implicit none
 
 
*     This segment of the GLBAK program
*     GLBAS -- Frame realization using the same methods as used in glorg
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'


* PASSED Variables
      integer*4 iout   ! Output unit number
 
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)

 
      integer*4 ema_data(1)

      integer*4 i, j, it   ! Loop counters
      integer*4 lnum,      ! Local site number for global (-1 if site not
                           ! used in this experiment).
     .          gtol_map   ! Function to map global site/svs number to local.
      integer*4 cnd_bits_copy   ! Copy of the [zyz]tran, [xyz]rot, scale bits
                           ! in case parameters turned off during iterations


*     pmu_parts(3,3,max_glb_sites)  - Polar motion UT1 partial derivatives

      real*8 pmu_parts(3,3,max_glb_sites)

      logical looks_ok
 
      common / globc_ema / ema_data

****  For this day, set up the options for apply_cond_full
       
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
* MOD TAH 190524: In GLBAK list of frame sites saved in cov_sites
*     to avoid overwriting the use_sites array
         use_pos(i) = cov_sites(i)
         use_rat(i) = cov_sites(i)
      end do
* MOD TAH 190711: Check that the sites are actually used in this
*     solution (needed if process noise is small and so that 
*     sigma on site is still resonable
      do i = 1, gnum_sites
         lnum = gtol_map( i, ltog_sites, gnum_sites )
         if( lnum.le.0 ) call sbit(use_pos,i,0)
         if( lnum.le.0 ) call sbit(use_rat,i,0)
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

****  Unlike glorg we will save the input covaraince matrix in 
*     memory so that it can be recovered when the systems is
*     iterated. (sol_vav_parm is next to cov_save_parm in 
*     memory so save one extra column).
*     (The 66*num_glb_parn allows for CA' matrix (cat) which 
*     is 7*num_glb_parn).
      call save_cov_sav(ema_data(icov_sav_parm), num_glb_parn,
     .     ema_data(isol_parm+2*num_glb_parn+66*num_glb_parn),'S') 

      do it = 1, num_stab_iter
          if( cnd_pos_bits.ne.0 ) then
* MOD TAH 190524: Make suer correct arrays are passed.  Final 
*     solution in glbak is aved in cov_sav_parm, sol_sav_parm.
*     These are arrays that the input solutions are read into.
* MOD TAH 200201: Use copy of parameters to be estimated incase
*     changed due to small number of sites.
              cnd_bits_copy = cnd_pos_bits  ! use copy in call
              call apply_cond_full( iout, parn_site(1,1,1), 2,
     .               pmu_parts, gnum_sites, ema_data(icov_sav_parm), 
     .               ema_data(isol_sav_parm), apr_val_site(1,1,1),
     .               cond_var(1,1),  num_glb_parn, 
     .               ema_data(isol_parm+2*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+14*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+28*num_glb_parn),
     .               use_pos, cnd_bits_copy, cnd_hgt_var(1),
     .               use_ratio(1), 'Position', gsite_names, 
     .               stab_site_var(1,1), stab_site_err,
     .               stab_nsig, stab_rel, stab_min_dh(1), 
     .               stab_min_rms(1), stab_min_dne(1), 
     .               it, org_opts, list_file)
          end if

*         Do the rates
          if( cnd_rat_bits.ne.0 ) then

              call apply_cond_full( iout, parn_site(1,2,1), 2, 
     .               pmu_parts, gnum_sites, ema_data(icov_sav_parm), 
     .               ema_data(isol_sav_parm), apr_val_site(1,1,1),
     .               cond_var(1,2),  num_glb_parn, 
     .               ema_data(isol_parm+2*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+14*num_glb_parn),
     .               ema_data(isol_parm+2*num_glb_parn+28*num_glb_parn),
     .               use_rat  , cnd_rat_bits, cnd_hgt_var(2),
     .               use_ratio(2), 'Velocity', gsite_names,
     .               stab_site_var(1,2), stab_site_err,
     .               stab_nsig, stab_rel, stab_min_dh(2), 
     .               stab_min_rms(2), stab_min_dne(2), 
     .               it, org_opts, list_file)
          end if

*         If we are not at the last iteration, re-read the loose 
*         covarinace matrix and re-update the parameter changes.
          if ( it.lt.num_stab_iter ) then
              call save_cov_sav(ema_data(icov_sav_parm), num_glb_parn,
     .          ema_data(isol_parm+2*num_glb_parn+66*num_glb_parn),'R') 
          end if
      end do
 
***** Now check that we covariance matrix is OK
 
      call check_covar(ema_data(icov_sav_parm), num_glb_parn, 
     .                 num_glb_parn)
   
***** Thats all
      return
      end
 
CTITLE SAVE_COV_SAVE

      subroutine save_cov_sav(cov_sav_parm, num_glb_parn,
     .     cov_copy,direct) 

      implicit none

*     Routine to copy cov_sav_parm plus sol_sav_parm to 
*     an array for later retrival, direct == S and for retrevial
*     direct == R.

      integer*4 num_glb_parn   ! number of global parameters
      real*8 cov_sav_parm(num_glb_parn,num_glb_parn+1),  ! Input+solution
     .       cov_copy(num_glb_parn,num_glb_parn+1)   ! saved copy

      character*(*) direct   ! 'S' to save, 'R' to recover

      integer*8 I8   ! Force vector size to I*8

      data I8 / 1 /

***** See which direction we need to go.
      if( direct(1:1).eq.'S' ) then
*        use VIS routine 
         call dwmov8(cov_sav_parm,1,cov_copy,1,
     .              num_glb_parn*(num_glb_parn+I8))
      else
*        Copy back
         call dwmov8(cov_copy,1,cov_sav_parm,1,
     .              num_glb_parn*(num_glb_parn+I8))
      endif

***** Thats all
      return
      end


      subroutine dump_cov_parm( String, cov_parm, num_glb_parn)

      implicit none

      character*(*) String
      integer*4 num_glb_parn
      real*8 cov_parm(num_glb_parn,num_glb_parn)

      print *,'DUMP ',string, cov_parm(1,1), cov_parm(1,2)

      end


