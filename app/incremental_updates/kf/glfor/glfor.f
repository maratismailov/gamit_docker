      subroutine glfor(ms_type)

      implicit none
 
 
*     This the forward global Kalman filter program.  It is made
*     up of three segments plus the main:
*     GLFOO -- Opens/reads/closes files
*     GLFOP -- Computes O-C and partial derivatives
*     GLFOK -- Actually does the Kalman filtering.
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED Variables

* ms_type   - Character string which indicates that glfor is a program
*          (MAIN) or a subroutine (SUBR)

      character*(4) ms_type

* LOCAL variables.
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Segment control used by GLFOO.
*               - pcontrol = 1, open the files we need
*               - pcontrol = 2, open new global solution
*               - pcontrol = 3, save the current solution
*               -               covariance matrix
*               - pcontrol = 4, write control block of GLOBK common
*               -               and covariance matrix of no back
*               -               solution
*   trimlen     - Length of string
 
      integer*4 cov_dcb(16), ema_data(1), i, j, ierr,
     .    pcontrol, trimlen

*   chi_too_big - Indicates that the chi**2 is too large. Run
*     will abort.  There will be problems in GLOUT if a back
*     solution is run.  But other than that stop should be clean
*   looks_OK    - Scan of prefit residuals, indcate that they are
*     resonabley small

      logical chi_too_big, looks_OK
 
*   timr        - LIBHS timer function
*   running_time  - Running time for filter
 
      real*4 timr, running_time

*   ephem_prev  - Previous satellite empheris JD 
*   secs        - Used in call to systime for current time
*   glb_var     - Variance multiplier for this global.  Used 
*                 make chi**2 closer to 1 and more uniform
*   glb_diag    - Diagonal scaling 

      real*8 ephem_prev, secs, glb_var, glb_diag 

*   apply_sd   - Indicates that we should apply SD corrections
*   apply_nut  - Indicates that we should apply nutation corrections
*   kbit       - Check bit set.

      logical apply_sd, apply_nut, kbit

* glb_used - Logical to denote if global file used
 
      logical glb_used 
 
*   cname       - Name of next solution equivalanced to iname

      character*(sort_recl) cname
 
*   cr          - Carriage return
 
      character*1 cr

*   outline     - Line to be output to the status report
*   progname    - Name of program

      character*256 outline, progname

*   lenprog     - Length of program name
*   rcpar       - Returns runstring

      integer*4 lenprog, rcpar, FmpClose
      
* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file.
*   for_chi, bak_chi -- Forward and backwarsd chi**2/f
      real*4 for_chi, bak_chi
 
* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data

*                        approximately 1900.
      data  ephem_prev / 2415020.d0 /
      data I8 / 1 /
 
***** Get the runstring, and open the GLOBK common, and the sort file

      cr = char ( 13 )
      running_time = timr(0)
      lenprog = rcpar(0,progname)

      pcontrol = 1

      call GLFOO (ms_type)
 
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

*     Set coseismic to show not applied
      do i = 1, num_eq
         eq_co_applied(i) = .false.
      end do

*     Get the apriori values of any log terms in the extended model
*     (needed specifically for log terms since these can be estimated
*     The other extended model components can not be estimated and 
*     therefore we do not need to save.  The apr_val_log terms are used
*     when log estimates are found in the input hfiles).
      do i = 1, gnum_sites
         call get_nonlog(i,apr_val_log(1,i))
      end do

* MOD TAH 030116: If user has asked, report the parameters to be estimated
      if( kbit(prt_opts,22) .or. kbit(org_opts,22) .or.
     .    kbit(crt_opts,22) ) then
          call rep_plst(log_unit)
          if( log_unit.ne.6 ) call rep_plst(6)
      end if

* MOD TAH 090930: Output the rename list to the log file or unit 6.
      if( kbit(prt_opts,30) .or. kbit(org_opts,30) .or.
     .    kbit(crt_opts,30) ) then
          call rep_renm(log_unit)
          if( log_unit.ne.6 ) call rep_renm(6)
      end if

***** set ptd_updated false (set true is changes made)
      ptd_updated = .false.

***** Now loop over the experiments in the global solution, adding
*     each one to the solution.
      i = 0 
      chi_too_big = .false.
      do while ( i.lt.num_glb_sol .and. .not.chi_too_big )
          i = i + 1
 
          read(100,iostat=ierr,rec=i) cname, glb_var,
     .                                for_chi, bak_chi
          
          call report_error('IOSTAT',ierr,'read',sort_file,0,
     .                      'GLFOR')

*         MOD TAH 970430: Break up the glb_var if it is < 0.
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
 
*                                     ! Open the next global solution
              pcontrol = 2
              call GLFOO(ms_type)
 
*****         Get epoch change between experiments (gepoch_epxt
              gepoch_prev = gepoch_expt
              gepoch_expt = cepoch_expt

* MOD TAH 961122: Check the start and end epoch based on the start
*             of the first experiment, and end of last.
* MOD TAH 970624: Ignore the direction and check for the eariler and
*             latest times. 
              if ( i.eq.1 ) then
                  gepoch_start = cepoch_start
                  gepoch_end = cepoch_end
                  call sbit(mul_pmu_opt,31,0)
              end if
              gepoch_start = min(gepoch_start,cepoch_start)
              gepoch_end = max(gepoch_end,cepoch_end)
*                                                          ! Years
              deltat = (gepoch_expt-gepoch_prev)/365.25d0
              if( i.eq.1 ) deltat = 0
*             We need a small deltat so that the earthquake modules
*             know which direction time is going.
              if( sort_direction.eq.-1.and.deltat.eq.0) deltat = -1.d-12 

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
                      if( log_unit.ne.6 ) then
                         write(log_unit,90) deltaephem
                      end if
                  end if
              end if
 
              running_time = timr(-1)
 
              write(*,100) i,(isol_obs+2*cnum_parn-icov_parm+128)*
     .                    4.0/1024./1024., running_time, glb_var, 
     .                     glb_diag
  100         format(' Global ',i4,' using ',f7.1,' Mb. Running time ',
     .                f8.2,' Scaling by ',F10.3,F11.8)
              if( log_unit.ne.6 ) then
                  write(log_unit,100) i,
     .                       (isol_obs+2*cnum_parn-icov_parm)*
     .                       4.0/1024./1024., running_time, glb_var,
     .                       glb_diag 
              end if
              write(outline,110) i,(isol_obs+2*cnum_parn-icov_parm+128)*
     .                    4.0/1024./1024., glb_var, glb_diag
  110         format('File ',i4,' using ',f8.2,
     .               ' Mb. Scaled by ',F10.3,F11.8)
              call report_stat('status',progname,'glfor',glb_inp_file,
     .                         outline,0)

*****         Now scale the covariance matrix

              if( glb_var.ne.1.0d0 ) then
                  call dwsmy8(glb_var, ema_data(icov_obs), 1, 
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
* MOD TAH 070926: Added explicit satellite epoch
              call glb_upd_apr( cepoch_expt, cepoch_expt, .true., 
     .                          ema_data(isol_obs), .true.)

*             Get the markov statistics for time dependent statellite
*             markov elements.
              call get_mar_svs( cepoch_expt )

*             Now compute nutation series corrections. 
              if( apply_nut ) then
                  call compute_std( cepoch_expt, vnut_ang_apr )
              else
                  vnut_ang_apr(1) = cnut_ang_apr(1)
                  vnut_ang_apr(2) = cnut_ang_apr(2)
              end if

              call get_sd_mid( cepoch_expt, apply_sd, gnum_sites)
 
*****         Now get the O-C and partials
              call GLFOP( looks_ok ) 
 
*****         Now actually do the Kalman filtering

              if( looks_ok ) then
 
                  call GLFOK(i, glb_used )

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
                  for_chi = dchi_save
                  
                  write(100,iostat=ierr,rec=i) cname, glb_var,
     .                                for_chi, bak_chi
                  
                  call report_error('IOSTAT',ierr,'writ','SRT file'
     .                                 ,0,'glfor') 

*****             If needed save the current solution. Only if the solutions
*                 is actually used. Write out the solution file even if not used
*                 since this is needed to for the bookeeping on size of the file. 
 
                  if( glb_bak_soln ) then
                      pcontrol = 3
                      call GLFOO (ms_type)
                  end if

*****             If estimated write out the pole positoin and
*                 UT1 estimates.

                  call print_pmu( 6 ,
     .                 ema_data(icov_parm), ema_data(isol_parm),
     .                 num_glb_parn)
                  if( log_unit.ne.6 ) call print_pmu(log_unit,
     .                 ema_data(icov_parm), ema_data(isol_parm),
     .                 num_glb_parn)
              end if
*                         ! No error getting next name
          end if
*                         ! Looping over the eperiments
      end do

***** Save the run time for this solution.

      call systime( grun_time, secs)
      grun_time(6) = secs
      grun_time(7) = 0

*     Update the gamit models if we have applied pole tide
* MDO TAH 200220: Only change flags if actually updated
      if( ptd_updated ) then 
         if( kbit(ptide_opt,1) ) then
            call sbit(ggamit_mod,19,1)
*           Now set bits based on mean pole
            if( mean_pole_def(1:6).eq.'IERS20' ) then
               call sbit(ggamit_mod,26,1)
               call sbit(ggamit_mod,23,0)
               call sbit(ggamit_mod,21,0)
            elseif( mean_pole_def(1:6).eq.'IERS10' ) then
               call sbit(ggamit_mod,26,0)
               call sbit(ggamit_mod,23,1)
               call sbit(ggamit_mod,21,0)
            elseif( mean_pole_def(1:6).eq.'IERS96' ) then
               call sbit(ggamit_mod,26,0)
               call sbit(ggamit_mod,23,0)
               call sbit(ggamit_mod,21,1)
            endif
         endif
*        See if ocean pole tide applied
         if( kbit(ptide_opt,2) ) then
*           Now set bits based on mean pole
            if( mean_pole_def(1:6).eq.'IERS20' ) then
               call sbit(ggamit_mod,27,1)
               call sbit(ggamit_mod,25,0)
               call sbit(ggamit_mod,24,0)
            elseif( mean_pole_def(1:6).eq.'IERS10' ) then
               call sbit(ggamit_mod,27,0)
               call sbit(ggamit_mod,25,1)
               call sbit(ggamit_mod,24,0)
            elseif( mean_pole_def(1:6).eq.'IERS96' ) then
               call sbit(ggamit_mod,27,0)
               call sbit(ggamit_mod,25,0)
               call sbit(ggamit_mod,24,1)
            endif
         endif
*****    Now aee if models removed. 
*        If SE pole tide removed, set bit 19 to 0 to show
*        not applied
         if( kbit(ptide_opt,3) ) then
            call sbit(ggamit_mod,19,0)
         endif

*        If set the ocean pole tide bits to 0 
         if( kbit(ptide_opt,4) ) then
             call sbit(ggamit_mod,27,0)
             call sbit(ggamit_mod,25,0)
             call sbit(ggamit_mod,24,0)
         endif
         write(*,200) trim(mean_pole_def), ptide_opt, ggamit_mod
 200     format('Updated pole-tide: Mean pole ',a,' PTIDE OPT ',o6,
     .          ' GAMIT_MOD o',o11.11) 
      endif


***** Save the ephemeris time for the last solution.
      do i = 1, cnum_svs
         svs_epoch(i) = csvs_epoch
      end do

***** Write control and if we did not save the solution at all save now
      pcontrol = 4
      call GLFOO (ms_type)

      if( log_unit.ne.6 ) close(log_unit)
      close(100)
       
* MOD TAH 961028: Close the sol_file dcb 
      ierr = FmpClose( cov_dcb )

      write(*,'(1x)')
 
***** Thats all
      end
 
      subroutine glfoo(ms_type)

      implicit none
 
 
*     This segment of the forward global Kalman filter program will:
*     GLFOO -- Opens/reads/closes files
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
* PASSED Variables

* ms_type   - Character string which indicates that glfor is a program
*          (MAIN) or a subroutine (SUBR)

      character*(4) ms_type

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
*               - pcontrol = 2, open new global solution
*               - pcontrol = 3, save the current solution
*               -               covariance matrix
*               - pcontrol = 4, write control block of GLOBK common
*               -               and covariance matrix if no back
*               -               solution
*   vma_bytes   - Number of bytes needed.
* MOD TAH 150524: Introduced integer*8 version for large 
*     solutions >20000 parameters
 
      integer*4 cov_dcb(16), ema_data(1), ierr,
     .    len_read, pcontrol, trimlen, fmpclose
      integer*8 vma_bytes

      logical load_upd  ! Function that returns true if loads are be 
                        ! updated.

* MOD TAH 190603: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
 
      common / progcon / pcontrol, cov_dcb
 
 
 
      common / globc_ema / ema_data
 
***** Depending on pcontrol do the desired operation
 
*                               ! Get the runstring, and
      if( pcontrol.eq.1 ) then
*                               ! open the GLOBK common, and the sort file

****      We only need do the following if this is a main program
          if( ms_type(1:4).eq.'MAIN' ) then  

***           Decode the runstring
              call decode_glfor_run
 
***           Open the globk common, and read the segments
              call open_globk( glb_com_file, ierr)
              call report_error('FmpOpen',ierr,'open',glb_com_file,1,
     .                      'Segment GLFOO1')
 
              call rw_globk_block('R','CONTROL',glb_control,ierr)
              call report_error('FmpRead',ierr,'read','control block',1,
     .                      'Segment GLFOO2')
              istart_vma = 0
              call rw_globk_block('R','MARKOV', glb_markov ,ierr)
              call report_error('FmpRead',ierr,'read','markov block',1,
     .                      'Segment GLFOO3')
 
              call rw_ema_block('R','EMA',glb_ema,ierr)
              call report_error('VREAD',ierr,'read','ema block',1,
     .                      'Segment GLFOO4')
          end if
 
 
***       Now open the sort file
* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file. Increaed reclength by 8 bytes
          open(100, file=sort_file, iostat=ierr, status='old',
     .         access='direct', recl=sort_recl+8+8)
          call report_error('IOSTAT',ierr,'open',sort_file,1,
     .                      'Segment GLFOO5')
 
***       Open the output file for the temporary covariance matrices
          call rw_glb_covar('C', cov_dcb, ema_data)
 
***       Initalize the covariance matrix for the solution
*****     Assign the memory needed for this run
          call get_max_glb_deriv
          call glfor_mem( ema_data )
          vma_bytes = 0
          write(*,*) 'Initial memory mapped'  
          call glb_for_map( cnum_parn, vma_bytes ) 
          call init_glb_cov ( ema_data(icov_parm), ema_data(isol_parm))
*         Initialize the initial multi-day pmu
          gmul_pmu_ep(1) = 0
      end if
 
 
***** See if global solution file operation.
      if( pcontrol.eq.2 ) then
 
*****     Open the current global solution file, and read the header, and
*         names
          call FmpOpen(cglb_dcb,ierr,glb_inp_file,'RO',1)
          call report_error('FmpOpen',ierr,'open',glb_inp_file,0,
     .                      'Segment GLFOO6')
 
***       Read the header
          call rw_glb_header('R',ierr)
          call report_error('FmpRead',ierr,'read',glb_inp_file,0,
     .                      'Segment GLFOO7')
 
*         MOD TAH 150825: Get the PRN to SVN information
          call read_svinf_rec
 
****      Get the names and local to global conversions
          call rw_names_block('R')

*         Apply any name changes that are needed
          call eq_name_change('NO') 
 
****      Get the parameter codes for this solution
 
          call get_param_codes( gpar_codes ) 

*****     See if we need to set up a new group of multi-day polar
*         motion/UT1 values

          call setup_mul_pmu


*****     Save wobble, UT1 and nutation angles
          call save_glb_apr

****      Now get the 
*                         ! See if nutation/PMU estimated
          call scan_parm

* MOD TAH 130418: See if we need to read in the load values from the 
*         input binary file 
          if( load_upd() ) then
              call read_sinf_recs
          end if

 
*****     Map the ema area usage and read the solution and covariance
*         matrix into ema
 
C          call get_max_glb_deriv
 
          call glb_for_map( cnum_parn, vma_bytes )
 
* MOD TAH 190603: Change to writd8 to allow for >32767 parameters. 
          if( cnum_parn.gt.32767 ) then 
             call readd8(cglb_dcb,ierr,ema_data(icov_obs),
     .              (I8*128)*cnum_par_vals, len_read, crec_par_vals)
          else
             call readd(cglb_dcb,ierr,ema_data(icov_obs),
     .                  128*cnum_par_vals, len_read, crec_par_vals)
          endif

C          call dump_covobs(ema_data(icov_obs), ema_data(isol_obs), 
C     .         cnum_parn)
 
          call report_error('VREAD',ierr,'read','global solution',0,
     .                      'Segment GLFOO8')
          ierr = fmpclose(cglb_dcb)
 
 
      end if
 
***** See if saving current solution
 
      if( pcontrol.eq.3 ) then
 
*         Save the current solution
          call rw_glb_covar('W', cov_dcb, ema_data(icov_parm))
 
      end if
 
***** Save the common and possible the covariance matrix
 
      if( pcontrol.eq.4 ) then
          if( ms_type(1:4).eq.'MAIN' ) then 
              call rw_globk_block('W','CONTROL',glb_control,ierr)
              call rw_globk_block('W','MARKOV' ,glb_markov ,ierr)
              call rw_ema_block('W','EMA',glb_ema,ierr)
          end if 

*         If this not MAIN, see we should the common
          if( ms_type(1:4).ne.'MAIN' .and.
     .        trimlen(glb_com_file).gt.0 .and.
     .        .not.glb_bak_soln ) then
             call rw_globk_common('W')
          end if
   
          if( .not.glb_bak_soln ) then
              call rw_glb_covar('W', cov_dcb, ema_data(icov_parm) )
          end if
          ierr = fmpclose(cov_dcb)
      end if
 
***** Thats all
      return
      end
 
CTITLE GLFOP
 
      subroutine glfop( looks_ok )

      implicit none
 
 
*     This segment of the forward global Kalman filter program:
*     GLFOP -- Computes O-C and partial derivatives
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
 
      integer*4 cov_dcb(16), ema_data(1), pcontrol

*   Looks_OK    - Indicates that prefit residuals look fine

      logical Looks_ok
 
 
      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
***** Compute the O-C values (Note the ema usage was mapped in GLFOO
*     when this global solution was first read)
      call glb_o_minus_c( ema_data(isol_obs), looks_ok )
 
***** Now compute the partial derivatives which we need
 
      call glb_partials( ema_data(ipart_pnt), ema_data(ia_part) )

***** Now check the sol_obs array and see if we should suppress
*     bad elements

      call check_param( ema_data(isol_obs) )
 
***** Thats all
      return
      end

CTITLE CHECK_PARAM 

      subroutine check_param( sol_obs )

      implicit none

*     Check o-minus-c size and remove parameter if too large

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

      real*8 sol_obs(cnum_parn)

* Local variables
*  type, indx -- Type and index foe parameter
*  orb_el, sv_num -- Orbital elements and svs number 

      integer*4 i, type, indx , orb_el, sv_num, save_indx, lpsec 

*  tol  -- User specfied max_prefit_diff converted to correct units.
      real*4 tol

*  eop_types(6) -- Descriptive term for EOP

      character*8 eop_types(6) 

*  converged -- Logical to indicate that all bad site coordinates have
*     been removed

      logical converged

      data eop_types / 'X-Pole','Y-Pole','X-Rate','Y-Rate',
     .                 'UT1-AT','LOD' / 

****  Pre-scan the site coordinates to make sure that after rotation
*     and translation that they are OK

      converged = .false.
      do while ( .not.converged ) 
          call scan_site( sol_obs, converged )
      end do


*     Loop over the parameters and check the size of o-minus-c
      do i = 1, cnum_parn

         call decode_code( gpar_codes(i), type, indx )

         call get_tol_pref( type, indx, max_prefit_diff, tol ) 

*        Do check sites here (type 7-9) since they are checked
*        above.  Check only those parameters that we will actually
*        use (indicated by indx_pnt(i) being greater than zero).
         if( abs(sol_obs(i)).gt.tol .and. type.gt.9 .and.
     .       indx_pnt(i).gt.0 ) then

*            Setting this to zero will cause the parameter to be
*            removed from solutionn.
             save_indx = indx_pnt(i)
             indx_pnt(i) = 0

*            Tell user that is is bad
             if( type.eq.13 .and. indx.le.4 ) then
                write(*,100) eop_types(indx), sol_obs(i)
                if( log_unit.ne.6 )
     .          write(log_unit,100) eop_types(indx), sol_obs(i) 
 100            format(' BAD PREFIT EOP parameter ',a8, 
     .                 ' Difference from apriori ',F18.5)
             else if ( type.eq.14 .and. indx.le.2 ) then
                write(*,100) eop_types(indx+4), sol_obs(i)/15.d0
                if( log_unit.ne.6 )
     .          write(log_unit,100) eop_types(indx+4), sol_obs(i) 
          
             else if( type.eq. 58 ) then

*               See if this is UT1, in which case we will try to sort
*               out leap second differences
                lpsec = nint(sol_obs(i)/15000.d0)
                if( abs(sol_obs(i)-lpsec*15000.d0).lt.tol ) then
*                   Value OK, just warn user
                    write(*,110) i, lpsec
 110                format('**WARNING** Adjusting leap-seconds for',
     .                     ' parameter ',i5,' by ',i5,' seconds')
                    sol_obs(i) = sol_obs(i) - lpsec*15000.d0
                    indx_pnt(i) = save_indx  
                else
                    write(*,115) i, sol_obs(i)/15.d0
 115                format(' BAD PREFIT for UT1 parameter ',i5,
     .                     ' Error ',f12.4,' ms')
                end if
 
             else if( type.ne. 51 ) then
                write(*,120) i, cnum_parn, type, indx, sol_obs(i)
                if( log_unit.ne.6 )
     .          write(log_unit,120) i, cnum_parn, type, indx, sol_obs(i)
 120            format(' BAD PREFIT for parameter ',i4,'/',i4,
     .                 ' Type/index ',2i5,' Value ',d12.5)
             else 
 
*               Get the orbital element and satellite number
                call decode_code( indx, orb_el, sv_num )
 
*               Check to see if being estimated
                if( parn_svs(orb_el, ltog_svs(sv_num)).ne.0 ) then
*                   Orbit being estimated.  We have a problem:
                    if( log_unit.ne.6 )
     .              write(log_unit,130) orb_el, sv_num,
     .                           gsvs_names(ltog_svs(sv_num)),
     .                           sol_obs(i) 
                    write(*,130) orb_el, sv_num,
     .                           gsvs_names(ltog_svs(sv_num)),
     .                           sol_obs(i)
 130                format(' BAD PREFIT for',
     .                     ' orbital element #',i2,' for local',
     .                     ' SV ',i3,1x,a,' Difference ',f22.6)
                end if
             end if

          end if
      end do

****  Thats all
      return
      end 

      subroutine debug_ut1( cov_parm, np, mes)

      implicit none

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

      real*8 cov_parm(num_glb_parn,num_glb_parn)
      integer*4 np
      character*(*) mes
      write(*,100) mes, np, cov_parm(np,np)
 100  format(a,1x,i4,1x,F12.3)
      return
      end

CTITLE SCAN_SITE

      subroutine scan_site( sol_obs, converged )

      implicit none

*     Routine to scan the site o-c vector, estimate a rotation and
*     translate and remove any "bad" sites

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED Variables

* sol_obs(cnum_parn)  -  Solution vector
* converged  - Logical set true when no stations deleted

      real*8 sol_obs(cnum_parn)
      logical converged

* LOCAL VARIABLES

*  i,j,k   - Loop counters
*  type, indx  -- Type of parameter and its index (see glb_hdr_def.h)
*  tt, it      -- Test version of type and indeX
*  ipivot(6)   -- Pivot elements for inversion 
*  ngood       -- Number of good sites left.
*  ned         -- Number of editied sites.
*  mesite      -- Maximum error site ,
*  meparn      -- Parameter number of maxium error site 

      integer*4 i,j,k, type, indx, ngood, ned,  ipivot(6), tt, it ,
     .          mesite, meparn 

*  eop_der(6)  - Derivatives for EOP and translation
*  neq(6,6)    - Normal eqautions for estimation
*  seq(6)      - Solution vector
*  scale(6)    - Scaling vector for inversion
*  res         - Site coordinate residual.
*  err(3)      - XYZ coordinate errors 
*  max_err     - Maximum error
*  eop_diff(3) - Differnece between eop estimate from Pole position
*                and from site corrridates.
*  tolvar      - Tolerance for variance of rotational variances (to
*                catch short baseline large sigmas)
      real*8 eop_der(6), neq(6,6), seq(6), scale(6), res, err(3) ,
     .       max_err , eop_diff(3) , tolvar

****  Clear the normal eqautions and soluion
      do i = 1, 6
         seq(i) = 0.d0
         do j = 1, 6
            neq(i,j) = 0.d0
         end do
      end do
      do i = 1,3
         eop_diff(i) = 0.d0
      end do 

****  Now loop over parameters getting the normal eqautions
      ngood = 0 
      do i = 1, cnum_parn
         call decode_code( gpar_codes(i), type, indx )
         if( type.ge.7 .and. type.le.9 .and. indx_pnt(i).gt.0 ) then
             ngood = ngood + 1 
             do j = 1,6
                eop_der(j) = 0.d0
             end do 
             call pmu_part(indx, eop_der, 
     .            apr_val_site(1,1,ltog_sites(indx)), type-6, cut1_apr)
             eop_der(type-6+3) = 1.d0

*** DEBUG:   
             if( glb_com_file(1:8).eq.'test.com' ) then
                  write(*,*) 'DEBUG: EOP Estimates ',
     .                 gsite_names(ltog_sites(indx)), ' Cooord ',indx,
     .                 ' Error ', sol_obs(i),' UT1 ', cut1_apr(1:2), i
             end if

*            Accumulate the normal equations
             do j = 1,6
                seq(j) = seq(j) + sol_obs(i)*eop_der(j)
                do k = 1,6
                   neq(j,k) = neq(j,k) + eop_der(j)*eop_der(k)
                end do
             end do
         end if
      end do

      if( ngood.gt.6 ) then

*         Enough data for solution
          call invert_vis( neq, seq, scale, ipivot, 6,6,1)

      else
          do j = 1,6
             seq(j) = 0.d0
          end do
      end if 

* MOD TAH 140106: Added test of sigma (for short baseline results)
      tolvar = max_eop_rot**2/1.d-5 ! Data variance of 3mm^2
      if ( neq(1,1)+neq(2,2)+neq(3,3).gt.tolvar ) then
          write(*,180) max_eop_rot, (neq(1,1)+neq(2,2)+neq(3,3))*1.d-5
          if( log_unit.ne.6 ) write(log_unit,200)  max_eop_rot, 
     .                              (neq(1,1)+neq(2,2)+neq(3,3))*1.d-5
 180      format('* Rotation variance large: Tolerance ',F10.2,
     .           ' mas, Sum of EOP varinances ',E15.4,' mas^2')
          converged = .true.
          RETURN
      end if

*     Now check the quality of the position estimate
      ned = 0
      max_err = 0
      if( abs(seq(1)).gt.max_eop_rot .or.  
     .    abs(seq(2)).gt.max_eop_rot .or.
     .    abs(seq(3)).gt.max_eop_rot ) then 
          write(*,200) max_eop_rot, seq
          if( log_unit.ne.6 ) write(log_unit,200) max_eop_rot, seq
 200      format(' Rotation TOO Large removing. Tolerance ',F10.2,
     .           ' mas',/,
     .           ' Rots (XYZ, mas) ',3f10.2,' Trans (m) ',3f10.6)
      end if 

      do i = 1, cnum_parn
         call decode_code( gpar_codes(i), type, indx )
         if( type.ge.7 .and. type.le.9 .and. indx_pnt(i).gt.0 ) then
             do j = 1,6
                eop_der(j) = 0.d0
              end do 
              call pmu_part(indx, eop_der, 
     .             apr_val_site(1,1,ltog_sites(indx)), type-6, cut1_apr)
              eop_der(type-6+3) = 1.d0

*             Now check residual
              res = sol_obs(i)
              do j = 1,6
                 res = res - seq(j)*eop_der(j)
              end do

*             If rotation too large remove
              if( abs(seq(1)).gt.max_eop_rot .or.  
     .            abs(seq(2)).gt.max_eop_rot .or.
     .            abs(seq(3)).gt.max_eop_rot ) then 
                   do j = 1,3
                      sol_obs(i) = sol_obs(i) - seq(j)*eop_der(j)
                   end do 
              end if 

*             See if passes tolerance
              k = 0
              if( abs(res).gt. max_prefit_diff ) then
                  ned = ned + 1 
                  if( abs(res).gt.max_err ) then
                      mesite = indx
                      meparn = i
                      max_err = abs(res)
                  end if
              end if
          end if

*         See if this is rotation paramter 
          if( type.eq.13 .or. type.eq.14 ) then 

*             If rotation too large remove
              if( abs(seq(1)).gt.max_eop_rot .or.  
     .            abs(seq(2)).gt.max_eop_rot .or.
     .            abs(seq(3)).gt.max_eop_rot ) then 
                   if( type.eq.13 .and.indx.eq.1 ) then
                       eop_diff(1) = sol_obs(i) 
                       sol_obs(i) = sol_obs(i) + seq(1)
                   else if( type.eq.13 .and.indx.eq.2 ) then
                       eop_diff(2) = sol_obs(i) 
                       sol_obs(i) = sol_obs(i) + seq(2)
                   else if( type.eq.14 .and. indx.eq.1 ) then
                       eop_diff(3) = sol_obs(i) 
                       sol_obs(i) = sol_obs(i) + seq(3)
                   end if 
              end if
           end if 
      End do

*     See if we have bad sites
      if( ned .gt.0 ) then
          i = meparn
          call decode_code( gpar_codes(i), type, indx )
          do j = i-(type-7), i-(type-9) 
              k = k + 1
              call decode_code(gpar_codes(i), tt, it )
              if( it.eq. indx 
     .           .and. tt.ge.7 .and. tt.le.9 ) then 
                  indx_pnt(j) = 0
                  err(k) = sol_obs(j)
              end if
          end do
          write(*,120) gsite_names(ltog_sites(indx)), err 
          if( log_unit.ne.6 )
     .    write(log_unit,120) gsite_names(ltog_sites(indx)), err 
  120     format(' BAD PREFIT coordinates for site ',a8,
     .           ' Diffs from apriori ',3F15.3,' m')
      end if

****  See if we edited anything
      if( ngood.lt.6 .or. ned.eq.0 ) then
          converged = .true.
      else
          converged = .false.
      end if

****  Report any large roations removed
      if( eop_diff(1).ne.0 .or. eop_diff(2).ne.0 .or. 
     .    eop_diff(3).ne.0 ) then 
          write(*,140) (eop_diff(i), seq(i), i = 1,3)
          if( log_unit.ne.6 )
     .    write(log_unit,140) (eop_diff(i), seq(i), i = 1,3)
 140      format('Large Rotation removed: Input EOP estimates ',
     .           3f12.3,' (Xp,Yp,UT mas) '/,
     .           '                        dPosition estimates ',
     .           3F12.3,' (Xp,Yp,UT mas) ')
      end if 

****  Thats all
      return
      end

CTITLE GLFOK
 
      subroutine glfok(ns, glb_used)

      implicit none
 
 
*     This segment of the forward global Kalman filter program:
*     GLFOK -- Does the Kalman filtering
*
*                                             07:56 PM SAT., 8 Aug., 1987
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

*   ns - Solution number.  If this is one we may do a direct
*        copy

      integer*4 ns

* glb_used - Logical to denote if global file used
 
      logical glb_used 
 
*   cov_dcb(16) - DCB buffer for the temporary storage of
*               - covariance matrces.
*   ema_data(max_vma_space) - Ema data area (used to hold all
*               - of the KF matrices.  Array is dynamically
*               - mapped during run)
*   i           - Loop counter (for looping over experiments)
*   ierr        - Fmpread error
*   pcontrol    - Controls segment GLFOO (not used here)
 
      integer*4 cov_dcb(16), ema_data(1), pcontrol

*   dir_copy   - Logical to indicate a direct of the initial
*                covariance matrix

      logical dir_copy 
 
 
      common / progcon / pcontrol, cov_dcb
 
 
      common / globc_ema / ema_data
 
***** Call the routines to add in the next solution into the
*     Kalman filer solution
*     First compress the parameters from the current SOLVK solution
*     which were not used in the solution
 
      call remove_params ( ema_data(icov_obs),   ema_data(isol_obs),
     .                     ema_data(ipart_pnt),  ema_data(ia_part)  )

* MOD TAH 961122: See if we can do a direct-copy of the covariance matrix
* MOD TAH 961220: Check if user said not to allow direct copying.
      if( ns.eq.1 .and. .not.no_direct_copy ) then
          call check_dc( ema_data(icov_parm), dir_copy)
      else
          dir_copy = .false.
      endif

      if( .not.dir_copy ) then 
          call glb_for_filter( ema_data(icov_parm),
     .                         ema_data(isol_parm),
     .                         ema_data(irow_copy),
     .                         ema_data(ipart_pnt),  ema_data(ia_part),
     .                         ema_data(itemp_gain), ema_data(ikgain),
     .                         ema_data(icov_obs),   ema_data(isol_obs),
     .                         glb_used)
      else
         write(*,300)
         if( log_unit.ne.6 ) write(log_unit,300)
 300     format(' Direct copy of initial covariance matrix invoked')
         call copy_dir( ema_data(icov_parm), ema_data(isol_parm),
     .                  ema_data(irow_copy),
     .                  ema_data(ipart_pnt),  ema_data(ia_part),
     .                  ema_data(icov_obs),   ema_data(isol_obs))
         glb_used = .true. 
      end if

 
***** Thats all
      return
      end
 
CTITLE REP_PLST

      subroutine rep_plst(iout)

      implicit none

*     Routine to report the estimated parameters in a solution
*

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES
      integer*4 iout  ! Output unit number

* LOCAL VARIABLES
      integer*4 i,j,k  
      integer*4 np, trimlen
      character*80 par_line

      character*4 pt(3,3)
      data pt / 'XPOS','YPOS','ZPOS','XDOT','YDOT','ZDOT',
     .          'NLOG','ELOG','ULOG' /

***** Output the stations parameters
      write(iout,110) num_glb_parn
 110  format('LIST OF ',i6,' PARAMETERS TO BE ESTIMATED',/,
     .       'Site         Par#  Parameter list')
      do i = 1, gnum_sites
         par_line = ' '
         np = 0
         do j = 1,2
            do k = 1,3
               if( parn_site(k,j,i).ne. 0 ) then
                   if( np.eq.0 ) np = parn_site(k,j,i)
                   par_line((k-1)*5+(j-1)*15+1:) = pt(k,j)
               end if
            end do
         end do
         do k = 1,3
            if( parn_log(k,i).ne.0 ) then
                par_line((k-1)*5+32:) = pt(k,3)
            end if
         end do
         if( np.gt.0 ) then
            write(iout,150) gsite_names(i),np, 
     .                       par_line(1:trimlen(par_line))
 150        format(A8,i6,2x, a)
         endif
      end do

****  Thats all for the momment
      return
      end

      subroutine dump_covobs(cov_obs, sol_obs, cnum_parn)

      implicit none

      integer*4 cnum_parn, i
      real*8 cov_obs(cnum_parn, cnum_parn), sol_obs(cnum_parn)
      write(*,100) cnum_parn
 100  format('Dumping diagonal for ',i5,' parameters')
      do i = 1, cnum_parn
         write(*,120) i, sol_obs(i), cov_obs(i,i)
 120     format(I4,1x,2d20.8)
      end do
      return
      end


CTITLE REP_RENM

      subroutine rep_renm(iout)

      implicit none

*     Routine to report the renames generated all the input eq-file
*     and specifically the names generated with the BREAK command.
*

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES
      integer*4 iout  ! Output unit number

* LOCAL VARIABLES
      integer*4 i,j , datee(5), dates(5) 
      real*8 sec
      character*8 oname   ! Original name with _GPS removed 
                          ! if it appears that way.

***** Output the stations parameters
      write(iout,110) num_renames
 110  format('RN* LIST OF ',i6,' RENAMES to be used',/,
     .       'RN*        SiteOld  SiteNew   [HFcode]        ',
     .       'Start  Date       End Date                 ',
     .       'dPOS (m)             Type     #')
      do i = 1, num_renames
         call mjd_to_ymdhms(rn_times(1,i),dates,sec)
         call mjd_to_ymdhms(rn_times(2,i),datee,sec)
         oname = rn_codes(1,i)
         if( oname(5:8).eq.'_GPS' ) oname(5:8) = ' ' 
         if( rn_dpos(1,i).eq.0.d0 .and. rn_dpos(2,i).eq.0.d0 .and.
     .       rn_dpos(3,i).eq.0.d0 ) then 
            write(iout,210) oname,rn_codes(2,i), rn_hfiles(i),
     .              dates, datee, i
 210        format('RN  RENAME ',a8,1x,a8,1x,a16,1x,I4,4(I3.2),
     .              2x,I4,4(I3.2),35x,'  ! # ',i6)
         else
            write(iout,230) oname,rn_codes(2,i), rn_hfiles(i),
     .              dates, datee, (rn_dpos(j,i),j=1,3),rn_types(i),i
 230        format('RN  RENAME ',a8,1x,a8,1x,a16,1x,I4,4(I3.2),
     .              2x,I4,4(I3.2),3F10.4,1x,a4, '  ! # ',i6)
         endif
      end do

****  Thats all
      return
      end


