CTITLE INIT_CTOGOBS
 
      subroutine init_ctogobs

      implicit none
 
*     This routine initializes the variables in the ctogobs
*     program.  Most variables are initialized locally, here
*     we set the default running conditions when no command
*     is passed.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*     None
 
* LOCAL VARIABLES

*  i   - Loop counter

      integer*4 i

*  status_default - Character string with the default output
*                   options
      character*80 status_default 

      data status_default / 
     .' ALL -CLK_JMP_P -PASS2 -DD_E -DD_S -DD2S -PC_BIAS -PRESCAN' /
 
*     For the moment set the total number of channels to max value.
*     (If future versions of the cfile contain this information we
*     can use actuall rather than max value)
 
      num_chan = max_gchannels
*     default the range estimate of clock jump to 100 times the
*     nominal frequency stability and ~1 microsecond (which ever
*     is larger).
      rng_jump_tol = 100.d0
      rng_jump_min = 1500.d0
 
*     Default the tolerance on setting a millisecond jump.  (set
*     to 100 microseconds.  
* MOD TAH 030127: Set to 10 us from 100 us
* MDD TAH 040917: Set to 0 us for the moment.  There is something
*     wrong with the code when this happens and this will stop the
*     error (pior to this the reset_tol was foiced to zero in prescan_clk)
* MOD RWK 150120: Use the scalar clock frequency since fL1 now an array
C     reset_tol   =  fL1*10.0d-6
      reset_tol   =  fClk*10.0d-6
 
*     Set the bad range residual tolorance to 10  sigma  and 1000
*     cycles (~200 m). 
* MOD TAH 990917: Added default for rng_res_max as ~400m
      rng_res_tol  = 10.d0
      rng_res_min  = 1000.d0
      rng_res_max  = 2000.d0

*     Ionospheric delay discontinity tolerances
C     dt_ion_tol = 240.0   ! Seconds
      dt_ion_tol =  30.0   ! Seconds
      do i = 1, max_cfiles
C        ion_rw_tol(1,i) = 4.00  ! multiplier from last change
C        ion_rw_tol(2,i) = 0.8   ! minimum change before multiply
*                                ! (any jump > 0.8 cycles caught)
*                                ! (NOTE: LG here is twice the GAMIT
*                                !        plot value.)
C        ion_rw_tol(3,i) = 5.0   ! Maximum amount
* MOD TAH 970812: New high-lat defaults
         ion_rw_tol(1,i) = 6.00  ! multiplier from last change
         ion_rw_tol(2,i) = 2.0   ! minimum change before multiply
         ion_rw_tol(3,i) = 5.0   ! Maximum amount

*        Clear the widelane bias values
         WL_bias_site(i) = 0.d0

*        Set the number of azimuth mask entries to zero
         num_azmask(i) = 0
      end do

      do i = 1, max_gsvs
         WL_bias_svs(i) = 0.d0
      end do

*     Set the phase fit tolerances for range clocks and phase clocks
*     Default values (may be updated with commands) *Defaults changed
*     from 1000, 500 for range fit tol.
      phs_fit_tol(1) = 2000.d0
      phs_fit_tol(2) = 1000.d0
      phs_fit_tol(3) =  200.d0
      phs_fit_tol(4) =  100.d0

*     Set the trim_oneways defaults for data to be removed
      min_dtl_bias = 120      ! seconds
      min_good_bias =  8
      min_dtr_end  = 0.1      ! delete if less than 10% of data after
*                             ! bias flag
      min_good_end =  24

*     Defaults for double difference ceeaning:  These we will change
*     once we see see the data spacing. (Old 100 25 5)
      max_wl_ret = 100        ! number of widelane data to return
      max_dd_ret = 50         ! initial number data to return for double 
*                             ! differences
*     New default for LG (old was 5)
      max_lg_use = 10
      tol_one_way_fix = 10.d0 ! Max length of time (in seconds) to allow
*                             ! one-way data fixing.

 
*     Set the maxiumum number of iterations in getting the range
*     clock estimates.
      max_rclk_iter = 3

*     Set the editing condition on max_number of ddscan bias flags
*     that can be added to zero so that no editing will be done
*     (default)
      max_scan_edit = 0

*     Minimim one way data to be retained (default 30)
      min_ow_data = 30

*     Set the normal point size to zero so that no normal pointing
*     will be done
      np_size  = 0
      np_start = 1

*     Set the relative weight of the clock process noise relative
*     to the range and phase data noise.  (Default is 10 times
*     weaker).  A larger value gives less weight to process noise.
      rel_clk_wght = 10.0

*     Set the roots to file names
      rng_clk_root  = ' '
      phs_clk_root  = ' '
      phs_res_root  = ' '
      sng_diff_root = ' '
      dd_outfile    = ' '
      gobs_outfile  = ' '
      summary_file  = 'autcln.sum'
      mf_name       = ' '
      igs_clk_file  = ' '

*     Data use values
      min_ctog_elev        = 10.0*pi/180.d0
      min_out_elev         = 15.0*pi/180.d0  
      use_gamit_elc        = .true.
      use_cview_edit       = .false.
      usr_ignore_gaps      = .false. 
      use_MM_ranges        = .false.
      do_one_bg            = .true.
      gaps_flagged         = .true.
      remove_ms_jump       = .true.
      remove_first_bias    = .false. 
      apply_phs_clk        = .true.
      pf_remove_bf         = .false.
      use_postfit          = .false.
      write_igs_clk        = .false.
      acbias_out           = .false.
      igs_clk_samp         = 10

* MOD TAH 091228: Make defaukt be no L1 only data to be processed.
*     Use L1ONLY command to allow L1 only data to be used.
      nol1only             = .true.  ! Set so that L1 only data will
                                     ! not be processed
* MOD TAH 160826: Added option to prefit station clocks.  Default
*     is true.  Use PREFIT_CLK N to turn off
      prefit_clk           = .true.
* MOD TAH 200509: Added as option to remove millisecond clock
*     jumps at start.  Use prefit_clk commmand to turn on.
      prescan_ms           = .false.

      do i = 1, max_cfiles
         gap_size(i) = 1
      end do

      do i = 1, (max_cfiles-1)/32+1
         phs_res_sites(i) = -1
         scan_sites(i)    = -1
      end do

      num_pre_edit = 0

*     Set the default output status with character string
*     containing the default options.
      call decode_option(status_default, status_rep_opts, 
     .                   num_ctog_status, status_rep, -1 )

*     Bias flag removal (Old defaults: 12.0 3.0 3600. 2.0)
      dchi2_ratio = 10.d0
      dchi2_min_val = 3.d0     
      dchi2_max_sep = 1800.d0   ! Seconds
      dchi2_gap_fact = 5.d0

*     Cycle slip detection in Double difference
*     cleaning loop.  Set WL default for AS data
*     (Non-AS data defaults -- 4.0  0.8  4.0 
*                           -- 4.0  0.5  1.0)
      dd_wl_tol(1) =  5.0
      dd_wl_tol(2) =  2.0
      dd_wl_tol(3) = 10.0
      dd_lc_tol(1) = 3.0
      dd_lc_tol(2) = 0.35
      dd_lc_tol(3) = 0.80

****  Set up iteratation values for postfit editing and
*     phase clock estimation iteration.
*     pc_max_iter must an even number. MAx number of iterations
*     Set default iterations to 2 (since apply phase clock is
*     now turned on as default (MOD TAH 000222)
C     pc_max_iter = 30
      pc_max_iter =  2
*     Converged when RMS changes by < 0.01 (1%)
      pc_convd    = 0.01
*     Allow non-integer biases after the 0'rd iteration.
      pc_non_int   = 0
*     Set to just 0 the number of full analysis iterations
      pc_full_anal  = 0
*     Set maximum correlation in bias solution to 0.1 (closer to
*     zero this is made the slower runs).
      pc_max_corr   = 0.1 
      
*     Start edting data after the 5th iteration
      pc_start_edit = 5
*     Amount to "over-shoot" when applying mean-one-way phase
*     coorections
      pc_over_shoot = 1.5d0

*     Check residuals on any satellite whose rms is 1.5 times site
*     rms
      edit_postfit = .false.
      pf_svs_ratio =  1.5
*     Edit 4-sigma residuals.
      pf_nsig      = 4.0
*     Max residual to restore
      pf_maxres    = 0.05
*     Discard all data in pass if >30% edited.
      pf_bad_ratio = 0.30

*     Max RMS for postfit residuals
      pf_max_rms   = 0.50d0

*     Set the L2 fact so that the data are used to determine the
*     L2 type.  When doing phase clocks with L1 only data this value
*     will be set to 0.
      fdd_L2_fact  = -1

*     Default max clk rms for a site to be used as reference clock 
*     (20 cylces = 3.8 m) and max allowed for output (100 cycles = 19m)
      rms_ref_clk = 20
      rms_max_clk = 100 

*     Defaults for resolving widelanes (LC_AUTCLN command)
* MOD TAH 080716: New defaults for DIR algorithm
      resolve_wl  = .false.
      lca_type = 'DIR'
      min_wl_tol  =   50     ! Minimum number of points
      dchi_wl_tol =   1.5d0  ! Ratio test tolerance
      msig_wl_tol =  1.d-2   ! Baseline length dependent ionospheric
                             ! delay constraint (make smaller when
                             ! ionosphere is very active).
C  Defaults for SEQ algorithm (get set when lc_autcln command encountered).
C      dchi_wl_tol =   10     ! Dchi change from best choice.
C      msig_wl_tol =  0.30d0  ! Maximum sigma allowed for wl to be
                             ! resolved.
      mdev_wl_tol =  0.25d0  ! Maximum deviation from integer allowed

* MOD TAH 130330: Set the default option for DPH output to satellites 
*     and sites (consistent with eariler versions)
      dph_output = 'MULTIPLE'

* MOD TAH 180307: Default not to apply iondelay to omc in cfile
      app_ion = .false.

* MOD TAH 180321: Set default to remap glonass ambiguities
      remap_glonass = .true.

* MOD TAH 200617: See trim segments treu.
      trim_seg = .true.

      return
      end
 
CTITLE OUT_CTOGOBS_HEAD
 
      subroutine out_ctogobs_head
 
      implicit none

*     This routine prints out a header message to tell
*     user what is going on.
 
* INCLUDES
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*     None
* LOCAL VARIABLES

      character*12 hsver
 
      write(*,100) caprog_name, hsver(ctogobs_version)
100   format(/,'++++++++++++++++++++++++++++++++',/,
     .         1x,a8,'           Version ',a,/,
     .         '++++++++++++++++++++++++++++++++',/)
 
*     Thats all
      return
      end
 
CTITLE GET_CTOGOBS_RUN
 
      subroutine get_ctogobs_run

      implicit none
 
*     Routine to read the ctogobs runstring:
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*     None
 
* LOCAL VARIABLES
 
*   rcpar   - Reads the runstring
*   len_run - Length of the runstring element
*   nc      - Incremental counter for number of cfiles
*   nr      - counter through runstring
*   ir      - Counter for position in runstring.  Can differ
*             depending on if autcln of ctogobs is run.
*   i       - Loop counter
*   trimlen - Length of string
 
      integer*4 rcpar, len_run, nc, nr, ir, i, trimlen
 
*     dummy   - Dummy string for getting excess cfiles
*     use_cseries - Cfile series to use (used if dfile name is
*                   passed
*     runstring - Runstring read (mainly for program name)
 
      character*8 dummy, use_cseries
      character*256 runstring  

***** Get the name of the program to see if ctogobs or autcln

      len_run = rcpar(0,runstring)
      call get_caprog( runstring, caprog_name)

      call out_ctogobs_head
 
*     Start decoding the runstring:
*     Get the name of the command file
 
      len_run = rcpar(1, CtoGobs_com_file )
      if( len_run.gt.0 ) then
          use_command_file = .true.
      else
          use_command_file = .false.
      end if

****  See if option passed to report defaults 
      dummy = ctogobs_com_file
      call casefold(dummy)
      if( dummy(1:8).eq.'DEFAULTS' ) then
          num_cfiles = 1
          cf_codes(1) = 'DEF '
          cfiles(1)   = 'Default settings'
          write(*,120)
 120      format(/,' CURRENTS DEFAULTS ARE:')
          call rep_run_params(6)
          stop 'Defaults reported'
      end if

*     If we are running ctogobs get the name of the Gobs file and
*     description. 
      if( caprog_name(1:4).eq.'CTOG' ) then
*         Get the name of the output file or code for its generattion.
          len_run = rcpar(2, Gobs_code)
 
*         Get the experiment descpription
          len_run = rcpar(3, ctogobs_title )
          ir = 3
      else
          Gobs_code = ' '
          ir = 1
      end if

*     Get the new cfiles series (if passed)
      len_run = rcpar(ir+1, new_cfile_series)
      if( len_run.le.0 ) new_cfile_series = ' '
 
****  Now loop getting all of the cfile names
      nc = 0
      nr = 0
      len_run = 1
      do while ( len_run.gt.0 .and. nc.lt. max_cfiles)
          len_run = rcpar(ir+2+nr, runstring)

*         See if dfile name passed instead of cfile name.  Based
*         on d as first character of name.
          if( runstring(1:1).eq.'d' .or. runstring(1:1).eq.'D' ) then
              df_name = runstring
              call read_dfile( df_name, nc, cfiles, max_cfiles )
*             see if we want to update the series to be used
              len_run = rcpar(ir+3+nr, use_cseries)
              if( len_run.gt.0 ) then
                  do i = 1, nc
                     cfiles(i)(6:6) = use_cseries(1:1)
                  end do
              end if
              nr = nr + 2
          else

*             Interpret the results as a cfile name
              cfiles(nc+1) = runstring
              if ( len_run.gt.0 ) nc = nc + 1
              if ( len_run.gt.0 ) nr = nr + 1
          end if
      end do
 
*     See if we have a list of cfiles
*                         ! Probably no runstring, print help
      if( nc.eq.0 ) then
          call proper_runstring(
     .         caprog_name(1:trimlen(caprog_name)) // '.hed',
     .         caprog_name,0 )
          call proper_runstring('autcln.hlp',caprog_name,-1)
      end if
 
*     If len_run is not equal zero, see if anything more
*     in the runstring.
      do while ( len_run.gt.0 )
          len_run = rcpar(ir+2+nr, dummy)
*                                     ! Print warning
          if( len_run.gt.0 ) then
              write(*,200) nr-3, dummy(1:len_run)
 200          format(' ** WARNING ** Too many c-files.  File ',
     .            I3,1x,a,' Ignored')
          end if
          nr = nr + 1
      end do
 
      num_cfiles = nc
 
      write(*,220) caprog_name, num_cfiles
 220  format(1x,a8,'processing ',i3,' cfiles')
 
***   Thats all
      return
      end
 
CTITLE READ_CFHEADS
 
      subroutine read_cfheads( merr )

      implicit none
 
*     This routinw will loop over the cfile list and read all the
*     header records.  If it has a problem with the c-file it
*     will remove it from the list.  The names of the sites will
*     be put in the Gobs_header at this time.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/gobs_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*  merr  -- Maxiumum IOSTAT error to occurr.

      integer*4 merr
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error on reading files.
*   i,j,k   - Loop counters
*   trimlen - Length of string
 
 
      integer*4 ierr, i, j, trimlen
 
*   series  - C-file series (extracted from name)
 
      character*4 series

      character*80 message   ! Message for report_stat
 
***** Start looping over the cfiles
 
      i = 0
      merr = 0

      do while ( i.lt.num_cfiles )
          i = i + 1
          call open_cf(100+i, cfiles(i), ierr)
*                                 ! OK open, read headers
          if( ierr.eq.0 ) then
              call read_cf1(100+i,'SH', ierr)
          end if   
*****     See if it error in reading the first record
*                                     ! Probably no a cfile.
          if( ierr.ne.0 ) then
              write(*,100)  cfiles(i)(1:trimlen(cfiles(i)))
 100          format(' Problem with ',a,' IGNORING this file')

*             Move the list of cfile names down
              do j = i+1, num_cfiles 
                 cfiles(j-1) = cfiles(j)
              end do

*             Decrement counters so that we get the next cfile.
              num_cfiles = num_cfiles - 1
              close(100+i)
              i = i - 1
              merr = max(merr, ierr)
           else
 
*             Proceed with reading.  Get the site code and
*             series number from the cfile name.
              call get_code_series( cfiles(i), cf_codes(i), series)
              curr_cfile_series = series
              call casefold(cf_codes(i))

*             Now read the remaining header records.
* MOD TAH 940516:V2.07:Read all cf2 record to get the extra
*             value with the cubic term in it.
              call read_cf2(100+i,'ALL', ierr )

* MOD TAH 950710:V2.11: Save the clock poly nominial information
              do j = 1, cf_nclock
                 apr_clk_poly(j,i) = cf_clock(j)
              end do

              do j =  cf_nclock+1, cf_maxclk
                 apr_clk_poly(j,i) = 0.d0
              end do

* MOD TAH 000227: Save the long name of the site
              long_names(i) = cf_sitnam

* MOD TAH 200511: Save the fL1 and fL2 now if the values are
*             non-zero.  This is in case some c-files have
*             zero frequencie because no data (happens for
*             GNSS systems more than GPS).
* MOD TAH 200617: Added reporting of changes in frequencies.
              do j = 1, cf_nsat
                 if( cf_fL1(j).ne.0 ) then
                    if(fL1(j).ne.0 .and. fL1(j).ne. cf_fL1(j) ) then
                       write(message,150) 'L1', fL1(j), cf_fL1(j),
     .                    cf_codes(i)
 150                   format(a,' frequency changed from ',F9.3,
     .                    'MHz to ',F9.3,' Mhz site ',a)
                       call report_stat('WARNING','AUTCLN',
     .                     'read_cfheads','',message,i)
                    endif
                    fL1(j) = cf_fL1(j)
                 endif 
*                Now check L2
                 if( cf_fL2(j).ne.0 ) then
                    if(fL2(j).ne.0 .and. fL2(j).ne. cf_fL2(j) ) then
                       write(message,150) 'L2', fL2(j), cf_fL2(j),
     .                    cf_codes(i)
                       call report_stat('WARNING','AUTCLN',
     .                     'read_cfheads','',message,i)
                    endif
                    fL2(j) = cf_fL2(j)
                 endif
              end do

****          Read next record type

              call read_cf3(100+i,'ALL', ierr )

* MOD TAH 950710: Code below not needed anymore
C             Saves the station clock polynomial coefficients
C             do j = 1,3
C                apr_clk_poly(j,i) = cf_preval(j+3)
C             end do

* MOD TAH 940516: If extra value is given, it is the cubic term so
*             save.
C             if( cf_nextra.eq.1 ) then
C                 apr_clk_poly(4,i) = cf_extra(1)
C             else
C                 apr_clk_poly(4,i) = 0.d0
C             end if

* MOD TAH 971204: Save the cfile prior values
              do j = 1, cf_nparam
                 cf_apr_save(j,i) = cf_preval(j)
              end do 
          end if
 
****      continue processing if no error
          if( ierr.eq.0 ) then
              call save_cfinfo(i)
*             Tell use about the data
              write(*,200) i, cf_codes(i), rcv_types(i), cf_nsat
 200          format(' Cfile ',i3,1x,a4,' receiver SW ',a3,', ',
     .                         i3,' Satellites')
          end if

      end do

* MOD TAH 200628: Save the cf_gnss system
      sv_gnss = cf_gnss
 
****  Save the number of parameters we need to estimate
 
      num_param = num_cfiles + num_sat

****  See if we have any data (i.e., did any of the cfiles open and
*     read OK.

      if( num_cfiles.le.0 ) then
          write(*,300) caprog_name
 300      format(' NO VALID Cfiles READ; ',a8,'Terminating')
          call report_stat('fatal','autcln','read_cfiles',' ',
     .         'No valid cfiles read',0)
      end if
 
****  Thats all
      return
      end
 
CTITLE SAVE_CFINFO
 
      subroutine save_cfinfo(cf)

      implicit none
 
*     This routine will save the necessary cfile header information
*     for processing the cfiles.  Some information is based on the
*     first cfile only (and check for subsequet ones).
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   cf      - Cfile number being checked
 
 
      integer*4 cf
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
 
      integer*4 i,j
 
****  Save/check the number of epoch.  (In the future we coudl
*     introduce a "series number here so that cfiles which differ
*     would be treated as a different series.  For the moment they
*     must be compatible with having passed through SOLVE)
 
*                         ! save common values
      if( cf.eq.1 ) then
          num_ep = cf_nepoch
          num_sat = cf_nsat
          sampling = cf_inter
          do i = 1, num_sat
*             Save the list of PRN numbers
              prn_list(i) = cf_ischan(i)  
*             For each PRN number save where it is in list
              prntol(prn_list(i)) = i    
          end do
      else
 
*         check common information the same
          if( num_ep.ne.cf_nepoch ) call cf_unmatch('EPOCHS',cf)
          if( num_sat.ne.cf_nsat  ) call cf_unmatch('NSAT',cf)
          if( sampling.ne.cf_inter) call cf_unmatch('INTER',cf)
          do i = 1, num_sat
              if( prn_list(i).ne.cf_ischan(i))
     .                     call cf_unmatch('PRN Number',i*100+cf)
          end do
      end if
 
*     save the original sampling interval for the receiver
 
      orig_sampling(cf) = cf_ircint
 
*     Save the number of data types for this cfiles
      num_gdata_types(cf) = cf_ndat

*     Save the datatype codes
      do i = 1, cf_ndat
          data_types_s(i,cf) = cf_dattyp(i)
      end do

      do i = cf_ndat +1, max_gdata_types
          data_types_s(i,cf) = 0
      end do
 
      do i = 1, num_sat
          do j = 1, cf_ndat
             lambda(i,j,cf) = cf_lambda(i,j)

* MOD TAH 970811: Check that none are zero
* MOD TAH 990511: Removed default setting so that we can process L1 only data
C            if( j.le.2 .and. cf_lambda(i,j).eq.0 ) then
C               lambda(i,j,cf) = -1.d0
C               write(*,300) cf_codes(cf), prn_list(i)
C300            format('** WARNING ** Site ',a,' PRN ',i2.2,
C    .                ' has zero wavelength factor. Set to -1')
C            end if
* MOD TAH 091228: See if wavelength factor is consistent with data_types.
*           (ie., if a data_type is zero, then factor should be zero).
               if( data_types_s(j,cf).eq.0 ) then
                  lambda(i,j,cf) = 0
               end if 
          end do
          do j = cf_ndat+1,4
              lambda(i,j,cf) = 0
          end do
      end do
 
      rcv_types(cf) = cf_rcvrsw
 
 
****  Thats all
      return
      end
 
CTITLE CF_UNMATCH
 
      subroutine cf_unmatch( type, file )

      implicit none
 
*     Routine to report the header values in the cfiles do not
*     match.
 
* INCLUDES
*     None
* PASSED VARIABLES
 
*   file        - C-file number with PRN number coded in upper part
 
 
      integer*4 file
 
*   type        - Type of parameter tested.
 
 
      character*(*) type
 
      write(*,100) file, type
 100  format(' ** WARNING ** for PRN+Cfile ',i5,1x,a,' does not match')
 
      return
      end
 
CTITLE GET_CODE_SERIES
 
      subroutine get_code_series( cfile, code, series )

      implicit none
 
*     General routine to analyze the cfile name and return the
*     station code and series name from the name.  The name
*     structure is assumed to be:
*     <directory>/cCODES.???<trailing>
*     where directory is is the directory
*         CODE is the four character code for the site
*         S    is the cfile series.
*     The routine works by finding the last / character and
*     decoding name after that.
 
* INCLUDES
*  None
* PASSED VARIABLES
 
*   cfile   - Name of the cfile
*   code    - code extracted from name
*   series  - series extracted from name
 
 
 
      character*(*) cfile, code, series
 
* LOCAL VARIABLES
 
*   len_name    - Length of the c-file name
*   i           - loop counter
*   trimlen     - Returnes length of string
 
 
 
      integer*4 len_name, i, trimlen
 
****  Start at the end of the string in try to find
*         last /
 
      len_name = trimlen(cfile)
      if ( Len_name.eq.0 ) then
          write(*,100)
 100  format(' **WARNING** Zero length cfile name passed to',
     .        ' GET_CODE_SERIES.  Returning blanks')
          code   = ' '
          series = ' '
          RETURN
      end if
      i = len_name
      do while ( i.gt.1 .and. cfile(i:i).ne.'/')
          i = i - 1
      end do
 
*     Adjust i so that it points to c in cfile name
      if( cfile(i:i).eq.'/' ) i = i + 1
 
*     Save the code and series
      code = cfile(i+1:i+4)
      series = cfile(i+5:i+5)
 
*     Thats all
      return
      end
 
CTITLE CTOG_MEM
 
      subroutine ctog_mem(vma_data)

      implicit none
 
*     This routine computes the amount of memory needed for
*     this run and then calls malloc to set aside the memory
*
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/gobs_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*   vma_data    - Area where arrays will be saved.
 
 
      integer*4 vma_data
 
* LOCAL VARIABLES
 
*   tot_memory  - Total memory required (I*4 words)
*   offset      - Offset between location and start of
*               - vma_data array
*   cse         - Product of channels by sites by epochs
*   cpe         - Product of parameters by epochs
*   trimlen     - Length of string
*   memassign   - Integer*8 Function to assign memory
 
      integer*4 tot_memory, cse, cpe, trimlen
      integer*8 offset, memassign

      character*256 outline
 
****  Compute the memory needed. (32 padding is just for
*     margin of error)
      cse = num_chan*num_cfiles*num_ep
      cpe = (num_sat+num_cfiles)*num_ep
 
* NOTE:  If more arrays are added be sure to update the total
*        number of I*4 words of memory to allocate.
      tot_memory = cse*15 + cpe*3 + 8
      
* MOD TAH 970828: See if we need to add space for azel values
      if( (trimlen(phs_res_root).gt.0) .or. 
     .    (apply_phs_clk .and. pc_max_iter.gt.1)  ) then
          tot_memory = tot_memory + cse*2
          write(*,'(a)') 'Allocating memory for azel values'
      end if    
      if( (np_size.gt.1 .and. use_postfit) ) then 
          tot_memory = tot_memory + cse*2
          write(*,'(a)') 'Allocating memory for Normal points'
      end if

* MOD TAH 000222: Add space for satellite clock corrections
*     if needed
      if(trimlen(igs_clk_file).gt.0 ) then
          tot_memory = tot_memory + num_sat*num_ep*2
      end if
 
      write(outline,100) tot_memory/(1024.0**2)*4
 100  format('Allocating ',F8.2,' Mbytes for run')
* MOD TAH 980417: removed extra 'autcln' argument added.
      call report_stat('status',caprog_name,'ctog_mem',' ',
     .       outline,0)
 
****  See if we can get memory
      offset = memassign(tot_memory,1,loc(vma_data))
 
****  See if we have enough.
 
      if ( offset.eq.0 ) then
          write(*,120) caprog_name
 120      format(' *** DISASTER *** Not enough memory to run',a8,/,
     .            ' Either run on a larger machine or reduce',
     .            ' number of cfiles.')
          call report_stat('fatal',caprog_name,'ctog_mem',' ',
     .         'Not enough memory to run program',0)
      end if
 
****  see where we need to start in vma_data 
      iL1r_phs_cse   = offset
      iL2r_phs_cse   = iL1r_phs_cse   + cse*2
      iL1r_rng_cse   = iL2r_phs_cse   + cse*2
      iL2r_rng_cse   = iL1r_rng_cse   + cse*2
      iL1_cyc_cse    = iL2r_rng_cse   + cse*2
      iL2_cyc_cse    = iL1_cyc_cse    + cse*2
      ictol_cse      = iL2_cyc_cse    + cse*2
      idata_flag_cse = ictol_cse      + cse
      ibf_type_cse   = idata_flag_cse + cse
      iparams_cse    = ibf_type_cse   + cse
      
      iazel_cse      =  iparams_cse   + cpe*2

* MOD TAH 970828: Add amounts for azel and mfile adjustements if needed
      if( (trimlen(phs_res_root).gt.0) .or. 
     .    (apply_phs_clk .and. pc_max_iter.gt.1)  ) then
          ipf_dphs_cse =  iazel_cse   + cse*2          
      else 
          ipf_dphs_cse =  iazel_cse 
      end if

* MOD TAH 000222: See if we need to save the satellite clock
*     "removed" term (We need num_sat*num_ep real*8 values)
      if( np_size.gt.1 .and. use_postfit ) then
          isvcL1_ce = ipf_dphs_cse + cse*2
      else
          isvcL1_ce = ipf_dphs_cse
      end if
      
*     pf_phs_cse will need cse*2 words of storage.      
      if( trimlen(igs_clk_file).gt.0 ) then
          ipar_flag_cse = isvcL1_ce + num_sat*num_ep*2
      else
          ipar_flag_cse = isvcL1_ce 
      end if

****  We now have all the memory allocated that we will
*     need.
      return
      end
 
CTITLE CERRFL_TO_DF
 
      subroutine cerrfl_to_df(cf_ierfl, cf_data_flag, data_flag)

      implicit none
 
*     This routine will map the cfile error flag to the appriopriate
*     bits in data_flag (see Gobs_data.h) for definition.
 
 
* INCLUDES
*   None
 
* PASSED VARIABLES
 
*   cf_ierrfl   - Cfile I*2 error flag
*   cf_data_flag - Cfile data flag (needed for bits 28 and 29) 
 
      integer*2 cf_ierfl, cf_data_flag
 
*   data_flag   - Gobs_file data flag (bit mapped)
 
      integer*4 data_flag

      logical kbit
 
* LOCAL VARIABLES
*     None
 
*     Test each of the possible cases
* MOD TAH 080622: Save data bits 28 and 29 which are the dcb flags
      data_flag = 0
      if( kbit(cf_data_flag,28) ) call sbit(data_flag,28,1)
      if( kbit(cf_data_flag,29) ) call sbit(data_flag,29,1)

 
      if( cf_ierfl.ne.0 ) then
*                                         ! Unweighted
          if( cf_ierfl.eq. -1 ) then
              call sbit(data_flag, 14, 1)
*                                         ! Deleted data
          else if( cf_ierfl.eq.  2 ) then
              call sbit(data_flag, 4,  1)
*                                         ! Low SNR
          else if( cf_ierfl.eq.  3 ) then
              call sbit(data_flag,  2, 1)
*                                         ! Low elev
          else if( cf_ierfl.eq.  4) then
              call sbit(data_flag, 15, 1)
*
          else if( cf_ierfl.eq.  7) then  ! Bad MODEL due phase center cutoff
              call sbit(data_flag, 21, 1)
*                                         ! Bias flag
          else if( cf_ierfl.eq. 10 ) then
              call sbit(data_flag, 32, 1)
*                                         ! Unknown Cfile flag
          else
              call sbit(data_flag, 21, 1)
          end if
      end if
 
****  Thats all
      return
      end
 
CTITLE INIT_RECEIVER  
 
      subroutine init_receiver  

      implicit none
 
*     This routine will set the receiver clock allan stadard
*     devaitions based on the receiver type.  (May be changed
*     later with the command file).  The site dependent elevation
*     angles and SNR's are also set.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*     None
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
 
 
      integer*4 i

****  Set the return values based on receiver sampling rate

      max_wl_ret = 100        ! number of widelane data to return
      max_dd_ret = (50*30)/sampling   ! initial number data to return for double 
*                             ! differences
* MOD TAH 960821: Check to make sure value is not too large

      if( max_dd_ret.gt. max_max_dd_ret ) max_dd_ret=max_max_dd_ret/2      
      if( max_dd_ret.lt.20 ) max_dd_ret = 20
      
*     New default for LG (old was 5)
      max_lg_use = (10*30)/sampling
      if( max_lg_use.lt.5 )  max_lg_use = 5
      if( max_lg_use.gt.max_dd_ret ) max_lg_use =max_max_dd_ret/4
      
****  Loop over the receivers
 
      do i = 1, num_cfiles

          if( rcv_types(i).eq.'ROG ' ) then
              rc_allan_sd(i) = 0.1d-9
              site_oelev(i) = 20.0*pi/180.d0
              site_celev(i) = 20.0*pi/180.d0
* If decoded with sxr then 5 and 4 are more appropriate
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'MIN ' ) then
              rc_allan_sd(i) = 1.d-12
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 15.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'COR ' ) then
              rc_allan_sd(i) = 1.d-9
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 15.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'GES ' ) then
              rc_allan_sd(i) = 1.d-9
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 15.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'TRM ' ) then
              rc_allan_sd(i) = 3.d-8
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 10.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'ASH ' ) then
              rc_allan_sd(i) = 3.d-8
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 10.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          else if( rcv_types(i).eq.'TRB'  ) then
              rc_allan_sd(i) = 3.d-8
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 10.0*pi/180.d0
*             Strictly these values are for SRX decoding
* MOD TAH 980513: Changed defaults to 2/2 from 4/5
              site_snr(1,i) = 2
              site_snr(2,i) = 2
*                                                 ! Unknow type set big
          else
              rc_allan_sd(i) = 3.d-8
              site_oelev(i) = 15.0*pi/180.d0
              site_celev(i) = 15.0*pi/180.d0
              site_snr(1,i) = 2
              site_snr(2,i) = 2
          end if
      end do

*     Initialize the satellite clock allan std dev.
      do i = 1, num_sat
          rc_allan_sd(i+num_cfiles) = 1.d-9
      end do
 
*     Now set the data noise based on data_type_s
      do i = 1,num_cfiles
*                                             ! P-code L1
          if( data_types_s(3,i).eq.3 ) then
              rng_noise(i) = 5.d0
*                                                 ! CA-code L1
          else if( data_types_s(3,i).eq.5) then
              rng_noise(i) =  25.d0
          else
              rng_noise(i) =  50.d0
 
*             Tell user we don't know data type
              write(*,100) cf_codes(i), data_types_s(3,i)
 100          format(' ***WARNING**** Unknown range type at ',a4,
     .                ' data_type is ',i3)
          end if

*         Set all of the phase noise values to 0.1 cycles
          phs_noise(i) = 0.1d0
      end do
 
****  Thats all
      return
      end
 
CTITLE UPDATE_INIT_CLK
 
      subroutine update_init_clk( params_cse, par_flag_cse )

      implicit none
 
*     This routine will scan the list of clock estimates at the
*     sites and satellites until it finds the first valid one.
*     These values are then sayed as initial values for the
*     next iteration.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   par_flag_cse(num_param,num_ep)  - Parameter quality flag.  Used here
*           on the second iteration to set the apriori clock values.
 
 
      integer*4 par_flag_cse(num_param,num_ep)
 
*   params_cse(num_param,num_ep)   - Estimates of the clock at the first
*                       - epoch of data.
 
 
      real*8 params_cse(num_param,num_ep)
 
* LOCAL VARIABLES
 
* j - Loop counter
 
      integer*4 i, j
 
* found - Logical to indicate that we have found a valid first
*         value
*  kbit - Function to test bit.
 
 
      logical found, kbit
 
***** Scan to find the first valid value.  Do as a while loop
*     in case we have a site with no data.  If we don't find a
*     value then probably no data and therefore we don't need
*     to do anything special.
      do i = 1,num_param
          found = .false.
          j = 0
          do while ( .not.found .and. j.lt.num_ep )
              j = j + 1
              if( .not.kbit(par_flag_cse(i,j),1) ) then
                  init_clk_val(i) = params_cse(i,j)
                  found = .true.
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE INIT_CLK
 
      subroutine init_clk( iter )

      implicit none
 
*     This routine will initialize the information about the clocks.
*     It sets the apriori values for the site clocks, and converts the
*     Allan std. dev. to a step variance appropriate for spacing of the
*     data,
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   iter  - current iteration
 
 
      integer*4 iter
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   ns      - Number of stations (short version)
 
 
      integer*4 i, ns
 
*     Loop over the ground sites setting apriori values and
*     variances
 
      do i = 1, num_cfiles
 
          apr_clk_val(i) = init_clk_val(i)
C         if( iter.gt.1 ) apr_clk_epoch(i) = 0
 
*         The 100.d0 (sec) here is the nominal time at which
*         the Allan Std.dev. is given.
          apr_clk_var(i) = (rc_allan_sd(i)*fClk)**2*100.d0*sampling
          apr_clk_sd(i) = sqrt(apr_clk_var(i))
 
      end do
 
*     Set the satelite clock epoch estimates to 0 (we have none) in
*     the first iteration.  In the second we will have the values
*     from the first run.
      ns = num_cfiles
      do i = 1, num_sat
          apr_clk_val(i+ns) = init_clk_val(i+ns)
          apr_clk_epoch(i+ns) = 0

*         The 100.d0 (sec) here is the nominal time at which
*         the Allan Std.dev. is given.
          apr_clk_var(i+ns) = (rc_allan_sd(i+ns)*fClk)**2*
     .                    100.d0*sampling
          apr_clk_sd(i+ns) = sqrt(apr_clk_var(i+ns))
      end do
 
***** Thats all
      return
      end
 
CTITLE INIT_NORM
 
      subroutine init_norm( ep, know_svs_clk )

      implicit none
 
*     This routine will clear the normal equations and solution
*     vector for each epoch of data processed.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ep      - Epoch number (used to see how much clock may have
*           - changed.
 
      integer*4 ep

*   know_svs_clk  -- Logical to say we know the satellite clocks

      logical know_svs_clk
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
*   ns          - NUmber of cfiles (used for convenience when
*                 satellites are used.)

      integer*4 i,j, ns

*   dep         - Change in epoch since last values

      real*8 dep 

****  If this is epoch 1, set all satellite clocks as unknown.
*     (This is because satellite clock information is not passed
*     in the cfiles).

      if( ep.le.1 ) then
          do i = 1, num_sat
             svs_clk_known(i) = know_svs_clk 
          end do
      end if

****  Firstly clear normal equations and solution vector
      do i = 1, num_param
          sol_eq(i) = 0.d0
          do j = 1, num_param
              norm_eq(i,j) = 0.d0
          end do
      end do
 
****  Now put on the apriori weights for each of the clocks.
      do i = 1,num_cfiles
          dep = ep-apr_clk_epoch(i)
*         First case: This is first data epoch from station and
*         so give it a high weight (using rclock as apriori)
          if( dep.eq.0 ) dep = .1
*         second case: data has not started yet.  Set to 1 so that
*         negative will not be put in normal equations. 
          if( dep.lt.0 ) dep =  1
          norm_eq(i,i) = 1.d0/(apr_clk_var(i)*dep)/rel_clk_wght
      end do

****  Set the apriori variances on those satellite clocks that we
*     know

      ns = num_cfiles
      do i = 1, num_sat
          if ( svs_clk_known(i) ) then
              dep = ep-apr_clk_epoch(i+ns)
              if( dep.le.0 ) dep = 1000
*             Initializse clock solution (MOD TAH 980606)
              sol_eq(i+ns) = 0.d0
              norm_eq(i+ns,i+ns) = 1.d0/(apr_clk_var(i+ns)*dep)/
     .                                   rel_clk_wght
          end if
      end do
     
 
***** Thats all
      return
      end
 
CTITLE good_bf
 
      logical function good_bf( data_flag, kine_flag, mask )

      implicit none
 
*     Function which returns true if the data flags and the kine_flag
*     say the data is OK, and their is a boas flag on the data point.
 
* PASSED VARIABLES
 
*   data_flag       - Data flag as defined in gobs_data.h
*   kine_flag       - Kinematic flags as defined in gobs_data.h
*   mask            - Mask to be applied to data_flag
*   loc_mask - Local mask with bias flags turned off.
 
      integer*4 data_flag, kine_flag, mask, loc_mask

* LOCAL VARIABLES

*   data_OK  - function that returns true if the data is OK
*   kbit     - function that returns true is bit is set

      logical data_OK, kbit
 
***** And the data flag with mask and see what we get.
*     Make sure that the bias flags are not in the mask.
      loc_mask = mask
      call sbit(loc_mask,31,0)
      call sbit(loc_mask,32,0)

 
      if( data_OK(data_flag, kine_flag, loc_mask) .and.
     .    (kbit(data_flag,31) .or. kbit(data_flag,32)) ) then
          good_bf = .true.
      else
          good_bf = .false.
      end if
 
****  Thats all
      return
      end

CTITLE  SET_RNG_MASK
 
      subroutine  set_rng_mask( rng_mask_loc )

      implicit none
 
*     Routine to set the range mask data editing conditions.  Here
*     we mask out low-elevation data flagged by GAMIT.
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
 
*   rng_mask_loc        - Mask to be and'd with data_flag to
*                   - to see if we use data
 
 
      integer*4 rng_mask_loc
 
***** Turn all bits on and then turn off those we don't want
      rng_mask_loc = -1
      call sbit(rng_mask_loc,31,0)
      call sbit(rng_mask_loc,32,0)
*                             ! GAMIT Low-elev bit.
      if( .not.use_gamit_elc ) call sbit(rng_mask_loc,15,0)
*                             ! GAMIT Marginal data bit
      if( .not.use_cview_edit ) call sbit(rng_mask_loc,14,0)

*     Set to ingnore bit showing below final elevation cutoff
      call sbit(rng_mask_loc,27,0)

****  If user has said to use minimac ranges, then clear this
*     bit from the rng_mask.
      if( use_MM_ranges) call sbit(rng_mask_loc,26,0)

*     Ignore the DCB bits 28 and 29
      call sbit(rng_mask_loc,28,0)
      call sbit(rng_mask_loc,29,0)

 
****  Thats all
      return
      end
 
CTITLE  SET_PHS_MASK
 
      subroutine  set_phs_mask( phs_mask_loc, phs_bias_mask_loc )

      implicit none
 
*     Routine to set the phase mask for data editing conditions.  Here
*     we mask out low-elevation data flagged by GAMIT.

      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   phs_mask_loc    - Mask to be and'd with data_flag to
*                   - to see if we use data.  This mask ignores the
*                     bias flags
*   phs_bias_mask_loc   - Mask to be and'd with data flag.  This maskes
*                     also checks the bias parameters.
 
      integer*4 phs_mask_loc, phs_bias_mask_loc
 
***** Turn all bits on and then turn off those we don't want
      phs_mask_loc = -1
      phs_bias_mask_loc = -1
 
*     For the phase mask turn off the checks on the bias flags
      call sbit(phs_mask_loc,31,0)
      call sbit(phs_mask_loc,32,0)
*                             ! GAMIT Low-elev bit.
      if( .not.use_gamit_elc ) call sbit(phs_mask_loc,15,0)
      if( .not.use_gamit_elc ) call sbit(phs_bias_mask_loc, 15, 0 )

*     Set to ingnore bit showing below final elevation cutoff
      call sbit(phs_mask_loc,27,0)
      call sbit(phs_bias_mask_loc,27,0)

*                             ! GAMIT Marginal data bit.
      if( .not.use_cview_edit ) call sbit(phs_mask_loc,14,0)
      if( .not.use_cview_edit ) call sbit(phs_bias_mask_loc, 14, 0 )

*     Set mask to ignore ranges marked as minimac
      call sbit(phs_mask_loc,26,0)
      call sbit(phs_bias_mask_loc,26,0)

*     Set mask to ignore no double difference data.  This bit is set
*     when phase clocks are computed after cleaning.  It can be ignored
*     except when we are computing statistics and restoring data
* MOD TAH 040513: Turned off ignoring no-double differences. Should
*     solve problem when bias fixing of getting bias flags on sequences
*     of data that can not be fixed.
C     call sbit(phs_mask_loc,23,0)
C     call sbit(phs_bias_mask_loc,23,0)

*     Ignore the DCB bits 28 and 29
      call sbit(phs_mask_loc,28,0)
      call sbit(phs_mask_loc,29,0)
      call sbit(phs_bias_mask_loc,28,0)
      call sbit(phs_bias_mask_loc,29,0)
 
****  Thats all
      return
      end
 
CTITLE RNG_OMC
 
      real*8 function rng_omc(L1r_rng,L2r_rng, rcv_clk, svs_clk, f1, f2)

      implicit none
 
*     Routine to compute the raange residual useing the current
*     clock values
 
* INCLUDES
 
      include '../includes/const_param.h'
 
* PASSED VARAIBLES
 
*   L1r_rng, L2r_rng    - L1 and L2 range residuals (if L2 rng is
*                       - is zero means just signal frquency range)
*   rcv_clk, svs_clk    - Current estimates of receiver and satellite
*                       - clocks (L1 cycles)
*   f1, f2              - L1 and L2 frequencies, named for consistency
*                         with naming in sb phs_omc, avoiding confusion
*                         with arrays fL1 and fL2 in ctogobs.h, used in
*                         phs_omc (but not here)
 
 
      real*8 L1r_rng, L2r_rng, rcv_clk, svs_clk, f1, f2
 
* LOCAl VARIABLES
 
*    res                - Intermediate residual calculation.
 
 
      real*8 res
 
***** See if we want dual frequency range correction
      if( L2r_rng.ne.0.d0 ) then
          res = (L1r_rng - (f2/f1)*L2r_rng)/
     .                (1.d0 - (f2/f1)**2)
      else
          res = L1r_rng
      end if
 
****  Now remove in the apriori clock information
 
      res = res - (rcv_clk - svs_clk )
 
      rng_omc = res
 
****  Thats all
      return
      end
 
CTITLE jmp_sort
 
      subroutine jmp_sort( num, mean, list )

      implicit none
 
*     This routine uses an exchande sort algormithm to sort
*     abs(mean) into descending order.  There are num values
*     in ilist.
 
*   num     - Number of values to be sorted
*   list(num)  - of the original indexing of the data
 
 
 
      integer*4 num, list(num)
 
 
      real*8 mean(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   biggest_one    - Smallest integer in current pass.
*   swap   - Value used to swap real*8
*   iswap  - the index to be swapped.
 
 
 
 
      integer*4 i,j, biggest_one , iswap
 
 
      real*8 swap
 
****  set up list for indexing
      do i = 1, num
         list(i) = i
      end do
 
****  Start loop using exchange sort
 
      do i = 1, num
          biggest_one = i
          do j = i+1, num
              if( abs(mean(j)).gt. abs(mean(biggest_one)) ) then
                  biggest_one = j
              end if
          end do
 
*****     See if we should swap
          if( biggest_one.gt. i ) then
              swap = mean(biggest_one)
              iswap = list(biggest_one)
              mean(biggest_one) = mean(i)
              list(biggest_one) = list(i)
              mean(i) = swap
              list(i) = iswap
          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end
 
CTITLE EST_CLK_RNG
 
      subroutine est_clk_rng(ep,  L1r_rng_cs, L2r_rng_cs,
     .                ctol_cs, data_flag_cs, params_cs,
     .                par_flag_cs )

      implicit none
 
*     Subroutine to estimate the clocks at this epoch using the
*     range data.  The statitics on the range data are also
*     accumulated in this pass.  If there is "large range" error
*     then it will edited here as well.  The rng_res_tol is set
*     by default and my be changed in the command file.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*    ctol_cs(num_chan, num_cfiles)  - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   par_flag_cs(num_param)         - parameter estimate flag
*   ep      - Epoch number being processed.
 
 
      integer*4 ctol_cs(num_chan, num_cfiles),
     .    data_flag_cs(num_chan, num_cfiles),
     .    par_flag_cs(num_param), ep
 
*   L1r_rng_cs(num_chan, num_cfiles) - L1 range measurements (L1 cylces)
*   L2r_rng_cs(num_chan, num_cfiles) - L2 range measurements (L2 cylces)
*   params_cs(num_param)    - Clock parameter estimates (L1 cycles)
 
 
      real*8 L1r_rng_cs(num_chan, num_cfiles),
     .    L2r_rng_cs(num_chan, num_cfiles), params_cs(num_param)
 
* LOCAL VARIABLES
 
 
*  i,j,k        - Loop counters
*   ns          - parameter number of satellite
*   ipivot(max_cfiles_max_gsvs)  - Pivot elmenent for inversion
*   site_with_worst - Site with the worst res/sigma value.
*   chan_with_worst - Satellite with the worst res/sigma value.
 
 
      integer*4 i,j, ns, ipivot(max_neq),
     .    site_with_worst, chan_with_worst
 
*   data_ok     - Logical function which is true if data OK
*   kbit        - Logical function to see if bit is set.
 
 
      logical data_ok, kbit
 
*   res         - Generic residual values
*   rng_omc     - Function to return range residual
*   data_var    - Data varaince at one site
*   scale(max_neq)  - Scale factors for invert_vis
*   worst_res   - Worst postfit range residual.  Scan the residuals
*                 first to find this value.  If it is outside of
*                 tolerance then delete the range.  In the next
*                 iteration, the range will be computed with out
*                 this value and the quality of the remaining data
*                 will be cheched.
 
 
 
      real*8 res, rng_omc, data_var, scale(max_neq),
     .    worst_res
 
***** Loop over the data incrementing the normal equations.
*     Use the previous iterations estimate of the clock offset
*     (Initially these are the rclock estimates (corrected for
*      the polynomial, and then later they are estimates from
*      the range data)
C     do i = 1, num_param
C        apr_clk_val(i) = params_cs(i)
C     end do
 
      do i = 1,num_cfiles
          data_var = rng_noise(i)**2
          do j = 1, actual_max_chan
              if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
 
*                 Compute OminusC
                  ns = num_cfiles + ctol_cs(j,i)
                  res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(ctol_cs(j,i)), fL2(ctol_cs(j,i))  )
 
*                 Increment normal equations based on the data quality
                  norm_eq(i,i)   = norm_eq(i,i)   + 1.d0/data_var
                  norm_eq(ns,ns) = norm_eq(ns,ns) + 1.d0/data_var
                  norm_eq(ns,i)  = norm_eq(ns,i)  - 1.d0/data_var
                  norm_eq(i,ns)  = norm_eq(i,ns)  - 1.d0/data_var
                  sol_eq(i)      = sol_eq(i)      + res/data_var
                  sol_eq(ns)     = sol_eq(ns)     - res/data_var
*                 Show that we have some data on this rcv/sat
                  call sbit(par_flag_cs(i),1,0)
                  call sbit(par_flag_cs(ns),1,0)
              end if
          end do
      end do
 
****  Check the parameter flag to see if we have data on a s
*     parameter.  If we no not set its diagonal to one.
      do i = 1, num_param
*                                         ! No data
          if( kbit(par_flag_cs(i),1) )  then
              norm_eq(i,i) = 1.d0
          end if
      end do
 
*     Now invert the system
      call invert_vis(norm_eq, sol_eq, scale, ipivot, num_param,
     .                max_neq,1)
 
*     Save as apriori and save values in parameter array
      do i = 1, num_param
          params_cs(i) = apr_clk_val(i) + sol_eq(i)
          apr_clk_val(i) = params_cs(i)
 
*         if we have data, save this epochs value as the
*         apriori for the next epoch
          if( .not.kbit(par_flag_cs(i),1) )  then
              apr_clk_epoch(i) = ep

*             if this is a satellite set that we know know the 
*             clock value approximately.
              if( i-num_cfiles.gt.0 ) then
                  svs_clk_known(i-num_cfiles) = .true.
              end if
          end if
      end do
 
***** Scan all the residuals to find the worst.  If it
*     is outside tolerance then edit the point and do not
*     include in RMS calculation.
 
      worst_res = 0.d0
      site_with_worst = 0
 
      do i = 1,num_cfiles
          do j = 1, actual_max_chan
              if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
 
*                 Compute OminusC
                  ns = num_cfiles + ctol_cs(j,i)
                  res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(ctol_cs(j,i)), fL2(ctol_cs(j,i) ) )
 
                  if( abs(res/rng_noise(i)).gt.worst_res .and.
     .                abs(res).gt.rng_res_min ) then
                      worst_res = abs(res/rng_noise(i))
                      site_with_worst = i
                      chan_with_worst = j
                  end if
               end if
           end do
       end do
 
*     Sum the residual into the rms accumulators if
*     it is in tolerance.  Else set the data_flag_cse
*     to indicate bad data. (If nobody exceeded the min limit
*     (rng_res_min) then worst_res will be zero.)
      if( worst_res .gt. rng_res_tol .and. site_with_worst.gt.0 ) then
 
*         Recompute the worst residual
          i = site_with_worst
          j = chan_with_worst
          ns = num_cfiles + ctol_cs(j,i)
          res = rng_omc( L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                   apr_clk_val(i),apr_clk_val(ns),
     .                   fL1(ctol_cs(j,i)), fL2(ctol_cs(j,i) ) )
*         Set the bad range bit (see Gobs_data.h)
          write(*,200) cf_codes(i), prn_list(ctol_cs(j,i)), ep,
     .           res/fL1(ctol_cs(j,i))*vel_light
 200      format(' Editing Range at ',a4,' on PRN ',
     .           I2.2,' Epoch ',I5,' Residual ',
     .           F20.3,' m')
          call sbit(data_flag_cs(j,i),13,1)
      end if

*     Check to see if we haveany flagged range data which is really OK
      do i = 1,num_cfiles
          do j = 1, actual_max_chan
             if( kbit(data_flag_cs(j,i),13) ) then  ! Flagged as bad
*                 see if OK.
*                 Compute OminusC
                  ns = num_cfiles + ctol_cs(j,i)
                  res = rng_omc( L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                           apr_clk_val(i),apr_clk_val(ns),
     .                           fL1(ctol_cs(j,i)), fL2(ctol_cs(j,i) ) )
                  if( abs(res/rng_noise(i)).lt.rng_res_tol .or.
     .                abs(res).lt. rng_res_min ) then 
                      write(*,300) cf_codes(i), 
     .                       prn_list(ctol_cs(j,i)), ep,
     .                       res/fL1(ctol_cs(j,i))*vel_light
 300                  format(' Restoring Range at ',a4,' on PRN ',
     .                       I2.2,' Epoch ',I5,' Residual ',
     .                       F20.3,' m')
                      call sbit(data_flag_cs(j,i),13,0)
                  end if
             end if
          end do
      end do
 
*     Now accumulate statistics on data
      do i = 1,num_cfiles
          do j = 1, actual_max_chan
              if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
 
*                 Compute OminusC
                  ns = num_cfiles + ctol_cs(j,i)
                  res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                    apr_clk_val(i),apr_clk_val(ns),
     .                    fL1(ctol_cs(j,i)), fL2(ctol_cs(j,i))  )
 
                  sum_rng_var(i) = sum_rng_var(i) + res**2
                  sum_rng_num(i) = sum_rng_num(i) + 1

              end if
          end do

      end do
 
****  Thats all
      return
      end
      
CTITLE DATA_OK
 
      logical function data_ok( data_flag, kine_flag, mask )

      implicit none
 
*     Function which returns true if the data flags and the kine_flag
*     say the data is OK.   The user should set mask such that bias
*     flags are not detected.  (Bits 31 and 32 set off)
 
* PASSED VARIABLES
 
*   data_flag       - Data flag as defined in gobs_data.h
*   kine_flag       - Kinematic flags as defined in gobs_data.h
*   mask            - Mask to be applied to data_flag
*   cand            - Common version of iand function in HP1000/libhp1000
 
      integer*4 data_flag, kine_flag, mask, cand
 
***** And the data flag with mask and see what we get
 
      if( cand(data_flag, mask).eq.0 .and. kine_flag.eq.0 ) then
          data_OK = .true.
      else
          data_ok = .false.
      end if
 
****  Thats all
      return
      end
        
CTITLE WRITE_CLK
 
      subroutine write_clk( root, params_cse, par_flag_cse )

      implicit none
 
*     This routine will write out the clock estimates into
*     separate files for each site and satellite.  It will then
*     generate a plot command file to plot the values.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   par_flag_cse(num_param,num_ep)         - parameter estimate flag
 
 
      integer*4 par_flag_cse(num_param,num_ep)
 
*   params_cse(num_param)    - Clock parameter estimates (L1 cycles)
 
 
      real*8 params_cse(num_param,num_ep)
 
*   root     - Root name for the output files
 
 
      character*(*) root
 
* LOCAL VARIABLES
 
*   i,j  - Loop counter over stations and  epochs
*   num_out  - Records number of values output
*   trimlen  - Length of string
*   ierr     - IOSTAT error
*   plot_on_page - Keeps track of number of plots on the
*                  page (currently 4 per page)
 
 
      integer*4 i,j, num_out, trimlen, ierr, plot_on_page
 
*   kbit     - Function to test if bit is set
 
 
      logical kbit
 
*   min_sc, max_sc, ds  - Min and max values and difference for
*                         scales
*   omin_sc,imax_sc     - Original min and max scales for working
*            where on the plot the labels should go.
*   jump  -  Value needed to adjust plot for jumps in the
*            clock
*   clk   - clock value in usec after accounting for jumps
*   lab_offset - Estimate of offset in label (vertically) to
*           try to stop the text overwriting the data.
 
 
      real*8 min_sc, max_sc, omin_sc, omax_sc, ds, jump,
     .    clk, lab_offset
 
*   root    - Root name for output file
*   outfile - Name of output file
 
 
      character*128 outfile

***** Loop over all the sites creating and writeing the output
*     file and generating the cplot command file.
      
      outfile =  root(1:trimlen(root)) // '.plt'
      open(200, file= outfile , iostat=ierr, status = 'unknown' )
c      open(200, file=root(1:trimlen(root)) // '.plt' , iostat=ierr )
      plot_on_page = 0
 
      do i = 1, num_param
 
          if( i.le.num_cfiles ) then
              outfile = root(1:trimlen(root)) // '.' // cf_codes(i)
              open(201,file=outfile, iostat=ierr, status='unknown')
              write(201,100,iostat=ierr) cf_codes(i),
     .                 rcv_types(i), init_clk_val(i)/fClk*1.d6
 100          format('* Clock information for site ',a,
     .               ' receiver ',a4,'. Initial offset ',F10.3,
     .               ' usec', /,
     .               '* Epoch   clk (us)  Flag')
          else
              write(outfile,110) root(1:trimlen(root)),
     .              prn_list(i-num_cfiles)
 110          format(a,'.PRN_',i2.2)
              open(201,file=outfile, iostat=ierr, status='unknown')
              write(201,120,iostat=ierr) prn_list(i-num_cfiles),
     .                init_clk_val(i)/fClk*1.d6
 120          format('* Clock information for Satellite PRN ',i2,
     .               '.  Initial offset ',F10.3,/,
     .               '* Epoch   clk (us)  Flag')
          end if
 
*****     Set the jump value to remove the initial clock offset
          jump = init_clk_val(i)/fClk*1.d6
 
*         write header information into data file
 
          min_sc = 100.d9
          max_sc = -100.d9
          num_out = 0
          do j = 1, num_ep
             if( .not.kbit(par_flag_cse(i,j),1) ) then
 
*                OK, estimate is good, see if there is a jump
                 if( kbit(par_flag_cse(i,j),2) .and. j.gt.1) then
 
*                     There has been a jump so set jump for
*                     continuity
                      jump = (params_cse(i,j)
     .                      - params_cse(i,j-1))/fClk*1.d6 + jump
                 end if
 
                 clk = params_cse(i,j)/fClk*1.d6 - jump
 
                 write(201,300) j, clk, jump, par_flag_cse(i,j)
 300             format(I5,1x,2F15.6,1x,o8)
                 min_sc = min(min_sc, clk )
                 max_sc = max(max_sc, clk )
                 num_out = num_out + 1
             end if
          end do
*         If we output data, make entry in cplot file
          if( num_out.gt.0 ) then
 
              plot_on_page = plot_on_page + 1
 
*             Clip to nearest microsecond. Save original scale
*             so we know if plot trends down or up
              omax_sc = max_sc
              omin_sc = min_sc
              max_sc = int(max_sc) + 1
              min_sc = int(min_sc) - 1
              ds = max_sc - min_sc
 
*             Set whether the labels are at the top or bottom
*             of the window.  (Since the plot starts at zero
*             we see if the averge of scales is positive or
*             negative.
              if( (omax_sc + omin_sc )/2.lt.0 ) then
                   lab_offset = 0
              else
                   lab_offset = (max_sc-min_sc) - 5*(ds/20)
              end if
 
*             Now do commands
              write(200,400) 0.05+(4-plot_on_page)*.225,
     .                       0.05+(5-plot_on_page)*.225,
     .                       outfile, num_ep+1, min_sc, max_sc,
     .                       lab_offset+(ds/20)*4,
     .                       lab_offset+(ds/20)*3,
     .                       lab_offset+(ds/20)*2,
     .                       lab_offset+(ds/20)*1
 400          format('  head 2 1 ',/,' font 5x7',/,
     .              '  view 0.1 0.95 ',2f8.4,/,
     .              '  file ',a,/,
     .              '  x_field 1 1 0',/,
     .              '  y_field 1 2 0 "Clock (usec)"',/,
     .              '  read ',/,
     .              '  x_scale  0 ',i6,/,
     .              '  y_scale  ',f10.3,1x,f10.3,/,
     .              '  point 1',/,
     .              '  char 1.8',/,
     .              '  line 0',/,
     .              '  draw ',/,
     .              '  xmn -1 0 ',/,
     .              '  xmx -1 0 ',/,
     .              '  ymx -1 0 ',/,
     .              '  ymn -1 1 "Clock (usec)"',/,
     .              '  poly ep usec 1 1', /,
     .              '  fit 0',/,
     .              '  pdr ',/,
     .              '  label 10 ',f10.3,' 1 0  :p1',/,
     .              '  label 10 ',f10.3,' 1 0  :p2',/,
     .              '  label 10 ',f10.3,' 1 0  :p3',/,
     .              '  label 10 ',f10.3,' 1 0  :h1',/ )
 
              if( plot_on_page.eq.1 ) then
                  write(200,'('' xmx -1 1 "Epoch number"'')')
              end if
              if( plot_on_page.eq.4 ) then
                  plot_on_page = 0
                  write(200,'('' ident'',/)')
                  write(200,'('' erase'',/)')
              end if
 
          end if
          close(201)
      end do
 
      close(200)
 
***** Thats all
      return
      end
 
CTITLE INIT_RMS
 
      subroutine init_rms( sum_var, sum_num, num )

      implicit none
 
*     Routine to clear the summation variables for accumulating
*     the statitics on the range residuals
 
* PASSED VARIABLES
 
* num  - Number of values to initialize
* sum_num(num) - Number of values in each sum
 
 
      integer*4 num, sum_num(num)
 
* sum_var(num) - Sum of residuals squared
 
 
      real*8 sum_var(num)
 
* LOCAL
 
* i    - Loop counter
 
      integer*4 i
 
***** Loop over values clearing them
      do i = 1, num
         sum_num(i) = 0
         sum_var(i) = 0.d0
      end do
 
*     Thats all
      return
      end
 
CTITLE UPDATE_STATS
 
      subroutine update_stats( params_cse, par_flag_cse, converged )

      implicit none
 
*     This routine will finish up the range noise rms calculation
*     and scan the parameter list to compute the rms of the clocks
*     at each site.
 
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   par_flag_cse(num_param,num_ep)         - parameter estimate flag
 
 
      integer*4 par_flag_cse(num_param,num_ep)
 
*   converged  - Logical to indicate that clock variance estimates
*                have converged
 
 
      logical converged
 
*   params_cse(num_param)    - Clock parameter estimates (L1 cycles)
 
 
      real*8 params_cse(num_param,num_ep)
 
* LOCAL VARIABLES
 
*   i,j,  - Loop counters
*   sum_clk_num - Number of values summed for each clock
*   start_ep, end_ep, step_ep  - The first epoch with valid data, the
*            end of the epoch range to scan, and the step in epochs for
*            receivers whose natural data spacing is different to the
*            that used in these cfiles.
 
 
      integer*4 i,j, sum_clk_num, start_ep, end_ep, step_ep
 
*   dclk  - Change in clock between epochs (if there is no jump)
*   sum_clk_var - Summation of clk differences squared
*   clk_var  - Estimate of the clock variance
*   clk_sd   - Estimate of the clock standard allan sd.  Used for
*              checking converge.
 
 
      real*8 dclk, sum_clk_var, clk_var, clk_sd
 
*   kbit  - Logical function to test if bit is on
*   start_found  -  Indicates that the first epoch of data has been found
*          for those sites which have a sampling interval greater than the
*          standard one.
 
 
      logical kbit, start_found
 
***** Finish up the range noise rms calculations
      do i = 1, num_cfiles
         if( sum_rng_num(i).gt.0 ) then
             rng_noise(i) = sqrt(sum_rng_var(i)/sum_rng_num(i))
 
*            Make sure that the range noise does not go to
*            zero.
             if( rng_noise(i).lt.1.d0 ) then
                 rng_noise(i) = 1.d0
             end if
         end if
      end do
 
***** Start looping over all the clock parameters
      converged = .true.
 
      do i = 1, num_param
          sum_clk_num = 0
          sum_clk_var = 0.d0
 
*         work out where we need to start if the orig_sampling of
*         data is greater than the current sampling, i.e., there will
*         gaps in the data which we need to account for.  We only
*         check this for stations.  (The satellites should have
*         results at the sampling interval of the data)
          start_ep = 2
          end_ep   = num_ep
          step_ep  = 1
 
*         Check the recievers to see if we need update the step
          if( i.le.num_cfiles ) then
              if( orig_sampling(i).gt.sampling ) then
 
*                 Scan up to the first existant observation to get the
*                 start epoch and set the step size
                  j = 0
                  start_found = .false.
                  do while ( j.lt.num_ep .and. .not.start_found )
                      j = j + 1
                      if( .not.kbit(par_flag_cse(i,j),1) ) then
*                         found the a good epoch of data so start here
                          step_ep  = orig_sampling(i)/sampling
                          start_ep = j+ step_ep
                          end_ep   = num_ep
                          start_found = .true.
                      end if
                  end do
              end if
          end if
 
          do j = start_ep, end_ep, step_ep
 
*             See if this estimate has no jump and that the
*             previous estimate is OK
              if( .not.kbit(par_flag_cse(i,j-step_ep),1) .and.
     .                 par_flag_cse(i,j).eq.0  ) then
                  dclk = params_cse(i,j) - params_cse(i,j-step_ep)
                  sum_clk_var = sum_clk_var + dclk**2
                  sum_clk_num = sum_clk_num + 1
C                 write(*,8000) i,j, dclk/fClk*1.d6,
C    .                 sqrt(sum_clk_var/max(1,sum_clk_num))/fClk*1d6
C8000             format(2i5, 2F15.3)
              end if
          end do
 
*****     Now update the Allan std deviation of this clock
          if( sum_clk_num.gt.0 ) then
*             halve the value because of the diffencing
              clk_var = sum_clk_var/sum_clk_num
              clk_sd  = sqrt(clk_var/(100.d0*sampling*step_ep))/fClk
 
****          Make sure that the value does not get too small;
*             limit the value to 0.01d-9.
              if( clk_sd.lt.0.01d-9 .and. i.le.num_cfiles ) then
                  write(*,100) cf_codes(i), clk_sd
 100              format(' Allan Std.Dev. at ',a,' below threshold.',
     .                   ' Value is ',D12.3,'. Limiting to 0.01d-9')
                  clk_sd = 0.01d-9
                  clk_var = (clk_sd*fClk)**2*(100.d0*sampling)
              end if
 
*****         Make sure value is not too large.  If too big then there
*             can be problems detecting milliseconf time jumps.
              if( clk_sd.gt. 100.d-9 .and. i.le.num_cfiles ) then
                  write(*,120) cf_codes(i), clk_sd
 120              format(' Allan Std.Dev. at ',a,' above threshold.',
     .                   ' Value is ',D12.3,'. Limiting to 100.d-9')
                  clk_sd = 100.d-9
                  clk_var = (clk_sd*fClk)**2*(100.d0*sampling)
              end if

*****         Limit the clock allan standard on the clocks
              if( i.gt. num_cfiles ) then
                  if( clk_sd.gt.1.d-9 ) then
                      write(*,140) prn_list(i-num_cfiles), clk_sd
 140                  format('  Allan Std.Dev. for PRN ',i2,
     .                       ' above threshold.',
     .                       ' Value is ',D12.3,'. Limiting to 1.d-9')
                      clk_sd = 1.d-9
                      clk_var = (clk_sd*fClk)**2*(100.d0*sampling)
                  end if
                  if( clk_sd.lt.0.01d-9 ) then
                      write(*,150) prn_list(i-num_cfiles), clk_sd
 150                  format('  Allan Std.Dev. for PRN ',i2,
     .                       ' below threshold.',
     .                       ' Value is ',D12.3,'. Limiting to 0.01d-9')
                      clk_sd = 0.01d-9
                      clk_var = (clk_sd*fClk)**2*(100.d0*sampling)
                  end if
              end if

*             If we have a prior estimate of the clock variance
*             check to see if have converged.
              if( rc_allan_sd(i).ne.0.d0 ) then
                  if( abs((rc_allan_sd(i)-clk_sd)/rc_allan_sd(i)).gt.
     .                0.2d0 ) converged = .false.
              end if
 
              rc_allan_sd(i) = sqrt(clk_var/(100.d0*sampling))/fClk
              rc_allan_num(i) = sum_clk_num
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE REPORT_CLK_RNG
 
      subroutine report_clk_rng ( iter, converged, uno )

      implicit none
 
*     This routine will report the quality of the clocks at all
*     all the stations and satellites, and the range rms quality
*     at the ground stations.
 
 
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
*  iter - Current iteration counter.
*  uno     - Unit number for output.  If this is the last iteration
*            Then the summary unit is written to, else unit 6 (stdout)
*            is used.
 
      integer*4 iter, uno

* converged - True if clock iteration has converged.

      logical converged 
 
* LOCAL VARIABLES
 
*  i  - loop counter
 
 
      integer*4 i

      character*12 hsver
 
***** Loop over all the parameter values reporting the clock
*     allan std deviation

      if( uno.eq.uns ) then
          write(uno,100) caprog_name, hsver(ctogobs_version)
 100      format(/,1x,a8,'SUMMARY FILE: Version ',a)
      end if

* MOD TAH 200628: Report GNSS used.
      write(uno,120) iter, sv_gnss
 120  format(/' Clock and Range noise statistics at iteration ',i2,
     .        ' SYS ',a1,/,
     .   ' Site/PRN    Allan SD@100  #     Range rms    #',/,
     .   '             sec  (ppb)            (mm)' )
 
      do i = 1, num_cfiles
           write(uno,200) cf_codes(i), rc_allan_sd(i)*1.d9,
     .                  rc_allan_num(i),
     .                  rng_noise(i)/fClk*vel_light*1000.d0,
     .                  sum_rng_num(i), rcv_types(i)
 200       format(1x,a4,4x,F12.6,1x,i6,3x,f12.1,1x,i6,1x,a4)
      end do
 
      do i = num_cfiles+1, num_param
           write(uno,220) prn_list(i-num_cfiles), rc_allan_sd(i)*1.d9,
     .                  rc_allan_num(i)
 220       format(1x,'PRN_',i2.2,2x,F12.6,1x,i6)
      end do
 
***** Thats all
      return
      end
 
CTITLE GET_RNG_RESETS
 
       subroutine get_rng_resets(L1r_rng_cse, L2r_rng_cse,
     .        data_flag_cse, params_cse, par_flag_cse,
     .        ctol_cse )

      implicit none
 
*     This routine will scan the clock adjusts at each site and if
*     there is a millisecond jump (with tolerance) it will adjust the
*     ranges to remove the jump.  (These will be recorded as "cycle
*     slips" in the range data when the gobs file is written)
 
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.    
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep), 
     .    ctol_cse(num_chan,num_cfiles,num_ep)
 
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*       i,j,k       - Loop counters
 
      integer*4 i,j,k
 
*       kbit        - Checks if bit is set
 
      logical kbit
 
*   last_clk    - Value of last clock offset.  Used to see if there
*               - is a millisecond clock jump (L1 cycles)
*   millisec_jumps  - Value which represents the number of millisecond
*                 jumps that we need to add. (L1 cycles)
*   dclk    - Change in clock when there is a jump (L1 cycles)
*   ms_dclk - dclk expressed in multiples of 1 millisecond
*           - jumps (L1 cycles)
 
      real*8 last_clk, dclk, ms_dclk, millisec_jumps
 
****  Loop over all the receivers
 
      do i = 1, num_cfiles
          last_clk = 0.d0
          millisec_jumps = 0.d0
          do j = 1, num_ep
 
*             See if jump in clock
              if( kbit(par_flag_cse(i,j),2) ) then
 
*                 See if the jump matches an whole number of millisecs.
                  dclk = params_cse(i,j) - last_clk        
                  ms_dclk = nint(dclk/(fClk*1.d-3)) * (fClk*1.d-3)
                  if( abs(dclk-ms_dclk).lt.reset_tol ) then
 
*                     We have a clock reset.  Set the last clk to
*                     to the even millisec number and adjust the ranges
                      millisec_jumps = millisec_jumps + ms_dclk
                  end if
              end if

*             Save the value of the last clk estimate
              last_clk = params_cse(i,j)
 
****          Now adjust the ranges by last_clk value if non zero.
              if( millisec_jumps.ne.0.d0 ) then
                  do k = 1, actual_max_chan
 
*                     See if we have data in this channel
                      if( .not.kbit(data_flag_cse(k,i,j),30) ) then
                          L1r_rng_cse(k,i,j) = L1r_rng_cse(k,i,j) -
     .                                        millisec_jumps
 
*                         if we have L2 data then update as well.
*                         (NOTE: Zero L2 range is used to indicate no
*                                 L2 range data, so do not update unless
*                                 we actually have range data)
                          if( L2r_rng_cse(k,i,j).ne.0.d0 ) then
                              L2r_rng_cse(k,i,j) = L2r_rng_cse(k,i,j)
     .                         - millisec_jumps*
     .                         fL2(ctol_cse(k,i,j))/fL1(ctol_cse(k,i,j))
                          end if
                      end if
                  end do
              end if
*                             ! Looping over epochs
          end do
*                             ! Looping over sites
      end do
 
***** Thats all
      return
      end
 
CTITLE LTOC

      integer*4 function ltoc( ctol, list, num_chan)

      implicit none

*     Function to return the channel number corresponding to the  
*     list satellite (list is the numbe in the list of PRN's)

* PASSED VARIABLES

* num_chan  - Number of channels used in this experiment
* ctol(num_chan)  - List numbers in each of the channels
* list      - number of satellite that we are looking for
*             (returns -1 if satellite not observed)

       integer*4  num_chan, ctol(num_chan), list

* LOCAL VARIABLES

* i         - Loop counter

       integer*4 i

*****  set the channel number to -1 and then search

       ltoc = -1
       do i = 1, num_chan
          if( list.eq.ctol(i) ) then
              ltoc = i
              RETURN
          end if
       end do

*****  Thats all
       return
       end


 
 
