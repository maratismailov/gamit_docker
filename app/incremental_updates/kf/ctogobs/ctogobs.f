 
      program ctogobs

      implicit none
 
*     This is the pre-processor program for the SOLVG GPS Kalman
*     filter program.  Its main task is merge a set of user specified
*     c-files and generate the Gobs_file for SOLVG.  The runstring for
*     the program is:
*     % ctogobs <command file> <Gobs_file> <descp> <new cfile series> 
*                    [....list of cfiles...]  or [dfile] [use series]
*     or
*     % autcln  <command file> <new cfile series> 
*                    [....list of cfiles...] or [dfile] [use series]
*
*     where <command file> is an optional command file (commands are
*             given the ctogobs.hlp file.   If no command file is
*             given (default file generation) then '' should be used
*             as a place holder in the runstring
*         <Gobs_file> is the optional name for the output (binary)
*             Gobs_file.  The options here are
*              (1) If no name is given (just '' in runstring)
*                 then the name will be generated as gYYMMDDS.gobs
*                 where yymmdd is the year, month, day of the data from
*                 c-files and S is the series character from the cfile.
*             (2) A single character is given, and it will used as the
*                 series designator
*             (3) A full arbitary file name (along with directory
*                 information is given.
*       <descp> is the experiment description enclosed in single
*               quotes e.g. 'TREX00 Day 1'
*       <new cfile series> new series for cfiles (if not given then no
*               new cfiles are written). '+' may be used to increment
*               series.
*         [....list of cfiles...] is the list of cfiles to be merged.
*       or [dfile] is the name of a dfile with the option of an additional
*               element which gives the correct series.
*
* UNIT NUMBERS:
*  100    - command file and dfile (in get_ctog_runstring).
*  101-100+num_cfiles - open cfiles.
*  200    - Plot files for range cloks, phase clocks and residuals
*  201    - Data file for the abobe types of plots
*  202    - dd.srt file
*  203    - Summary file
*  204    - Used to read mfile (only open for short time).
*  400+prn - Single difference files

*
* INCLUDE FILES
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/gobs_def.h'
      include 'ctogobs_com.h' 
c      external ctog_cmds_bd 
 
* MAIN PROGRAM VARIABLES
 
 
* Large size dynamically mapped variables
 
* vma_data(1)   - This the array that will ultimately
*               - contain all the data.  It is dynamically
*               - mapped using malloc.
* vma_real(1)   - Double precision vma_data.  Equivalanced to vma_data
*               - to ensure that the real*8 values fall on 32 bit
*               - boundaries
 
      integer*4 vma_data(1), ierr
 
      real*8 vma_real(1)
 
      equivalence ( vma_data, vma_real)

* LOCAL VARIABLES

* trimlen  - returns the length of the used portion of a string
* date(5)  - Run-time date (used for timing)
* sectag   - Current second (used for timing)
* cpass    - Cleaning pass number (counter for 4 cleaning iterations.)


      integer*4 trimlen, date(5), cpass
      real*8 sectag

* outline - Line to output to report status

      character*128 outline

* scan_edits  - Logical to indicate that there were scan edits and
*     we need to keep iterating
* converged   - True when phase clock algorithm converged.
* fatal_exist - Set true if fatal file exists.
* cref_upd    - Set true if the reference clock list is updated
*               results in additional iteration.

      logical scan_edits, converged, fatal_exist, cref_upd
      
***** Initial the program and print out a header message

      call init_ctogobs
 
***** Read the ctogobs runstring:
 
      call get_ctogobs_run

****  Check to see if GAMIT.fatal exists
      inquire( file = 'GAMIT.fatal', exist = fatal_exist )
      if( fatal_exist ) then
          call report_stat('FATAL','AUTCLN','autcln',' '
     .                  ,'GAMIT.fatal exists: AUTCLN not executed',0)
      end if
 
***** Now loop over all of the cfiles and read the headers
*     (We will use this to set the memory size needed by the
*      program).
 
      call read_cfheads( ierr )
      if( ierr.ne.0 ) then
          call report_stat('fatal','autcln','main',' ',
     .           'Error reading one or more cfiles', ierr)
      end if

***** Set the frequency quantities based on fL1 and fL2 from the c-files
      call set_freqs 

*     Set the clock statistics by receiver type before reading
*     commands (maybe updated)
 
      call init_receiver 
 
*     Now read the ctogobs command file (if given).
 
      if( use_command_file ) then
          call read_ctog_cmds
      end if

C MOD TAH 180331: Moved set of frequencies until after command file
C     is read due to remap_glonass commmand affects frequencies
***** Set the frequency quantities based on fL1 and fL2 from the c-files
C      call set_freqs 

*     Open up the summary file:
      if( trimlen(summary_file).gt.0 ) then
          uns = 203
          call open_lu(uns, summary_file, ierr, 'unknown')
          call report_error('IOSTAT',ierr,'open',summary_file,
     .                      0,'CtoGobs/Main')
          if( ierr.ne.0 ) uns = 6
      end if
 
****  Get the memory allocation needed to process data
      call ctog_mem(vma_data(1))

****  Now read all of the data from the c-files
      call systime( date,sectag)
      write(*,100) '+Start cfile read    ', date, sectag
 100  format(a,i5,4i3,1x,f6.2,:,' Pass ',i2)
      call report_stat('STATUS','autcln','main',' ',
     .                 'Start: Reading cfiles',0)

          
      call read_cfdata(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .    vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .    vma_data(iL1_cyc_cse ), vma_data(iL2_cyc_cse ),
     .    vma_data(ictol_cse ),   vma_data(idata_flag_cse),
     .    vma_data(ibf_type_cse), vma_data(iparams_cse), 
     .    vma_data(ipar_flag_cse), 
     .    vma_data(iazel_cse),    vma_data(ipf_dphs_cse), 
     .    vma_data(isvcL1_ce),     ierr )
      if( ierr.ne.0 ) then
          call report_stat('fatal','autcln','main',' ',
     .                     'Error reading cfiles',ierr)
          stop 'AUTCLN: Error reading some cfiles'
      end if

***** Now process the data to get epoch by epoch clocks
*     and close estimates of the cycle slips

*     Now do from here to scanning loop iteratively incase
*     we need to discard a site/sv combination. cpass count
*     is to force the loop exit (hopefully will not be needed).
      cpass = 0
      scan_edits = .true.
      cref_upd   = .true.

      call init_scan_edits( num_cfiles, num_sat, num_dd_flags,
     .                      max_cfiles )

***   For GLONASS re-fit the clocks 
      if( cf_gnss.eq.'R' ) then
          call fit_glonass_clks( 
     .           vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .           vma_data(ictol_cse ), vma_data(idata_flag_cse),
     .           vma_data(iparams_cse), vma_data(ipar_flag_cse) )
      endif

* MOD TAH 200509: Added option to scan for and remove milli-
*     second clock jumps.  This option may cause problems with
*     time tags if cleaned c-files are re-modeled with model
*     because on pseudo-range residuals are updated not the 
*     observed pseudo-ranges themselves.
      if ( prescan_ms ) then
         call prescan_msj(vma_data(iL1r_phs_cse), 
     .      vma_data(iL2r_phs_cse),
     .      vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .      vma_data(iL1_cyc_cse ), vma_data(iL2_cyc_cse ),
     .      vma_data(ictol_cse ), vma_data(idata_flag_cse),
     .      vma_data(iparams_cse), vma_data(ipar_flag_cse) )
      endif

*     Make mutiple passes through this loop so that phase clocks
*     computed with all double difference slips detected.
      do while ( ((scan_edits .or. cref_upd) .and. cpass.lt.5)
     .            .or. cpass.lt.3 )

          cpass = cpass + 1

* MOD TAH 160826: See if we will prefit clock polynomials.  We
*         limit this to first 2 passes (to catch breaks and offsets)
          if( cpass.le.5 .and. prefit_clk ) then
*             Determine new clock polynomials fits based on R range
*             data
              call systime( date,sectag)
              write(*,100) '+Start prefit clocks ', date, sectag, cpass
              write(outline, 105) cpass
 105          format('Prefit clocks from range data. Pass ',i2)
              call report_stat('status','autcln','main',' ',
     .                     outline, 0)

*             Now make linear fit to clocks  
              call fit_igs_clk( 6 , vma_data(iparams_cse), 
     .             vma_data(ipar_flag_cse),
     .             vma_data(iL1_cyc_cse), 
     .             vma_data(iL2_cyc_cse), 
     .             vma_data(ictol_cse), 'R')

*             The estimates are returned in save_clk(2,..) and we
*             apply these to the apr_clk_poly values that come from
*             the c-files and apply the corrections to the range and
*             phase data (so it looks like the apr_clk_poly values
*             were used in model).
              call app_prefit_clk( vma_data(iL1r_phs_cse), 
     .           vma_data(iL2r_phs_cse),
     .           vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .           vma_data(iL1_cyc_cse ), vma_data(iL2_cyc_cse ),
     .           vma_data(ictol_cse ), vma_data(idata_flag_cse),
     .           vma_data(iparams_cse), vma_data(ipar_flag_cse) )
         endif 
 
          call systime( date,sectag)
          write(*,100) '+Start Range clocks  ', date, sectag, cpass
          write(outline, 110) cpass
 110      format('Estimating clocks from range data. Pass ',i2)
          call report_stat('status','autcln','main',' ',
     .                     outline, 0)

          call process_rng(vma_data(iL1r_phs_cse), 
     .        vma_data(iL2r_phs_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(iL1_cyc_cse ), vma_data(iL2_cyc_cse ),
     .        vma_data(ictol_cse ), vma_data(idata_flag_cse),
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse) )

*         Routine to write out the Range estimates of the clocks
*         in IGS clock format.  Currently not activated.
          if( write_igs_clk .and. cpass.ge.2  ) then
              call fit_igs_clk( 6 , vma_data(iparams_cse), 
     .                        vma_data(ipar_flag_cse),
     .                        vma_data(iL1_cyc_cse), 
     .                        vma_data(iL2_cyc_cse), 
     .                        vma_data(ictol_cse), 'R')
              call wr_igs_clk( 6, vma_data(iparams_cse), 
     .                         vma_data(ipar_flag_cse),
     .                         vma_data(isvcL1_ce),
     .                         vma_data(idata_flag_cse),
     .                         vma_data(ictol_cse ), cref_upd, 'R' )
          end if
 
*         See if we are to flag gaps with bias flags
          if( gaps_flagged ) then
              call set_phs_mask( phs_mask, phs_bias_mask ) 
              call flag_gaps(cpass, vma_data(ictol_cse), 
     .                       vma_data(idata_flag_cse),
     .                       vma_data(ibf_type_cse) ) 
* MOD TAH 200614: Added call to trim_shortseg to remove short blocks
*             of data that are impossible to patch
              if( cpass.eq.1 .and. trim_seg )
     .        call trim_shortseg(vma_data(ictol_cse), 
     .                          vma_data(idata_flag_cse),'Y',6)
          end if

*****     If we are on second iteration, then add the prescan
*         phase routine to align the range and phase clocks.
          if( cpass.ge.2 ) then
              call systime( date,sectag)
              write(*,100) '+Start align phase', date, sectag, cpass
              write(outline, 115) cpass
 115          format('Prealigning phase data. Pass ',i2)
              call report_stat('status','autcln','main',' ',
     .                         outline, 0) 

              call align_phs(cpass, vma_data(iL1r_phs_cse), 
     .            vma_data(iL2r_phs_cse),
     .            vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .            vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .            vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .            vma_data(ibf_type_cse), 
     .            vma_data(iparams_cse),  vma_data(ipar_flag_cse)) 

          end if
              
*         Now determine all of the cycles values which need to be added to
*         the one-way data.  (We also recompute the clocks at this time
*         based on the phase data.)
          call systime( date,sectag)
          write(*,100) '+Start phase clocks  ', date, sectag, cpass
          write(outline, 120) cpass
 120      format('Estimating clocks from phase data. Pass ',i2)
          call report_stat('status','autcln','main',' ',
     .                     outline, 0)
          call process_phs(vma_data(iL1r_phs_cse), 
     .            vma_data(iL2r_phs_cse),
     .            vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .            vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .            vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .            vma_data(ibf_type_cse), 
     .            vma_data(iparams_cse),  vma_data(ipar_flag_cse),
     .            vma_data(iazel_cse)  ,  cpass  )

*         Scan the double differences to make sure we have no unflagged
*         slips.
          call systime( date,sectag)
          write(*,100) '+Start Scan 1        ', date, sectag, cpass
          write(outline, 130) cpass
 130      format('Scanning Double difference for slips. Pass ',i2)
          call report_stat('status','autcln','main',' ',
     .                     outline, 0)
* MOD TAH 990519: Changed to pass in the iteration.  After first
*         iteration we ignore gaps in data.
          call scan_dd(cpass,vma_data(iL1r_phs_cse), 
     .            vma_data(iL2r_phs_cse),
     .            vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .            vma_data(ictol_cse), vma_data(idata_flag_cse),
     .            vma_data(ibf_type_cse) )

*         Now check to see if we have exceeded the maxiumum number
*         of dd_scan biased flags added.  If we have then edit the
*         offending site/sv combinations, reset the number of
*         cycles and try again.
          call chk_scan_edit( cpass, scan_edits, 
     .            vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .            vma_data(ictol_cse), vma_data(idata_flag_cse),
     .            vma_data(ibf_type_cse) )

      end do

      call report_clk_rng ( cpass, converged, uns )
      call rep_num_dd_flag( uns, cpass )

*     Now clean the double differences
      call systime( date,sectag)
      write(*,100) '+Start DD clean 1    ', date, sectag 
      call report_stat('status','autcln','main',' ',
     .                 'Cleaning data. First pass', 0)
      call clean_dd(1, vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(ibf_type_cse), 
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse)  )

* MOD TAH 940506: Iterate the cleaning loop three more times
      cpass = 1
      do while ( cpass.lt.3 )
         cpass = cpass + 1
         call systime( date,sectag)
         write(*,100) '+Start DD clean iter ', date, sectag, cpass
         write(outline, 150) cpass
 150     format('Cleaning data. Iteration ',i2)
         call report_stat('status','autcln','main',' ',outline,0)
         call clean_dd(cpass, vma_data(iL1r_phs_cse), 
     .        vma_data(iL2r_phs_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(ibf_type_cse), 
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse)  )

*        Get the wide-lane site dependent bias.
         call get_mean_owwl(vma_data(iL1r_phs_cse), 
     .        vma_data(iL2r_phs_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .        vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .        vma_data(ibf_type_cse), 
     .        vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )

      end do
      
* MOD TAH 990523: Call routine to flatten the double differences
      call systime( date,sectag)
      write(*,100) '+Start Flat DD  ', date, sectag 
      call report_stat('status','autcln','main',' ',
     .                     'Start Flat DD',0)

      call flat_dd(1, vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(ibf_type_cse),
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse)  )


* MOD TAH 970109: Do final fit to the phase clock values.  Only
*     needed if we are computing one-way phase residuals.
      
      if ( apply_phs_clk ) then

          call systime( date,sectag)
          write(*,100) '+Start Final phase clocks  ', date, sectag
          call report_stat('status','autcln','main',' ',
     .                     'Final Phase clock fit',0)
*         MOD TAH 970112: Iterate the fitting of the one-way phase day
*         for clock estimates.
* MOD TAH 990524: Start pass at 1 so that clocks are directly
*         estimated
*          cpass = 0
          cpass = 1
          converged = .false.

          do while ( cpass.le.pc_max_iter-1 .and. .not.converged )
             cpass = cpass + 1
             call systime( date,sectag)
             write(*,100) '+Start OW clean iter ', date, sectag, cpass
             write(outline, 220) cpass
 220         format('+Phase clock and bias estimation pass ',i3)
             call report_stat('status','autcln','main',' ',outline,0)

*            Now try to remove offsets in the one-way data.  This
*            call only needed at start, after that we apply the 
*            mean LC residuals to the data.
C            if( int(cpass/2)*2-cpass.eq.0 .and. cpass.lt.3 ) then
             if( int(cpass/2)*2-cpass.eq.0 .and. cpass.lt.0 ) then
                 write(*,240) cpass
 240             format('+CLEAN Oneway absolute pass ',i2)
                 call clean_ow(cpass, vma_data(iL1r_phs_cse), 
     .                vma_data(iL2r_phs_cse),
     .                vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .                vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .                vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .                vma_data(ibf_type_cse), 
     .                vma_data(iparams_cse), vma_data(ipar_flag_cse)  )

             end if

*            Now estimate phase clock values.  On even iterations no
*            biases are estimated.  On Odd passes, the estimates after
*            a bias flag are not used until we have sufficient data to
*            make an estimate of the bias.

*            On the odd passes (1st/3rd) estimate the biases explicitly
*  MOD TAH 011229: Only allow the full PC solution is there are not
*            more than the max allowed channels.
             if( num_chan.gt.max_gchannels ) pc_full_anal = 0
             if( int(cpass/2)*2-cpass.ne. 0 .and. 
     .               cpass.le.pc_full_anal ) then
                 call est_bps(cpass, vma_data(iL1r_phs_cse), 
     .                vma_data(iL2r_phs_cse),
     .                vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .                vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .                vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .                vma_data(ibf_type_cse), 
     .                vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )
             else

*                Estimate the clock terms.
                 call proc_phsfin(cpass, vma_data(iL1r_phs_cse), 
     .                vma_data(iL2r_phs_cse),
     .                vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .                vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .                vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .                vma_data(ibf_type_cse), 
     .                vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )

*                Now remove the mean offsets in each segment of data.
                 call rm_ow_mean(cpass, vma_data(iL1r_phs_cse), 
     .                vma_data(iL2r_phs_cse),
     .                vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .                vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .                vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )
             end if

*            Report the RMS scatter if an even pass number.
             if( int(cpass/2)*2-cpass.eq. 0 ) then

                 write(*,260) cpass
 260             format('+Postfit RMS of residuals after pass ',i3) 
                 call comp_pf_rms(6, cpass, vma_data(iL1r_phs_cse), 
     .              vma_data(iL2r_phs_cse),
     .              vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .              vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .              vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )
                 call check_convd(cpass, converged )

*                Only on the even passes do we check the editing.
                 if( edit_postfit .and. cpass.ge.pc_start_edit ) then

* MOD TAH 051227: Added call to trim_oneways
                     call trim_oneways(vma_data(ictol_cse), 
     .                             vma_data(idata_flag_cse),'N')

                     call comp_elev_rms(uns, cpass, 
     .                 vma_data(iL1r_phs_cse), 
     .                 vma_data(iL2r_phs_cse),
     .                 vma_data(iL1_cyc_cse),   vma_data(iL2_cyc_cse),
     .                 vma_data(ictol_cse),  vma_data(idata_flag_cse),
     .                 vma_data(iparams_cse), vma_data(ipar_flag_cse),
     .                 vma_data(iazel_cse), 'NOR' )
                     call pf_edit(cpass,
     .                  vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse), 
     .                  vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .                  vma_data(ictol_cse), vma_data(idata_flag_cse),
     .                  vma_data(iparams_cse), vma_data(ipar_flag_cse),
     .                  vma_data(iazel_cse) )
                 end if
             end if


          end do

*****     Save number of iterations needed
          pc_num_iter = cpass

*         Report the rms to the summary file
* MOD TAH 970812: Report RMS values if more than 1 iteration.
          if( uns.ne.6 .and. cpass.ge.2 ) 
     .                       call report_rms(uns, pc_num_iter)

*         Compute an elevation fit to the phase residuals
          if( cpass.gt.2 )             
     .    call comp_elev_rms(uns, cpass, vma_data(iL1r_phs_cse), 
     .              vma_data(iL2r_phs_cse),
     .              vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .              vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .              vma_data(iparams_cse),  vma_data(ipar_flag_cse),
     .              vma_data(iazel_cse), 'REP' )
     
*         Write postfit resiudals.        
          if( trimlen(phs_res_root).gt.0 ) then
              call write_phs(phs_res_root, 
     .                 vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .                 vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .                 vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .                 vma_data(ictol_cse), vma_data(idata_flag_cse),
     .                 vma_data(iparams_cse), vma_data(ipar_flag_cse),
     .                 vma_data(iazel_cse)   )
          end if
        
****      Now we need to reset the cycles back to integer values.
          call reset_cyc( vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse))       

*         Try to remove bias flags in the one ways.
          if( pf_remove_bf ) then
             cpass = 0 
             call systime( date,sectag)
             call report_stat('status','autcln','main',' ',
     .                 'One-way bias flag removal',0)
             write(*,100) '+Start Postfit remove BF ', date, sectag
             call pf_ow_sf(cpass, vma_data(iL1r_phs_cse), 
     .            vma_data(iL2r_phs_cse),
     .            vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .            vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .            vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .            vma_data(ibf_type_cse), 
     .            vma_data(iparams_cse), vma_data(ipar_flag_cse)  )
          end if

****      If user asked, now output the phase clock estimates or the
*         phase postfit resiudals.
*         Write phase clock estimates.
          if( trimlen(phs_clk_root).gt.0 ) then
              call write_clk(phs_clk_root, 
     .                 vma_data(iparams_cse), vma_data(ipar_flag_cse) )
          end if
        
****      See if we going write out an igs_clock file
          if( write_igs_clk ) then
              call fit_igs_clk( uns , vma_data(iparams_cse), 
     .                        vma_data(ipar_flag_cse),
     .                        vma_data(iL1_cyc_cse), 
     .                        vma_data(iL2_cyc_cse), 
     .                        vma_data(ictol_cse), 'P')
              call wr_igs_clk( uns, vma_data(iparams_cse), 
     .                         vma_data(ipar_flag_cse),
     .                         vma_data(isvcL1_ce),
     .                         vma_data(idata_flag_cse),
     .                         vma_data(ictol_cse ), cref_upd, 'P' )
              if( cref_upd ) then
                 print *,'CLOCKS: Updating Phase Reference list'
                 call fit_igs_clk( uns , vma_data(iparams_cse), 
     .                           vma_data(ipar_flag_cse),
     .                           vma_data(iL1_cyc_cse), 
     .                           vma_data(iL2_cyc_cse), 
     .                           vma_data(ictol_cse), 'P')
                 call wr_igs_clk( uns, vma_data(iparams_cse), 
     .                           vma_data(ipar_flag_cse),
     .                           vma_data(isvcL1_ce),
     .                           vma_data(idata_flag_cse),
     .                           vma_data(ictol_cse ), cref_upd, 'P' )
              endif   

          end if
      end if


*     Scan the one-way data and remove data with closely spaced 
*     bias flags and biases near the tail

      call trim_oneways(vma_data(ictol_cse), vma_data(idata_flag_cse),
     .                  'Y')

*     Recompute the post-fit RMS so that we have the correct number of measurements
      cpass = pc_max_iter
      call comp_pf_rms(6, cpass, vma_data(iL1r_phs_cse), 
     .              vma_data(iL2r_phs_cse),
     .              vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .              vma_data(ictol_cse),    vma_data(idata_flag_cse),
     .              vma_data(iparams_cse),  vma_data(ipar_flag_cse)  )

*     See if LC_AUTCLN command used and Wide lanes will be resolved
*     (and biases passed to solve).
      if( resolve_wl ) then

         call report_wlstat(uns)
         if( lca_type(1:3).eq.'SEQ' )
     .   call zero_dd_wl(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(ibf_type_cse), vma_data(iazel_cse))
         if( lca_type(1:3).eq.'DIR' ) then
         call est_dd_wl(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(ibf_type_cse), vma_data(iazel_cse) )
         endif

C        call est_wlfull(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
C    .        vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
C    .        vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse),
C    .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
C    .        vma_data(ibf_type_cse), vma_data(iazel_cse))

      endif


*     See if normal points are to be formed: 
      if( np_size.gt.1 ) then
          call systime( date,sectag)
          call report_stat('status','autcln','main',' ',
     .                     'Forming normal points',0)
          write(*,100) '+Start flag Normal Pt', date, sectag
          call flag_np(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .        vma_data(iL1r_rng_cse),  vma_data(iL2r_rng_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse)  )
      end if

*     Update the cfiles
      call systime( date,sectag)
      call report_stat('status','autcln','main',' ',
     .                 'Outputing clean c-files',0)
      write(*,100) '+Start update Cfiles ', date, sectag
      call update_cfiles(vma_data(iL1r_phs_cse), vma_data(iL2r_phs_cse),
     .        vma_data(iL1_cyc_cse),  vma_data(iL2_cyc_cse),
     .        vma_data(iL1r_rng_cse),  vma_data(iL2r_rng_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse), 
     .        vma_data(iparams_cse), vma_data(ipar_flag_cse),
     .        vma_data(ipf_dphs_cse)   )

*     Report the bias flag status
      call report_bf( vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse),
     .        vma_data(ictol_cse), vma_data(idata_flag_cse),
     .        vma_data(ibf_type_cse) )

*     If desired dump the run parameters for this analysis
      call rep_run_params(6)
      if( uns.ne.6 ) call rep_run_params(uns)
 
*     Now generate the Gobs file.
      if( caprog_name(1:4).eq.'CTOG' ) then
          call systime( date,sectag)
          write(*,100) '+Start Gen Gobs       ', date, sectag
          call gen_gobs(vma_data(iL1_cyc_cse), vma_data(iL2_cyc_cse), 
     .                vma_data(idata_flag_cse),
     .                vma_data(iL1r_rng_cse), vma_data(iL2r_rng_cse), 
     .                vma_data(iparams_cse), vma_data(ipar_flag_cse) )
      end if
      call systime( date,sectag)
      write(*,100) '+Finished run         ', date, sectag
      call report_stat('status','autcln','main',' ','Finished',0)
      close(uns)
      stop 'Normal finish of autcln'
 
****  Thats all
      end
      
*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include 'ctog_cmds_bd.f' 
