CTITLE proc_ctog_cmd
 
      subroutine proc_ctog_cmd(inline, indx, iel, finished, 
     .                         first_word )

      implicit none
 
*     This routine will process the command line read from the input
*     file

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/cfile_def.h'      
      include 'ctogobs_com.h'
 
 
* PASSED VARIABLES
 
*   iel     - command number
*   indx    - Position in string.
 
      integer*4 iel, indx
 
*   inline  - Line read from command file
*   first_word  - First word of command
 
      character*(*) inline, first_word
 
*   finished    - Indicates we are finished
 
      logical finished
 
* LOCAL VARIABLES
 
*   trimlen     - Returns length of non-blank portion of string
*   ierr        - IOSTAT errors on decode
*   i           - Loop counter
*   jel         - Entruy number when looking up site names
*   jndx        - Pointer to position in temporary strings
*   prn_read    - PRN number in pre-edit command
*   start_ep_read, stop_ep_read - Start and stop epoch numbers in
*                 pre_edit command
*   ivalues(10) - Integer values to be read from command line
 
      integer*4 trimlen, ierr, i, jel, jndx, prn_read,
     .          start_ep_read, stop_ep_read, ivalues(10),
     .          jerr

*   cf_code_read - Cfile code read in commands for specific cfiles
*   ans          - Generic answer variable (yes or no)

      character*4 cf_code_read, ans

 
*   values(10)  - Array for reading values which are
*               - then converted to internal units
 
      real*8 values(10)
 
*   cd          - Dummy character
 
      character*80 cd
 
****  Check iel to see if OK
      if( iel.le.0 ) then
 
          write(*, 50) iel, (inline(1:max(trimlen(inline),1)))
  50      format(' ERROR ',i3,' decoding ',a)
          RETURN
      end if
 
****  Now process command
 
      goto( 100,  200,  300,  400,  500,  600,  700,  800,  900, 1000,
     .     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     .     2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
     .     3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     .     4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 
     .     5100      ) iel
 
*     END command
 100  continue
          finished = .true.
          RETURN
 
*     RNG_JUMP_TOL: Tolerance for the number of times the clock
*     stability will be considered at jump and a minimum which will be
*     considered a jump .i.e., both numbers must be exceeded for a jump
 200  continue
          call multiread(inline, indx, 'R8', ierr, values, cd, 2)
          if( values(1).gt.0 ) rng_jump_tol = values(1)
*         Convert value from usec to cycles
          if( values(2).gt.0 ) rng_jump_min = values(2)*1.d-6*fClk
          RETURN
 
*     CLK_RESET_TOL: Tolerance for jump to be taken as a msec reset
*     in the clock (usec, default 10 usec)
 300  continue
          call read_line(inline, indx, 'R8', ierr, values, cd)
          reset_tol = values(1)*1.d-6*fClk
          RETURN
 
*     RNG_RESID_TOL: tolerance for bad range residuals.  The first is
*     multiplier on the sigma on the data and second is a minimum value 
*     in meters.
* MOD TAH 990917: Added new argument for max range error allowed.
 400  continue
          call multiread(inline, indx, 'R8', ierr, values, cd, 2)
          if( values(1).gt.0 ) rng_res_tol = values(1)
          if( values(2).gt.0 ) rng_res_min = values(2)/vel_light*fClk
          call read_line(inline, indx, 'R8', ierr, values(3), cd)
          if( values(3).gt.0 .and. ierr.eq.0 ) 
     .                         rng_res_max = values(3)/vel_light*fClk

          RETURN

*     RNG_NOISE: Lets user specficy the noise in the range measurements
*     by site (input units, mm (converted to cycles, and will be updated
*     during ctogobs run.)
 500  continue
          call GetWord( inline, cf_code_read, indx)
          call read_line(inline, indx, 'R8', ierr, values,cd)
*         If the cf_code_read is ALL then copy the range noise
*         to all sites, other wise find the actual site
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              do i = 1, num_cfiles
                 rng_noise(i) = values(1)
              end do
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.gt.0 ) rng_noise(jel) = values(1)
          end if
          RETURN

*     MAX_RCLK_ITERATIONS: Maximum nnumber of range clock iterations
 600  continue
          call read_line(inline, indx, 'R8', ierr, values, cd)
          if( values(1).gt.0 ) max_rclk_iter = nint(values(1))
          RETURN
 
*     REL_CLK_WEIGHT : Weight to be given to the clock noise model in
*     estimating the clocks
 700  continue
          call read_line(inline, indx, 'R8', ierr, values, cd)
          rel_clk_wght = values(1)
          RETURN
 
*     RNG_CLK_ROOT: Lead part of the name to be given to the rclk
*     solution output.  If no root is given then the range clock
*     solution will not be output
 800  continue
          call read_line(inline, indx, 'CH', ierr, values, rng_clk_root)
          RETURN

*     PHS_CLK_ROOT: Lead part of the name to be given to the rphs
*     solution output.  If no root is given then the phase clock
*     solution will not be output
 900  continue
          call read_line(inline, indx, 'CH', ierr, values, phs_clk_root)
          RETURN
 
*     PHS_RES_ROOT: Lead part of the name to be given to the phase 
*     residual output.  If no root is given then the phase residuals
*     will not be output
1000  continue
          call read_line(inline, indx, 'CH', ierr, values, phs_res_root)

* MOD TAH 130329: See if single file or multiple file option has been 
*         passed.
          call read_line(inline, indx, 'CH', ierr, values, cd)
          if( ierr.eq.0 ) then 
             if( cd(1:1).eq.'N' .or.cd(1:1).eq.'n' ) then
                dph_output = 'NONE'
             elseif( cd(1:1).eq.'S' .or.cd(1:1).eq.'s' ) then
                dph_output = 'SINGLE'
             elseif( cd(1:1).eq.'M' .or.cd(1:1).eq.'m' ) then
                dph_output = 'MULTIPLE'
             else
                call report_stat('warning','autcln',
     .             'proc_ctog_cmd/phs_res_root unknown file type ',
     .                   ' ',cd,0)
             end if
          end if 

          RETURN

*     SNG_DIFF_ROOT: Lead part of the name for single difference files
*     (These may be used with the program mon_data to processing single
*     difference kinematic data).  The single differences are between
*     first site and subsequent sites.  The residual site list is used 
*     select which sites to write to single difference files.
1100  continue
          call read_line(inline, indx, 'CH', ierr, values, 
     .                   sng_diff_root)
          RETURN

*     RESIDUAL_SITE: Lets the user specify which sites should be output
*     to the residual files
1200  continue
          call decode_option(inline(indx:), cf_codes, num_cfiles, 
     .                       phs_res_sites, -1)
          RETURN

*     RCV_ALLAN_SD: Allows the user to specify the standard deviations
*     of the clocks at each site.
1300  continue
          call GetWord( inline, cf_code_read, indx)
          call read_line(inline, indx, 'R8', ierr, values,cd)
*         If the cf_code_read is ALL then copy the allan standard
*         deviation to all sites, other wise find the actual site
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              do i = 1, num_cfiles
                 rc_allan_sd(i) = values(1)*1.d-9
              end do
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.gt.0 ) rc_allan_sd(jel) = values(1)*1.d-9
          end if
          RETURN

*     REMOVE_BIAS_COND -- Chi**2 constrasts for removing biases
*         The three values specified are:
*         (1) Chi**2 ratio
*         (2) minimum value of smallest chi**2 and
*         (3) Gap over which double difference baiases will not be
*             attempted (seconds)
1400  continue
          call multiread(inline, indx, 'R8', ierr, values,cd, 3)
          dchi2_ratio = values(1)
          dchi2_min_val = values(2)
          dchi2_max_sep = values(3)
*         Get the optional gap scaling factor (ignore error since
*         optional)
          call read_line(inline,indx,'R8',ierr, values(1), cd)
          if( ierr.eq.0 ) then
              dchi2_gap_fact = values(1)
          end if
          RETURN

*     USE_GAMIT_ELEV - Allows user to specify if to use GAMIT Elevation
*     cut-off angle or not. If Yes is specifiied then the gamit cutoff
*     is used
1500  continue
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              use_gamit_elc = .true.
          else
              use_gamit_elc = .false.
          end if
          RETURN

*     USE_CVIEW_EDIT - Allows user to specify if to use CVIEW editing
*     flag (-1) or not. If Yes is specifiied then cview flags -1 will
*     not be used.
1600  continue
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              use_cview_edit = .true.
          else
              use_cview_edit = .false.
          end if
          RETURN

*     IGNORE_GAPS - Lets user specify that gaps should be ignored
*     when forming acceptable double difference during cleaning.  
*     This option should only be used for cleaned data when the
*     GAMIT elevation cutoff and cview edits are used.
1700  continue
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              usr_ignore_gaps = .true.
          else
              usr_ignore_gaps = .false.
          end if
          RETURN

*     EDIT_SITE_SV - Allows user to specify site/satellite combinations
*     over spcific epoch ranges to be edited and not used in determining
*     clock behavior and double difference editing
1800  continue
          call GetWord( inline, cf_code_read, indx)
          call read_line(inline, indx, 'I4', ierr, prn_read,cd)
          call read_line(inline, indx, 'I4', ierr, start_ep_read,cd)
          call read_line(inline, indx, 'I4', ierr, stop_ep_read,cd)

*         Increment the number of pre-edits.  If we have too many then
*         warn user and skip this one; else add to the list
          num_pre_edit = num_pre_edit + 1
          if( num_pre_edit.gt.max_pre_edit ) then
              write(*,1820) max_pre_edit, inline(1:trimlen(inline))
1820          format('** WARNING ** Too many EDIT_SITE_SV commands,',
     .               ' Max allowed is ',i4,/,
     .                a,' being ignored **')
              num_pre_edit = num_pre_edit -1
              RETURN
          end if

*         If the cf_code_read is ALL then set the pre_edit site number
*         zero and this will edit all sites
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              pre_edit(1,num_pre_edit) = 0
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.le.0 ) then
*                 Site not found so ignore this line
                  num_pre_edit = num_pre_edit -1
                  RETURN
              end if
              pre_edit(1,num_pre_edit) = jel
          end if

*         Save the rest of the information
          pre_edit(2,num_pre_edit) = prn_read
          pre_edit(3,num_pre_edit) = start_ep_read
          pre_edit(4,num_pre_edit) = stop_ep_read
          RETURN

*     PHS_FIT_TOL - Tolerances in deciding if a cycle slip has
*     occurred in pre-fit clock fit data
1900  CONTINUE
          call multiread(inline,indx,'R8',ierr, values,cd, 4)
          phs_fit_tol(1) = values(1)
          phs_fit_tol(2) = values(2)
          phs_fit_tol(3) = values(3)
          phs_fit_tol(4) = values(4)
          RETURN

*     STATUS_REPORT - Allows user to tailor the output of the program
*     by selecting which quanities will be output.
2000  CONTINUE
          call decode_option(inline(indx:), status_rep_opts, 
     .                       num_ctog_status, status_rep, -1)
          RETURN

*     DD_REPORT - Allows user to specify a file of the format readable
*     by cview and to specifiy which types of double differences should
*     by output to this file (i.e., ALL, FIXED, NOT_FIXED')
2100  CONTINUE
          call read_line(inline, indx, 'CH',ierr, values, dd_outfile)
          call read_line(inline, indx, 'CH',ierr, values, dd_out_opts)
          call casefold(dd_out_opts)
          RETURN

*     MIN_ELEVATION - Minimum elevation to which data will be cleaned.
*     Once this value has been set in ctogobs, the data below this 
*     elevation will not be useable later without furthur cleaning.
*     WARNING: Use of this command will overwrite the default setting
*          for site dependent elevation cutoffs.
2200  continue
          call read_line(inline, indx, 'R8',ierr, values,  cd)
          min_ctog_elev = values(1)*pi/180.d0

*         See if output minimum elevation passed.  If not then 
*         set value to same as min_ctog_elev.
          call read_line(inline, indx, 'R8',ierr, values,  cd)
          if( ierr.eq.0 ) then
               min_out_elev = values(1)*pi/180.d0
          else
               min_out_elev = min_ctog_elev
          end if
          if( min_out_elev.lt.min_ctog_elev ) then
               write(*,2210) 
 2210          format('**NOTE** MIN_ELEVATIONS around wrong way:',
     .                ' Switching to make correct')
               values(1)    = min_out_elev
               min_out_elev = min_ctog_elev
               min_ctog_elev = values(1) 
          end if
          do i = 1, num_cfiles
             site_celev(i) = min_ctog_elev
             site_oelev(i) = min_out_elev
          end do 
*         Indicate that we will not be using the GAMIT cutoff.            
          use_gamit_elc = .false.
          
          RETURN

*     TRIM_ONEWAY_TOL - Lets user set the tolerances in triming the
*     oneway data to remove small segments of data.  The four values
*     input are:
*     min_dtl_bias - mininum time in seconds between bias flags
*     min_good_bias - minunum number of epochs between bias flags and
*                     maximum time over which a bias flag will be 
*                     removed in one-way data.
*     min_dtr_end   - Fraction of total duration of data allowed for
*                     a bias flag at the end of one-way sequence. 
*     min_good_end  - Number of epochs of data allowed before last bias
* MOD TAH 200617: Added trim_seg and negative value for default)
*     Y/N         - trim_seg option. 
*                     flag.
*     If these conditions are not meet then the one-way data are flagged.
* MOD TAH 200617: Added negative default option
2300  continue
          call multiread(inline,indx,'R8',ierr, values, cd, 4)
          if( values(1).gt.0 )
     .    min_dtl_bias = values(1)    ! seconds
          if( values(2).gt.0 )
     .    min_good_bias = nint(values(2))   ! Number of epochs
          if( values(3).gt.0 )
     .    min_dtr_end   = values(3)         ! Fraction remaining at end
          if( values(3).gt.0 )
     .    min_good_end  = nint(values(4))   ! Number of epochs
*         see if trim_seg option given (optional)
          call read_line(inline,indx,'CH',ierr, ivalues(1),cd)
          call casefold(cd)
          if( cd(1:1).eq.'Y' ) trim_seg = .true.
          if( cd(1:1).eq.'N' ) trim_seg = .false.
          RETURN

*     DD_RETURN_SIZE - Lets user sets the number of data to try to use
*     when double difference cleaning.  The arguments are:
*     max_wl_ret  - Maximum number of values to return in one-ways
*                   for widelanes (when available)
*     max_dd_ret  - Maximum number of double differences to return
*     max_lg_use  - Maximum number of values to use in the LG fit
*     tol_one_way_fix - Maximum number of seconds in gap and still
*                   allow one-way bias removal (for WL data) (seconds)
*     (negative values will keep default)
 
2400  continue
          call multiread(inline,indx,'I4',ierr, ivalues, cd, 3)
          if( ivalues(1).gt.0 )
     .    call set_dd_ret('WL Returns', ivalues(1), max_wl_ret, 
     .                     max_max_wl_ret)
          if( ivalues(2).gt.0 )
     .    call set_dd_ret('DD Returns', ivalues(2), max_dd_ret, 
     .                     max_max_dd_ret)
          if( ivalues(3).gt.0 )
     .    call set_dd_ret('LG Use', ivalues(3), max_lg_use, 
     .                     max_max_dd_ret)

*         Get the one-way fix tolerance (seconds) (Optional)
          call read_line(inline,indx,'I4',ierr, ivalues(1),cd)
          if( ivalues(1).gt.0 ) tol_one_way_fix = ivalues(1)


          RETURN

*     DD_FIT_TOL  - Lets user set tolerances for flagging cycle slips
*     in the double difference data.  The two values entered are for
*     WL lane observables and one for LC double differences.  The WL is given
*     as number of WL sigmas (based on range noise estimates from clock
*     solution) and a max amount allowed (cycles).  The DD are given in cycles.
2500  continue
          call multiread(inline,indx,'R8',ierr, values, cd, 6)
          dd_wl_tol(1) = values(1)
          dd_wl_tol(2) = values(2)
          dd_wl_tol(3) = values(3)
          dd_lc_tol(1) = values(4)
          dd_lc_tol(2) = values(5)
          dd_lc_tol(3) = values(6)
          RETURN

*     USE_MM_RANGE - Lets user say whether MiniMac ranges should be used.
2600  continue 
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              use_mm_ranges = .true.
          else
              use_mm_ranges = .false.
          end if
          RETURN

*     ALLOW_ONE_BG - Lets user say whether to a second pass should be
*     made if no double differences can be found because all satellites
*     at a station have a bias or gap at the same time.
2700  continue 
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              do_one_bg     = .true.
          else
              do_one_bg     = .false.
          end if

*         Force flagging of gaps if this feature is turned on
          if( do_one_bg ) gaps_flagged = .true.
          RETURN

*     ION_JUMP_TOL: Lets user specify the time cap and random walk
*     parameters for detecting jumps in the ionospheric delay.
*     NOTE: the time limit given by the last use of this command
*     will be used.
2800  continue
          call GetWord( inline, cf_code_read, indx)
          call multiread(inline, indx, 'R8', ierr, values,cd,4)
          dt_ion_tol = values(1)
          if( values(2).le.0 ) then
              write(*,2820) 
2820          format('**INVALID ion_rw_tol given in',
     .               ' ION_JUMP_TOL command')
              RETURN
          end if

*         If the cf_code_read is ALL then set the ion_rw_tol for      
*         all sites
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              do i = 1,num_cfiles
                 ion_rw_tol(1,i) = values(2)
                 ion_rw_tol(2,i) = values(3)
                 ion_rw_tol(3,i) = values(4)
              end do
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.gt.0 ) then
                  ion_rw_tol(1,jel) = values(2) 
                  ion_rw_tol(2,jel) = values(3) 
                  ion_rw_tol(3,jel) = values(4)
              end if
          end if
          RETURN

*     SCAN_SITE: Lets the user specify which sites should be scanned in
*     double diffrences before cleaning in double differences
2900  continue
          call decode_option(inline(indx:), cf_codes, num_cfiles, 
     .                       scan_sites, -1)
          RETURN

*     REMOVE_MS_JUMP - Lets user decide if millisecond jumps should be
*     removed from the ranges as the new cfiles are written.
3000  continue 
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              remove_ms_jump = .true.
          else
              remove_ms_jump = .false.
          end if
          RETURN

*     FLAG_GAPS - Lets user decide if gaps are to be flagged rather
*     than treated as gaps.  (When this option is used, ignore_gaps
*     will be set true.  Also allow_one_bg will automatically turm
*     this feature on.
3100  continue 
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              gaps_flagged = .true.
          else
              gaps_flagged = .false.
          end if
          RETURN

*     SUMMARY_FILE - Name of a file for output of the summary information
*     rather than it just being sent to the screen.
 3200 continue
          call read_line(inline, indx, 'CH',ierr, values, summary_file)
          if( ierr.ne.0 ) summary_file = '6'
          RETURN

*     REMOVE_FIRST_BIAS - Allows user to say that the first bias 
*     flag on the one-way data should be removed.  (Mainly for gamit
*     processing)
3300  continue 
          call read_line(inline, indx, 'CH',ierr, values, ans)
          call casefold(ans)
          if( ans(1:1).eq.'Y' ) then
              remove_first_bias = .true.
          else
              remove_first_bias = .false.
          end if
          RETURN

*     GAP_SIZE:  Size of gap to be allowed at each sites.  Automatically
*     turns on the flag_gaps options.
3400  continue
          call GetWord( inline, cf_code_read, indx)
          call read_line(inline, indx, 'I4', ierr, ivalues,cd)

*         If the cf_code_read is ALL then set the ion_rw_tol for      
*         all sites
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              do i = 1,num_cfiles
                 gap_size(i) = ivalues(1)
              end do
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.gt.0 ) then
                  gap_size(jel) = ivalues(1)
              end if
          end if
          gaps_flagged = .true.
          RETURN

*     MAX_SCAN_EDIT - Set the maximum number of dd bias flags that
*     can added during scanning before all the site/sv data is edited.
 3500 continue
          call read_line(inline, indx, 'I4', ierr, max_scan_edit,cd)
* MOD TAH 050217: See if second argument passed with min_ow_data
          call read_line(inline, indx, 'I4', ierr, ivalues,cd)
          if( ierr.eq.0 .and. ivalues(1).ne.-1 ) then
              min_ow_data = ivalues(1)
          endif
          RETURN

*      NP_SET      - Set the normal pointing parameters (size and start)
 3600  continue
          call read_line(inline, indx, 'I4', ierr, np_size,cd)
          call read_line(inline, indx, 'I4', ierr, np_start,cd)

*         Make sure np_size is an odd number
          if( np_size-int(np_size/2)*2.eq.0 ) then
              write(*,3605) np_size + 1
 3605         format('**WARNING** NP size must be odd, changing value',
     .               ' to ',i3)
              np_size = np_size + 1
          end if                        
          
*         If second argument not passed then default to starting at
*         first epoch.
          if( ierr.eq.-1 ) np_start = 1
          RETURN

*     SITE_PARAMS  - Sets elevation cutoff and SNR limits by site
 3700 continue
          call GetWord( inline, cf_code_read, indx)

*         Get values for cleaning min elevation, output min. elev
          call read_line(inline, indx, 'R8', ierr, values(1),cd)
          call read_line(inline, indx, 'R8', ierr, values(2),cd)
*         Get values for min L1 and L2 SNR values.
          call read_line(inline, indx, 'I4', ierr, ivalues(1), cd)
          call read_line(inline, indx, 'I4', ierr, ivalues(2), cd)

*         If the cf_code_read is ALL then set the ion_rw_tol for      
*         all sites
          call casefold(cf_code_read)
          if( cf_code_read.eq.'ALL ' ) then
              do i = 1,num_cfiles
                 site_celev(i) = values(1)*pi/180.d0
                 site_oelev(i) = values(2)*pi/180.d0
*                If values around the wrong switch so that celev
*                is smaller than oelev.
                 if( site_oelev(i).lt.site_celev(i) ) then
                     write(*,3710)  cf_code_read
 3710                format('**NOTE** SITE_PARAM elevations wrong order'
     .                      ,' for site ',a4,'.  Switching')
                     values(1)  = site_oelev(i)
                     site_oelev(i) = site_celev(i)
                     site_celev(i) = values(1)
                 end if
                 site_snr(1,i) = ivalues(1)
                 site_snr(2,i) = ivalues(2)
              end do
          else
              jndx = 1
              call get_cmd(cf_code_read, cf_codes, num_cfiles,
     .                         jel, jndx)
              if( jel.gt.0 ) then
                 site_oelev(jel) = values(1)*pi/180.d0
                 site_celev(jel) = values(2)*pi/180.d0
                 if( site_oelev(jel).lt.site_celev(jel) ) then
                     values(1)  = site_oelev(jel)
                     site_oelev(jel) = site_celev(jel)
                     site_celev(jel) = values(1)
                 end if
                 site_snr(1,jel) = ivalues(1)
                 site_snr(2,jel) = ivalues(2)
              end if
          end if
          
*         Indicate that we will not be using the GAMIT cutoff.            
          use_gamit_elc = .false.

          RETURN
          
*     APPLY_PHS_CLK  - no arguments, sets true that phase clocks should
*     applied when writing out the cfiles.
*     APPLY_PHS_CLK [MAX Iter] [Non-int] [Converged %] [Over shoot]
*
 3800 continue
 
          apply_phs_clk = .true.
          call read_line(inline, indx, 'I4', ierr, ivalues(1), cd)
          if( ierr.eq.0 ) pc_max_iter = ivalues(1)
*         Make sure value is even.          
          if( int(pc_max_iter/2)*2-pc_max_iter.eq.-1 ) 
     .                             pc_max_iter = pc_max_iter+1  
*         Get number of iterations before non-integer cycles used
          call read_line(inline, indx, 'I4', ierr, ivalues(1), cd)
          if( ierr.eq.0 .and. ivalues(1).ge.0 ) pc_non_int = ivalues(1)
*         Get percentage chnage fro convergeence
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 .and. values(1).gt.0 ) pc_convd = values(1)
*         Get over shoot when means removed.
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 .and. values(1).gt.0 ) pc_over_shoot =values(1)
*         Get maximum iterations with full analysis code         
          call read_line(inline, indx, 'I4', ierr, ivalues(1), cd)
          if( ierr.eq.0 .and. ivalues(1).ge.0 ) pc_full_anal =ivalues(1)
*         Get the maximum correlation to be used in full code.
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 .and. values(1).gt.0 ) pc_max_corr  = values(1)

          RETURN
          
*      USE_POSTFIT <mfile name>
 3900  continue
          call gen_mfile_name(inline, indx, df_name, mf_name)
          call open_mf(204, mf_name, ierr)
          if( ierr.eq.0 ) call read_mf(204, ierr)
          close(204)
          if( ierr.ne.0 ) then
              mf_name = ' '
              use_postfit = .false.
          else
              use_postfit = .true.
          end if
          RETURN
              
 
*      POSTFIT_EDIT: Edit postfit editing.
*      POSTFIT_EDIT [Start Iter] [Nsigma] [Max Restote] [Max RMS]
 4000  continue

*         Get the itertation at which we should start editing
          edit_postfit = .true.
          apply_phs_clk = .true. 
          call read_line(inline, indx, 'I4', ierr, ivalues(1), cd)
          if( ierr.eq.0 ) pc_start_edit = ivalues(1)

*         Get the n-sigma limit
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 ) pf_nsig = abs(values(1))

*         Maximum Residual to restore
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 ) pf_maxres = abs(values(1))

*         Maximum RMS for residuals (cycles)
          call read_line(inline, indx, 'R8', ierr, values(1), cd)
          if( ierr.eq.0 ) pf_max_rms = abs(values(1))
            
          RETURN

*      PF_REMOVE_BF: Removes bias flags in one-ways after the phase
*      clocks have been applied.
 4100  continue
          pf_remove_bf = .true.
          RETURN 

*      AZ_MASK : Lets user specify azimuth dependent elevation angle
*      mask.  Form is az_mask <site> <az1> <el1> <az2> <el2> <az3> ...
 4200  continue
          write(*,*) inline 
          call read_az_mask( inline, indx )
          RETURN

*      IGS_CLK_FILE: Specifies name of IGS clock file and the set of
*      sites to be used for reference
 4300  continue
          call read_line(inline, indx, 'CH', ierr, values, igs_clk_file)
          call report_error('IOSTAT',ierr,'decode','igs_clk_file',0,
     .                      'PROC_CTOG_CMDS')
          if( ierr.eq.0 ) then
              write_igs_clk = .true.
          else
              write_igs_clk = .false.
              igs_clk_file = ' '
              RETURN
          end if

*         Now get the sampling interval
          call read_line(inline,indx,'I4',ierr,igs_clk_samp,cd)

*         Now get the RMS for reference clocks and the max rms allowed in
*         output
          call read_line(inline,indx,'R8',ierr,rms_ref_clk,cd)
          call read_line(inline,indx,'R8',ierr,rms_max_clk,cd)

*         Now see if list of reference stations passed
          num_ref_clk = 0     
          ierr = 0
          do while ( ierr.eq.0 )
             call read_line(inline,indx,'CH',ierr,values,
     .                      cf_code_read)
             if( ierr.eq.0 ) then
                 call casefold(cf_code_read)
                 jndx = 1
                 call get_cmd(cf_code_read, cf_codes, 
     .                            num_cfiles, jel, jndx)
                 if( jel.gt.0 ) then
                     num_ref_clk = num_ref_clk + 1
                     ref_clk_code(num_ref_clk) = cf_code_read
                 end if
             end if
          end do

          RETURN

*      USE_ORIG_BF used to control use of original 
*      bias flags in data set.  With version 3.15, the default
*      is no.
 4400  continue
          use_orig_bf = .true.
          RETURN
	  
*      MAX_CHAN allows user to specify max channels in receivers
*      override the normal max.  The full (very slow) clock 
*      estimation can not be used in this if num_chan > max_gchannels
 4500  continue
         call read_line(inline, indx, 'I4', ierr, num_chan, cd)	  
	 RETURN 
                     
*      LC_AUTCLN: Set parameters to have autcln resolve the widelanes
*      for solve
 4600  continue
          resolve_wl = .true.
* MOD TAH 080523: See if algorith type past
          jndx = indx
          call GetWord(inline, lca_type, jndx)
          call casefold(lca_type)
          if( lca_type(1:3).eq.'DIR' ) then
              indx = jndx
          elseif( lca_type(1:3).eq.'SEQ' ) then
              indx = jndx
              dchi_wl_tol =   10     ! Dchi change from best choice.
              msig_wl_tol =  0.30d0  ! Maximum sigma allowed for wl to 
                                     ! resolved.
          else
              call check_num(lca_type,jerr) 
              if( jerr.eq.0 ) then
*                 Argument is a number, so process and set
*                 default lca_type: 
* MOD TAH 080714: Make default type DIR
                  lca_type = 'DIR'
              else   ! Some unknown string
                  call report_stat('warning','autcln',
     .                 'proc_ctog_cmd/lc_autcln',' ',lca_type,0)
                  lca_type = 'DIR' 
                  indx = jndx
              end if
          endif
              
*         See if arguments passed
          call read_line(inline, indx, 'I4', ierr, ivalues, cd)
          if( ierr.eq.0 ) then
              if( ivalues(1).gt.0 ) min_wl_tol = ivalues(1)
          endif
          if( ierr.eq.0 ) then
              call read_line(inline, indx, 'R8', ierr, values, cd)
              if( values(1).gt.0 ) dchi_wl_tol = values(1)
          end if
          if( ierr.eq.0 ) then
              call read_line(inline, indx, 'R8', ierr, values, cd)
              if( values(1).gt.0 ) msig_wl_tol = values(1)
          end if
          if( ierr.eq.0 ) then
              call read_line(inline, indx, 'R8', ierr, values, cd)
              if( values(1).gt.0 ) mdev_wl_tol = values(1)
          end if
* MOD TAH 050217: Added option to output acbias.dat file.
          if( ierr.eq.0 ) then
              call read_line(inline, indx, 'CH', ierr, values, cd)
              if( cd(1:1).eq.'Y' .or.cd(1:1).eq.'y' ) 
     .                 acbias_out = .true.
         end if
         RETURN

*     L1ONLY: Sets nol1only so that L1-only data will be 
*     processed (reversed from original nol1only command)
* MOD AH 091228: Reverse the nol1only commmand to L1only and
*     this then sets it false. 
 4700 continue
         nol1only = .false.
         RETURN

* MOD TAH 200508: Added second option to prescan clocks for milli-
*     second jumps.
*     PREFIT_CLK <Y/N> <Y/N>: Tell autcln to prefit the linear clocks to the
*     range residuals to remove any drifts in the clocks.  
*     Should make autcln independent of I-file used in model
 4800 continue
         call read_line(inline, indx, 'CH', ierr, values, cd)
         prefit_clk = .true.
         if( ierr.eq.0 ) then
            if( cd(1:1).eq.'N' .or.cd(1:1).eq.'n' ) 
     .           prefit_clk = .false.
         endif
* MOD TAH 200509: See if prescan_ms option given. 
         call read_line(inline, indx, 'CH', ierr, values, cd)
         prescan_ms = .false.
         if( ierr.eq.0 ) then
            if( cd(1:1).eq.'Y' .or.cd(1:1).eq.'y' ) 
     .           prescan_ms = .true.
         endif 
         RETURN 

*     APP_ION command to apply ionospheric delay to omc-values
*     when c-files are written.  First part of code to process
*     GLONASS.
 4900 continue
         app_ion = .true.
         RETURN

*     REMAP_GLONASS <Y/N>
 5000 continue
*        See what option was passes
         call GetWord(inline, cd, indx) 
         if( cd(1:1).eq.'Y' .or. cd(1:1).eq.'y' ) then
            remap_glonass = .true.
            do i = 1,cf_nsat
               fL1u(i) = cf_fL1(i)
               fL2u(i) = cf_fL2(i)
            end do
         else
            remap_glonass = .false.
            do i = 1,cf_nsat
               fL1u(i) = fL1(i)
               fL2u(i) = fL2(i)
            end do
         endif
         RETURN

*     Next command
 5100 continue
       

***** Thats all for the moment
      END

CTITLE READ_AZ_MASK

      subroutine read_az_mask( inline, indx)

      implicit none

*     Routine to read the az_mask entries from the line passed
*     by the user.

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARAIBLES
* indx  -- Current position in string
* inline -- Command buffer read from user file

      integer*4 indx
      character*(*) inline

* LOCAL VARIABLES
* ierr -- IOSTAT error reading line
* jel  -- Station number based on cf_code 
* na   -- Running sum of number azimuth entries

      integer*4 ierr, jel, na, jndx

* word -- Word read from line (contains station code)
* cdum -- Dummay entry for readline

      character*4 word, cdum

* rval -- Real*4 number read from line
      real*4 rval


***** Get the name of the site:
      call GetWord( inline, word, indx)

*     Get the site code:
      call casefold(word)
      jndx = 1
      call get_cmd( word, cf_codes, num_cfiles, jel, jndx )

*     If site is found, then process the rest of the line.
      write(*,*) 'Site ', word, ' Number ',jel
      if( jel.gt. 0 ) then
         na = 0
         ierr = 0
         do while ( ierr.eq.0 .and. na.lt.max_azmask )
            call read_line( inline, indx, 'R4', ierr, rval, cdum)
            if( ierr.eq.0 ) then
                na = na + 1
                az_mask(1,na,jel) = rval*pi/180.d0
*               Get the elevation cutoff
                call read_line( inline, indx, 'R4', ierr, rval, cdum) 
                az_mask(2,na,jel) = 0
                if( ierr.eq.0 ) az_mask(2,na,jel) = rval*pi/180.d0
            end if
            write(*,*) jel, na, az_mask(1,na,jel)*180/pi, 
     .                          az_mask(2,na,jel)*180/pi
         end do

*        Save number of values
         num_azmask(jel) = na
      end if

****  Thats all
      return
      end


 
