  
*     This is the common block include for the CTOGOBS and
*     AUTCLN programs.

* Set the version for ctogobs/autcln
* 2.07 - Changed the dimensioning of maxprm and maxlab to 
*        handle new cfile format (Only Ghandlers needed 
*        recompiling.
* 2.08 - Modified methods to treat clock jumps in the range data
*        Modified trim_one_ways to iterate incase there remain
*        short peices of data at the ends of one-ways.
* 2.09 - Added more conditions on allowing one-bias-gap to be
*        invoked.  Now it will not be used on the first iteration
*        cleaning or if the gap in one-ways is greater than 5% of
*        maximum gap over which bias flags will be removed.
*        Also added a second scanning loop with satellites in
*        reverse order after the first cleaning pass.
* 2.10 - Made the small gap condition for one-bias-gap grow with
*        cleaning pass number (5%,10%,20% of the max separation
*        over which biases will be removed).
*        Removed the second scanning loop (did not seem to do 
*        anything useful)
*        Increased the number of cleaning passes since these are
*        quite quick.
* 2.11 - Updated cfile format to add in new variables.  No changes
*        make to cleaner programs.
* 2.12 - Added the following new commands:
*        max_scan_edit -  Sets the number of epochs that can have
*           bias flags added during scanning before the all of the
*           data for that site/sv combination are removed from the
*           solution.  Used to remove satellites with "burns" and
*           receivers that have failed.  Default is 0 which will
*           stop this feature being evoked.
*        np_set       - Number of epochs of data to be used in
*           normal points and the offset of from the beginning to
*           of the data to start the normal pointing.  All epochs
*           of data (without bias flags) must be present to have the
*           normal point formed.
*        min_elev command extended to all additional argument for
*           miniumum elevation to be used in solved.  (Allows
*           cleaning to a lower elevation angle than will actually
*           be processed).  The additional argument is min_out_elev
*           and is set equal to min_ctog_elev if second argument
*           not given.
*  2.13 - Introduced new data checking commands and more site dependent
*         data flag.  The new command SITE_PARAMS allows the user
*         to specify by site min. cleaning elevation, output min.
*         elev, and min L1 and L2 SNR values.  Defaults are set
*         based on the type of receiver.  I new summary feature
*         giving number of data edited by reason is also added.
*  2.14 - Refined the patching algorithm. By checking the extimates
*         of the jumps from different observable types (only useful
*         when range data available); and by forcing the polynomial
*         coefficients to be the same on each side of the jump
*  2.15 - Increased the dimensioning of the number of pre-edits
*         and the number of cfiles.  (100->1000 and 50->100)
*  2.16 - Switched the order of the cleaning and output elevations
*         angles to be consistent with internal code.  (Cleaning
*         elevation angle comes first, but code will switch if wrong).
*  2.17 - Fixed a dimensioning problem in max_max_dd_ret.  Needs to be
*         same as max_max_wl_ret.  This error can result in a bounds
*         check error in get_dd_data.
*  2.18 - Added a <100 meter difference between L1 and L2 range as
*         a check on range quality.  (ctog_utm.f modified).  Bit 4
*         of data_flag used to show data bad (same bit as used equal
*         L1 and L2 range, but since makex now traps this condition,
*         this bit is no longer used.)
*  2.19 - Implemented report_stat output to monitor the process of the
*         program.
*  2.20 - Fixed small bug setting the defaults for dd_ret when small
*         sampling interval used (960821)
*  2.21 - Changed Ghandler create_cf routine so that cfiles will now
*         be overwritten even if they exist.
*-------------------------------------------------------------------
*  3.00 - Philosphical changes that allow autcln to read an mfile
*         and clean using postfit residuals rather than prefit (new
*         command USE_POSTFIT <mfile>.  Other new commands are 
*         APPLY_PHS_CLK -- The estimates of phase clocks are applied
*         to the one-way data.
*         EDIT_POSTFIT  -- Allows editing on one-way postfit residuals
*  3.01 - Adding sorting of the cfile list in the dfile when only the
*         dfile is used to gve the cfile names.
*  3.02 - Fixed bug associated with truncating half cycles when phase
*         clocks applied to half wavelength data.  Added check for zero
*         lambda factors. 
*  3.03 - Added saving of azimuth and elevation (R*4) and mfile adjustments
*         (both if needed).  Az El needed when phs_res_root used; and
*         mfile adjusments needed when np_set used with postfit residuals.
*  3.04 - Added pf_remove_bf which will remove bias flags from one-ways after
*         phase clocks have been applied.
*  3.05 - Added feature which allow apriori station coordinates and orbital
*         elements to be different between cfile and mfile used for 
*         postfit residuals.
*  3.06 - Modified the restore code in pf_edit mode to be more robust in
*         allowing data to be restored.  Fixed a problem with lost bias
*         flags when the rms of residuals was large. (980305)
*  3.07 - Introduced explicit estimation of bias parameters in phase clock
*         estimation mode.  Should help in global networks.  Bacause this
*         option takes a long time and in many cases is not needed, we
*         limit the number of iterations done this way.
*  3.08 - Updated to read new-cfile and m-file formats.  Some small changes
*         to the way post-fit residuals are computed.
*  3.09 - Two major changes made to the program: (a) Code reworked to allow
*         cleaning of L1 only data; (b) Introduced command for AZ_MASK which
*         lets users give elevation angle cutoff as a function of Azimuth.
*         An impact on normal processing from the (a) changes is that double
*         differences are now scanned allowing gaps upto the gap_size for the
*         site.  
*  3.10 - Refined the phase clock algorithm to do a better job when range data
*         is corrupt.  The new algorithm calls flat_dd to remove mean biases
*         in the double differences before attempting the final phase clock fit.
*  3.11 - Added computation of averged phase residuals to better estimate the
*         amount of correlated data noise.
*  3.12 - Improved the handling of bad clock and range data so that less data
*         is discarded.  Introduce rng_res_max as the maximum allowable range
*         error.  Also added check on pf_max_rms so that large RMS sites are
*         removed. 
*  3.13 - Fixed dimensioning of mfile arrays that cause problems with many stations
*         and atmospheres.
*  3.14 - Added feature to allow output of IGS clock files.  Needed to save 
*         cf_svcL1 (removed portion of satellite clock) so that it can be added
*         back later.  Also introdced new command and options.
*         Other changes: Defaulted allow_onebg to true, set
*  3.15 - Refined the igs clock module to better line up the phase and range 
*         clocks.  Also added default to not use existing bias flags.  The
*         use_orig_bf command can be used to use the original bias flags.
*  3.16 - Small non-printed statement added to proc_phasfin.f which seems
*         to fix at HPUX f90 bug.  Small similar change to ctogobs_dd.f
*  3.17 - Increased number of channels allowed to 14 and added max_chan
*         command that allows this to be increased further
*  3.18 - Cleaned up code associates with making igs clock files. 030222.
*  3.19 - Added LC_AUTCLN command which works with solve LC_AUTCLN
*         bias fixing mode so that WL ambiguities will be resolved in
*         autcln and the results passed to solve. 031213
*  3.20 - Modified to handler new solve M-files that have average and PWL atm
*         delays.  Reworked the LC_AUTCLN code in this version to generate
*         and use non-redundant double-difference biases. 040702
*  3.21 - Modified memory allocation to br 64-bit compatable.
*  3.22 - Modified to handle 1020 version of cfiles (050206)
*       - Added min_ow_data which sets minimum length of one way before site/sv
*         is deleted completed
*       - acbias_out -- Logical set in hidden lc_autcln command to write the
*         old acbias.dat file.
*  3.23 - Fixed problem with satellite clock being zero for some downweighted
*         data in model
*       - Added checks on missing L2 range data for zero_dd_wl code.
*       - Added P1 and P2 residuals to the DPH phase files (051116)
*  3.24 - Added additional calls to trim_oneway to eliminate short segments
*         of data during postfit editing (can remove high rms residuals)
*         Added better detection of need for bias flags when data restored during
*         postfit editing.
*         Fixed P1/P2 residual output. (060109)
*  3.25 - Small changes to add get_cmd call and treatment of short spans of data
*         new the ends of data.  May reduce singular matrix problems with LC_AUTCLN
*         (061017) 
*  3.26 - Small changes to setting whether WL-ambiguity can be fixed. (070312)  
*  3.27 - Added code to flag data with no double differences before resolving
*         ambiguities.  Fix should solve problems with singular matrices in solve
*         and bad dependent biases (070403)
*  3.28 - Added NOL1ONLY command to stop processing of L1 only data (needed in
*         large network processing).  Also changed some prefit range 
*         error detection (080512)
*  3.29 - Added SEQ and DIR options to lc_autcln command. The SEQ option is
*         original autcln sequence ambiguity resolution code.  The DIR option
*         is the new code that directly estimates the ambiguities.
*  3.30 - Replaced nol1only command with L1only command which now must
*         be given if L1 only data is to be processed.  All L2 data is
*         deleted with this command even for dual frequency receivers.
*  3.31 - Added ELMEAN, averaged phase residuals by site in 1-deg bins to
*         output summary file.  sh_plot_elmean can be used to plot the values.
*         Output status AZNAMEAN will also output azimuth means (<35 deg elev)
*         and satellite nadir means
*  3.32 - Updated to read and write 10.41 cfile version (dry/wet zen names
*         added (130119)
*  3.33 - Updated to write phase residuals to single file rather than by PRN
*         Argument added phs_res_root command.
*  3.34 - Updated to handle multi-GNSS by computing frequencies quantities
*         at execution rather than having them stored in libraries/includes/
*         const_param.h.
*  3.35 - Added removal of linear clocks based on psueorange before processing
*         Feature added with PREFIT_CLK command
*  3.36 - GNSS processing for GLONASS. Implemented by mapping to a common
*         frequency and applying ionospheric delay correction to omc values
*         output by autcln.  New APP_ION command added to apply ion correction.
*         Add lif1 and lif2 to compute ion delay at L1 in cycles.  Value set in 
*         set_freqs.f
*  3.37 - Added option to control whether Glonass ambiguities are integer or
*         scaled to frequency of satellite (default, non-integer, remap_glonass Y)
*  3.38 - Mod's for handle ECOMC model and mapping of GLONASS frequencies during
*         c-file writing.  190711.
*  3.39 - Mods to prefit clock and removing millisecond jumps in clocks.  Fix
*         Beidou clocks when L7 is mapped to L6 for processing.
*  3.40 - Improved GLONASS processing, added trim_shortseg option to 
*         trim_oneway command.  Added verify_scan to check UNFLAGGED slips.  
*         Added reporting of GNSS in outut files.
*  

*    ctogobs_version - Version number for ctogobs
      character*(*) ctogobs_version

      parameter ( ctogobs_version= '3.40' )

*---------------------------------------------------------- 
*     Parameter statements for ctoglb
 
*   max_cfiles      - Maximum number of cfiles which can
*                   - processed at one time (Not much use
*                   - to have it larger than the number of
*                   - recievers (max_grcv)
*   max_gprn        - Maximum PRN number expected.  (For
*                   - looking up list number (used PRN's) from
*                   - the PRN itself
*   max_max_wl_ret  - Maximum number of one-way data returns
*                     for WL calculations
*   max_max_dd_ret  - Maximum number of double difference
*                     returns
*   max_patch_poly  - Maximum order of the polyonmial to
*                     fir for patching
*   max_unflg       - Maximum number of unflagged slips
*                     allowed during patching
*   max_dd_edit     - Maximum number of data points to be
*                     edited during double difference cleaning
*                     at a sinle epoch (i.e., as data is scaned
*                     around this point.  This value must be at
*                     twice the max_max_dd_ret.
*   max_pre_edit    - Maxumim number of pre-editing entries that
*                     can be specified (These are for editing bad
*                     ranges and satellites)
*   max_scan        - Maximum number of points to be scanned at
*                     one time in the double differences.
*   num_ctog_cmds   - number of ctgobs commands (set larger than
*                     needed for easy addition of commands)
*   ctog_cmd_len    - Length of the command strings
*   num_ctog_status - Number of status indicators which can be
*                     set
*   ct_max...       - Maximum dimensions needed for mfile decoding
*   max_azmask      - Maximum number of entries in the azimuth
*                     mask per station.
 
 
      integer*4 max_cfiles, max_gprn,  max_max_wl_ret,
     .          max_max_dd_ret, max_patch_poly, max_unflg,
     .          max_dd_edit, max_pre_edit, num_ctog_cmds, 
     .          ctog_cmd_len, num_ctog_status, max_scan,
     .          ct_maxsat, ct_maxorb, ct_maxprm, ct_maxatm,
     .          max_azmask
 
      parameter ( max_cfiles      = 100 )
      parameter ( max_gprn        =  50 )
      parameter ( max_max_wl_ret  = 2000 )
      parameter ( max_max_dd_ret  = 2000 )
      parameter ( max_patch_poly  =   3 )
      parameter ( max_unflg       =  50 )
      parameter ( max_dd_edit     = 2*max_max_dd_ret )
      parameter ( max_pre_edit    = 1000 )
      parameter ( max_scan        =   3 )
      parameter ( num_ctog_cmds   =  55 )
      parameter ( ctog_cmd_len    =  16 )
      parameter ( num_ctog_status =  16 )
      parameter ( max_azmask      =  36 )
* MOD TAH 190623: Updated values to be consistent with ECOMC model.
*     THESE VALYES MUST BE KEPT CONSISTENT WITh CF_.... equivalanents.
* MOD TAH 200210: Increased from 32 to 35
* MOD TAH 201031: Increased from 35 to 45 for Beidou
      parameter ( ct_maxsat       =  45 )
      parameter ( ct_maxorb       =  22 )  ! Increased 15->22 
      parameter ( ct_maxatm       =  49 ) 
       
c    parameter ( ct_maxprm       =  13 + ct_maxsat*ct_maxorb)
* MODTAH 000214: Revised algorithm to get max number of cfile parameters
      parameter ( ct_maxprm       = 3 + 3*ct_maxatm + 
     .                          (ct_maxorb+2)*ct_maxsat + 6 )

*         max_neq       - Maximum size of normal equations to be used
*                   - includes bias parameters
 
      integer*4 max_neq
 
      parameter ( max_neq = max( (max_cfiles*max_gchannels + 
     .                      max_cfiles+max_gsvs),
     .                      max_cfiles*max_gprn ) )
* NOTE to g77 users: The above line uses a feature that will be included
* in g77 in the future (despite the development of g77 being ended).  
* We strongly recommend that you upgrade to the currently supported 
* compiler gfortran.  If g77 is still used,  replace the above lines with
C      parameter ( max_neq = max_cfiles*max_gprn )
* If parameter declarations are changed, then the above max calculation
* should be manually checked.  All other compilers will do this automatically.
 
*     Pointers to the dynamically mapped memory areas.  These are the
*     array entry points in the vma_data array.  We use malloc to
*     to get the address in memory of the first address after which
*     there is enough memory storage.  The offset between this
*     address and vma_data(1) is use to get the first array element
*     in vma_data which can be used.  The array index to store the
*     rest of the data are then computed.
 
*   iL1r_phs_cse    - Entry to store the L1 phase residuals
*                   - for all channels, stations and epochs.
*                   - (cycles @ L1) (R*8)
*   iL2r_phs_cse    - Entry to store the L2 phase residuals
*                   - for all channels, stations and epochs.
*                   - (cycles @ L2) (R*8)
 
*   iL1r_rng_cse    - Entry to store the L1 range residuals
*                   - for all channels, stations and epochs.
*                   - (cycles @ L1) (R*8)
*   iL2r_rng_cse    - Entry to store the L2 range residuals
*                   - for all channels, stations and epochs.
*                   - (cycles @ L2) (R*8)
 
*   iL1_cyc_cse     - Entry for the number of cycles to added
*                   - to the c-file L1 phase residual to "small"
*                   - one-way residuals (R*8).  For the 
*                   - proprocessor we do not try to resolve

*   iL2_cyc_cse     - Entry for the number of cycles to added
*                   - to the c-file L2 phase residual to "small"
*                   - one-way residuals (R*8)
 
*   ictol_cse       - Entry for the mapping of channels to PRN
*                   - PRN number for each site and epoch.  This
*                   - procedure saves about half the memory
*                   - of the program provided receivers have
*                   - 8 of less channels.
 
*   idata_flag_cse  - Data flags for all the data (Meaning same
*                   - as data_flag in gobs_data.h definition.
*                   - (I*4)

*   ibf_type_cse    - Bias flag type by channel, site and epoch.
*                     The definition of the bits in the bias flag
*                     type is:
*                     1  - Exisiting bias flag in cfile (32)
*                     2  - Added in process_phs due to large jump
*                     3  - Added in process_phs due ION jump
*                     4  - Flagged cap (when flag_gap option used)
*                     5  - LCDD flagged bias while scanning
*                     6  - WLOW flagged bias
*                     7  - LCDD flagged bias while DD cleaning   
 
*   iparams_cse     - These are estimated clock offsets for
*                   - all sites and satellites for all epochs
*                   - There are initally computed using the
*                   - range maeasurements, and then updated using
*                   - phase measurements.  Units on all clocks is
*                   - L1 cycles.
*   ipar_flag_cse   - Bit mapped quality indicator for the
*                   - parameter estimates.  The current settings are:
*                   - BIT  Meaning
*                   -   1  Estimate of clock is not valid (usally no
*                   -      range data available for estimate)
*                   -   2  There is a jump at this epoch.
*                   -   3  There is a millisecond jump at the this epoch
*   ipf_dphs_cse    - Start of the entries for the post-fit changes
*                     in L1 phase.  (Saved when c-files are read and
*                     used to un-correct phase data when normal points
*                     are formed.

*   iazel_cse       - Start of the azimuth/elevation storage

*   isvcL1_ce       - Satellite clock contribution removed by model.
*                     Needs to be added back to output clock value

 
      integer*8 iL1r_phs_cse, iL2r_phs_cse, iL1r_rng_cse, iL2r_rng_cse,
     .    iL1_cyc_cse, iL2_cyc_cse, ictol_cse, idata_flag_cse, 
     .    ibf_type_cse, iparams_cse, ipar_flag_cse, ipf_dphs_cse, 
     .    iazel_cse, isvcL1_ce
 
*     Variable needed for the ctogobs run.
 
*   rclock_s(max_cfiles)    - Initial estimates of the clocks
*                   - at the receivers from cfiles (L1 cycles)
*   rc_allan_sd(max_cfiles+max_gsvs) - Allan std.dev. at ~100 seconds
*                   - for each of the receivers.  Defaults
*                   - based on receiver type and may be
*                   - changed in the command file. (sec/sec)
*   apr_clk_val(max_cfiles+max_gsvs)    - Apriori estimate of
*                   - clock at each epoch (based on previous
*                   - epoch unless there is a jump) (L1 cycles)
*   init_clk_val(max_cfiles+max_gsvs)   - Clock values at the
*                     start of data span.  Initially based on
*                     rclock, but after first iteration is based
*                     on estimates. (L1 cycles)
*   apr_clk_var(max_cfiles+max_gsvs)    - Step variance of
*                   - the clocks i.e. variance of change in
*                   - clocks 1 epoch apart. (Only used for
*                   - ground stations currently)
*   apr_clk_sd(max_cfiles+max_gsvs)     - Step standard deviation
*                   - (used when sqrt of variance needed)
*   apr_clk_poly(4,max_cfiles)          - Polynomial coefficients
*                     for clock at each station.  (Units: sec,
*                     sec/sec and sec/sec**2).  Read initially
*                     from cfile type 3 record. Cubic term allowed
*                     read from extra in cfile record 2 (sec/sec**3)
*   rel_clk_wght    - Relative weight of the clock process noise
*                     to the range data.  (A larger value means
*                     the range data have more weight.)
*   rng_noise(max_cfiles) - Estimates of the range noise.  The
*                   - apriori is based on data type in reciever
*                   - and may be updated with commands
*   phs_noise(max_cfiles) - Phase noise to be assumed for each station
*                    (L1 cycles)
*   sum_rng_var(max_cfiles) - Sum of the one-range residuals
*                     squared (for computing rms)
*   norm_eq(max_neq,max_neq)    - Normal
*                   - equations for estimation of clocks
*   sol_eq(max_neq) - Solution vector for clocks
*   rng_jump_tol    - Tolerance when testing for range clock jumps
*                   - (Default 100).  May be updated with command file.
*                     Value given as multiple of step variance of
*                     the clock at the site.
*   rng_jump_min    - Minumim jump to be considered a jump (cycles)
*   rng_res_tol     - Tolerance on the quality of the range residuals
*                     (default is 100) (times the rng_noise)
*   rng_res_min     - Minimum range residual to be considered a bad
*                     data point (cycles)
*   rng_res_max     - Maximum range residual allowed (cycles)
*   reset_tol       - Tolerance on match od jump size to a 1 milli
*                   - second reset (L1 cycles)

*   dt_ion_tol      - Maximum time in seconds over which continuity
*                     of the ionospheric delay will be checked.
*   ion_rw_tol(3,max_cfiles) - Multiplier from last dion to see if 
*                    slip has occurred (1) and minium values for
*                    this site (cycles) (2) and (3) Maxiumum value
*                    allowed (cyc).
*   curr_ion(max_gsvs, max_cfiles) - Current estimates of the L1
*                 phase ionospheric delay.  With P1/P2 data, the
*                 ranges are used to get this estimate.  For C1
*                 data it is propagated from last good observation.
*   curr_dion(max_gsvs, max_cfiles) - Change in ionospheric delay
*                 from the epoch before.  Used to see how quickly
*                 ionospheric delay is changing.

*   phs_fit_tol(4)  - Tolerance on checking phase fit to clocks:
*                     1 - Mean offset based on range data
*                     2 - deviation from mean based on range data
*                     3 - Mean offset based on phase data
*                     4 - devaition from mean based on phase data
*                    (All units L1 cycles)

*   min_dtl_bias    - Minimum time duration (secs) between bias flags
*                     flags
*   min_dtr_end     - Minimum ratio of time difference between last
*                     bias flag and end of data, and the total duration
*                     of data (total duration in computed from the 
*                     number of good data).  (Ratio 0-1).

*   data_lft(4,max_max_wl_ret) - Data for L1, L2, P1 and P2 in the
*                     left segment for patching
*   data_rgh((4,max_max_wl_ret) - Data for L1, L2, P1 and P2 in the
*                     right segment for patching
*   data_wl(4,max_max_wl_ret) - Data for L1, L2, P1 and P2 in the
*                     left segment for patching: working arrays
*   data_wr((4,max_max_wl_ret) - Data for L1, L2, P1 and P2 in the
*                     right segment for patching: working arrays

*   unflg_mag(max_unflg) - Magnitudes of the unflagged slips with
*                     in one patch interval

*   curr_L1_slip, curr_L2_slip - Current estimates of the additional
*					  cycles to be added due to cycle slips

*   dchi2_ratio - Ratio that the best (modified by dchi2_min_val) chi**2
*                 to next best must satisfy in order for a bias flag
*                 be removed.
*   dchi2_min_val - Min value of the best chi**2 to use.  This is 
*                 implemented as 
*                 min_chi2+dchi2_min_val*exp(-min_chi2/dchi2_min_val)
*   dchi2_max_sep - Max separation allowed to have the bias flag
*                 removed (seconds) 
*   dchi2_gap_fact - The multiplier to be used to reduce the chi**2
*                 quality value by to account for the size of gap in
*                 the data and the amount of data being used for the 
*                 patch.

*   dd_wl_tol(3)  - Size of jump in WL which will be flagged based 
*                   on number of WL sigma. (1) Multiplier on scatter;
*                   (2) Miniumum value and (3) Maxiumum value (cyc)
*   dd_lc_tol(3)  - Size of jump on LC DD which will be flagged 
*                   as a slip (cycles).	(1) Multiplier on scatter,
*                   (2) minimum to be flagged; (3) Maximum to be
*                   flagged (cyc).
*   tol_scale     - A scaling factor applied to the jump tolerance
*                   with each iteration so that we don't keep flagging
*                   and patching the same data.

*   cf_apr_save(ct_maxprm,max_cfiles) -- Save of cf_preval arrays
*   rms_ref_clk   - RMS clock scatter needed to be considered a reference
*                   clock for satellite output (cycles)
*   rms_max_clk   - Largest RMS allowed for a clock to be output to clock
*                   file (cycles)
 
      real*8 rclock_s(max_cfiles), rc_allan_sd(max_cfiles+max_gsvs),
     .    apr_clk_val(max_cfiles+max_gsvs),
     .    init_clk_val(max_cfiles+max_gsvs),
     .    apr_clk_var(max_cfiles+max_gsvs), 
     .    apr_clk_sd(max_cfiles+max_gsvs) , apr_clk_poly(4,max_cfiles),
     .    rel_clk_wght, rng_noise(max_cfiles), phs_noise(max_cfiles), 
     .    sum_rng_var(max_cfiles),
     .    norm_eq(max_neq,max_neq),
     .    sol_eq(max_neq),
     .    rng_jump_tol, rng_jump_min, rng_res_tol, rng_res_min, 
     .    rng_res_max, reset_tol, dt_ion_tol, ion_rw_tol(3,max_cfiles),
     .    curr_ion(max_gsvs, max_cfiles), 
     .    curr_dion(max_gsvs, max_cfiles),
     .    phs_fit_tol(4), min_dtl_bias, min_dtr_end,
     .    data_lft(4,max_max_wl_ret), data_rgh(4,max_max_wl_ret),
     .    data_wl(4,max_max_wl_ret), data_wr(4,max_max_wl_ret),
     .    unflg_mag(max_unflg), curr_L1_slip, curr_L2_slip,
     .    dchi2_ratio, dchi2_min_val, dchi2_max_sep, dchi2_gap_fact,
     .    dd_wl_tol(3), dd_lc_tol(3), tol_scale, 
     .    cf_apr_save(ct_maxprm,max_cfiles), rms_ref_clk, rms_max_clk

* SAVED CLOCK INFORMATION
* save_clk(2,max_cfiles+max_gprn) -- Saved values of the clock fit (cycles
*         and cycles/sec
* save_epc(max_cfiles+max_gprn)   -- Saved epoch offset for clock (sec)
* save_rms(max_cfiles+max_gprn)   -- Save phase rms fit to linear clock (cyc)
* save_drms(max_cfiles)           -- Save rms of differences in clock
*                                    values.  Used to select new reference
*                                    sites if needed.

      real*8 save_clk(2,max_cfiles+max_gprn), 
     .       save_epc(max_cfiles+max_gprn), 
     .       save_rms(max_cfiles+max_gprn),
     .       save_drms(max_cfiles) 

*   min_ctog_elev  - Minumum elevation angle for cleaning (default 0)
*                    (rads).
*   min_cfile_elev - Minumum elevation angle in cfiles (with zero error
*                    flag (rads)
*   min_out_elev   - Minumim elevation angle to be set when writing
*                    cfiles (rads)
*   site_celev(max_cfiles) - Elevation cutoffs by site for cleaning
*   site_oelev(max_cfiles) - Elevation cutoffs by site for output
*   az_mask(2, max_azmask, max_cfiles) -- Azimuth mask pairs given
*          in form az1, el1, az2, el2, ... EL1 is the elevation 
*          cutoff between az1 and az2; el2 is the cutoff between
*          az2 and az3 etc.  


      real*4 min_ctog_elev, min_cfile_elev, min_out_elev,
     .       site_celev(max_cfiles), site_oelev(max_cfiles),
     .       az_mask(2, max_azmask, max_cfiles)
 
*   num_cfiles      - Number of cfiles to be processed.
*   num_ep          - Number of epochs for preprocessing
*   num_obs         - Number of one-way observations in these
*                     cfiles.
*   num_sat         - Number of satellites in cfile.
*   num_param       - Number of clock parameters to estimate
*                   - (num_cfiles + num_sat)
*   data_types_s(max_gdata_types, max_cfiles)   - Data types by
*                   - statoin and data tyype index (see Gobs_header.h
*                   - for definition).
*   lambda(max_gsvs, max_gdata_types, max_cfiles)  - Wavelength factors
*                   - for each channel at each station
*   prn_list(max_gsvs)  - numerical list of PRN numbers
*                   - used. 6 23 14 ..
*   prntol(max_gprn)    - A list giving the local list number
*                   - of a prn indexed by the prn. i.e.,
*                   - prn 23 = list entry 2
 
*   num_chan        - Number of channels needed
*   actual_max_chan - Actual maximum number of channels used in
*                   - any receiver in network.
*   apr_clk_epoch(max_cfiles+max_gsvs)  - The epoch at which
*                   - the last clock value was known.
*   sampling        - Sampling interval (sec)
*   orig_sampling(max_cfiles)   - Original sampling interval for each
*                   - reciever.
*   rng_mask        - Editing mask to be used with range data
*                   - (and'd with data_flag)
*   phs_mask        - Editing mask to be used with phase data.  (Bias
*                   - parameters are not checked with this mask i.e.,
*                   - bits 31 and 32 are not checked.
*   phs_bias_mask   - Editing mask for phase data with the bias
*                   - parameters checked.
*   sum_rng_num(max_cfiles)  - Number of values in computing
*                     range noise. 
*   rc_allan_num(max_cfiles+max_gsvs) - Number of values in the
*                     estimates of allan standard deviation

*   max_rclk_iter   - Maxiumum iterations to get range clock 
*                     estimates and statistics.

*   num_gdata_types(max_cfiles) - Number of data types for each cfile.
*                    (Need to save so that we write the correct value
*                     out when the cfiles are written).

*   min_good_bias  - Minumim number of good data need between bias 
*                    flags
*   min_good_end   - Minimum number of good data needed after a bias
*                    flag at the end of data.
*   num_good(max_gsvs,max_cfiles) - Number of good data in the 
*                    in the one-ways after  trimming.
*   num_deleted(max_gsvs,max_cfiles)   - Number of data deleted 
*                    during the trimming operation.
*   num_bias_flag(max_gsvs,max_cfiles) - Number of bias flags for each
*                    one way sequence.

*   dd_site_list(max_cfiles) - Ordered list of site for cleaning
*                    starting with best (lowest rnage residuals)
*   dd_svs_list(max_gsvs)    - Ordered list of satellites for
*                    cleaning starting with best (lowest clock rms
*                    first)
*   dd_site_search(max_cfiles-1) - Ordered list of sites to
*                    search when looking for double diff 
*                    combinations
*   dd_svs_search(max_gsvs-1)   -  Ordered list of satellies to 
*                    search when looking for double diff
*                    combinations

*   tol_one_way_fix  - Maximum number of epochs in cap and
*                    one way data cleaning still allowed.
*   max_wl_ret, max_dd_ret  - User selected max values to 
*                    return in left and right segments
*   max_lg_use       - Number of values to use for LG combination
*                    when patching.
*   ep_lft(max_max_wl_ret), ep_rgh(max_max_wl_ret) - Epoch numbers
*                    of the data in the left and right blocks.
*   ep_wl(max_max_wl_ret), ep_wr(max_max_wl_ret)   - Epoch numbers
*                    for the working array
*   num_lft, num_rgh  - Number of data in left and right segments
*   num_wl, num_wr    - Number of data in working left and right
*                    segments.
*   dd_ep(max_max_dd_ret,2) -Epoch numbers in left and right segments
*                    for double differences
*   dd_num(2)        - Number of epochs in dd left and right segments.

*   unflg_num         - Number of unflagged slips in one patching
*                       cycle
*   unflg_ep(2,max_unflg) - Epoch numbers with suspected unflagged
*                       slips and there companion from which the jump
*                       was computed.
*   unflg_sv(2,max_unflg) - Saves the second site and satellite when
*                       an unflagged slip is found in the double 
*                       differnces.  I needed all one-ways will be
*                       searched to find the jump.  

*   curr_ion_ep(max_gsvs, max_cfiles) - The epoch number to which the
*               curr_ion values apply (Range ion will be used differnce
*               from current epoch is >~ 1hour.


*   dd_edit_num       - Number of data pointed edited
*   dd_edit_ep(max_dd_edit) - Epochs of the edited data.

*   status_rep     - Bit mapped word for status bits to be
*                       set for reporting program status

*   num_pre_edit      - Number of pre-editing entries given
*   pre_edit(4, max_pre_edit) - List of pre-edit entries.  The
*                       ordering is:
*                     (1) Site number (0 for all sites)
*                     (2) SVS PRN number (0 for all Satellites)
*                     (3) Start Epoch number
*                     (4) End epoch number

*   scan_sites((max_cfiles-1)/32+1) - Bit mapped list of sites
*                      to be scanned.
*   num_scan           - Number of observations in current scan
*                        list
*   scan_ep(max_scan)  - Epochs in the current scan list

*   phs_res_sites((max_cfiles-1)/32+1) - Bit mapped word giving
*                      the sites to output in the residual files
*                      and the single differece files.

*   dd_scan_sitess((max_cfiles-1)/32+1) - Bit mapped word giving
*                      the sites to scanned in double diffs. before
*                      cleaning.

*   gap_size(max_cfiles) - Size of the gap allowed when gaps are
*           flagged before procesing.  (Not used if the gaps are
*           not flagged.)

*   uns   - Unit number for summary file (If it is a file the unit
*           number is 203; if it is the screen then the unit number
*           is 6.

*   max_scan_edit - Maximum number of bias flags that can be added
*           during scanning before all data are edited.
*   num_dd_flags(max_cfiles, max_gprn)  - Number of dd bias flags
*           added during scanning.  If too many are added then
*           complete site/svs sequence is removed. 
*   min_ow_data - Minimum amount of ow data needed for site/sv to
*           be retained.  (Set with max_scan_edit command)

*   np_size       - Size of normal points in epochs (all data must
*           be present with no bias flags)
*   np_start      - First epoch to start normal pointing from.

*   site_snr(2,max_cfiles) - Minimum SNR to allow at L1 and L2.
*   edit_counts(32,max_cfiles) - Counts of the number of edited
*                   data by site.

*   lcg_pol(2)  -  Order of polynomial to fit to LC data and LG data
*   num_azmask(max_cfiles) -- Number of entries in the azimuth mask
*         for each station.
*   num_ref_clk - Number of clocks used to define the reference
*                 clock series (either computed or specified by user)
*   namb          ! Entry in ambiquity tables (computed from site and satellite 
*                 ! number)
*   amb_tab(5,max_cfiles*max_gprn)  ! Ambiquity table.  Entries in table are:
*                 ! 1 -- Site 1 (starts as ref_site)
*                 ! 2 -- Site 2 (Second site in baseline)
*                 ! 3 -- SVS  1 (starts as Ref_svs)
*                 ! 4 -- SVS  2 (second satellite in DDiff)
*                 ! 5 -- Status: 1 if amb OK, -1 if reference site/svs
*                 !      -2 No data on satellite

 
      integer*4 num_cfiles, num_ep, num_obs, num_sat, num_param,
     .    data_types_s(max_gdata_types, max_cfiles),
     .    lambda(max_gsvs, max_gdata_types, max_cfiles), 
     .    prn_list(max_gsvs),
     .    prntol(max_gprn), num_chan, actual_max_chan,
     .    apr_clk_epoch(max_cfiles+max_gsvs), sampling,
     .    orig_sampling(max_cfiles), rng_mask, phs_mask, 
     .    phs_bias_mask, sum_rng_num(max_cfiles),
     .    rc_allan_num(max_cfiles+max_gsvs), max_rclk_iter, 
     .    num_gdata_types(max_cfiles),
     .    min_good_bias, min_good_end,
     .    num_good(max_gsvs,max_cfiles), 
     .    num_deleted(max_gsvs,max_cfiles),
     .    num_bias_flag(max_gsvs,max_cfiles),
     .    dd_site_list(max_cfiles), dd_svs_list(max_gsvs),
     .    dd_site_search(max_cfiles-1), 
     .    dd_svs_search(max_gsvs-1),
     .    tol_one_way_fix, max_wl_ret, max_dd_ret, max_lg_use,
     .    ep_lft(max_max_wl_ret), ep_rgh(max_max_wl_ret),
     .    ep_wl(max_max_wl_ret), ep_wr(max_max_wl_ret),
     .    num_lft, num_rgh, num_wl, num_wr,
     .    dd_ep(max_max_dd_ret,2), dd_num(2),
     .    unflg_num, unflg_ep(2,max_unflg), unflg_sv(2,max_unflg),
     .    curr_ion_ep(max_gsvs, max_cfiles),
     .    dd_edit_num, dd_edit_ep(max_dd_edit), status_rep,
     .    num_pre_edit, pre_edit(4, max_pre_edit),
     .    phs_res_sites((max_cfiles-1)/32+1),
     .    scan_sites((max_cfiles-1)/32+1), num_scan, scan_ep(max_scan),
     .    gap_size(max_cfiles), uns,max_scan_edit, min_ow_data,
     .    num_dd_flags(max_cfiles, max_gprn),
     .    np_size, np_start, site_snr(2,max_cfiles),
     .    edit_counts(32,max_cfiles) , lcg_po(2), 
     .    num_azmask(max_cfiles), num_ref_clk, 
     .    namb,amb_tab(5,max_cfiles*max_gprn) 

 
*   use_command_file    - Indicates that we should use the command
*                   - file.
*   out_dd_file         - Logical to indicate that we should output
*                     double difference status file
*   out_sd_file         - Logical to indicate that we should output
*                     single difference files.
*   use_gamit_elc       - Indicates that GAMIT elevation cutoff should
*                     be retained during cleaning.  This will limit the
*                     amount of data that can be used later.
*   use_cview_edit      - Indicates that the cview edit flag -1 should
*                     be used (i.e., adopt previously cleaned data)
*   use_MM_ranges   - Logical set true (by user) if they want to use
*                     the MiniMac range data (Generally this data is
*                     pretty bad and can be ignored since the clocks 
*                     are usually good).
*   ignore_gaps         - Indicates that gaps should be ignored when
*                     scanning the double differences.  Should only be
*                     set for clean data when the gamit elevation limit
*                     and cview editing flags are used. 
*   usr_ignore_gaps - User selection of ignore gaps.  May not be used if
*                     user selects gaps to be flagged (flag_gaps)
*   allow_one_bg    - Set to allow one bias or gap in the second satellite
*                     of the single differences.  This should allow gap
*                     jump on all satellites to be repaired.
*   do_one_bg       - User control to say that one bias or gap should be
*                     be allowed (command ALLOW_ONE_BG)
*   remove_ms_jump  - Option to remove millisec jumps in the clock before
*                     writing new cfiles
*   gaps_flagged    - Indicates that gaps have had bias flags added and 
*                     therefore can be ignored.

*   remove_first_bias - Logical to indicate that the first bias should
*                     be removed when cfiles are written

*   dd_bias_added   - Set true if bias flag added during double difference
*                     cleaning (causes cleaning loop to be repeated.)
*   unflg_fnd       - Logical set true to indicate that jump was found
*                     in fit_obs.  (initialzed in inc_cyc_es and should
*                     checked after each call so that the second site
*                     and satellite can be saved.

*   svs_clk_known(max_gsvs) - Indicates that the satellite clock
*                   - values are approximately known.  (Based on
*                   - using satellite once.)

*   kin_as_sep_sites - Set true if we wish to treat kinematic sites as
*                     separate sites (instead of as stochastic site)

*   lcg_force(2)    - Set true if the polynomial coefficients of LC and
*                     LG are to be forced equal.

*   apply_phs_clk   - Set true if the phase clocks are to be used to
*                     update the one-way phase data.

*   pf_remove_bf    - Set true is bias flags are tested and removed after
*                     phase clocks have been computed. (Default is false).
*   write_igs_clk   - Logical set true if we are to write the igs clock file
*   igs_clk_samp    - Sampling interval for IGS clock output
*   use_orig_bf     - Set true to use any original bias flags in the data.
*                     With version 3.15, the default is not use the original
*                     ones since many of these seem bogus.
*   acbias_out      - Logical set if acbias.dat file is to be written
*   nol1only        - set true with NOL1ONLY command and stop L1 data only being
*                     processed.
*   prefit_clk      - Set true to prefit linear trends to the station clocks
*                     making autcln independent of I-file
*   prescan_ms      - Set true to scan for millisecond clock jumps before 
*                     processing. New version 3.39. TAH 2020509.

* MOD TAH 180307: Ver 3.35 Added app_ion logical
*   app_ion         - Set to apply ionospheric delay to output omc values
* MOH TAH 200617: Ver 3.40 Added trim_seg
*   trim_seg        - Set true to remove short segments of data between gaps
*                     on first cleaning pass (Default true, set false with
*                     trim_oneway_tol command.
 
      logical use_command_file, out_dd_file, out_sd_file,
     .        use_gamit_elc, use_cview_edit, ignore_gaps,
     .        usr_ignore_gaps, gaps_flagged,  remove_first_bias, 
     .        dd_bias_added, unflg_fnd,
     .        svs_clk_known(max_gsvs), use_MM_ranges, allow_one_bg,
     .        do_one_bg, remove_ms_jump, kin_as_sep_sites ,
     .        lcg_force(2), apply_phs_clk, pf_remove_bf, write_igs_clk,
     .        use_orig_bf, acbias_out, nol1only,  prefit_clk, 
     .        prescan_ms, app_ion, trim_seg

      logical remap_glonass ! Set true (default) to map glonass FDMA channel
              ! ambiguities to non-integer values (ser with remap_glonass <Y/N>
              ! command.

      integer*4 igs_clk_samp

*   caprog_name   - Name of the program being run.  Either ctogobs or
*                   autcln

      character*16 caprog_name
 
*   cf_codes(max_cfiles)    - Site codes from the cfile names
*   rcv_types(max_cfiles)   - The receiver types at each of the
*                   - stations.
*   new_cfile_series        - New series letter for updated cfiles.
*                   - passed in the runstring. '+' will cause
*                   - the series to incremented by one letter or
*                   - goto a 'a' if current series is numeric.
*                    (Four characeters are allocated even though
*                     only one is used because of bugs in Sun OS with
*                     single character strings.
*   curr_cfile_series - Series for the cfiles being read in.
*   dd_out_opts     - Options for output double difference combinations
*                     during cleaning.  Options are:
*                     ALL  -- All combinations output as they are tried.
*                     FIXE -- Bias flag removed combinations only
*                     NOTF -- Bias flag not removed combinations only.
*   unflg_type(max_unflg) - Type of observable used to deduce an unflaged
*                     bias

*   edit_names(32)        - Codes for the types of edits

*   ref_clk_code(max_cfiles) - List of sites used to define to the 
*                     reference clock

      character*4 cf_codes(max_cfiles), rcv_types(max_cfiles),
     .            new_cfile_series, curr_cfile_series, dd_out_opts,
     .            unflg_type(max_unflg), edit_names(32),
     .            ref_clk_code(max_cfiles)

*   long_names(max_cfiles) -- Long names for the sites (should contain
*                 DOMES number in first 20 characters)

      character*32 long_names(max_cfiles)
 
*   ctogobs_title   - Description of expermiment passed in gtoobs
*                   - rustring.  This may be modified before being
*                   - saved in the Gobs_file.
      character*64 ctogobs_title

*   rng_clk_root    - Root to file name for output of clock errors
*                     derived from range data 
*   phs_clk_root    - Root to file name for output of clock errors
*                     derived from phase data
*   phs_res_root    - Root to file name for output of phase residuals
*   sng_diff_root   - Root to file name for output of single difference
*                     files.  (Single differences are from first site
*                     is list of c-files)

      character*32 rng_clk_root, phs_clk_root, phs_res_root, 
     .             sng_diff_root

 
*   CtoGobs_com_file - Command file for ctoobs.
*   Gobs_code       - Either the name or a single character passed
*                     in the runstring which is used to generate the
*                     gobs_outfile name
*   Gobs_outfile    - Name of the Gobs file.  May be modfiied from
*                   - name passed if one of the short forms of the
*                   - name is used.
*   cfiles(max_cfiles)  - List of cfiles to be used.
*   dd_outfile      - Name of file to save the double difference
*                     combinations used in fixing biases (same format
*                     as dd.srt file).  Output to this file is controlled
*                     by dd_out_opts.
*   summary_file    - Name of the file for the output of a summary of
*                     the results.  May be 6 for output to screen.
*   igs_clk_file    - Name of the igs clock file to output
 
 
      character*128 CtoGobs_com_file, Gobs_code, Gobs_outfile, 
     .              cfiles(max_cfiles), dd_outfile, summary_file,
     .              igs_clk_file

*   lca_type - Algorithm for lc_autcln. SEQ is original sequential
*     method, DIR is ver 3.29 direct estimation method.

      character*4   lca_type  
 
*    dph_output     - Type of DPH file to output
*                     NONE - no DPH file written
*                     SINGLE - DPH file output by site only
*                     MULTIPLE - DPH out by site and by satellite
 
      character*8   dph_output 

* MOD TAH 200628: Added saved GNSS from cf_gnss values (assumed all the 
*     same)
      character*1 sv_gnss  ! Saved GNSS system

*...................................................................

      common / ctogobs_mem / iL1r_phs_cse, iL2r_phs_cse, iL1r_rng_cse,
     .    iL2r_rng_cse, iL1_cyc_cse, iL2_cyc_cse, ictol_cse,
     .    idata_flag_cse, ibf_type_cse, iparams_cse, ipar_flag_cse,
     .    ipf_dphs_cse, iazel_cse, isvcL1_ce

      common / ctog_mat / norm_eq, sol_eq  

      common / ctogobs_com / rclock_s, rc_allan_sd,
     .    apr_clk_val, init_clk_val, apr_clk_var, apr_clk_sd, 
     .    apr_clk_poly, rel_clk_wght, sum_rng_var,
     .    rng_noise, phs_noise, 
     .    rng_jump_tol, rng_jump_min, rng_res_tol, rng_res_min, 
     .    rng_res_max,
     .    reset_tol, dt_ion_tol, ion_rw_tol, curr_ion, curr_dion,
     .    phs_fit_tol, min_dtl_bias, min_dtr_end, 
     .    data_lft, data_rgh, data_wl, data_wr,
     .    unflg_mag, curr_L1_slip, curr_L2_slip,
     .    dchi2_ratio, dchi2_min_val, dchi2_max_sep, dchi2_gap_fact,
     .    dd_wl_tol, dd_lc_tol, tol_scale, cf_apr_save, rms_ref_clk,
     .    rms_max_clk, save_clk,  save_epc,  save_rms, save_drms,
     .    min_ctog_elev, min_cfile_elev,
     .    min_out_elev, site_celev, site_oelev, az_mask,
     .    num_cfiles, num_ep, num_obs, num_sat,
     .    num_param, data_types_s, lambda, prn_list, prntol, num_chan,
     .    actual_max_chan, apr_clk_epoch, sampling, 
     .    orig_sampling, rng_mask, phs_mask, phs_bias_mask, 
     .    sum_rng_num, rc_allan_num, max_rclk_iter, num_gdata_types,
     .    min_good_bias, min_good_end, num_good, num_deleted, 
     .    num_bias_flag,
     .    dd_site_list, dd_svs_list, dd_site_search, dd_svs_search,
     .    tol_one_way_fix, 
     .    max_wl_ret, max_dd_ret, max_lg_use, ep_lft, ep_rgh,
     .    ep_wl, ep_wr, num_lft, num_rgh, num_wl, num_wr,
     .    dd_ep, dd_num, unflg_num, unflg_ep, unflg_sv, curr_ion_ep, 
     .    dd_edit_num, dd_edit_ep, status_rep,
     .    num_pre_edit, pre_edit, phs_res_sites,
     .    scan_sites, num_scan, scan_ep, gap_size, uns,
     .    max_scan_edit, min_ow_data, num_dd_flags, np_size, np_start,
     .    site_snr, edit_counts, lcg_po, num_azmask, num_ref_clk,
     .    namb, amb_tab, 
     .    use_command_file, out_dd_file, out_sd_file, 
     .    use_gamit_elc, use_cview_edit, ignore_gaps,
     .    usr_ignore_gaps, gaps_flagged,  remove_first_bias,  
     .    dd_bias_added, unflg_fnd, use_MM_ranges, allow_one_bg, 
     .    do_one_bg, remove_ms_jump, kin_as_sep_sites, lcg_force,
     .    apply_phs_clk, pf_remove_bf, write_igs_clk, use_orig_bf,
     .    acbias_out, nol1only, prefit_clk, prescan_ms, 
     .    app_ion, trim_seg, remap_glonass, igs_clk_samp,  
     .    svs_clk_known, caprog_name, cf_codes, rcv_types, 
     .    new_cfile_series, curr_cfile_series, dd_out_opts, unflg_type,
     .    edit_names, ref_clk_code,
     .    rng_clk_root, phs_clk_root, phs_res_root, sng_diff_root,
     .    long_names,
     .    ctogobs_title, CtoGobs_com_file, Gobs_code, Gobs_outfile, 
     .    cfiles, dd_outfile, summary_file, igs_clk_file, lca_type,
     .    dph_output, sv_gnss 
 
*...................................................................

*   ctog_commands(num_ctog_cmds) - List of commands for ctogobs
*   status_rep_opts(num_ctog_status) - Status report options.  Allows
*         user to set the output detail for the program. Options are
*         given in autcln.hlp file and ctog_cmds_bd.f
*          

      character*(ctog_cmd_len) ctog_commands(num_ctog_cmds),
     .                         status_rep_opts(num_ctog_status)

      common / ctog_cmds_com / ctog_commands, status_rep_opts

*-------------------------------------------------------------------

*     Block of include file for processing postfit residuals
*

*   WL_bias_site(max_cfiles) - Average value of one-way widelane by
*                 site
*   WL_bias_svs(max_gsvs) - Average value of one-way widelane by
*                 satellite
*   WL_RMS(max_cfiles)  - RMS scatters of widelanes

*   LC_bias(max_cfiles) - Average LC one-way phase residuals by
*                 site
*   LC_RMS(max_cfiles)  - RMS scatter of LC residuals by site.
*   ALL_RMS, ALL_PRV    - Over all RMS value of residuals, and
*                         all_prv RMS
*   PR_RMS(max_cfiles)  - RMS scatter from previous iteration

*   LC_svs_bias(max_gsvs, max_cfiles) - Average LC one-way phase 
*                 residuals by site and satellite
*   LC_svs_RMS(max_gsvs, max_cfiles)  - RMS scatter of LC residuals 
*                 by site and satellite 

*   LC_svs_ams(max_gsvs, max_cfiles)  - RMS scatter of 25-point averages
*                 of the LC residuals by site and satellite
*   LC_AMS(max_cfiles)  -- RMS scatter of the averages over all sites

*   LC_elv_bias(18, max_cfiles) - Average LC one-way phase 
*                 residuals by site and 5 degree elevation bins
*   LC_elv_RMS(18, max_cfiles)  - RMS scatter of LC residuals 
*                 by site and 5 degreee elevation bins
*   pc_convd   - change in postfit rms to use to decide if
*                iteration has converged.  (Used is old/new-1 <
*                pc_convd
*   pf_svs_ratio - Ratio of satellite RMS to overall RMS at a station
*                to decide if data should be edited.
*   pf_nsig      - Postfit ratio of residual to sigma for point
*                to be edited.
*   pf_maxres    - Maximum size of residual to restore (cycles).
*   pf_bad_ratio - Ratio of bad-to-good data on a LC sequence of
*                data before the whole sequence is deleted (sequence
*                is data between bias flags).
*   pc_over_shoot     - Multiplier for adjusting the means in phase
*                       clock iteration
*   pc_max_corr  - Maximum correlation to consider when solving for 
*                bias parameters in phase clock code
*   pf_max_rms   - Largest RMS allowed in postfit residuals (cyc)
*   pf_zdep(2,max_cfiles) ! Elevation angle dependent model coefficients
*                  for noise.  The squared values are saved in cycles**2.

      real*8 WL_bias_site(max_cfiles), WL_bias_svs(max_gsvs),
     .       WL_RMS(max_cfiles),
     .       LC_bias(max_cfiles), LC_RMS(max_cfiles),
     .       ALL_RMS, ALL_PRV,
     .       PR_RMS(max_cfiles),
     .       LC_svs_bias(max_gsvs,max_cfiles), 
     .       LC_svs_RMS(max_gsvs,max_cfiles),
     .       LC_svs_AMS(max_gsvs,max_cfiles), LC_AMS(max_cfiles),
     .       LC_elv_bias(18,max_cfiles), LC_elv_RMS(18,max_cfiles),
     .       pc_convd, pf_svs_ratio, pf_nsig, pf_maxres, pf_bad_ratio,
     .       pc_over_shoot, pc_max_corr, pf_max_rms,
     .       pf_zdep(2,max_cfiles)

*   WL_num(max_cfiles) - Number of data in WL statistics
*   LC_num(max_cfiles) - Number of data in LC statistics
*   PR_num(max_cfiles) - Number of data in LC statistics from previous
*                        iteration.
*   LC_svs_num(max_gsvs, max_cfiles) - Number of data in site
*            and satellite statistics

*   LC_svs_anm((max_gsvs, max_cfiles) - Number of averaged values
*            in RMS calculation by site and satellite
*   LC_anm(max_cfiles) -- Number of values in averaged RMS by site
*   LC_elv_num(18, max_cfiles) - Number of data in site
*            and elevation bin statistics 
*   pc_max_iter       - Maximum of iterations to be used
*            in generating postfit residuals.
*   pc_num_iter       - Actual number of iterations needed.
*   pc_non_int        - Number of phase clock iterations before
*                       we start using non-integer adjustments to
*                       biases
*   pc_full_anal      - Number of (slower) interations to make with
*                       full solution.
*   pc_start_edit     - Number of iteration of phase clock calcs
*                       before editing starts.  (needs a few iterations
*                       to settle down).

      integer*4 WL_num(max_cfiles), LC_num(max_cfiles),
     .          PR_NUM(max_cfiles),
     .          LC_svs_num(max_gsvs, max_cfiles),
     .          LC_svs_anm(max_gsvs, max_cfiles), LC_anm(max_cfiles), 
     .          LC_elv_num(18, max_cfiles), pc_max_iter, pc_num_iter,
     .          pc_non_int , pc_full_anal, pc_start_edit

 
*   site_np     - Parameter number for first site coordinate
*               - at current site
*   atm_np      - Parameter number for first atmospheric
*               - zenith delay at current site
*   acon_np     - Parmaeter for constant part of atm delay
*   grad_np(2)  - Parameter number for first gradient
*               - parameter delay at current site
*   site_clk_np - Parameter number for first site clock at
*               - current site
*   svs_np      - Parameter number of first orbital element
*               - (offset to current PRN is computed)
*   sof_np      - Parameter number for svs antenna offset
*   svs_clk_np  - Parameter number for first satellite clock (no used
*                 after version 940 of m-file).
*   eop_np      - Parameter number of EOP parameters.
  
      integer*4 site_np, atm_np, acon_np, grad_np(2), site_clk_np, 
     .    svs_np,sof_np, svs_clk_np, eop_np

*   use_postfit  - Logical set true to indicate that postfit
*                  residuals should be used
*   edit_postfit - Set true if we are to edit the residuals.

      logical use_postfit, edit_postfit 
      
*   mf_name      - Name of the mfile
*   df_name      - dfile name (passed in runstring)
*   mf_root      - Root name to mfile read from command line 

      character*128 mf_name, mf_root, df_name    

      common / comp_com / WL_bias_site, WL_bias_svs,
     .       WL_RMS, LC_bias, LC_RMS, PR_RMS, ALL_RMS, ALL_PRV,
     .       LC_svs_bias, LC_svs_RMS, LC_svs_AMS, LC_AMS,
     .       LC_elv_bias, LC_elv_RMS,
     .       pc_convd, pf_svs_ratio, pf_nsig, pf_maxres, pf_bad_ratio,
     .       pc_over_shoot, pc_max_corr, pf_max_rms, pf_zdep, 
     .       WL_num, LC_num, PR_num, LC_svs_num, LC_svs_anm, LC_anm,
     .       LC_elv_num, 
     .       pc_max_iter, pc_num_iter, pc_non_int, 
     .       pc_full_anal, pc_start_edit,
     .       site_np, atm_np, acon_np, grad_np, site_clk_np, 
     .       svs_np, sof_np, svs_clk_np, eop_np, 
     .       use_postfit, edit_postfit,
     .       mf_name, mf_root, df_name

 
*     Additional variables  needed to use the new bias parameter  estimation
*     in autcln.
 
*          num_allp - Total number of parameters  being
*                   - estimated (includes  clocks and bias
*                   -parameters)
*   sp_bp(2, max_cfiles*max_gsvs)        - Site and satellite  list
*                   - entry by parameter number
*   bp_sp(max_cfiles, max_gsvs)      - parameter number assigned
*                   - to each possible bias flag.
*   ep_bp(2, max_cfiles, max_gsvs)   - Start and last epoch
*                   - number used for each bias flag, saved by
*                   - cfile and satellite.
 
*   elim_bp(4,max_cfiles)   - Bit mapped array that markes which bias
*                   - parameters are to be eliminated.
*   elim_np(max_cfiles, max_gsvs)    - Parameter numbers of the
*                   - bias parameters to eliminated
 
      integer*4 num_allp, sp_bp(2, max_cfiles*max_gsvs),
     .    bp_sp(max_cfiles, max_gsvs), ep_bp(2, max_cfiles, max_gsvs),
     .    elim_bp(4,max_cfiles), elim_np(max_cfiles, max_gsvs)
 
*       bp_elim_needed  - Logical to denote that we need to eliminate
*                   - some bias parameters
 
      logical bp_elim_needed
 
      common / est_bpcom / num_allp, sp_bp, bp_sp, ep_bp, 
     .     elim_bp, elim_np, bp_elim_needed
 
*...................................................................
 
*
***** Entries needed for the flat_dd sequence
*
* max_bf_table  -- Length of bias flag table 

      integer*4 max_bf_table

      parameter ( max_bf_table = 10*max_cfiles*max_gprn)

****  Entries needed

* fdd_max_mean(2)  -- Mean value of double differences (L1 or LC cycles depending
*                  on L2 wavelength factor) and MWWL
* fdd_max_rms(2)   -- RMS values of double differences

      real*8 fdd_max_mean(2), fdd_max_rms(2)

* bf_table(5, max_bf_table)  -- Table of all bias flags.  The four entries
*    are:
*    (1,*) -- Epoch of the bias flag
*    (2,*) -- Site number
*    (3,*) -- SVS numbe
*    (4,*) -- Duration of bias flag (i.e., number of epochs between bias flag
*             and last good data point before next bias flag).  Once the
*             bias has been resolved; this duration is set to a negative value.
*    (5,*) -- Reports status of bias.  Bit mapped
*             1 -- WL is resolved to integer value
*             2 -- DD Combination found with fixed biases that can be tested
*             3 -- DD combination found with fixable biases (i.e., combination
*                  did involve biases that could be fixed although they were
*                  not. 
*             4 -- Autcln thinks it can resolve the NL bias as well.
*             5 -- Bias flag has been used in a double difference (needed to 
*                  stop redundant double differences being output.
* bf_index( max_cfiles ) -- Pointer to the first bias flag in the bf_table for
*             each station. 
* tot_bf   -- Total number of bias flags

* ref_ent  -- Entry number in the bf_table
* ref_site -- Reference site for the longest sequence
* ref_svs  -- Reference satellite for the longest sequence
* ref_ep   -- Epoch of the reference bias flag.

* fdd_L2_fact -- L2 wavelength factor used in flat double differences.  If this
*                value is set to zero then only L1 will be in fitting phase
*                clocks (even for dual frequency data).  This is the only way
*                getting flat L1 residuals.
* fdd_max_dd  -- Maxiumum number of double differences found

      integer*4 bf_table(5, max_bf_table), bf_index(max_cfiles), tot_bf,
     .          ref_ent, ref_site, ref_svs, ref_ep, fdd_L2_fact, 
     .          fdd_max_dd


      common / fdd_com / fdd_max_mean, fdd_max_rms,
     .                   bf_table, bf_index, tot_bf,
     .                   ref_ent, ref_site, ref_svs, ref_ep, 
     .                   fdd_L2_fact, fdd_max_dd

*     Additional variables needed for widelane ambiquitity
*     resolutions

      integer*4 wl_ref_site  ! Reference site number for DD WL
                       ! ambiguities
     .,   wl_ref_svs   ! reference satellite list number for DD WL
                       ! ambiguities
     .,   wl_ref_dur   ! Duration in epochs of WL reference
     .,   wl_ref_bf    ! bf_table entry for WL Reference

     .,   min_wl_tol   ! Minimum number of values in wide lane for
                       ! it to be resolved.
     .,   ambnum_mw(max_cfiles*max_gprn) ! Number of values in MW-WL
                       ! estimates in the DIR algorithm.
     .,   site_dcb(max_cfiles) ! DCB type by station (0, 1 and 2) (set in EST_COMB)
                               ! Type of DCB by station based
                               ! on bit 28 and 29 of data flag

      logical resolve_wl   ! Logical set true if LC_AUTCLN command used
                       ! resolve the widelane ambiguities
c
c  orderlen((max_cfiles*(max_cfiles+1))/2)  ! Ordered list of
c                        ! lengths which is modified by the sigma of the
c                        ! widelanes (larger sigma wide lanes are moved down
c                        ! the list).


      real*8 baselens((max_cfiles*(max_cfiles+1))/2)
     .,      orderlen((max_cfiles*(max_cfiles+1))/2) 
     .,   wl_conf(max_bf_table)  ! WL confidence for each bias
     .,   nl_conf(max_bf_table)  ! NL confidence for each bias
     .,   dchi_wl_tol            ! Tolerance for fixing wide lane
     .,   msig_wl_tol            ! Maximum sigma for mean widelane
     .,   mdev_wl_tol            ! Maximum deviation for mean widelane

     .,   ambest_lc(max_cfiles*max_gprn) ! LC estimate of ambiguity
     .,   ambest_ex(max_cfiles*max_gprn) ! EX-WL estimate of ambiguity (L2-L1)
     .,   ambest_mw(max_cfiles*max_gprn) ! MW-WL estimate of ambiguity (L2-L1)
     .,   ambsig_lc(max_cfiles*max_gprn) ! LC sigma of estimate of ambiguity
     .,   ambsig_ex(max_cfiles*max_gprn) ! EX-WL sigma of estimate of ambiguity (L2-L1)
     .,   ambsig_mw(max_cfiles*max_gprn) ! MW-WL sigma of estimate of ambiguity (L2-L1)

     .,   mwbias_site(max_cfiles)   ! Computed DCB biases in MW-WL by station
     .,   mwbias_svs(max_gprn,3)    ! Computed DCB biases by satellite for the three
                                    ! types of DCB biases (0, 1 and 2)
      real*8 Cov_owex(max_cfiles*max_gprn,max_cfiles*max_gprn)  ! Covariance matrix 
                                    ! for one-way EW-WL (needed for double diff sigma)

      common / rwl_i4 / wl_ref_site,  wl_ref_svs, wl_ref_dur
     .,      wl_ref_bf, min_wl_tol, ambnum_mw, site_dcb, resolve_wl
     .,      Cov_owex

      common / rwl_r8 / baselens, orderlen, wl_conf, nl_conf 
     .,      dchi_wl_tol, msig_wl_tol, mdev_wl_tol  
     .,      ambest_lc, ambest_ex, ambest_mw, ambsig_lc, ambsig_ex 
     .,      ambsig_mw, mwbias_site, mwbias_svs 


*------------------------------------------------------------------------------

* Frequency quantities now SV dependent (for GLONASS) and set dynamically in ctogobs/set_freqs.f
   
* fl1(ct_maxsat) -- L1 frequency for constellation from C-file (SV-specific only if GLONASS)
* fl2(ct_maxsat) -- L2 frequency for constellation from C-file (SV-specific only if GLONASS)
*    additional frequencies to be added when autcln and solve can handle more than 2
* dfsf(ct_maxsat) -- Diff/sum ratio for WL and NL   dfsf = (fL1-fL2)/(fL1+fL2)
* sfdf(ct_maxsat) -- Sum/Diff ratio for WL and NL   sfdf = (fL1+fL2)/(fL1-fL2)
* lcf1(ct_maxsat) -- Multiplier for LC  from L1     lcf1 = 1.d0/(1.d0 - (fL2/fL1)**2) )
* lcf2(ct_maxsat) -- Multiplier for LC  from L2     lcf2 = -(fL2/fL1)/(1.d0 - (fL2/fL1)**2)     
* lgf1(ct_maxsat)- - Multiplier for LG from L1      lgf1 = -fL2/fL1
* lgf2(ct_maxsat)- - Multiplier for LG from L21     lgf2 = 1.d0 
* exf1(ct_maxsat) -- Multiplier for EX-WL from L1   exf1 = 1.d0 
* exf1(ct_maxsat) -- Multiplier for EX-WL from L2   exf2 = -fL1/fL2
* pcf1(ct_maxsat) -- Multiplier for PC from P1      fL1**2/(fL1**2-fL2**2)
* pcf2(ct_maxsat) -- Multiplier for PC from P2     -fL2**2/(fL1**2-fL2**2)   
* fClk            -- Scalar L1 frequency used for clock, = center freq for Glonass, fL1(1) for others
* fL1u(ct_maxsat) -- L1 frequency as originally in the c-files (unmapped)
* fL2u(ct_maxsat) -- L2 frequency as originally in the c-files (unmapped)

      real*8 fL1(ct_maxsat) ,fL2(ct_maxsat)
     .     , fL1u(ct_maxsat),fL2u(ct_maxsat) 
     .     , dfsf(ct_maxsat),sfdf(ct_maxsat)
     .     , lcf1(ct_maxsat),lcf2(ct_maxsat)
     .     , lgf1(ct_maxsat),lgf2(ct_maxsat)
     .     , exf1(ct_maxsat),exf2(ct_maxsat)
     .     , pcf1(ct_maxsat),pcf2(ct_maxsat)                                                         
     .     , fClk 
                                              
* ION conversions
*  l1tecu - Change in range (m) per TECU at L1   l1tecu = 40.3d0/fL1**2*1.d16
*  l2tecu - Change in range (m) per TECU at L2   l2tecu = 40.3d0/fL2**2*1.d16
      real*8 l1tecu(ct_maxsat), l2tecu(ct_maxsat)

* ADDED TAH 180308: Values to compute ionospheric delay (I1 L1 ion delay in cycles)
* lif1(ct_maxsat) -- Multiplier for I1  from L1     lif1 =  1.0/(1.0 - (fL1/fL2)**2)
* lif2(ct_maxsat) -- Multiplier for I1  from L2     lif1 =  -(fL1/fL2)/(1-(fL1/fL2)**2)
* fL1_R0, fL2_R0, fL3_R0 -- Refeence frequencies for GLONASS to map values to a common
*                    frequency 
      real*8 lif1(ct_maxsat), lif2(ct_maxsat)
      real*8 fL1_R0, fL2_R0, fL3_R0
                                                                         
      common/ freq_com / fl1,fl2, fL1u, fL2u, dfsf,sfdf,lcf1,lcf2
     .   , lgf1,lgf2, exf1,exf2,pcf1,pcf2,fClk,l1tecu,l2tecu, lif1, lif2
     .   , fL1_R0, fL2_R0, fL3_R0


