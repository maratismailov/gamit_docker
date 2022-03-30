*
*  Include for declaration for variables needed by the track GPS
*  kinematic GPS analysis program.
*
      include '../includes/sp3_def.h'

* PARAMETERS
*
* max_site -- Maximum number of sites that can be processed.
* max_kine -- Maximum number of kinematic sites that can be processed at
*             at one time
* max_prn  -- Maximum number of prns allowed
* max_chan -- Maximum number of satellites that can be observed at 
*             one time.
* max_epochs -- Maximum number of epochs that can be processed
* max_ep_wrds -- Maximum number of 32-bit words needed for 1 bit per epoch
* max_ambs   -- Maximum number of ambiquities allowed (each measurement
*               points to a particular ambiquity value
* max_parm   -- Maximum number of parameters that can be estimated
*               simultaneously.  There is 1 zenith delay per-site
*               and 3 coordinates per kinematic site
* max_dtype  -- Maximum number of obervable types that can be used
*               at one time 
* max_obs    -- Maximum number of observations at one time
* max_wlrng  -- Maximum WL range allowed (used to compute the maximum
*               search range
* max_search -- Largest number of ambiquity values we can search over
* max_amb_save -- Maxumim number of ambiquity choices to save while
*               searching for ambiquities
* max_anal_type -- Maximum number of analysis types that can be run
*     (For example, analysis maybe with L1, L2 and LC).
* max_edits     -- Maximum number of edits that can be specified
* max_atm_ref   -- Maxumum number of atmospheric delay values that
*     can be read
* max_dabuf     -- Maximum size needed for the buffer used to write
*               direct access file with covariance matrix, solution 
*               vector and ambiquity numbers during the smoothing run.
* max_stopgo    -- Maximum number of stop-go entries that can be saved
*               based on ID types 2 and 3 in rinex files

      integer*4 max_site, max_kine, max_prn, max_chan, max_epochs,
     .          max_ambs, max_parm, max_dtype, max_obs,
     .          max_wlrng, max_search, max_amb_save, max_anal_type,
     .          max_edits, max_atm_ref,
     .          max_dabuf, max_ep_wrds, max_stopgo 

      parameter ( max_site =  10 )
      parameter ( max_kine =   9 )
      parameter ( max_prn  = 100 )
* MOD TAH Increased max_chan for GNSS solutions.
      parameter ( max_chan = 40 )
      parameter ( max_epochs = 86400 )   ! 12-hr flight with 2Hz data
      parameter ( max_ep_wrds = max_epochs/32 + 1 )
      parameter ( max_ambs = max_site*max_prn*10 )  ! allow 10 ambiquities per
*                                          site per satellite

      parameter ( max_dtype = 3 )        ! Assume 3 observables at one time
                                         ! e.g. L1+L2+P1
      parameter ( max_obs = max_dtype*max_chan*max_site )
*     Allocate parameters for atmospheric delay, XYZ per kinematic sites
*     and bias parameters for each observation type.
* MOD TAH 071029: Increased max_parm to allow for all data. Upto 5 breaks on
*     each satellite at each station.
c      parameter ( max_parm = max_site + 3*max_kine + max_obs)
      parameter ( max_parm = max_site + 3*max_kine + 
     .                       5*(2*max_sat*max_site))

      parameter ( max_dabuf = (max_parm+1)*max_parm + 3*max_parm )

      parameter ( max_wlrng  = 12 )      ! 12 cycles in MW widelane
      parameter ( max_search = 6*max_wlrng*max_wlrng)  ! 6 is approximately
                                         ! 77/17 used to get NL1 from MW-WL
               
      parameter ( max_amb_save = 1000 )
      parameter ( max_anal_type = 10  )
      parameter ( max_edits    = 1000 )
      parameter ( max_atm_ref  =  970 )
      parameter ( max_stopgo   = 1024 )

      integer*4 max_wl   ! Max number of saved widelanes for CSLIP repair
      parameter ( max_wl = 100 )

* track_version -- Version of the track program:
* 1.00 -- Initial version 990920
* 1.01 -- Version 991029: Added float biases for biases that are not resolved
* 1.02 -- Version 991228: Added backward running filter and the ability of not
*         doing any SEARCH type ambiquity resolution.
* 1.03 -- Version 000617: Added smoothing filter feature.
* 1.04 -- Version 000819: Better reading of TEQC nav files.  Extra line sometimes
*         has two arguments which caused problems.
* 1.05 -- Version 010616: Added Ion variance factor to floating ambiquity resolution
*         Corrected some problems with WL book keeping and use in bf_float algorithms
* 1.06 -- Version 020205: Better cycle slip detection with SD observable, fixed bias  
*         fixing when reference bias not fixed already, added diagnostic for why
*         bias not fixed.
* 1.07 -- Bug fixes and work-arounds for HPUX 11/f90.
* 1.08 -- Added new features:
*         1. Editing of phase data if epoch RMS exceeded tolerance.
*         2. Automatically remove millisecond slips in phase from uZ data.
* 1.09 -- Added new features
*         1. Better checking of the MW-widelane for detecting cycle slips
*         2. Added EXCLUDE_SVS command to remove a satellite completely
*         3. Fixed residuals output code
*         4. Added WLS_ROOT and RWL_ROOT to output the EX and MW widelanes from
*            the final ambiquiuities and raw (no cycle slips corrected)
*         5. Added passing day number string and week number string through the
*            runstring
* 1.10 -- Added new features
*         1. Detection of ID records 2 and 3 in rinex files that denote
*            ID 2 Start KINEMATIC mode; ID 3 start STATIC mode.  Option added
*            that will reduce process noise to zero during STATIS mode.
* 1.11 -- Some small bug fixes and output of atmopsheric delays at fixed sites
* 1.12 -- Introduced mwwl_jump command to allow control of jump in MW-WideLane 
*         that will be flagged as cycle slip (used to be 5 cycles hardwired)
*         Changed calculation of WL sigmas so it can never be greater then RMS.
* 1.13 -- Updated to handle sp3c format; updates to handle IIRM satellites (blk 5)
*         changes to theory to had bad range data. 
* 1.14 -- Added phase center models for both stations and satellites.  New commnd
*         antmod_file specifies ANTEX format file to used, Ante_off command
*         modified to all antenna specification
* 1.15 -- Cleaned up double difference bias selection to give more weight to the
*         reference PRN combinations, Added command to add ambiquity.
* 1.16 -- Fixed problems with using Gamit Zenith delasy estimates and with epoch
*         alignemt for 10-Hz data (reduced time match tolerance from 0.1 to 0.01 s)
* 1.20 -- Major updates to ambiguity resolution code.  Introduced elevation 
*         angle depependent data noise (in data_noise command); Make clock 
*         error tolerance depend on position uncertainty.
* 1.21 -- Fixed bug that caused jumps in smooth solutions.  Added correlation between
*         height and atmospheric delay estimates to NEU output.  Updated ambiguity
*         selection. Added REF_NEU and USE_GPTGMF.  Added addition of WL averaging.
*         Add command TIME_UNIT for time unit on process noise
* 1.22 -- Added set_msec_bf command to allow bias flag when the MW-WL shows a jump
*         due to millisec range jumps.
* 1.23 -- Fixed problems with ambiguities being found to be not-dependent, resorted
*         ambiguities when O status during fixing.  Tuned the post-fit editing.
* 1.24 -- Added features from trackRT. Rinex files simply concatenated, velocities in
*         site_pos command and new commnd TIMEDEP_PROCNS
* 1.25 -- Introduced a dump file out and fixed bugs with multistation residual 
*         rms calculations.
* 1.26 -- Added reading of ionex files (IGS format) which are interpolated to yeild
*         ionospheric delay estimates, added new LOS ion file that can directly read
*         fixed a bug in EXWL chi-squared contribution in ambiguity fixing (should 
*         improve long-baseline ambiguity resolution, improved bias fixing iterations,
*         Added atm_scale option for kinematic site atmospheric delay estimation
*         (Scale factor for atm-delay/height difference partial).  New MIN_TOLS 
*         command which sets minimum values (had been hard wired in prior versions)
* 1.27 -- Added Relative humidity argument to USE_GPTGMF command for compatability 
*         with GAMIT use of GPT. Default left at 0.00 (used 0.50 for GAMIT compatability).
* 1.28 -- Fix error in detecting bad range data based on clock estimates.  Program
*         now more robust to detecting bad data
* 1.29 -- Modified the calculation of widelanes with o-minus-c values to account for model
*         differences between L1 and L2; WLS_root outputs with models applied; RWL_root just
*         uses the raw data with no model effects accounted for.  Fixed a bug in reading
*         the last epoch in the IONEX files which contain RMS values as well as TEQ maps.
* 1.30 -- Updated read_all_rinex etc to read 2.11 rinex with multiple constelations (GPS
*         only processed.  Add new RM_CSLIP command to remove cycle slips before ambiquity
*         resolution.  
* 1.31 -- (180321) Updated to GNSS processing. Introduced new command tr_gnss <RGEC> to
*         which GNSS to process,
* 1.40 -- (200227) Updates for ocean tides and VMF1 from UTas group.  Added use_blq
*         command to apply ocean tide 
* 1.41 -- (200519) Updates to solid-Earth model and behavior of use_blq command.

      character*4 track_version
      parameter ( track_version = '1.41' )

* RINEX DATA saved array information (used so that rinex files only
* ----------
* need to be read once.)
* L1o_all_cse(max_chan,max_site,max_epochs) -- L1 phase measurements (cycles)
* L2o_all_cse(max_chan,max_site,max_epochs) -- L2 phase measurements (cycles)
* P1o_all_cse(max_chan,max_site,max_epochs) -- P1 range measurements (m)
* P2o_all_cse(max_chan,max_site,max_epochs) -- P2 range measurements (m)
* kine_xyz(3,max_kine,max_epochs)           -- XYZ position of the kinematic
*     platforms at each epoch

      real*8  L1o_all_cse(max_chan,max_site,max_epochs),
     .        L2o_all_cse(max_chan,max_site,max_epochs), 
     .        P1o_all_cse(max_chan,max_site,max_epochs), 
     .        P2o_all_cse(max_chan,max_site,max_epochs), 
     .        kine_xyz(3,max_kine,max_epochs)

* ambiq_all(2,max_ambs)  -- Ambiquities at L1 and L2 (cycles).  The amb_point_cse
*     points to the particular entry needed for each measurement.
* ambiq_float(2,max_ambs)  -- Floating estimates of the ambiquity adjustments
* ambiq_flsig(2,max_ambs)  -- Sigmas for the floating ambiquities.
* wls_all(2,max_ambs)    -- Estimates of residual MW and Extra-wide lanes after
*                           removing the ambiguities in ambiq_all.
* wls_sig(2,max_ambs)    -- Standard deviation of the residual widelanes.
*                           (Estimates based on RMS and number of data).
* data_start(max_site)   -- Start of data for each file (MJD)
* data_end(max_site)     -- End of data at each site (MJD)
* ref_start -- Reference start time.  It is from this time that the epochs
*              are counted.  It is the time when at least two receivers have
*              started
* ref_sec   -- Reference start expressed as seconds from first SP3 file
*              time.
* usr_start -- Start time speficied by the user (If not specified then
*     eariliest possible start time is used. (MJD)
* usr_interval -- Data interval selected by the user.  If not specified 
*     then the interval for the rinex file is used.
* kine_out(3,max_kine) -- Kinematic position to output
* wl_tau       -- Correlation time for WL-sigma calculation (sec)
* lambda(1), lambda(2) -- L1/L2 wavelength
 

      real*8 ambiq_all(2,max_ambs), wls_all(2,max_ambs),
     .       ambiq_float(2,max_ambs), ambiq_flsig(2,max_ambs),
     .       wls_sig(2,max_ambs),
     .       data_start(max_site), data_end(max_site), 
     .       ref_start, ref_sec, usr_start, usr_interval,
     .       kine_out(3,max_kine), wl_tau, lambda(2)

* elev_cse(max_chan,max_site,max_epochs)    -- Elevation angle for each
*     measurement (deg)
* az_cse(max_chan,max_site,max_epochs)    -- Azimuth of each measurement (deg)
* sec_offset(max_site,max_epochs)  -- Time tag offset from the nominal
*     time tag of the measurement

      real*4 elev_cse(max_chan,max_site,max_epochs),
     .       az_cse(max_chan,max_site,max_epochs),
     .       sec_offset(max_site,max_epochs) 

* data_flag_cse(max_chan,max_site,max_epochs) -- Bit mapped data flag for
*     each measurement.  The mapping of the bits is:
*     Bit  Meaning
*       1  Unless set, there is no data at this epoch
*       2  if set, data is bad (deleted as an outlier)
*       3  L1 cycle slip at this epoch
*       4  L2 cycle slip at this epoch
*       5  Data below elevation cutoff (this needs to be dynamically computed
*          as the kinematic position estimates change)
* MOD TAH 040316: Added additional bits to the data flag to record KINEMATIC
*          and STATIC MODES
*       6  When set receiver is in KINEMATIC mode (as set by ID 2 in rinex
*          records).
* MOD TAH 160414: Add bit to show ambiguity added by user so that cslip_repair
*          will not remove it
*       7  Set when BF added by user. (ss_abf).
* amb_point_cse(max_chan,max_site,max_epochs)  -- Pointer to ambiquity value
*      to be used for each measurement.  If the entry is zero then the ambiquity
*      is not known.
* ctop_cse(max_chan,max_site,max_epochs) -- Mapping from channel number to
*      PRN for each measurement.
* num_chan_se(max_site,max_epochs)       -- Number of channels at each site
*      for each epoch


      integer*4 data_flag_cse(max_chan,max_site,max_epochs),
     .          amb_point_cse(max_chan,max_site,max_epochs),
     .          ctop_cse(max_chan,max_site,max_epochs),
     .          num_chan_se(max_site,max_epochs)

* prn_used(max_prn)  -- List of PRNS actually used in the data.
* num_prn            -- Total number of PRN's in the data set.
* num_site           -- Number of sites being processed.
* num_kine           -- Number of kinematic sites
* num_epochs         -- Total number of epochs of data in the data set. 
*                       (Epochs are counted from the first overlapping data
*                        points)
* num_ambs           -- Number of ambiquities used
* bf_ents(5,max_ambs) -- Information about the bias flags
*     Index  Meaning
*         1  site number
*         2  PRN
*         3  Starting Epoch
*         4  Ending Epoch
*         5  Bit 1 set once the bias is resolved with Wide-lanes
*            Bit 2 set once the bias is finally fixed with search
*            Bit 3 reserved for user fixed in amb_tab
*            Bit 4 set if not resolved due to "other" bias not fixed
* wls_ref(3,max_ambs) -- Other ambiquities used in forming the wide
*            lanes. (These are the numbers from bf_ents)


      integer*4 prn_used(max_prn), num_prn, num_site, num_kine,
     .          num_epochs, num_ambs, bf_ents(5,max_ambs),
     .          wls_ref(3,max_ambs)

* max_gap  -- Largest gap allowed
* min_good -- Minium number of good data needed for a bias flag to
*             be keep.
* data_mask  -- Bit pattern set flag bad data
* usr_nepochs -- Number of epochs specified by the user.  All possible
*     epochs are used if not given
* MOD AZ 190305: pin_holder for testing
* pin_holder -- for testing only - common
 
      integer*4 max_gap, min_good, data_mask, usr_nepochs, pin_holder

* EDITING TABLES
* num_edits  -- Number of edits given by user
* num_abf    -- NUmebr of user set bias flags
* num_rbf    -- Number of user biases to be removed.
* ss_exclude(max_prn) -- List of satellites to exclude (used when
*               when satellite not in SP3 file or is missing for
*               an interval of time
* num_exclude -- Number of excluded satellites
* ss_edit(2,max_edits) -- Site and satellite to be edited
* ss_abf(2,max_edits) -- Site and satellite for user bias flag
* ss_rbf(2,max_edits) -- Site and satellite bias to be removed

      integer*4 num_edits, ss_edit(2,max_edits), ss_abf(2,max_edits),
     .          ss_exclude(max_prn), num_exclude, num_abf, num_rbf,
     .          ss_rbf(2,max_edits)

* tt_edit(2,max_edits) -- Start and stop times to be edited (MJD)
* tt_abf(max_edits) -- time of user bf (next obs after this time) (MJD)
* tt_rbf(max_edits) -- time of bf to be removed (MJD)

      real*8 tt_edit(2,max_edits), tt_abf(max_edits), tt_rbf(max_edits)

* RM_CSLIP command variables/NLx_est implemenation
      integer*4 num_mwwl, num_exwl   ! Number of values before and after
                ! slip to average the MW and EX widelanes
      real*4 cslip_tex, cslip_tmw, 
     .       cslip_sex, cslip_smw    ! Max devation and sigma limit for 
                ! EX and MW fits in cycle-slip repair

      real*4 NLx_RR   ! Relative rank for the NLx cycle slip repair
     .,      NLx_addChi   ! Additive minimum chi**2 to add
     .,      NLx_minChi   ! Maximum value of least chi**2 vale 

* Multi-GNSS: Define variables for the reference L1 and L2 frequencies
      real*8 fR1, fR2   ! Reference frequencies for L1 and L2 
                        ! Systems such as Glonass will be mapped to these
                        ! frequencies (SP3 file common contains actual 
                        ! frequencies.
* Factors that will be computed from the Reference frequencies and
* used for the MW-WL, LC, LG, PC and EX combinations.

      real*8 dfsf, sfdf, lcf1, lcf2, lgf1, lgf2, pcf1, pcf2
     .     , exf1, exf2, l1tecu, l2tecu

*----------------------------------------------------------------------------

      common / track_com_b8 /  L1o_all_cse,  L2o_all_cse, 
     .         P1o_all_cse,   P2o_all_cse,                        
     .         kine_xyz, ambiq_all, ambiq_float, ambiq_flsig,
     .         wls_all, wls_sig, 
     .         data_start, data_end, ref_start, ref_sec,
     .         usr_start, usr_interval, tt_edit, tt_abf,  tt_rbf,
     .         kine_out, wl_tau, lambda, 
     .         fR1, fR2,  dfsf, sfdf, lcf1, lcf2, lgf1, lgf2, 
     .         pcf1, pcf2 , exf1, exf2, l1tecu, l2tecu

      common / track_com_b4 / elev_cse, az_cse, sec_offset,
     .         data_flag_cse, amb_point_cse, ctop_cse, num_chan_se,
     .         prn_used, num_prn, num_site, num_kine, num_epochs,
     .         num_ambs, data_mask, bf_ents, wls_ref,
     .         max_gap, min_good, usr_nepochs, num_edits,
     .         ss_edit, ss_abf, ss_exclude, num_exclude, 
     ,         num_abf, ss_rbf, num_rbf, num_mwwl, num_exwl,
     .         cslip_tex, cslip_tmw, cslip_sex, cslip_smw
     .,        NLx_RR, NLx_addChi, NLx_minChi, pin_holder

*----------------------------------------------------------------------------

* ESTIMATION DECLARATIONS:
* -----------------------
*

* APRIORI VARIANCES AND PROCESS NOISE

* cov_parm(max_parm, max_parm) -- Covariance matrix for the parameter
*      estimates
* sol_vec(max_parm) -- Solution vector for parameter estimates.
* cov_parm_sav(max_parm, max_parm) -- Saved covariance matrix.  Needed 
*      for the smoothing solution
* sol_vec_sav(max_parm) -- Saved solution vector needed for smoothing 
*      analysis.
* cov_obs(max_obs, max_obs)    -- Covariance matrix for the observations
* sol_obs(max_obs) -- Final o-minus-c values for the observations to be
*      used in the Kalman Filter.
* apart(max_parm,max_obs) -- Partial derivative array. NOTE: Store this
*      in transposed form to ensure multication down columns  
* kgain(max_parm,max_obs) -- Kalman gain matrix.  Again note in transposed
*      form
* obs_to_res(max_obs,max_obs) -- Matrix which multiplies a
*      set of double difference changes to the ambiquitues and generates a
*      new set of residuals.  Adjustment to the parameters are obtained
*      from the kalman gain matrix
* amb_to_res(max_obs,max_obs) -- Matrix which multiplies a
*      set of L1 and L2 changes to the ambiquitues and generates a
*      new set of residuals.
      real*8 cov_parm(max_parm, max_parm), sol_vec(max_parm),
     .       cov_parm_sav(max_parm, max_parm), sol_vec_sav(max_parm),
     .       cov_obs(max_obs, max_obs), sol_obs(max_obs),
     .       apart(max_parm,max_obs), kgain(max_parm,max_obs),
     .       obs_to_res(max_obs,max_obs), amb_to_res(max_obs,max_obs)

* ion_parm(max_site,max_site) -- Covariance matrix for ion estimate
* ion_vec(max_site) -- Estimate of the ionospheric delay at each site
* ion_cov(max_obs,max_obs)    -- OW ion estimate covariance matrix
* ion_obs(max_obs,max_obs)    -- Double difference ion delay
* ow_ion(max_obs)   -- Ionospheric delay data
* ow_iov(max_obs)   -- Variances for the ionospheric data
* ow_ionp(max_site,max_obs) -- Ion delay partials
* kg_ion(max_site,max_obs)  -- Ion Kalman gain matrix
* apion(max_site,max_obs) -- ion double difference partials
* ion_to_res(max_obs, max_obs) -- Mapping from changes in ion data
*      to residuals
* amb_to_ion(max_obs, max_obs) -- Mapping from changes in ambiquities
*      to changes in residuals

* ion_var  -- Constant part of ionospheric delay variance (cycles**2)
* ion_corr -- Correlated portion (multiplied by exp(-(b/tai)**2)
* ion_tau  -- Correlation length of ionosphere (m)
* ion_wght -- Weight to be given to ionosphere in computing RMS
* max_ion_jmp -- Largest jump allowed in the ionospheric observable
* mwwl_jmp -- Size of MW_WL jump that will be flagged

      real*8 ion_parm(max_site,max_site), ion_vec(max_site),
     .       ion_cov(max_obs,max_obs), ion_obs(max_obs),
     .       ow_ion(max_obs), ow_iov(max_obs),
     .       ow_ionp(max_site,max_obs), kg_ion(max_site,max_obs),
     .       apion(max_site,max_obs),
     .       ion_to_res(max_obs, max_obs),  
     .       amb_to_ion(max_obs, max_obs),
     .       ion_var, ion_corr, ion_tau, ion_wght, max_ion_jmp,
     .       mwwl_jmp

* curr_site_xyz(3,max_site) -- Position estimate at the current
*      epoch for all of the sites.  The kinematic sites have 
*      their values stored back in to kine_xyz array.
* ow_part(max_parm,max_obs) -- One way partial derivatives (m/m for
*      position and atm partials)
* ow_vec(max_obs) -- One way vector of o-minus-c values
* ow_var(max_obs) -- Variance of the one-way measurements
* rms_dd  -- RMS of the final double differences for each epoch
* rms_dd_avg -- Average RMS of the double differences (used for scaling
*     sigma)
* rms_edtol -- Tolerance on the RMS before we start editing data 
* min_lvar  -- Minumum phase variance for editing
* wrms_dd_site(max_site)  -- WRMS scatter of residials by site


      real*8 curr_site_xyz(3,max_site), ow_part(max_parm,max_obs),
     .       ow_vec(max_obs), ow_var(max_obs), rms_dd, rms_dd_avg,
     .       rms_edtol, min_lvar, wrms_dd_site(max_site)

* apr_site(3,max_site) -- Apriori variances for the site positions
*      (currently only the kinematic sites are assumed to have
*       non-zero values)  (m**2)
* pst_site(3,max_site) -- Aposterori variances for the site positions
*      (currently only the kinematic sites are assumed to have
*       non-zero values)  (m**2)

* mar_site(3,max_site) -- Markov process noise on the kinematic
*       site positions (m**2/epoch)
* apr_atm(max_site)    -- Apriori variance for atmpospheric 
*       delays.  (m**2)
* mar_atm(max_site)    -- Process noise for the atmospheric 
*       delays   (m**2/epoch)
* mar_atm_hgt(max_site) -- Process noise for atmosphere that depends on
*       height rate of change (m**2/(dhdt)**2/epoch (ie. value multiplied
*       by dh/dt in m/s to to process noise

* data_var(5,max_sat)  -- Data sigmas for L1, L2, P1, P2 data types
*     Plus elevation angle weighting (default 0.5)
*     Sigmas for LC and PC are computed from these
* out_sig_limit        -- Maximum error bar for outout of results
*     Default 1 m
* stopgo_dvar   -- Variance reduction for stop go mode

* tu_to_ep      -- Conversion factor from time unit to epoch
* min_lcsig     -- Minimum sigma for LC/L1/L2 estimates
* dynamic_tol   -- RMS scatter of PC solution to define moving site

      real*4 apr_site(3,max_site), pst_site(3,max_site), 
     .       mar_site(3,max_site),
     .       apr_atm(max_site), mar_atm(max_site),mar_atm_hgt(max_site),
     .       data_var(5,max_sat), out_sig_limit, stopgo_dvar, 
     .       tu_to_ep, min_lcsig, dynamic_tol

      real*4 mar_tdep(3,5*max_site)   ! Time dependence process noise
      real*8 mjd_tdep(2,5*max_site)   ! Epoch ranges for time dependence
      integer*4 nst_tdep(5*max_site)  ! Site numbers (-1 means all)
      integer*4 num_tdep              ! Number of time dependent entries.

* num_anal_type -- Number of analysis types to be run
* num_parm    -- Number of parameters to be estimated
* num_dtype   -- Number of data types used in current estimation 
* pos_parn(3,max_kine) -- Parameter numbers for the kinematic site
*                positions
* atm_parn(max_site)   -- Parameter numbers for the number of atmospheric
*                delay parameters to be estimated
* amb_parn(2,max_ambs) -- Ambiguity parameter numbers for each of bias parameters
*                         that need estimation.
* amb_save(2,max_ambs) -- Saved number Ambiguity parameter numbers for each of bias parameters
*                         that need estimation.
* parn_to_amb(max_parm)  -- Mapping from parameter number to ambiquity
*                number (NOTE: Same ambiquity number is used for L1 and L2).
* non_amb_parm    -- Number of non-ambiquity parameter numbers
* kine_to_site(max_kine) -- Mapping of kinematic site number to overall
*                site number.
* ow_des(2,max_obs)  -- Observation vector descriptor. The entries
*       1 -- Site  2 -- PRN for one-ways
* sd_des(2,max_obs)  -- Single difference pointer.  Points to the two
*      oneway measurements in the SD
* dd_des(4,max_obs)  -- Double difference pointer.  Points to the four
*      oneway measurements in the DD 
* dd_dsi(4,max_obs)  -- Double dofference pointers for the ionospheric
*      delay obervable.  Differs from dd_des when more than one obervable
*      is used.
* num_ow  -- Number of oneways
* num_sd  -- Number of SDiffs
* num_dd  -- Number of DDiffs.
* num_ion -- Number of ionospheric delay estimates
* num_ddion -- Number of ionospheric delay double differences
* num_ow_by_site(max_site)  -- Number of oneways at each site.
* num_dd_avg -- Number of double differences in computation of double
*     difference RMS
* num_dd_site(max_site) -- NUmber of double diffs by site

* tot_parm   -- Maxiumum number of parameters seen at any time during the 
*     filter run
* tac_parm   -- Actual total number of parameters (accounts for changes
*               in the number of observables between the float and data_type
*               solutions.
* darecl     -- Record length (I*4 words) needed to save the covariance matrix
*               and solution vector and ambiquity parameter numbers needed for
*               smooting filter run
* kine_OK(max_ep_wrds,max_kine) -- Bit set to 1 for each kinematic position
*               with a sigma < 1 m.
* ow_dt(max_obs) -- Bit mapped flag for data type in one-way array.  Bit 1 L1
*               bit 2 L2, bit 3 P1, bit 4 P3.  LC and PC set bits 1/2 and 3/4
* dd_dt(max_obs) -- Same quantity for double differences. 
* static_obs(max_site) -- Number of static measurements at a site currently
* stopgo_point(max_site)   -- Point number for stopgo mode

      integer*4 num_parm, num_dtype, pos_parn(3,max_kine),
     .          atm_parn(max_site), kine_to_site(max_kine),
     .          ow_des(2,max_obs), sd_des(2,max_obs),
     .          dd_des(4,max_obs), dd_dsi(4,max_obs), 
     .          num_ow, num_sd, num_dd,  num_dd_site(max_site), 
     .          num_dd_avg,num_ow_by_site(max_site), num_ion,
     .          num_ddion, num_anal_type, amb_parn(2,max_ambs),
     .          amb_save(2,max_ambs),
     .          parn_to_amb(max_parm), non_amb_parm, tot_parm, 
     .          tac_parm, darecl, kine_OK(max_ep_wrds,max_kine),
     .          ow_dt(max_obs), dd_dt(max_obs), static_obs(max_site),
     .          stopgo_point(max_site)

* kine_known  -- Logical that is set true once we have initial
*      kinematic site position estimates from P1 data.
* stopgo_mode -- When true, track will use the KINEMATIC/STATIC flags
*      in the rinex file to change the kinematic process noise.
* ante_off    -- Antenna information read
* atm_mtt     -- Set true to use the MTT mapping function and seasonal
*                model.  Set false to use GPT/GMF or VMF
* atm_gmf     -- Set true to use the GMF mapping function
*                Default value:false
* MOD AZ 190305: new variables for VMF and Tide implemenation
* atm_vmf     -- Set true to use the VMF mapping function - common
*                Default value:false
* MOD TAH 200225: Make use of ocean-tide BLQ files optional
* use_blq     -- Set true if ocean tide to be applied.
* leapsec_path -- leap second table
* nbody_path -- Sun&Moon position table
* site_otl --   Site name for BLQ files
* noetide     -- Set true to stop the application of the Earth tide 
* etide_2010  -- Use IERS 2010 E-tide model versus older analytic model
*                (Added TAH 200519).  Option in use 
* set_msec_bf -- Logical set true with the SET_MSEC_BF command if
*                bias flags are to be added to msec jumps in the MW-WL
* atm_scale(max_site)  -- Set true to atm dH scaling rather than zenith
*                delay estimate (set in atm_stats command).   Only
*                possible with kinematic sites).

      logical kine_known, stopgo_mode, ante_off, atm_mtt, noetide,
     .        set_msec_bf,  atm_scale(max_site), atm_gmf, atm_vmf,
     .        use_blq, etide_2010

* anal_types(max_anal_type)  -- List of analysis types to be done
* search_type  -- Type of search to be preformed
* float_type -- Type of analysis for float ambiquity analysis
* back_type  -- Type of back solution to run (BACK commmand)
* time_unit  -- String containing time unit for process noise models

      character*256 leapsec_path, nbody_path
      character*8  site_otl(max_site)
      character*16 anal_types(max_anal_type), search_type, float_type,
     .             back_type
      character*8  time_unit

      character*8  tr_gnss   ! String with GNSS systems to process (GREC
                   ! supported).  Set with TR_GNSS command

*----------------------------------------------------------------------------
      common / track_est_b8 / cov_parm, sol_vec,
     .        cov_obs, sol_obs,  apart, kgain, obs_to_res, 
     .        amb_to_res, curr_site_xyz, ow_part,  ow_vec, ow_var,
     .        rms_dd, rms_edtol, rms_dd_avg, min_lvar, wrms_dd_site,
     .        ion_parm, ion_vec, ion_cov, ion_obs, ow_ion, ow_iov,
     .        ow_ionp, kg_ion, apion, ion_to_res,  amb_to_ion,
     .        ion_var, ion_corr, ion_tau, ion_wght, max_ion_jmp,
     .        mwwl_jmp, cov_parm_sav, sol_vec_sav, 
     .        mjd_tdep

      common / track_est_b4 / apr_site, pst_site, mar_site,
     .       apr_atm, mar_atm, mar_atm_hgt, data_var, 
     .       out_sig_limit, stopgo_dvar, tu_to_ep, 
     .       min_lcsig, dynamic_tol,
     .       num_parm, num_dtype, pos_parn, parn_to_amb, non_amb_parm,
     .       atm_parn, amb_parn, amb_save, kine_to_site, 
     .       ow_des, sd_des,
     .       dd_des, dd_dsi, num_ow, num_sd, num_dd, num_dd_site,  
     .       num_dd_avg, num_ow_by_site,
     .       num_ion, num_ddion, num_anal_type,
     .       kine_known, stopgo_mode, ante_off, atm_mtt, noetide,  
     .       set_msec_bf, atm_scale, tot_parm, tac_parm, darecl, 
     .       kine_OK, ow_dt, dd_dt, static_obs, stopgo_point,
     .       mar_tdep, nst_tdep, num_tdep, atm_gmf, atm_vmf, 
     .       use_blq, etide_2010

      common / track_est_ch / anal_types, search_type, float_type,
     .       back_type, time_unit, tr_gnss, leapsec_path,nbody_path

*----------------------------------------------------------------------------

* BIAS FIXING VARIABLES BASED ON SEARCH
* ion_lim(2,max_chan)  -- lower and upper limits on the extra-wide
*     lane (cycles) based on a model ionspheric magnitude
* save_rms(max_amb_save) -- Saved RMS values of redisiduals for
*     ambiquity search
* save_rmsi(max_amb_save) -- Saved ION RMS Values
* rank_min_rms  -- Minimum RMS scatter to use when ranking
*     ambiquitity sets by inverse RMS**2
* comb_ivar(max_amb_save) -- Sum of inverse of RMS**2 sorted by decreasing
*    value
* comb_ion(max_amb_save)  -- Sum of inverse of RMS**2 for the ionospheric
*    delay residuals
* max_dd_elev(max_obs)  -- Maximum difference allowed in residual
*     while searching.  Elevation angle dependent.  If value
*     exceeded for one measurement then search is stopped.
* stratm_dcorr -- Size of structure function coefficient for
*     atmospheric delay (Units: m/sqrt(D km) when D is baseline
*     length
* stratm_base  -- Zero baseline length limit on the size of 
*     double difference residual allowed while searching,
* ion_ppm   -- Parts-per-million size of ionospheric 
* ion_ht    -- Height of ionospheric layer (m)
* saved_rank(max_ambs) -- Saved value of best rank for each ambiquity.
*     If ambiquity can not be reliably resolved, ambiquities set at
*     value with highest rank
* relrank_limit -- Relative rank ratio needed to fix biases.
* relrank_iter  -- Iteration count when relrank = 1.584 user tolerance
* float_limit(2) -- Limits on sigma for fixing float biases and max sigma
*     for which resolution will be attempted. 
*      Relrank_limit is used as fit quality.
* wl_fact  -- Scaling factor mixing the (res/sig)^2 + wl_fact*(wl/rms)**2 
*     in determining the fit quality
* lg_fact  -- Scaling factor mixing the ion delay variance into the fit
*     criterion (for short baselines should be 1)
* max_fit  -- Maximimum value allowed for fit parameter to have bias
*     fixed.

      real*4 ion_lim(2,max_chan), save_rms(max_amb_save),
     .       save_rmsi(max_amb_save),
     .       rank_min_rms, comb_ivar(max_amb_save), 
     .       comb_ion(max_amb_save),
     .       max_dd_elev(max_obs), stratm_dcorr, stratm_base,
     .       ion_ppm, ion_ht, saved_rank(max_ambs),
     .       relrank_limit, float_limit(2), wl_fact, lg_fact, max_fit,
     .       relrank_iter


* num_bf -- Number of bias flags at this epoch
* bf_flags(max_chan*max_site) -- Entries for the bias flags in the ambiquity
*           table
* num_amb_samp -- Number of samples (evenly spaced over the duration
*     of a bias flag to use in resolving ambiquities)
* max_tot_search -- Maximum number of values to try and search before
*     deferring to a later search (hopefully when more more ambiquities
*     have been fixed)
* num_tot_resolved -- Number of newly resolved bias flags with each
*     SEARCH iteration.


* num_sent -- Number of L1, L2 pairs of searchs in the table of
*     biases to be searched
* num_search(max_chan) -- Number of entries in each search line
* bf_search(2,max_chan,max_search) -- Biases to be search at L1 and
*     L2
* amb_search(max_chan) -- Ambiquity entry pointed to by each of the
*     entries to be searched.
* num_amb_save  -- Number of save ambiquities while searching
* save_amb(2,max_chan,max_amb_save) -- The saved set of ambiqiuities
* save_na(max_chan, max_amb_save)   -- Ambiquity numbers for the
*     saved ambiquities (incase numbers change in region searched)
* save_ep(max_amb_save)  -- Epoch numbers for each of save enties.
* save_nu(max_amb_save)  -- Number of ambiquities saved for entry
* rank_na(max_chan) -- List of ambiquity numbers to be ranked
* rank_nu -- Number of entries in rank_na
* num_ranked -- Total number of values ranked
* comb_amb(2,max_chan,max_amb_save) -- dL1 and dL2 values that are
*    associate with comb_ivar.
* comb_tot(max_amb_save)  -- Total number of values used in ranking
* float_iter  -- First iteration to use float feature.
* float_sample -- Sampling iterval to use in float solutions.

      integer*4 num_bf, bf_flags(max_chan*max_site),
     .          num_sent, num_search(max_chan),
     .          bf_search(2,max_chan,max_search),
     .          amb_search(max_chan), num_amb_save,
     .          save_amb(2,max_chan, max_amb_save),
     .          save_na(max_chan, max_amb_save), 
     .          save_ep(max_amb_save), save_nu(max_amb_save),
     .          rank_na(max_chan), rank_nu, num_ranked,
     .          comb_amb(2,max_chan,max_amb_save),
     .          comb_tot(max_amb_save), num_amb_samp,
     .          max_tot_search, num_tot_resolved, 
     .          float_iter, float_sample

*----------------------------------------------------------------------------
      common / track_sea_b4 / ion_lim, save_rms, save_rmsi,
     .          rank_min_rms, comb_ivar, comb_ion,
     .          max_dd_elev, stratm_dcorr, stratm_base,
     .          ion_ppm, ion_ht, saved_rank, relrank_limit, 
     .          relrank_iter, float_limit, wl_fact, lg_fact, max_fit, 
     .          num_bf, bf_flags, num_sent,
     .          num_search, bf_search, amb_search, num_amb_save,
     .          save_amb, save_na, save_ep, save_nu,
     .          rank_na, rank_nu, num_ranked,  
     .          comb_amb, comb_tot, num_amb_samp,
     .          max_tot_search, num_tot_resolved, float_iter,
     .          float_sample


*----------------------------------------------------------------------------

* COMMANDS and MISCELLANEOUS parameters declarations
*

* site_apr(3,max_site)   -- Apriori coordinates for the sites
*      (for kinematic sites, this position is approximately
*       the location at the start)
* site_geod(3,max_site)  -- Geodetic corrinates of sites computed 
*       from site_apr (co-lat, long (rad), ellip hgt (m)).
* site_int(3,max_site)   -- Initial site position at MJD
* site_vel(3,max_site)   -- Apriori site velocity (m/yr)
* baselens(max_site*(max_site-1)/2)  -- Baseline lengths (m). Used
*                           for ex-wl constraints
* site_epoch(max_site)   -- Epoch for site postion (MJD).
* site_offarp(3,max_site) -- NEU offsets for the antenna Reference
*       point for each site (m)
* site_offsL1(3,max_site) -- Phase center offset in NE and U for
*       the L1 phase center: Removed 0612128
* site_offsL2(3,max_site) -- Phase center offset in NEU for L2 
*       Removed 0612128 (sit_L12 replaces variable)
* elev_cutoff  -- Elevation angle cuttoff to use (deg)
* atm_offset(max_site) -- Correction to be applied to the atmospheric
*     data by site 
* ref_xyz(3)  -- reference XYZ coordinates to refer NEU values to
*     (if not specified using the REF_NEU command then site 1 is used) 

      real*8 site_apr(3,max_site), site_geod(3,max_site), 
     .       site_int(3,max_site), site_vel(3,max_site),
     .       baselens(max_site*(max_site-1)/2),
     .       site_ep(max_site),site_offarp(3,max_site),
     .       elev_cutoff, atm_offset(max_site), ref_xyz(3) 

* site_type(max_site)    -- Type of site for each file.  Options are
*     Type 0  -- Fixed site     F
*     Type 1  -- Kinematic site K
* rcv_type(max_site)  -- type of DCB offsets for receiver
*     P -- L1 range correction 
*     C -- L1 and L2 range correction
*     N -- No correction needed (P-code ranges)
*     Types can be found rcvant.dat 
* dcbs(max_sat) -- DCBs for each satellite (meters)
* dcb_mjd   -- Epoch for dcb values
* debug_start, debug_end -- Range of epochs over which detailed information
*     is printed
* lus  -- Logical number for summary file

      integer*4 site_type(max_site), debug_start, debug_end, lus
      character*2 rcv_type(max_site)
      real*8 dcbs(max_sat), dcb_mjd

* read_site_apr(max_site) -- Logical set true when we have read in the
*     apriori coordinates for a site
* write_res   -- Set true to indicate that residuals should be output
* static(max_kine)  -- Set true if site appears to be static

      logical read_site_apr(max_site), write_res, static(max_kine)

* site_names(max_site) -- 4-character names of the sites
* site_NUC(max_site)   -- Upper case version of site name
* swvers(max_site)     -- software versions for each receiver

      character*4 site_names(max_site), site_NUC(max_site)
      character*16 swvers(max_site)

* obs_file(max_site) -- Names of the observation files
* obs_file_type(max_site) -- Types of files R -- Rinex, X -- xfiles
*                       In Version 1.00, only Rinex is supported
* nav_file  -- Rinex navigation file.  Needed either for orbit or
*              for satellite clock parameters
* nav_file_type -- Type of ephemeris file given.  Only SP3 supported
*              in Version 1.00
* sp3_file  -- Ephemeris file in SP3 format.  
* log_file  -- Name of file for detailed log'ing of information
* bat_file  -- Name of file with commands in it.  Basic control file.
* ambin_file -- Name of file containing the values for the ambiquities.
*      (Format must coincide with the output from a track run with
*      same setup).
* sv_clk_file -- File containing the satellite clock information (not needed)
* dcb_file    -- Name of file with DCB values (needed for mixed receivers)

      character*256 obs_file(max_site), nav_file, sp3_file, 
     .              log_file, bat_file, ambin_file, sv_clk_file,
     .              dcb_file

* resid_root  -- Root part of the name for the residual file
* posit_root  -- Root part of the name for the geodetic position file
* posit_type  -- Type of position output (Either GEOD or NEU)
* wls_root    -- Root for widelanes corrected with ambiquities
* rwl_root    -- Root for widelanes not corrected with ambiquities
* sum_file    -- Name of summary file

      character*256 resid_root, posit_root, posit_type,
     .              wls_root, rwl_root, sum_file

* runday  -- Day number for files: Invoked with <day> in file names
* runweek -- Week number for files: Invoked with <week> in file names
*     (Both arguments are passed in the runstring).
* runstr(max_runstr) -- Update three strings that can be substituted with form
*      <S01> <S02>..  <S10> in commands.

      integer*4 max_runstr   ! Max number of <Snn> strings allows 
*     NOTE: max_snn also needs to be updated in read_bat/subsrdw
      parameter ( max_runstr = 20 ) 

      character*16 runday, runweek
      character*64 runstr(max_runstr)

      character*4 nav_file_type

      character*1 obs_file_type(max_site) 

      common / track_misc_b8 /  site_apr, site_geod, site_int,  
     .       site_vel, site_ep, site_offarp, elev_cutoff,
     .        atm_offset, ref_xyz, dcbs, dcb_mjd

      common / track_misc_b4 / site_type, read_site_apr,
     .       debug_start, debug_end, write_res, lus, static
      common / track_misc_ch / site_names, site_NUC, rcv_type,
     .       swvers, obs_file, nav_file, sp3_file,  log_file, dcb_file,
     .       bat_file, ambin_file, resid_root, posit_root, posit_type,
     .       nav_file_type, obs_file_type, sv_clk_file, wls_root,
     .       rwl_root, runday, runweek, sum_file, runstr

*----------------------------------------------------------------------------

* SP3 File definitions
*

* Nav file information (if need)
*  nav_clk(3,max_sat, max_eph_eps) - Clock poly nominal information from
*                                    the navigation file

      real*8  nav_clk(4,max_eph_eps, max_sat) 

* Satellite position information

*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_clk(max_sat)       - Clock error
*  svs_loc(3,max_sat)    - Phi, Lambda, and height of sattellite

      real*8 svs_clk(max_sat), svs_xyz(3,max_sat),
     .       svs_loc(3,max_sat)

* Ephemeris information
* start_nav, end_nav -- Start and stop times of the entries
*     found in the rinex nav file (used to compute duration
*     of orbit)
      real*8 toe_jd(max_sat), toe(max_sat), aode(max_sat),
     .    af0(max_sat), af1(max_sat), af2(max_sat),
     .    crs(max_sat), crc(max_sat), dn(max_sat), m0(max_sat),
     .    cuc(max_sat), cus(max_sat), ecc(max_sat), art(max_sat),
     .    cic(max_sat), cis(max_sat), om0(max_sat), i0(max_sat),
     .    w(max_sat), omd(max_sat), idt(max_sat),
     .    cflg12(max_sat), pflg12(max_sat), 
     .    svacc, svhealth, tgd, aodc(max_sat), start_nav,
     .    end_nav

* HEADER Information from individual rinex files
* sv_ndat  -- Number of data types in each file
* sv_mdat      
      real*8  sv_ndat(max_site),sv_mdat(max_site),sv_offarp(3,max_site)

* NOTE: xf_maxat is limited to 6 xf_mdat (L1,L2,P1,P2,C1,C2)
      integer*4 sv_dattyp(6,7,max_site) ! Saved xf_dattyp for
            ! each rinex file and by GNSS system.

      real*4 rxver(max_site)   ! Rinex version for each site.


* num_nav(max_sat) -- Number of entries in the NAV file for each PRN
* weekno        -- GPS Week number

      integer*4 num_nav(max_sat), prns(max_sat),
     .          weekno
     
      common / sv_com_b8 /  svs_clk, 
     .       svs_xyz, svs_loc, nav_clk,  sv_ndat, sv_mdat, sv_offarp 

      common / sv_com_b4 / num_nav, weekno, sv_dattyp, rxver

      common / sv_com_ep /  toe_jd, toe, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12,  svacc, svhealth, tgd, aodc,
     .    start_nav, end_nav 

*----------------------------------------------------------------------------
* ATMOSPHERE Tables
* -----------------
* atm_ref_num  -- Number of tabular points in tables
* use_atm_ref -- Logical to indicate that we should use tabular values

      integer*4 atm_ref_num
      logical use_atm_ref

* atm_ref_sep  -- Separation of points (days)
* atm_ref_mjd  -- Start of table (MJD)
* atm_ref_zen(max_atm_ref, max_site) -- Tabular values for the atmospheric
*                 zenith delay (m)
* ref_rel_hum  -- Reference relative humdity for use with GPT
      
      real*8 atm_ref_sep, atm_ref_mjd, 
     .       atm_ref_zen(max_atm_ref, max_site), ref_rel_hum

* atm_file -- Name of file with atmosphere values

      character*256 atm_file

      common / track_atm_b8 / atm_ref_sep, atm_ref_mjd, atm_ref_zen,
     .                        ref_rel_hum
      common / track_atm_b4 / atm_ref_num, use_atm_ref
      common / track_atm_ch / atm_file

*----------------------------------------------------------------------------
* Residual statistics entries
* sum_res(max_dtype,max_kine,max_sat)  - Sum of residual
* sum_var(max_dtype,max_kine,max_sat)  - Sum of residual squared
* sum_num(max_dtype,max_kine,max_sat)  - number of residual
* elv_res(max_dtype,max_kine,18)  - Sum of residual by elev bin
* elv_var(max_dtype,max_kine,18)  - Sum of residual squaresd by elev bin
* elv_num(max_dtype,max_kine,18)  - number of residual by elev bin

      real*8 sum_res(max_dtype,max_kine,max_sat),
     .       sum_var(max_dtype,max_kine,max_sat),
     .       elv_res(max_dtype,max_kine,18),
     .       elv_var(max_dtype,max_kine,18)
      integer*4 sum_num(max_dtype,max_kine,max_sat),
     .       elv_num(max_dtype,max_kine,18)

      common / track_stat_b8 / sum_res, sum_var, elv_res, elv_var
      common / track_stat_b4 / sum_num, elv_num

*****************************************************************************
* Phase center model for track.  Data format is antex but we allow 1x1 deg
* spacing in table. Specific string can be used for "antenna" to identify 
* specific site.
*
* antmod_file  -- Name of antenna phase center file.  (command may be used 
*      multiple times to read different files  Later results over write eariler 
*      ones
* ant_name(max_site) -- Antenna names for each site.  May be non-standard name
*      to model mutlipath at specific site (20-characters)
*
      integer*4 max_antmods   ! Maximum number of antmod_files allowed
      integer*4 max_zen, max_az ! Maximum number allowed to zen and azimuth table

      parameter ( max_antmods = 10 ) 
      parameter ( max_zen = 91  )
      parameter ( max_az  = 361 )

      character*256 antmod_file(max_antmods)
      character*20  ant_name(max_site)

      integer*4 num_antmods   ! Number of antmod files specified

* Site information
* num_sit_dph(2,max_site) -- Number of phase entries in zenith angle and azimuth 
*      directions for each site
      integer*4 num_sit_dph(2,max_site)

* sit_dphs(max_zen,max_az, 2,max_site) -- Phase center adjustments at function of zenith, 
*      azimuth, L1/L2 by site (mm)
* sit_dzn(3,2,max_site) -- zen1, zen2, dzen, az1, az2, daz for each site (deg)
* sit_L12(3,2,max_site) -- "North, East, Up" offsets at L1/L2 for each site (mm)

      real*4 sit_dphs(max_zen,max_az, 2,max_site)
      real*8 sit_dzn(3,2,max_site), sit_L12(3,2,max_site)

      common / track_phs_ch / antmod_file, ant_name
      common / track_phs_b4 / num_antmods, num_sit_dph, sit_dphs
      common / track_phs_r8 / sit_dzn, sit_L12 

*****************************************************************************

* COMMON BLOCK ENTRIES FOR IONEX FILES

      integer*4 max_lat, max_lng, max_tim    ! Dimensions for model
      parameter ( max_lat = 72 )    !  2.5 deg
      parameter ( max_lng = 73 )    !  5.0 deg
*     parameter ( max_tim = 13 )    !  2-hr interval for whole 24-hrs
      parameter ( max_tim = 25 )    !  1-hr interval for whole 24-hrs
                                    ! (changed 11/11/14; 309 2014 files)


      integer*2 ionex_tec(max_lat, max_lng, max_tim)  ! 0.1 TECU Values from ionex
      integer*4 ionex_num          ! Number of times in file.
      logical   use_ionex          ! Set true when ionex files read
      logical   use_ionlos(max_site)  ! Set true when ion LOS values to be read.
      logical   ion_known          ! Set true once the kinematic trajectory have been
                                   ! determined and ion computed.  Also set true if
                                   ! LOS ion values are read from file.

      character*256 ionex_file     ! Name of file with IONEX map
      character*256 ionlos_file(max_site)   ! Name of file with LOS ionospheric TEC values.
                                   ! One file per site.

      real*4 tec_los_cse(max_chan,max_site,max_epochs)  ! LOS TEC values for each
                                   ! measurement.  Save once the kine_postion is
                                   ! know
      real*4 tec_subI_cse(2,max_chan,max_site,max_epochs) ! Long/Lat of sub-ion point

      real*8 ionex_ht              ! Heught of ion-layer (assumed same for all maps)
      real*8 ionex_times(max_tim)  ! MJD of ionex map times
      real*8 ionex_dlat, ionex_dlng  ! Grid spacing in lat/long (deg).

      common / ionex_b2 / ionex_tec
      common / ionex_b4 / ionex_num, use_ionex, use_ionlos, ion_known, 
     .                    tec_los_cse, tec_subI_cse
      common / ionex_b8 / ionex_ht, ionex_times, ionex_dlat, ionex_dlng
      common / ionex_ch / ionex_file, ionlos_file 


*****************************************************************************






