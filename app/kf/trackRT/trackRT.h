*
*     trackRT common block that contains processing results
*     options and data.  File is derived from track_comm.h but
*     differs due to data sources and approachs

* PARAMETER BLOCK
*
* max_site -- Maximum number of sites that can be processed.
* max_sat  -- Maxumim number of satellites allowed
* max_eph_eps -- Maximum number of entries in SP3 file
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
     .          max_anal_type,
     .          max_sat, max_eph_eps, max_edits, max_atm_ref,
     .          max_dabuf, max_ep_wrds, max_stopgo 

      parameter ( max_site =  15 )
      parameter ( max_sat  = 100 )
      parameter ( max_eph_eps = 193 ) ! Not many needed because always updated
      parameter ( max_kine =  15 )
      parameter ( max_prn  = 100 )
      parameter ( max_chan = 100 )
      parameter ( max_epochs = 4 )   ! This is for back solution (maybe)
      parameter ( max_ep_wrds = max_epochs/32 + 1 )
      parameter ( max_ambs = max_site*max_prn*2 )  ! At most only two 
                                   ! ambiguities kept active at one time
      parameter ( max_dtype = 3 )  ! Assume 3 observables at one time
                                   ! e.g. L1+L2+P1
      parameter ( max_obs = max_dtype*max_chan*max_site )
*     Allocate parameters for atmospheric delay, XYZ per kinematic sites
*     and bias parameters for each observation type.

      parameter ( max_parm = max_site + 3*max_kine + 
     .                       (max_sat*max_site)/2)

      parameter ( max_dabuf = (max_parm+1)*max_parm + 3*max_parm )

      parameter ( max_anal_type = 10  )
      parameter ( max_edits    = 1000 )
      parameter ( max_atm_ref  =  970 )
      parameter ( max_stopgo   = 1024 )

* trackRT_version -- Versions of trackRT and trackRTr program
*  1.00  -- Initial version 100301 (2010 March 1)
*  1.10  -- Version to handler BNC 2.5 and greater that is now outputing
*           ascii rather than binary data (trackRTB handles BNC2.0) 110501
*  1.11  -- Version to handle missing satellites in SP3 files 110505
*  1.12  -- Handle bad PRN from BNC where G23 == 203.
*  1.13  -- Fixed time-equate iin SaveObsA so that a differnce of more
*           than 0.1 seconds in needed for a new epoch (trackRT only)
*           Added reference relative humidity for GPT.
*  1.14  -- Fixed problems with possibly missed cycle slips, and corrected
*           a bug in calculating receiver clock errors.
*  1.15  -- Updated to handle BNC 2.6 feed format.
*  1.16  -- Updated to handle Rinex 3 formats being included in recent version
*           of BNC ver 2.6-2.8.  Added ~ for home directory in file names (must
*           be first character).  Improved handing of missing satellites from
*           sp3 files. (121231)

      character*4 trackRT_version
      parameter ( trackRT_version = '1.16' )
*
************************************************************************
*     Unit numbers
*      80  -- Sp3 files
*      81  -- Antmod.dat files
*      82  -- svnav.dat file
*      83  -- DCB file
*      98  -- Summary file
*     100  -- Command file (also used for the upd_file)
*     101-100+max_site  -- rinex files in trackRTr
*     201-200+5*max_site  -- Output position files
*     

************************************************************************

* BASIC Information

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
* wls_sum(2,max_ambs)    -- Sum of MW and EX WL in double differnce from (current
*                           ambiqiuity estimates are removed)
* wls_sqr(2,max_ambs)    -- Sum of squares of MW/EX DD form
* wls_all(2,max_ambs)    -- Estimates of residual MW and Extra-wide lanes after
*                           removing the ambiguities in ambiq_all. (double diff)
*                           (Values are relative to current ambiq_all estimates).
* wls_sig(2,max_ambs)    -- Standard deviation of the residual widelanes.
*                           (Estimates based on RMS and number of data).
* curr_sdmw(max_ambs), curr_sdex(max_ambs) -- Current site-site single difference
*                           MW and EX widelanes (used for slip detection).
* curr_sdob(4,max_ambs)  -- Single difference residusls at current epoch
* curr_dd(2,max_ambs)    -- Current DD values (used for looking for slips) (LC or L1/L2)
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
* kf_prev_mjd, kf_curr_mjd -- MJD of previous epoch and of current epoch
*      used to to compute time difference for process noise
* kf_step  -- Time step in seconds 
 

      real*8 ambiq_all(2,max_ambs), wls_all(2,max_ambs),
     .       ambiq_float(2,max_ambs), ambiq_flsig(2,max_ambs),
     .       wls_sum(2,max_ambs), wls_sqr(2,max_ambs), 
     .       wls_sig(2,max_ambs),
     .       data_start(max_site), data_end(max_site), 
     .       ref_start, ref_sec, usr_start, usr_interval,
     .       kine_out(3,max_kine)
     .,      curr_sdmw(max_ambs), curr_sdex(max_ambs)
     .,      curr_sdob(4,max_ambs), curr_dd(2,max_ambs)
     .,      kf_prev_mjd, kf_curr_mjd, kf_step

* wls_num(max_ambs) -- Number of values in dWL estimate sums 

      integer*4 wls_num(max_ambs)

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
* bf_ents(6,max_ambs) -- Information about the bias flags
*     Index  Meaning
*         1  site number
*         2  PRN
*         3  Starting Epoch
*         4  Ending Epoch
*         5  Bit 1 set when we can use this bias 
*            Bit 2 set once the bias is finally fixed with
*            Bit 3 set if jump in wide lanes
*         6  Last ending edit (replces 4 when data edited during filter).
* wls_ref(3,max_ambs) -- Other ambiquities used in forming the wide
*            lanes. (These are the numbers from bf_ents)
* wls_obn(4,max_ambs) -- Observation number for current epoch of the
*            data in the double difference. 1,2 are SV1, site 1 and 2
*            3,4 are SV2 site 1,2
* ref_sv(max_site) -- Reference satellite for double differences.  
*     Each site may have its own reference to base station
* ssblk(max_site), seblk(max_site) -- Start and stop observation numbers
*     for each site


      integer*4 prn_used(max_prn), num_prn, num_site, num_kine,
     .          num_epochs, num_ambs, bf_ents(6,max_ambs),
     .          wls_ref(3,max_ambs), wls_obn(4,max_ambs), 
     .          ref_sv(max_site), ssblk(max_site), seblk(max_site)

* amb_dd(max_ambs)  -- number of dd for each ambiguity value (can be different
*     due to editing) 
* dd_amb(max_obs)   -- Ambiquity for each double difference (reverse lookup of 
*     amb_dd): NOTE: These pointers change with each observation type so use with
*     caution since after formdd analysis types (upto 6 L1,L2,LC,P1,P2,PC) values
*     point to last value (Other can be computed from number of DD's)
      integer*4 amb_dd(max_ambs), dd_amb(max_obs)

* max_gap  -- Largest gap allowed
* min_good -- Minium number of good data needed for a bias flag to
*             be keep.
* data_mask  -- Bit pattern set flag bad data
* usr_nepochs -- Number of epochs specified by the user.  All possible
*     epochs are used if not given

 
      integer*4 max_gap, min_good, data_mask, usr_nepochs

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

* clk_error -- Clock error at current epoch (r*4 version save by relative
*     epoch in sec_offset

      real*8 clk_error(max_site) 

*----------------------------------------------------------------------------

      common / track_com_b8 /  L1o_all_cse,  L2o_all_cse, 
     .         P1o_all_cse,   P2o_all_cse,                        
     .         kine_xyz, ambiq_all, ambiq_float, ambiq_flsig,
     .         wls_all, wls_sig, wls_sum, wls_sqr, 
     .         data_start, data_end, ref_start, ref_sec,
     .         usr_start, usr_interval, tt_edit, tt_abf,  tt_rbf,
     .         kine_out, clk_error, curr_sdmw, curr_sdex, curr_sdob
     .,        curr_dd, kf_prev_mjd, kf_curr_mjd, kf_step

      common / track_com_b4 / elev_cse, az_cse, sec_offset,
     .         data_flag_cse, amb_point_cse, ctop_cse, num_chan_se,
     .         prn_used, num_prn, num_site, num_kine, num_epochs,
     .         num_ambs, data_mask, bf_ents, wls_ref, wls_obn, ref_sv,
     .         ssblk, seblk, max_gap, min_good, usr_nepochs, num_edits,
     .         ss_edit, ss_abf, ss_exclude, num_exclude, 
     ,         num_abf, ss_rbf, num_rbf, wls_num, amb_dd, dd_amb

*----------------------------------------------------------------------------

* ESTIMATION DECLARATIONS:
* -----------------------
*

* APRIORI VARIANCES AND PROCESS NOISE

* cov_parmm(max_parm, max_parm) -- Covariance matrix for the parameter
*      estimates before update
* cov_parmp(max_parm, max_parm) -- Covariance matrix for the parameter
*      estimates after data update
* sol_vecm(max_parm) -- Solution vector for parameter estimate before update
* sol_vecp(max_parm) -- Solution vector for parameter estimates after update
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
      real*8 cov_parmm(max_parm, max_parm), sol_vecm(max_parm),
     .       cov_parmp(max_parm, max_parm), sol_vecp(max_parm),
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
* exwl_jmp -- Size of EW_WL jump that will be flagged (0.28 change for
*             1-1 slip).
* dd_jmp   -- Size of largest bias fixed residual before new cycle slip
* exwl_minsig -- Minimum sigma for mean EX-WL
* exwl_scale  -- Additional baseline length dependent sigma (10^(-6))
* exwl_elev   -- Elevation angle factor (1+fact/sin(el))
* mwwl_minsig -- Minimum sigma for mean MW-WL

      real*8 mwwl_jmp, exwl_jmp, exwl_minsig, exwl_scale, exwl_elev
      real*8 mwwl_minsig, dd_jmp
 
* curr_site_xyz(3,max_site) -- Position estimate at the current
*      epoch for all of the sites.  The kinematic sites have 
*      their values stored back in to kine_xyz array.
* curr_dxyz(3,max_site) -- Current offset in coordinates from 
*      aproiri locations.
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

      real*8 curr_site_xyz(3,max_site), curr_dxyz(3,max_site),
     .       ow_part(max_parm,max_obs),
     .       ow_vec(max_obs), ow_var(max_obs), rms_dd, rms_dd_avg,
     .       rms_edtol, min_lvar, wrms_dd_site(max_site)

* apr_site(3,max_site) -- Apriori variances for the site positions
*      (currently only the kinematic sites are assumed to have
*       non-zero values)  (m**2)
* mar_site(3,max_site) -- Markov process noise on the kinematic
*       site positions (m**2/sec)
* apr_atm(max_site)    -- Apriori variance for atmpospheric 
*       delays.  (m**2)
* mar_atm(max_site)    -- Process noise for the atmospheric 
*       delays   (m**2/sec)
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

      real*4 apr_site(3,max_site), mar_site(3,max_site),
     .       apr_atm(max_site), mar_atm(max_site),mar_atm_hgt(max_site),
     .       data_var(5,max_sat), out_sig_limit, stopgo_dvar, 
     .       tu_to_ep

      integer*4 num_edtol  ! Number of bad Lx residuals before cycle 
                           ! slip added
      integer*4 min_ldd    ! Minimum of phase double differences for
                           ! estimates to be made (default 4 so that
                           ! there is some redunancy).

* num_anal_type -- Number of analysis types to be run
* num_parm    -- Number of parameters to be estimated
* num_psav    -- Saved version of number of parameters for when data
*                is edited and update iterated.
* num_dtype   -- Number of data types used in current estimation 
* neam        -- Number of ambiguities per ambiquity (1 for L1/LC and 2 for L1+L2)
* site_parn(3,max_kine) -- Parameter numbers for the kinematic site
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
* num_dd_site(max_site)  -- Number of double differences per site

      integer*4 num_parm, num_psav, num_dtype, neam, 
     .          site_parn(3,max_kine),
     .          atm_parn(max_site), kine_to_site(max_kine),
     .          ow_des(2,max_obs), sd_des(2,max_obs),
     .          dd_des(4,max_obs), dd_dsi(4,max_obs), 
     .          num_ow, num_sd, num_dd,  num_dd_avg,
     .          num_ow_by_site(max_site), num_ion, num_ddion,
     .          num_anal_type, amb_parn(2,max_ambs),
     .          amb_save(2,max_ambs),
     .          parn_to_amb(max_parm), non_amb_parm, tot_parm, 
     .          tac_parm, darecl, kine_OK(max_ep_wrds,max_kine),
     .          ow_dt(max_obs), dd_dt(max_obs), static_obs(max_site),
     .          stopgo_point(max_site), num_dd_site(max_site)

* kine_known  -- Logical that is set true once we have initial
*      kinematic site position estimates from P1 data.
* stopgo_mode -- When true, track will use the KINEMATIC/STATIC flags
*      in the rinex file to change the kinematic process noise.
* ante_off    -- Antenna information read
* atm_mtt     -- Set true to use the MTT mapping function and seasonal
*                model.  Set false to use GPT and GMF.
* noetide     -- Set true to stop the application of the Earth tide 
* set_msec_bf -- Logical set true with the SET_MSEC_BF command if
*                bias flags are to be added to msec jumps in the MW-WL
* L1_only     -- Set true if only L1 data is available
* write_pos(max_site)  -- Set true when file for site is to be written
*               (if reference site and no atmospheric delay estimate then
*                file is not written)

      logical kine_known, stopgo_mode, ante_off, atm_mtt, noetide,
     .        set_msec_bf, L1_only, write_pos(max_site)

* anal_types(max_anal_type)  -- List of analysis types to be done
* search_type  -- Type of search to be preformed
* float_type -- Type of analysis for float ambiquity analysis
*  back_type  -- Type of back solution to run (BACK commmand)
* time_unit -- String containing time unit for process noise models

      character*16 anal_types(max_anal_type), search_type, float_type,
     .             back_type
      character*8  time_unit

*----------------------------------------------------------------------
      common / track_est_b8 / cov_parmm, sol_vecm, cov_parmp, sol_vecp,
     .        cov_obs, sol_obs,  apart, kgain, obs_to_res, 
     .        amb_to_res, curr_site_xyz, curr_dxyz,
     .        ow_part,  ow_vec, ow_var,
     .        rms_dd, rms_edtol, min_lvar, rms_dd_avg, wrms_dd_site, 
     .        mwwl_jmp, exwl_jmp, dd_jmp, exwl_minsig, 
     .        exwl_scale, exwl_elev,mwwl_minsig, 
     .        cov_parm_sav, sol_vec_sav

      common / track_est_b4 / apr_site, mar_site,
     .       apr_atm, mar_atm, mar_atm_hgt, data_var, 
     .       out_sig_limit, stopgo_dvar, tu_to_ep, num_parm, num_psav,
     .       num_dtype, neam, site_parn, parn_to_amb,
     .       non_amb_parm, atm_parn, amb_parn, amb_save, kine_to_site, 
     .       ow_des, sd_des,
     .       dd_des, dd_dsi, num_ow, num_sd, num_dd, num_dd_avg, 
     .       num_ow_by_site,
     .       num_ion, num_ddion, num_anal_type,
     .       kine_known, stopgo_mode, ante_off, atm_mtt, noetide,  
     .       set_msec_bf, L1_only, write_pos, tot_parm, tac_parm, 
     .       darecl, kine_OK, ow_dt, dd_dt, static_obs, stopgo_point,
     .       num_dd_site, num_edtol, min_ldd

      common / track_est_ch / anal_types, search_type, float_type,
     .       back_type, time_unit

*----------------------------------------------------------------------------

* BIAS FIXING VARIABLES BASED ON SEARCH
* ion_lim(2,max_chan)  -- lower and upper limits on the extra-wide
*     lane (cycles) based on a model ionspheric magnitude
* rank_min_rms  -- Minimum RMS scatter to use when ranking
*     ambiquitity sets by inverse RMS**2
* comb_ivar(max_amb_save) -- Sum of inverse of RMS**2 sorted by decreasing
*    value
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
* float_limit(2) -- Limits on sigma for fixing float biases and max sigma
*     for which resolution will be attempted. 
*      Relrank_limit is used as fit quality.
* wl_fact  -- Scaling factor mixing the (res/sig)^2 + wl_fact*(wl/rms)**2 
*     in determining the fit quality
* lg_fact  -- Scaling factor mixing the ion delay variance into the fit
*     criterion (for short baselines should be 1)
* max_fit  -- Maximimum value allowed for fit parameter to have bias
*     fixed.
* min_ambsig -- Minimum value for ambiguity sigma
* min_exsig  -- Minimum value for ex-wl sigma


      real*4 ion_lim(2,max_chan), 
     .       rank_min_rms, 
     .       max_dd_elev(max_obs), stratm_dcorr, stratm_base,
     .       ion_ppm, ion_ht, saved_rank(max_ambs),
     .       relrank_limit, float_limit(2), wl_fact, lg_fact, 
     .       max_fit, min_ambsig, min_exsig

* WL_avnum -- Maximum number to be used in computing average WL sigma
* WL_mnnum -- Minimum number need to allow MW-WL fixing.
      integer*4 WL_avnum, WL_mnnum 


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
     .          num_tot_resolved, 
     .          float_iter, float_sample

*----------------------------------------------------------------------------
      common / track_sea_b4 / ion_lim, 
     .          max_dd_elev, stratm_dcorr, stratm_base,
     .          ion_ppm, ion_ht, saved_rank, relrank_limit,
     .          float_limit, wl_fact, lg_fact, max_fit, min_ambsig,
     .          min_exsig, 
     .          num_bf, bf_flags, num_tot_resolved, float_iter,
     .          float_sample, WL_avnum, WL_mnnum

*------------------------------------------------------------
* SAVED Ambiguity resolution information for status reports
      real*4 asv_res(5,max_ambs) ! Saved residuals of best estmates
     .,      asv_sig(5,max_ambs) ! Saved sigmas
     .,      asv_chi(5,max_ambs) ! chi**2 conts for each
     .,      asv_rbn(3,max_ambs) ! RelRank, Best and Next best chi**2

      real*4 psv_res(max_obs)  ! Postfit residuals for upto 6 observables
     .,      psv_sig(max_obs)  ! Sigmas for residuals
     .,      psv_ae(2,max_obs) ! Azimuth and elevation (deg) of obs

      integer*4 avbad_num(max_obs)  ! number of bad residuals
      integer*4 avbad_sse(3,max_obs) ! Site/SV and last epoch of 
                               ! bad residuals
      integer*4 avbad_tot           ! Total number of bad residuals saved

      integer*4 psv_ss(4,max_obs)   ! Site and Satellite for DD
     .,         psv_amb(max_obs)    ! Ambiguity number

      character*2 psv_dtype(max_obs)   ! Data type 


      integer*4 asv_dl12(2,max_ambs)   ! Saved L1/L2 offsets
     .,         asv_resep(max_ambs)    ! Epoch at which ambiguity resolved
                                       ! (zero if not resolved).

      logical asv_used(5,max_ambs)  ! Saved Used

      character*8 asv_fcode(max_ambs)  ! Saved F codes

      common / asv_b4 / asv_res, asv_sig, asv_chi, asv_rbn, 
     .       psv_res, psv_sig, psv_ae, 
     .       psv_ss, psv_amb,
     .       avbad_num, avbad_sse, avbad_tot,
     .       asv_dl12, asv_resep, asv_used
      common / asv_ch / asv_fcode, psv_dtype

      
*----------------------------------------------------------------------------

* COMMANDS and MISCELLANEOUS parameters declarations

    

* site_apr(3,max_site)   -- Apriori coordinates for the sites
*      (for kinematic sites, this position is approximately
*       the location at the start)
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
* init_geod(3) -- Initial geodetic coordinates for plotting
*     lat long height (deg/deg/m)

      real*8 site_apr(3,max_site), 
     .       site_int(3,max_site), site_vel(3,max_site),
     .       baselens(max_site*(max_site-1)/2),
     .       site_ep(max_site),   site_offarp(3,max_site),
     .       elev_cutoff, atm_offset(max_site), ref_xyz(3), 
     .       init_geod(3,max_site) 

* site_type(max_site)    -- Type of site for each file.  Options are
*     Type 0  -- Fixed site     F
*     Type 1  -- Kinematic site K
* rcv_type(max_site)  -- type of DCB offsets for receiver
*     P -- L1 range correction 
*     C -- L1 and L2 range correction
*     N -- No correction needed (P-code ranges)
*     Types can be found rcvant.dat 
* dcbs(max_prn) -- DCBs for each satellite (meters)
* dcb_mjd   -- Epoch for dcb values
* debug(10) -- Range of epochs over which detailed information
*     is printed based in different types:  Values are paired
*     1,2 -- Original track start and end
*     3,4 -- TrackRT values

      integer*4 status_int !  Number of epochs at which status info is 
*                   written to the summary file.

      character*8 status_type  ! Characters denoting status output types;
                 ! P -- paramters, A -- Ambiguity resolution, 
                 ! W -- widelanes, R -- residuals

      logical RESET(max_site)   ! set true with RESET command to reset filter
      logical updread ! set true when the update file is read
                      ! (get sets false when the file can't be opened)
     .,       needupd  ! Set true when batch update has been specified.

* lus  -- Logical number for summary file

      integer*4 site_type(max_site), debug(10), lus
      character*2 rcv_type(max_site)
      real*8 dcbs(max_prn), dcb_mjd


* read_site_apr(max_site) -- Logical set true when we have read in the
*     apriori coordinates for a site
* write_res   -- Set true to indicate that residuals should be output
* static(max_kine)  -- Set true if site appears to be static

      logical read_site_apr(max_site), write_res, static(max_kine)

* site_names(max_site) -- 4-character names of the sites
* swvers(max_site)     -- software versions for each receiver

      character*4 site_names(max_site)
      character*16 swvers(max_site)

* rx_file(max_site) -- Names of the observation files
* sp3_dir      -- Directory for SP3 files
* sp3_root     -- Root for sp3 names (ig for igs default)
* sp3_file  -- Ephemeris file in SP3 format.  
* log_file  -- Name of file for detailed log'ing of information
* bat_file  -- Name of file with commands in it.  Basic control file.
* ambin_file -- Name of file containing the values for the ambiquities.
*      (Format must coincide with the output from a track run with
*      same setup).
* sv_clk_file -- File containing the satellite clock information (not needed)
* dcb_file    -- Name of file with DCB values (needed for mixed receivers)
* prt_root    -- Root part of print file name so that stdout goes to a file
*               (passed in the runstring with the -n option)

      character*256 sp3_dir, sp3_file, dcb_file,log_file, 
     .              bat_file, prt_root, ambin_file, sv_clk_file

* resid_root  -- Root part of the name for the residual file
* posit_root  -- Root part of the name for the geodetic position file
* posit_type  -- Type of position output (Either GEOD NEU XYZ DHU)
* csv_root    -- Root paty of then name for comma separated values for use 
*                AMCharts.
* wls_root    -- Root for widelanes corrected with ambiquities
* rwl_root    -- Root for widelanes not corrected with ambiquities
* sum_file    -- Name of summary file
* curr_timetag -- Current string being used in file names (new file
*     opened when name changes)
* upd_file    -- Name for a batch type command file that is read
*     during the trackRT run (file is read once created and then
*     must be deleted to be re-read).

      character*256 resid_root, posit_root, posit_type, csv_root, 
     .              wls_root, rwl_root, sum_file, curr_timetag,
     .              upd_file

      real*8 file_updint   ! Interval in days between file updates
                ! i.e, file closed and new one created (maybe fractional)

* runday  -- Day number for files: Invoked with <day> in file names
* runweek -- Week number for files: Invoked with <week> in file names
*     (Both arguments are passed in the runstring).
* runstr(10) -- Update three strings that can be substituted with form
*      <S01> <S02>..  <S10> in commands.

      character*16 runday, runweek
      character*64 runstr(10)

      character*4 sp3_root

      character*1 obs_file_type(max_site) 

      common / track_misc_b8 /  site_apr, site_int, site_vel, site_ep, 
     .       site_offarp, elev_cutoff, atm_offset, ref_xyz, init_geod,
     .       file_updint, dcbs, dcb_mjd, baselens
      common / track_misc_b4 / site_type, read_site_apr,
     .       debug, write_res, lus, static, status_int, RESET,
     .       updread, needupd

      common / track_misc_ch / site_names, rcv_type,
     .       swvers, sp3_dir, sp3_file,  log_file, dcb_file, bat_file, 
     .       prt_root, ambin_file, resid_root, posit_root, posit_type,
     .       csv_root, sp3_root, obs_file_type, sv_clk_file, wls_root,
     .       rwl_root, runday, runweek, sum_file, curr_timetag, runstr,
     .       upd_file, status_type

*----------------------------------------------------------------------------

* SP3 File definitions
*
*  sp3_xyz(3,max_sat, max_eph_eps) - Earth fixed coordinates of satellites
*       (XYZ, m)
*  sp3_sig(3,max_sat, max_eph_eps) - Sigmas for SP3 file coordinates
*       (XYZ, m)
*  sp3_time(max_ephs_eps) - MJD's for the ephemeris entries
*  sp3_sec(max_ephs_eps)  - Seconds of day part of the SP3 time tags running
*       from the first day of the file.
*  sp3_clk(max_sat, max_eph_eps) - The clock corrections for the satellites


* Nav file information (if need)
*  nav_clk(3,max_sat, max_eph_eps) - Clock poly nominal information from
*                                    the navigation file
*  nav_time(max_eph_eps)           - MJD's for the times given in the nav
*                                    file
      real*8  sp3_xyz(3,max_sat, max_eph_eps), 
     .    sp3_sig(3,max_sat, max_eph_eps), 
     .    sp3_time(max_eph_eps), sp3_sec(max_eph_eps),
     .    sp3_clk(max_sat, max_eph_eps),
     .    nav_clk(4,max_eph_eps, max_sat) 

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
      real*8  sv_ndat(max_site),sv_mdat(max_site),
     .     sv_dattyp(10,max_site),sv_offarp(3,max_site)

* num_sat -- Number of satellites
* num_sp3 -- Number of entries in SP3 file
* num_nav(max_sat) -- Number of entries in the NAV file for each PRN
* prns(max_sat) --  List of PRN's
* weekno        -- GPS Week number
* sp3_refmjd    -- Reference MJD for the SP3 file (i.e., MJD of first
*     day in file)

* prn_blk(max_sat)  -- Block numbers for satellites
 
      integer*4 num_sat, num_sp3, num_nav(max_sat), prns(max_sat),
     .          weekno, sp3_refmjd, prn_blk(max_sat)
     
      common / sv_com_b8 /  sp3_xyz,  sp3_sig, sp3_time, sp3_sec,
     .       sp3_clk, svs_clk, 
     .       svs_xyz, svs_loc, nav_clk,  sv_ndat, sv_mdat, 
     .       sv_dattyp, sv_offarp 

      common / sv_com_b4 / num_sat, num_sp3, num_nav, prns, weekno,
     .       sp3_refmjd, prn_blk

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
* ref_rel_hum  -- Reference relative humidity (default 0.0)
      
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
      integer*4 max_dna         ! Max SVS dnadir values

      parameter ( max_antmods = 10 ) 
      parameter ( max_zen = 91  )
      parameter ( max_az  = 361 )
      parameter ( max_dna = 15  )

      character*256 antmod_file(max_antmods)
      character*20  ant_name(max_site)

      integer*4 num_antmods   ! Number of antmod files specified

* Satellite information
* num_svs_dph(2,max_sat) -- Number of phase entries in dnadir and azimuth directions
*      for satellites
      integer*4 num_svs_dph(2,max_sat)
* svs_dphs(max_dna,1, 2,max_sat) -- Phase center adjustments at function of Nadir, azimuth,
*      L1/L2 by satellite (mm)
* svs_dna(3,2,max_sat) -- zen1, zen2, dzen, az1, az2, daz for each satellite (deg)
* svs_L12(3,2,max_sat) -- "North, East, Up" offsets at L1/L2 for each satellite (mm)

      real*4 svs_dphs(max_dna,1, 2,max_sat)
      real*8 svs_dna(3,2,max_sat), svs_L12(3,2,max_sat)

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
      common / track_phs_b4 / num_antmods, num_svs_dph, num_sit_dph,  
     .         svs_dphs, sit_dphs
      common / track_phs_r8 / svs_dna, svs_L12, sit_dzn, sit_L12 

*****************************************************************************

* Variable for trackRTr which reads rinex files

      integer*4 max_rxtype  ! Max rinex observable types
   
      parameter ( max_rxtype = 20 ) 

      integer*4 rx_ndat(max_site)  ! Number of data types in each rinex file
     .,         rx_dt(max_rxtype, max_site) ! Data codes for each RX file

      logical rx_not_open    ! Indicates the RX files are not open yet

      character*256 rx_file(max_site) !  Names of the observation files

      common / RTr_I4 / rx_ndat, rx_dt, rx_not_open
      common / RTr_CH / rx_file

*****************************************************************************

* Variables for getting command line from C into fortran
      character*4 ref_code       ! Code for reference site
      character*4 prc_code(max_site)  ! Codes for sites to be processed
      integer*4 num_prc          ! number of sites to process
      integer*4 rt_port          ! Port for comms (not needed in fortran)
      character*80 rt_machine    ! Machine with comms port (not needed in fortran)


      common / RTComl_CH / ref_code, prc_code, rt_machine
      common / RTComl_I4 / num_prc, rt_port



*******************************************************************************


