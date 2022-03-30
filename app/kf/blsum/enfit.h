c
c
c Start with parameter definitions
c --------------------------------
c max_ent -- the maximum number of baselines we can extract
c num_field -- the with of the ema field (15 integer*4 words)
c max_site -- the maximum number of sites that can be encountered
c max_files -- maximum number of solutions files which can be processed
c max_runstring -- maximum number of file names that can be passed in the
c     runstring
c max_ts  -- Maximum length of single timeseries
      integer*4 max_ent, max_ts

c
      integer*4 num_field, max_site, max_files, max_runstring
 
c
c tol -- the tolerance for the inversion of normal equations
c
      real*8 tol
 
c
*                                   ! number of baselines to be
      parameter ( max_ent = 5000000)
c                                   ! extracted
      parameter ( max_ts = 10000 )  ! Over 20 years of daily values.
c
*                                   ! width of ema field
      parameter ( num_field = 15 )
c
*                                   ! maximum number of sites allowed
      parameter ( max_site = 8192 )
c
*                                   ! maximum number of solution files
      parameter ( max_files = 5000 )
c
*                                                 ! maximum number of
      parameter ( max_runstring = max_files + 2 )
c                                     file names in runstring
c
*                                   ! tolerance on inverting normal
      parameter ( tol = 1.d-12   )
c                                   ! equations
c
c Output devices
c --------------
c solution_file -- the list of the solution files
c summary_file  -- file for summary of results
c values_file   -- file for saving baseline values
c
      character*64 solution_file(max_files), summary_file, values_file
 
c
c Solution information
c --------------------
c num_files -- actual number of files to be processed
c header    -- header records from the solution files
c batch     -- mode of processing
c sort      -- Indicates if we should sort data (SORT in runstring)
c num_soln  -- solution number form SOLVK
c max_soln  -- maxiumum number of solutions found
c min_out   -- Minumum number of values needed for output to
c              to the summary file
c
      integer*4 num_files, num_soln, max_soln, min_out    
c
c realsigma -- Set true with -RS option for realistic sigmas
      logical batch, sort, realsigma 
 
c
      character*160 header(max_files)
 
c
c site information
c ----------------
c site -- the site numbers for a baseline entry
c sol_epoch  -- epoch of the baselines being read
c base_data  -- the baseline data (length and sigma)
c names -- names of the sites
c num_site -- number of sites
c num_ent -- number of baseline entries found
c
      integer*4 num_site, site(2)
 
c
      real*8 sol_epoch, base_data(2)
 
c* Times(max_ts), Res(max_ts), Rerr(max_ts) -- Single component 
*     times, residuals and sigmas
      real*8  Times(max_ts), Res(max_ts), Rerr(max_ts)

      integer*4 num_ent
 
c
      character*8 names(max_site)
 
c
c information for estimating mean and slope
c -----------------------------------------
c ref_epoch -- the reference epoch for fitting slope
c ref_length -- the reference length
c min_epoch  -- julian date of first experiment
c max_epoch  -- julian date of last experiment
c mean_epoch -- julian date of mean epoch of experiments
c wgh_mean   -- weighted mean of baseline lengths minus reference
c               length (meters)
c sig_mean   -- sigma of mean (formal, meters)
c intercept  -- intercept at reference time minus ref_length (m)
c slope      -- slope (m/yr)
c sig_slope  -- sigma of slope (m/yr)
c
      real*8 ref_epoch, ref_length, min_epoch, max_epoch, mean_epoch,
     .    wgh_mean, sig_mean, intercept, slope, sig_slope
 
c
c Statistics variables
c --------------------
c wrms_mean -- weighted rms scatter of baseline lengths about the mean (m)
c nrms_mean -- normalized rms scatter
c wrms_slope -- weighted rms scatter about slope (m)
c nrms_slope -- normalized rms scatter about slope
* sig_mul    -- Multiplier for sigmas of parameters. Equals nrms if nrms>1
*               and equals 1 if nrms < 1.
c
      real*8 wrms_mean, nrms_mean, wrms_slope, nrms_slope, sig_mul
 
* num_ts -- number of values in time series
      integer*4 num_ts

c
c Variables used in getting solution
c ----------------------------------
c stat    -- array for accumulating statistics
c
      real*8  stat(4,2)
 
c
c line_st -- start line number for tplot
c line_en -- end line number for tplot
c
      integer*4 line_st, line_en

* MOD TAH 981229: New variables added to allow for mapping of
*     coordinates from one site to another

* max_match -- Maximum number of match values to average
      integer*4 max_match

      parameter ( max_match = 128 )

* match_sep -- Separation of sites in meters for sites to be
*     alligned

* coords(3,max_site) -- Coordinates of sites (based on first entries
*     found)

      real*8 match_sep, coords(3,max_site)

* si_ep(max_match)    -- Epochs of site i for match
* si_cd(3, max_match) -- Coordinates for site i
* si_sg(3, max_match) -- Sigmas for coordinates of site i

* sj_cd(3, max_match) -- Coordinates for site j
* sj_sg(3, max_match) -- Sigmas for coordinates of site j
* sa_ep(3, max_match)    -- Epochs of site j for match
* sa_cd(3, max_match) -- Coordinates for site j
* sa_sg(3, max_match) -- Sigmas for coordinates of site j


* ni, nj(3)  -- Number of entries in match for each site
* na(3)      -- Numebr of all site j entries

      real*8  si_ep(max_match), si_cd(3,max_match), 
     .        si_sg(3,max_match),
     .        sj_cd(3,max_match), sj_sg(3,max_match),
     .        sa_ep(3,max_match), sa_cd(3,max_match), 
     .        sa_sg(3,max_match) 

      integer*4 ni, nj(3), na(3)

c
c Now the common block definition
c -------------------------------
c
      common batch, sort, min_out, num_files, header,
     .    solution_file, summary_file, values_file,
     .    num_soln, max_soln, site, sol_epoch, base_data,
     .    names, num_site, num_ent,
     .    ref_epoch, ref_length, min_epoch, max_epoch, mean_epoch,
     .    wgh_mean, sig_mean, intercept, slope, sig_slope,
     .    nrms_mean, wrms_mean, nrms_slope, wrms_slope, 
     .    sig_mul, Times, Res, Rerr,
     .    stat, line_st, line_en, num_ts, realsigma

      common / match / match_sep, coords, si_ep, si_cd, si_sg,
     .        sj_cd, sj_sg , sa_ep, sa_cd, sa_sg,
     .        ni, nj, na 
cc
* MODS for ENFIT to allow more parameters to be estimates
* MOD TAH 030116: Added log fit as well as exp

* max_exp  -- Max number of exponential terms
* max_log  -- Max number of logarithm terms
* max_out  -- Max number of times to output estimates
* max_per  -- Max number of periodic terms
* max_parm -- Max number of parameters that can be estimated. 

      integer*4 max_exp, max_log, max_out, max_per, max_parm
      parameter ( max_exp  = 16 )
      parameter ( max_log  = 16 )
      parameter ( max_out  = 1000 )
      parameter ( max_per  = 16 )
      parameter ( max_parm = 2 + max_exp + max_log + 2*max_per )

* num_exp, num_log, num_per, num_parm  -- Number of exponentials, log, 
*     periodic and parameters
* num_out  -- Number of output times

      integer*4 num_exp, num_log, num_per, num_parm, num_out

* ep_exp(max_exp)  -- Epoch for exponentials
* ep_log(max_exp)  -- Epoch for logarithmic fit
* out_exp(max_out) -- Days to output exponentials (from start time
*                     of first exponential)
* out_log(max_out) -- Days to output logarithm (from start time of first log)
* tau_exp(max_exp) -- Decay times for exponentials (days)
* tau_est(max_exp) -- Current estimate of tau for current data set (days)
* exp_est(max_exp) -- Current estimate of the amplitude of the exponent
* tau_sig(max_exp) -- Uncertainity in tau (days).  Setting value will
*                     cause iterated solution for tau
* dtl_log(max_log) -- time offset for log (days)
* dtl_est(max_log) -- Current estimate of time offset for current data set (days)
* log_est(max_log) -- Current estimate of the amplitude of the log
* dlt_sig(max_log) -- Uncertainity in time offset(days).  Setting value will
*                     cause iterated solution for dlt
      real*8 ep_exp(max_exp), out_exp(max_out), tau_exp(max_exp),
     .       tau_est(max_exp), exp_est(max_exp), tau_sig(max_exp)
      real*8 ep_log(max_log), out_log(max_out), dtl_log(max_log),
     .       dtl_est(max_log), log_est(max_log), dtl_sig(max_log)

* con_exp(max_exp) -- Apriori sigmas for exponentials (mm)
* con_log(max_log) -- Apriori sigmas for logarithms (mm)
* con_per(max_per) -- Apriori sigmas for periodic terms (mm)
* per_per(max_per) -- Period of periodic terms (days)

      real*8 con_exp(max_exp), con_log(max_log), 
     .       con_per(max_per), per_per(max_per)

*
* norm_eq(max_parm,max_parm) -- Normal equations
* bvec(max_parm)             -- Solution vector
* finsol(max_parm)           -- Final solution
* apart(max_parm)            -- Partial derivatives
* sig_est(max_parm)          -- Sigmas of estimates
 
      real*8 norm_eq(max_parm,max_parm), bvec(max_parm),
     .       apart(max_parm), sig_est(max_parm), finsol(max_parm)

* est_tau(max_exp) -- Logical set true if tau is to be estimated
* est_dtl(max_log) -- Logical set true if dtl is to be estimated

      logical est_tau(max_exp), est_dtl(max_exp)

* fit_cmd_file -- Name of fit command file
* fit_val_file -- Name of output values files

* Name selection command definitions
* max_site_sel -- Max number of site selectoins allowed
      integer*4 max_site_sel
      parameter ( max_site_sel = 1024 )
* num_ss -- Number of site selections
      integer*4 num_ss
* site_sel(max_site_sel) -- Site selection code
      character*8 sel_site(max_site_sel)

      character*128 fit_cmd_file, fit_val_file 

      common / fit_com / ep_exp, out_exp, tau_exp, tau_est, 
     .      exp_est, tau_sig,
     .      con_exp, con_log, con_per, per_per,
     .      ep_log, out_log, dtl_log, dtl_est, 
     .      log_est, dtl_sig,
     .      norm_eq, bvec,  apart, sig_est, finsol,
     .      num_exp, num_log, num_per, num_parm, num_out, 
     .      est_tau, est_dtl, num_ss ,
     .      fit_cmd_file, fit_val_file, sel_site


