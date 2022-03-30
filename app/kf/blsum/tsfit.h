*     Include file for definitions in tsfit
      include 'tssum.h'
*
*     Parameters
      integer*4 max_par   ! Maxumum of paramaters to be estimaated
      integer*4 max_per   ! Maximum number of periodic terms
      integer*4 max_pat   ! Maximum number of site name patterns
      integer*4 max_eq    ! Maximum number of earthquakes
      integer*4 max_rn    ! Maximum number of renames
      integer*4 max_off   ! Maximum number of offsets
      integer*4 max_ptyp  ! Maximum number of parameter types to be estimated


      parameter( max_par =    100 )
      parameter( max_per =     10 )
      parameter( max_pat =    256 )
      parameter( max_eq  =    256 )
      parameter( max_rn  = 100000 )
      parameter( max_off = 100000 )
      parameter( max_ptyp  = 9 )


      integer*4 unr          ! Unit number for reading .pos files (may 
                             ! be 100 or 101 depending on use of list file)
     .,         num_per      ! Number of periodic terms
     .,         num_pat      ! Number of site patterns
     .,         num_eq       ! Number of earhquakes
     .,         num_rn       ! Number of renames
     .,         num_off      ! Number of offsets
     .,         num_par      ! Number of parameters to be estimated
     .,         pcode(max_par) ! Parameter codes
                             !  1 -- Offset
                             !  2 -- Rate
                             !  3 -- Periodic cos
                             !  4 -- Periodic sin
                             !  5 -- Break
                             !  6 -- Earthquake break
                             !  7 -- log coefficient ln(1+dt/tau)
                             !  8 -- exponent (1-exp(dt/tau))
                             !  9 -- Drate (when rate change after EQ)
     .,         rcode(max_ptyp) ! user defined list of parametes to be removed from output residuals	     
     .,         pindx(max_par) ! Index parameter that points to specific
                             ! Periodic, break, log or exponent term
     .,         ts_used(3)    ! Used number of data
     .,         min_sigscale  ! Minimum number of degrees of freedom to allow
                             ! NRMS rescaling
     .,         min_rsig      ! Minimum number needed for real sigma (30 default)


      integer*4 uns   ! Summary file unit number (200)
     .,         unv   ! Velocity file unit (202)
     .,         una   ! Output apriori file unit number (203)
     .,         unp   ! Output PBO-format velocity file unit number (206)
     .,         unt   ! Unit number for report edits as renames (204)
     .,         unf   ! Unit numebr for position offsets summary (207)
     .,         une(4,max_eq)  ! Unit numbers of earthquake output;
                      ! 1 for offset, 2 log and 3 exponent, 
                      ! 4 total accumulated motion


      logical  log_eq(max_eq)  ! Set true for log earthquake
     .,        exp_eq(max_eq)  ! Set true for exp earthquake
     .,        dtr_eq(max_eq)  ! Set true if rate change at time of eq
     .,        eq_out          ! Set true if logs to output
     .,        mean_only       ! Set true to estimate mean and no linear trend
     .,        resid_sum       ! Set true to write residual to SUM file (used with diff files).
     .,        pbov_header     ! Set true after pbovel header written


      real*8    per_per(max_per)  ! Periods of periodic terms (days)
     .,         per_sig(max_per)  ! Apriori sigma for periodic terms (mm)
     .,         av_res            ! Length of time to average residuals overs to
                                  ! better statistics
     .,         max_persig, max_eqsig  ! Maximum values for sigmas of periodic
                                  ! and post-seismic terms for them to be estimated
     .,         max_rmchi         ! Max chi allowed for periodic and log terms
                                  ! divided by their sigma for the parameter to
                                  ! be removed, ie. if sigma is large, but the 
                                  ! estimate significant, the parameter will not
                                  ! be removed. 
     .,         htscale           ! scale for height sigma for above testes

*   eq_pos(3,max_eq) - Cartesian position of each earthquake (m)
*   eq_rad(max_eq)   - Radius of influence of each earthquake (m)
*   eq_depth(max_eq) - Depth of Earthquake (used to scale the spatial
*                      dependent part of co-, pre- and post-seismic
*                      process noise (m).
*   eq_epoch(max_eq) - Epoch of each earthquake (Julian date)

*   eq_apr_coseismic(6,max_eq) - Apriori sigma of coseismic displacement
*                      in North, East and Up (m), and spatial dependent
*                      apriori sigma (NEU) (m at epicenter and formed
*                      as (d/l)**2 where l is distance and d is depth
*   eq_log_sig(6,max_eq)  -- Apriori sigma on log estimates (m) constant 
*                      term (NEU) and (D/R)**2 dependent terms 
*   eq_log_tau(max_eq)    -- Time constant for log estimates (days)
*   eq_exp_sig(6,max_eq)  -- Apriori sigma on log estimates (m) constant 
*                      term (NEU) and (D/R)**2 dependent terms 
*   eq_exp_tau(max_eq)    -- Time constant for log estimates (days)

      real*8  eq_pos(3,max_eq), eq_rad(max_eq), eq_depth(max_eq)
     .,       eq_epoch(max_eq) 
     .,       eq_dur(2, max_eq)
     .,       eq_mar_pre(6,max_eq), eq_mar_post(6,max_eq)
     .,       eq_apr_coseismic(6,max_eq) 
     .,       eq_log_sig(6,max_eq),eq_log_tau(max_eq)
     .,       eq_exp_sig(6,max_eq),eq_exp_tau(max_eq)

      logical eq_rename(max_eq)

* MOD TAH 170406: Added scale factor for eq_rad value.  Default is 1 and
*      value optional with eq_file command
      real*8 eq_rad_scale   ! Scale factor of eq_radius values

* min_rwvar  -- Mininum process noise when Kalman filter used (default
*     0.05 mm^2/yr same as sh_gen_stats
      real*8 min_rwvar
      real*8 RWvar(3)   ! Random walk noise (m^2/yr)
      real*8 prechi(3)  ! Prefit chi*2/f for Kalman filter.

* Fitting parameters
      real*8 norm_eq(max_par,max_par) ! Normal equations
     .,      bvec(max_par)     ! Solution vector
     .,      soln(max_par,3), solsig(max_par,3) ! Solution and sigma in mm
     .,      smet(max_par,3), smtsig(max_par,3) ! Solution and sigma in meters
     .,      solcov(max_par,max_par,3)  ! Full covariance of solution by comp.
     .,      apr_con(max_par,3)  ! Aprori variance constrains on parameters.
     .,      stats(4)          ! Statistics array (res*wgh, wgh, res^wgh, #) 
     .,      apart(max_par)    ! Parials for observation
     .,      cen_mjd           ! Center time for time series
     .,      outlog_days       ! Days after earthquake log output

       integer*4 pflag(max_par)  ! Flag for parameters. 
                                 ! 0 -- no data
                                 ! 1 -- data available to estimate
       integer*4 pnm_indx(max_par)  ! For each site name this pionts to the
                                 ! offset parameter for this site.
     .,          pnm_rn(max_par) ! Rename number that goes with each site name
                                 ! (Allows break names to be saved for output)

* Editing
      real*8 max_sigma(3)      ! Maximum sigma allowed
     .,      nsigma            ! If greater than 0, n-sigma delete condition.
     .,      sig_scale(3)      ! NRMS or scaled NRMS for real-sigma
     .,      taufin(3)         ! Final time constants with realsigma
     .,      wrms(3)           ! WRMS (mm) for each component.
     .,      wn_nrms(3)        ! NRMS from White noise model
     .,      mean_mjd(3)       ! Mean date for each component estimate
     .,      time_range(2)     ! MJD for start and stop of data to be used in fits


* Site renaming arrays
*   rn_times(2, max_rn) - Epochs over which rename should occur (start and
*                      stop Julian dates)
*   rn_dpos(3, max_rn)  - dX,dY,dZ postion change for renamed site (added
*                      to the estimate (m)).

      real*8 rn_times(2, max_rn), rn_dpos(3, max_rn),
     .       off_times(max_off), off_dpos(3, max_off)


      logical   use_constraints   ! Set true to use eq_file eq constraints
     .,         report_edits      ! Set true to report edits on time series
     .,         report_renames    ! Reports re-names from editing
     .,         report_constraints ! Set true to reoprt eq_constraints
     .,         real_sigma        ! Use realistic sigma algorithm

*   eq_codes(max_eq) - Two letter codes for each Earthquake (saved in a
*                      character*8 array for 32Bit manipulation)
*   rn_codes(2,max_rn) - Old and new station names.
*   rn_types(max_rn)   - Type of position change (XYZ/NEU)
* MOD TAH 971112:
*   rn_hfiles(max_rn)  - String to checked to see if glb_inp_file contains
*                        this string.  Restricted to 16 characters.

      integer*4 lrn_hf(max_rn)  ! Lengths of hfile renames.
      integer*4 num_pnm         ! Number of site names needed for 
                                ! parameter estimates.
      real*8 time_pnm(2,max_par)  ! Times over which parameter site name
                                ! applies
      real*8 time_act(2,max_par)  ! Actual start and stop times for break (missing data
                                ! can result in end not matching start of next segment).
      real*8 soff_pnm(3,max_par)  ! Sum of offsets for breaks which go with each
                                ! of the num_pnm site names.  These are the sums
                                ! that are added to original aprori.  The EXTENDED
                                ! OFFSET lines are changes at each epoch

      character*8 pnm_site(max_par)  ! Name of the sites needed.  Names
                         ! get associated with  each parameter 
                         ! (renames and eqs).
      character*8 pnm_type(max_par)   ! Type of break for each site name.
      

      character*8 eq_codes(max_eq), rn_codes(2,max_rn), rn_types(max_rn)
      
      character*8 off_codes(max_off), off_types(max_off)

      character*4 kfopt   ! Kalman filter option (Blank of none, RW for
                          ! random walk, FOGM First order Gauss Markov
                          ! (correlation time from RealSigma), WH white 
                          ! noise (no correlated noise)

      character*16 rn_hfiles(max_rn)

      character*24 site_pat(max_pat) ! Patterns for using a site

      character*256 cmdfile  ! Name of command file (maybe NONE for no comamnd
                             ! file).
     .,             eqfile   ! Name of eq file being read
     .,             sumfile  ! Name of summary file.
     .,             infile   ! name of position file being processed
     .,             listfile ! name of file with list
     .,             detroot  ! Root for output for each site
     .,             velfile  ! Name of velocity file to be output
     .,             outapr_file ! Name of output apriori file
     .,             eqo_root ! Root for output earthquake file
     .,             rename_file  ! File with edit renames
     .,             resroot  ! Root for output residual files.
     .,             out_pbovel ! Name of output PBO-format velocity file
     .,             out_posf   ! Output position adjustment in .vel format

*-----------------------------------------------------------------------
* Common blocks

      common / tf_i4 / unr, num_per, num_pat, num_eq, num_rn, num_off,
     .       lrn_hf, eq_rename,
     .       use_constraints,  report_edits, report_renames,
     .       report_constraints,real_sigma, num_par, pcode,
     .       pindx, ts_used, min_sigscale, min_rsig, log_eq, exp_eq,  
     .       dtr_eq, eq_out, mean_only, resid_sum, 
     .       uns, unv, unf, une, unt, una, unp, 
     .       num_pnm, rcode, pflag, pnm_indx, pnm_rn, pbov_header



      common / tf_r8 / per_per, per_sig, av_res, 
     .        max_persig, max_eqsig, htscale, max_rmchi,
     .        eq_pos, eq_rad, eq_depth, eq_epoch, eq_apr_coseismic, 
     .        eq_log_sig,eq_log_tau, eq_exp_sig, eq_exp_tau,
     .        eq_dur, eq_mar_pre, eq_mar_post, eq_rad_scale, 
     .        rn_times, off_times, rn_dpos, off_dpos, solcov, 
     .        norm_eq, bvec, soln, solsig, smet, smtsig, apr_con, 
     .        stats, apart, cen_mjd, outlog_days, 
     .        max_sigma, nsigma, sig_scale, taufin, wrms, wn_nrms,
     .        mean_mjd, time_range, time_pnm, time_act, soff_pnm, 
     .        min_rwvar, RWvar, prechi

      common / tf_ch / eq_codes, rn_codes, rn_types, rn_hfiles,
     .      off_codes, off_types, 
     .      site_pat, cmdfile ,eqfile, sumfile , infile ,listfile,
     .      detroot, velfile, outapr_file, pnm_site, pnm_type, 
     .      eqo_root, rename_file, resroot, out_pbovel, out_posf, 
     .      kfopt

*-----------------------------------------------------------------------

