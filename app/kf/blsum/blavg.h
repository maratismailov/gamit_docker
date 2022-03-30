c
c@NEWCM -- the control common block for NEWBS
c
c MOD JLD 870505 To increase size of file descriptors to 64 bytes for
c                full CI file names
c
c Start with parameter definitions
c --------------------------------
c max_ent -- the maximum number of baselines we can extract
c num_field -- the with of the ema field (15 integer*4 words)
c max_site -- the maximum number of sites that can be encountered
c max_files -- maximum number of solutions files which can be processed
c max_runstring -- maximum number of file names that can be passed in the
c     runstring
c
      integer*4 max_ent
 
c
      integer*4 num_field, max_site, max_files, max_runstring
 
c
c tol -- the tolerance for the inversion of normal equations
c
      real*8 tol
 
c
*                                   ! number of baselines to be
      parameter ( max_ent = 3000000 )
c                                   ! extracted
c
*                                   ! width of ema field
      parameter ( num_field = 10 )
c
*                                   ! maximum number of sites allowed
      parameter ( max_site = 1000 )
c
*                                   ! maximum number of solution files
      parameter ( max_files = 1  )
c
*                                                 ! maximum number of
      parameter ( max_runstring = max_files + 2 )
c                                     file names in runstring
c
*                                   ! tolerance on inverting normal
      parameter ( tol = 1.d-12   )
c                                   ! equations
c
c Now the common block definition
c -------------------------------
c
      common icrt, iprnt, batch, num_files, header,
     .    solution_file, summary_file, values_file,
     .    num_soln, max_soln, site, sol_epoch, base_data,
     .    names, num_site, num_ent,
     .    ref_epoch, ref_length, min_epoch, max_epoch, mean_epoch,
     .    wgh_mean, sig_mean, intercept, slope, sig_slope,
     .    nrms_mean, wrms_mean, nrms_slope, wrms_slope,
     .    norm_eq, b_vec, stat, line_st, line_en
 
c
c Output devices
c --------------
c icrt -- user terminal
c iprnt -- the print device
c solution_file -- the list of the solution files
c summary_file  -- file for summary of results
c values_file   -- file for saving baseline values
c
      integer*4 icrt, iprnt
 
c
      character*64 solution_file(max_files), summary_file, values_file
 
c
c Solution information
c --------------------
c num_files -- actual number of files to be processed
c header    -- header records from the solution files
c batch     -- mode of processing
c num_soln  -- solution number form SOLVK
c max_soln  -- maxiumum number of solutions found
c
      integer*4 num_files, num_soln, max_soln
 
c
      logical batch
 
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
 
c
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
c
      real*8 wrms_mean, nrms_mean, wrms_slope, nrms_slope
 
c
c Variables used in getting solution
c ----------------------------------
c norm_eq -- normal equations
c b_vec   -- b_vector for solutions
c stat    -- array for accumulating statistics
c
      real*8 norm_eq(3), b_vec(2), stat(4,2)
 
c
c line_st -- start line number for tplot
c line_en -- end line number for tplot
c
      integer*4 line_st, line_en
 
c
*     Values needed for averaging data
*      av_norm(3)       - Normal equations
*      av_b(2)          - B values
*      av_epoch         - Epoch of average
*      av_length        - Average length
*      av_slope         - Average slope
*      av_lsig          - Length sigma
*      av_ssig          - Slope sigma
*      start_ep         - Start epoch for averaging
*      stepep           - Length of averaging interval.
*      curr_ep          - Epoch of start of current block
*      ref_av_ep        - reference epoch for averaging
*      ref_av_bl        - reference baseline length
 
      real*8 av_norm(3), av_b(2), av_epoch, av_length, av_slope,
     .    av_lsig, av_ssig, start_ep, stepep, curr_ep, ref_av_ep,
     .    ref_av_bl

*      av_rms    - Sum of residuals**2 in the duration of the
*                  averaging
*      rms       - Sum of resiudals**2 in the total span
*      av_bm     - Bvec for finding mean after removeal of total
*                  offsettand slope.
*      num_tot   - Total number in findinmg the slope through
*                  the averages.
*      av_wrms   - WRMS of avergaes
*      av_nrms   - Normalized RMS of averaes
      real*8 av_rms, rms, av_bm, av_wrms, av_nrms
      integer*4 num_tot
 
*         av_num        - Number of values in current average
*         min_av_num    - Minimum number of values need for
*                         average to output.
 
      integer*4 av_num, min_av_num
 
*     Common block needed
 
      common / BLAVG_COM / av_norm, av_b, av_epoch, av_length,
     .    av_slope, av_lsig, av_ssig, start_ep, stepep, curr_ep,
     .    ref_av_ep, ref_av_bl, av_rms, rms, av_bm, av_wrms,
     .    av_nrms, av_num, min_av_num, num_tot
 
