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
      parameter ( max_ent = 1500000)
c                                   ! extracted
c
*                                   ! width of ema field
      parameter ( num_field = 15 )
c
*                                   ! maximum number of sites allowed
      parameter ( max_site = 2000 )
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
      logical batch, sort
 
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
     .    norm_eq, b_vec, stat, line_st, line_en

      common / match / match_sep, coords, si_ep, si_cd, si_sg,
     .        sj_cd, sj_sg , sa_ep, sa_cd, sa_sg,
     .        ni, nj, na 
c
