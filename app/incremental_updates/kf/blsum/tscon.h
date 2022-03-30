*     Include for frame realization in tscon.  Program to convert timeseries
*     files and re-realize the reference frame.

      integer*4 max_np       ! Maximum number of parameters that can be
                             ! estimated in frame transformation (up 10
                             ! when CM offsets added)
      parameter ( max_np = 10 ) 

      integer*4 max_stinf    ! Maximum number of station.info entries that
                             ! be saved
      parameter ( max_stinf = 128 )
*
      integer*4 max_nonsec
      integer*4 max_hierach  ! Maximum depth for hieractal site list
      integer*4 max_stab     ! Maximum stabilzaiton slots

      parameter ( max_nonsec =  6*4096 )
      parameter ( max_hierach = 10 )   ! Allow up to 10 depth of sites
      parameter ( max_stab = 2048  )   ! Max stabilzaiton slots

      integer*4 max_repscl   ! Maximum number of scale estimates allow
                             ! (~40-years of daily values)
      parameter ( max_repscl = 15000 ) 

      character*8 stab_site(max_stab,max_hierach)  ! Contains site names
                            ! for sites stabilzation site (Changed with
                            ! MOD TAH 130331.
      integer*4 stab_len(max_stab,max_hierach)  ! Length of names (for speed).
      integer*4 nin_stab    ! Number of input stabilization sites

* num_nonsec -- Number of non-secular terms
* param_nonsec(2,max_nonsec) -- Gives the site number (1) and the
*     type of non-secular term.  The type defined are:
*     1 -- Offset and Rate change (applied after reference epoch)
*     2 -- Periodic (cosine and sine terms)
*     3 -- Exponential (exp(-t/tau))
*     4 -- Logarithmic (ln(1+t/tau))
*    >5    Not defined.

      integer*4 num_nonsec, param_nonsec(2,max_nonsec)


*   cnd_parts_bits - Bit mapped word which sets which partials are to be
*                    estimated.  Bits 1-3 -- rotations, bits 4-6 translations
*                    bit 7 - scale, bits 8-10 - CM load offsets
*   stab_it - Number of stabilization iteration
*   used_index(max_site) -- Indices to entries being used at current evaluation
*     mjd
*   num_used  -- Number of sites in list
*   stab_index(max_site) -- Indiced to stabilzaiton sites (can be changed as
*     sites are edited)
*   num_stab  -- Number of stab sites for specific day
*   nps  -- Parameter number associated with scale estimate (set whether
*           mean height or scale estimate output.

      integer*4 cnd_parts_bits, stab_it
     .,         used_index(max_site),stab_index(max_site)
     .,         num_used, num_stab, nps

*    CMest -- Logical set true is CM estimates are to be made
      logical CMest    

*   first_mjd, last_mjd -- First and last time seen
*   apr_val_site(3,2,max_site)    - Apriori values for
*               - site XYZ and XYZ rate values (m and m/year).
*   site_epoch(max_site)  -- MJD for site position epoch
*   apr_val_nonsec(8,max_nonsec)  - Non-secular parameters for station
*                  positions.  The order is:
*                  JD -- Julian date either for periodic 0 phase
*                        of start of exponential/logarithm/steps
*                  parameter -- Either period or delay times (days)
*                  Paired for X Y and Z (read in as NEU and converted)
*                  coef 1    -- Exponential/logarithm constant;
*                               cosine term or offset (m)
*                  coef 2    -- skipped for exponential and log,
*                               sine term (m) or rate change (m/yr)
*                  (see also I*4 param_nonsec(2,max_nonsec)

*   cnd_hgt_var -- Height variance ratio
*   use_ratio   -- Ratio of difference between median and best to allow
*                  the worst to be used.
*   stab_nsig   -- N-sigma editing condition
*   stab_rel    -- Relative weight between constant and site dependent
*   stab_min_dh -- Minimum height dsigma
*   stab_min_dne -- Minumum NE dsigma
*   ref_coords(3,max_site) -- Reference site coordimates at ejmd
*   sigma_scale  -- Scale factor for sigmas; either -s=<value> in command
*                  line, or SIGSCALE command in command file.
*   units_scale(max_np)  -  Scaling for output units in out_org file

      real*8 apr_val_site(3,2,max_site)
     .,      site_epoch(max_site)  
     .,      apr_val_nonsec(8,max_nonsec)
     .,      first_mjd, last_mjd
     .,      cnd_hgt_var, use_ratio, stab_nsig
     .,      stab_rel, stab_min_dh, stab_min_dne
     .,      ref_coords(3,max_site), sigma_scale
     .,      units_scale(max_np) 

      real*8 stab_times(2,max_site)  ! Limitation on times over which
                      ! a site should be used in stabalization.

      logical debug   ! Set true for debug output
      logical jpl_xyz ! Set true if JPL XYZ format
      logical mea_xyz ! Set true if new Measures XYZ format with correlations
* MOD TAH 200720: Added new mea_raw type
      logical mea_raw ! Set true if 20/07 Measures TrendXYZ file (XYZ with no
                      ! offsets or other terms removed).
      logical mea_neu ! Set true if Measures NEU format
      logical rep_scale  ! Set true to reapply the scale corrections based 
                      ! on repscl file derived from GIPSY xfiles.
      logical fix_unre   ! Set true to fix the UNR East error.
      logical worg    ! Set true if the origin parameter estimates are to be
                      ! writen out

      character*8 gsite_names(max_site)  ! All site names
     .,   site_name   ! Current site name extracted from file name

      character*16 stab_hfs(max_stab) ! Restricted name of for stab_site list
                      ! this string needs to appear the .pos product name.
      character*16 runstring_refframe  ! Refrence frame from runstring

* Added values to fix UNR east errror
      integer*4 num_stinf    ! Number station.info lines found
      integer*4 stinf_sta(max_stinf)  ! Status of getting antenne NEU values
                      !  0 - not found yet, -1 antenna with NONE radome found
                      ! +1 - Antenna and rsdome found.

      real*8 stinf_mjd(max_stinf)         ! MJD of antenna change
      real*8 stinf_neu(3,2,max_stinf)     ! Antenna NEU offsets for G01 and G02 
                                          ! for each antenna entry

* Added values to add back UNR scale estimates
      real*8 repscl(3,max_repscl)  ! Scale estimates from GIPSY x-files
                      ! Values are MJD, scale scale+- (sigma not used)
      integer*4 num_repscl   ! Number of scale estimates read.

      character*16 stinf_ant(max_stinf)   ! Antenna name
      character*4  stinf_rad(max_stinf)   ! Radome name



      character*256 comfile  ! Command file
     .,    apr_file   ! Names of apriori coordinate file
                      ! (Included EXTENDED ENTRIES).  Can used multiple times
     .,    out_org    ! Name of output file with estimated pos_org parameters.
                      ! (Added TAH 170720).
     .,    repscl_file ! Name of file containing scale corrections (should be 
                      ! be consistent with input time series)    
     .,    stinf_file    ! Station.info file name
     .,    antmod_file   ! Antmod.data file name

      common / xyzI4 / stab_len, num_stab, num_nonsec, param_nonsec
     .,      cnd_parts_bits
     .,      stab_it, used_index, num_used, stab_index, nin_stab
     .,      debug, jpl_xyz, mea_xyz, nps, CMest, worg, mea_neu, mea_raw
     .,      rep_scale, fix_unre, num_stinf, stinf_sta, num_repscl

      common / xyzR8 / apr_val_site,site_epoch, apr_val_nonsec
     .,      first_mjd, last_mjd, cnd_hgt_var, use_ratio, stab_nsig
     .,      stab_rel, stab_min_dh, stab_min_dne, ref_coords
     .,      stab_times, sigma_scale, units_scale 
     .,      stinf_mjd, stinf_neu, repscl

      common / xyzCH / gsite_names, site_name, comfile, apr_file 
     .,       stab_hfs, stab_site, runstring_refframe, out_org
     .,       repscl_file, stinf_file,antmod_file 
     .,       stinf_ant, stinf_rad


