 
*---------------------------------------------------------------------
*                                                   SOLVK_MARKOV.FTNI
*     This is the block of the solvk common which contains the information
*     about the markov paramters and parmeter numbers.  This information
*     is not needed to run the FORSL or the BAKSL programs.
*
*     T.Herring                   10:59 PM SAT., 17 Jan., 1987
*--------------------------------------------------------------------
 
*   markov      - Header word for start of block
 
      integer*4 markov
 
*   atm_parn(2,max_sites)   - the parameter numbers of the
*               - atmospheric zenith delay and zenith rate
*               - parameters for each site.  Value is set to
*               - zero if parameter not estimated.
*   atm_az_parn(2,max_sites)    - the parameter numbers for
*               - atmospheric asymetric delay parameters.
*               - These parameters represent NS and EW values.
*               -(There is not rate associated with these values)
*   atm_map_parn(max_sites)     - the parameter number for mapping
*                 function adjustments
*   bclk_parn(max_bclk) - Parameter numbers for the baseline 
*               - dependent clocks.
*   bclk_bas_num(2,max_bclk) - Baselines (given as site 1 and site 2)
*               - for baseline dependent clocks
*   num_bclk    - Number of baseline dependent clocks in solution.
*   axo_parn(max_sites)     - the parameter numbers for the axis
*               - offset corrections.
*   crt_outopts - The output options for the solution to be
*               - printed to the users terminal.  This word is
*               - bit mapped with the following meanings:
*               - Bit   Meaning
*               -   1   Print the correlation matrix.
*               -   2   Print the Markov statistics used in this
*               -       solution.
*               -   3   Print the availability of propagation
*               -       medium calibrations for this solution.
*               -   4   Print the availability of contributions
*               -       for this solution
*               -   5   Print the soucre structure information
*               -       for this solution
*               -   6   Print the DELETE COUNTS for this solution.
*               - > 7   NOT USED
 
*   default_base_cont     - Default setting for baseline contributions
*   default_data_mask     - Default value for data masking for edits
*   default_ion_mask      - Default value for ion masking for edits
*   default_site_cont     - Default setting for site contributions
 
*   etd_ext_parn(12,max_sites)   - The parameter numbers for the
*               - extended earth tide partials
 
*   gam_parn    - parameter number for gamma.
 
*   nut_parn(2) - Parameter numbers for nutation in longitude
*               - obliquity.
 
*   plx_parn(max_sources)   - Parralax partials numbers for each
*               - source.
 
*   pmu_parn(2,3) - Parameter numbers for X,Y pole position
*               - and UT1. Value and rate of change (mas, mas/1000 sec)
*               - NOTE FORM: Value and rate and then down by x,y, and UT1
*   prt_outopts - The output options for the printer (if different
*               - to the users crt.) See CRT_OUTOPTS for meaning
*               - of bits.
 
*   severity    - Bit mapped word which determines the severity of
*               - the error messages which are output.  Most of
*               - error messages occurr when site and source names
*               - cannot be matched. If the bit is set an error
*               - message will be printed when:
*               - Bit   Meaning
*               -   1   Station name in Class 1 command (station
*               -       command) can not be matched with list
*               -       of sites in the experiment.
*               -   2   Source name  in Class 2 command (source
*               -       command) can not be matched with list
*               -       of sources for experiment.
*               -   3   Station name in apriori
*               -       file does not match names in experiment.
*               -       (Since the apriori files contain many
*               -       sites, setting this bit will in general
*               -       cause lots of error messages.)
*               -   4   Source name in apriori file does not match
*               -       in experiment (see Bit 3 warning)
*               -   5   OBSELETE at the moment
*               -   6   Station name in eccentricity file cannot
*               -       be matched to stations in this experiment.
*               - > 7   NOT USED
 
*   sit_parn(3,max_sites)   - Parameter numbers for each component
*               -(global XYZ) of the site position for each site
*   sou_parn(2,max_sources) - Parameter numbers for RA and Dec of
*               - each source.
*   sit_used(wrd_site)  - Bit set if a site is actually used in this
*               - solution.
*   sou_used(wrd_source)- Bit set is a source is actually used in this
*               - solution
 
*   tid_parn(3,max_sites)   - parameter numbers for h,l,and Lag
*               - at each site or for all the sites (in index 1)
*               - if global tides are to be estimated.
 
      integer*4 atm_parn(2,max_sites), atm_az_parn(2,max_sites),
     .    atm_map_parn(max_sites), 
     .    bclk_parn(max_bclk), bclk_bas_num(2,max_bclk), num_bclk,
     .    axo_parn(max_sites), crt_outopts, default_base_cont,
     .    default_data_mask, default_ion_mask, default_site_cont,
     .    etd_ext_parn(12,max_sites), gam_parn, nut_parn(2),
     .    plx_parn(max_sources), pmu_parn(2,3), prt_outopts, severity,
     .    sit_parn(3,max_sites), sou_parn(2,max_sources),
     .    sit_used(wrd_site), sou_used(wrd_source),
     .    tid_parn(3,max_sites)

* New pointers needed for new partials and command
*   eor_ut1_parn(4)  - Parameter numbers for the 4 UT1 diurnal ans
*       semidirunal signals
*   eor_xy_parn(6)   - Paramter numbers for the 6 PM terms, diurnal
*       retrograde, and semidiurnal Prograde and retrograde.

      integer*4 eor_ut1_parn(4), eor_xy_parn(6)
 
*   default_medium_cont   - Default setting for medium contributions
 
      integer*4 default_medium_cont
 
*   atm_az_mar(2,max_sites) - random walk statistics for the
*               - NS and EW components of the asymetric part
*               - of the atmospheric delay.
*               - Units: ps**2/sec
*   atm_az_apr(2,max_sites) - apriori variances for the azimuthal
*               - components of the atmospheric delay.
*               - Units: ps**2
*   atm_map_mar(max_sites)  - Random walk statistic for mapping 
*                 function Units: C**2/sec
*   atm_map_apr(max_sites)  - Apriori varianace for mapping 
*                 function Units: C**2
*   
*   bclk_apr(max_bclk) - apriori variances for the baseline 
*               - dependent clocks.
*   axo_mar(max_sites)  - random walk statistic for the axis
*               - offset parameter.
*               - Unit: m**2/s
*   axo_apr(max_sites)  - apriori variance for axis offset
*               - parameter.
*               - Unit: m**2
*   etd_ext_mar(12,max_sites)    - the random walk statistics for
*               - the 8 extended earth tide partials.
*   etd_ext_apr(12,max_sites)    - the apriori variances for the
*               - eight extended earth tide parameters.
*   eor_ut1_apr(4)  - apriori sigms for UT1 diurnal and semidiurnal
*               - signals (mas)
*   eor_ut1_mar(4)  - Markov sigmas for UT1 diurnal and semi (mas**2/sec)
*   eor_xy_apr(6)   - apriorio sigmas for diurnal and semidiurna;
*               - polar motion (mas)
*   eor_xy_mar(6)   - markov statistics for diurnal and semidiurnal
*               - polar motion (mas**2/sec)
*   gam_mar     - random walk statistics for gamma.
*               - Unit: /sec
*   gam_apr     - apriori variance for gamma.
*               - Unit: dimensionless
*   nut_mar(2)  - random walk statistics for the
*               - nutation in longtiude and obliquity
*               - parameters.
*               - Units: mas**2/sec
*   nut_apr(2)  - apriori variance for the nutation angles.
*               - Units: mas**2 for both
*   plate_motion(3,max_sites)   - Rotation vector for each site
*               - due to plate motion.  These values are given
*               - as the three Euler angle rates (rad/year).
*   plx_apr(max_sources)    - the aprioris varinaces for the
*               - paraxalax partials
*   plx_mar(max_sources)    - random walk statistics for the
*               - parralax parameters.
*   pmu_mar(2,3)  - random walk and IRW statistics for the
*               - X,Y pole position and UT1-AT earth orientation
*               - parameters.  NOTE: Value and rate and then by x,y,UT1
*               - Units: mas**2/sec, mas**2/(1000sec)*2/sec
*   pmu_apr(2,3)  - Apiori variances for the pole position and ut1.
*               - Units: mas**2 for all three, (mas/1000 sec)**2
 
*   sit_mar(3,max_sites)    - random walk statistics for the
*               - XYZ positions of site.
*               - Units: m**2/sec for all three.
*   sit_apr(3,max_sites)    - Apriori variance of the components
*               - of the site positions.
*               - Units: m**2
*   site_change(3,max_sites)    - changes to the XYZ cordinates
*               - of each site due to change in apriori positions,
*               - plate motions, and the eccentricity corrections.
*               - Units: m
*   site_pos_epoch(max_sites)   - the epochs of the site positions
*               - used in this solution. (Julian date)
*   sou_mar(2,max_sources)  - random walk statistics for the
*               - RA and Dec of each radio source.
*               - Units: mas**2/sec
*   sou_apr(2,max_sources)  - aproiri variance of the components
*               - of the source position.
*               - Units: mas**2
*   source_change(2,max_sources)    - changes to the RA and DEC
*               - of each source due to changes in the aprioris.
*               - Units: both mas.
*   tid_mar(3,max_sites)    - random walk statistics for the
*               - earth tide h,l and lag angle parameters at
*               - each site.  (If Global tides then index 1
*               - contains the values)
*               - Units: /sec, /sec, and deg**2/sec
*   tid_apr(3,max_sites)    - Apriori variance for the tidal
*               - parameters.
*               - Units: dimensionless for h and l, deg**2 for
*               -        lag)
 
      real*4 atm_az_mar(2,max_sites), atm_az_apr(2,max_sites),
     .    atm_map_mar(max_sites), atm_map_apr(max_sites),
     .    bclk_apr(max_bclk),
     .    axo_mar(max_sites), axo_apr(max_sites),
     .    etd_ext_mar(12,max_sites), etd_ext_apr(12,max_sites), 
     .    eor_ut1_apr(4), eor_ut1_mar(4), eor_xy_apr(6),
     .    eor_xy_mar(6),  gam_mar,
     .    gam_apr, nut_mar(2), nut_apr(2), plate_motion(3,max_sites),
     .    plx_apr(max_sources), plx_mar(max_sources), pmu_mar(2,3),
     .    pmu_apr(2,3), sit_mar(3,max_sites), sit_apr(3,max_sites),
     .    site_change(3,max_sites), site_pos_epoch(max_sites),
     .    sou_mar(2,max_sources), sou_apr(2,max_sources),
     .    source_change(2,max_sources), tid_mar(3,max_sites),
     .    tid_apr(3,max_sites)
 
*   batch_mode  - indicates solution is running in batch mode
*               - (set by passing markov control file through
*               - the runstring)
*   ecc_applied - Indicates that the eccentricity corrections have
*               - already been applied.  (This is used in multiple
*               - solutions)
*   eof_mar     - Indicates that the end of the Markov file
*               - has been found or that -2 was encountered.
*   gen_par     - Indicates that KAL_FILE will need to be
*               - generated (needed for multiple solutions)
*   glb_soln    - Indicate global solution file is to be saved.
*               - Set by giving GLB_FILE name (unless name of
*               - file NONE)
*   global_tides    - Indicates that global rather than site
*               - dependent tides should be estimated.
*   redo_clocks - Indicates that the apriroi clock model in the
*               - DATA_FILE should not be used.  The new clock
*               - parameters will be saved in the data file.
*   use_apr     - Indicates that site and source positions should
*               - updated to values in APR_FILE.  Set by giving
*               - file name.
*   use_ecc     - Indicate that the eccentricities in ECC_FILE
*               - should be used.  Set by giving the ECC_FILE
*               - name.
*   use_plate_motion    - Indicates that the plate motion model
*               - should be used.  Set by giving the plate motion
*               - model parameters.
 
 
      logical batch_mode, ecc_applied, eof_mar, gen_par, glb_soln,
     .    global_tides, redo_clocks, use_apr, use_ecc,
     .    use_plate_motion
 
*   last_markov_word    - Last word in the markov block
*   dummy_markov(127)   - 127 word padding to ensure that we do
*               - not overwrite anything when file is read.
 
* MOD TAH 920604 - Removed 3 words to allow for pmu_batch interval
*                  values.
      integer*4 last_markov_word, dummy_markov(124)
 
      equivalence ( markov, atm_parn )

* MOD TAH 920604

*   pmu_batch  - Time between fixed estimates of the polar motion
*                UT1 estimates (minutes).  If this value is set
*                to zero then standard RW process will be used.
*                These values are input in the pmu_mar command.

      real*4 pmu_batch(3)
 
*-----------------------------------------------------------------------
*     Start common declaration
 
      common / solvk_markov / atm_parn, atm_az_parn, atm_map_parn,
     .    bclk_parn, bclk_bas_num, num_bclk, axo_parn,
     .    crt_outopts, default_base_cont, default_data_mask,
     .    default_ion_mask, default_site_cont, etd_ext_parn, 
     .    eor_ut1_parn, eor_xy_parn, gam_parn,
     .    nut_parn, plx_parn, pmu_parn, prt_outopts, severity,
     .    sit_parn, sou_parn, sit_used, sou_used, tid_parn,
     .    default_medium_cont, atm_az_mar, atm_az_apr, 
     .    atm_map_mar, atm_map_apr, 
     .    bclk_apr, axo_mar, axo_apr, etd_ext_mar, etd_ext_apr, 
     .    eor_ut1_apr, eor_ut1_mar, eor_xy_apr,
     .    eor_xy_mar,  gam_mar, gam_apr, nut_mar, 
     .    nut_apr, plate_motion, plx_apr, plx_mar, pmu_mar, pmu_apr,
     .    sit_mar, sit_apr, site_change, site_pos_epoch, sou_mar,
     .    sou_apr, source_change, tid_mar, tid_apr, batch_mode,
     .    ecc_applied, eof_mar, gen_par, glb_soln, global_tides,
     .    redo_clocks, use_apr, use_ecc, use_plate_motion,
     .    last_markov_word, pmu_batch, dummy_markov
 
