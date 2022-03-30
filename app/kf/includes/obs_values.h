 
*--------------------------------------------------------------------
*                                                     OBS_VALUES.FTNI
*     The values part of the header records of the new Kalman filter
*     files.
*     These records contain the information needed by both the readin
*     software and the Kalman filter software.
*
*     The common block itself is mostly laidout in the same order as the
*     declaration.  The exceptions to this are the variables which give
*     the size and start block numbers of the various sections.  These
*     are positioned (and must remain) at the beginning of the first
*     block of the file.  The data base, version and experiment title
*     are also put at the font for convenience.
*
* MOD TAH 870730: Added atm_max_temp(max_sites) to the end of the values
*     block.  The spaced used is extracted from the dummy array area
*     so that the file size will not change.  Also defined bits 6 & 7 of
*     data_notes.
*
* MOD TAH 880202: Added usual_num_channels to give the number of channels
*     in this data.  This value is set in KALUP and in READIN for newer
*     KalObs files.  If FRNGE_QUAL is OK and yet one or more channels
*     is missing, then the FRNGE_QUAL is set to Z.
 
* MOD TAH 900131: Added new atmospheric quantities which can be entered
*     by the user at run time.  These value also have corresponding
*     time dependent values in the data records.  OBS_APR also changed
*     to add extended earth tide aprioris and UT1 and PM diurnal and
*     semi-diurnal corrections. Partials for these also added to data
*     record.

* MOD TAH 901118: Added check on the status of the tidal contribution
*     added to UT1. (Used BIT 11 of data_notes, and added ut1_tide to
*     obs_data)
*
*     T.Herring                   03:02 PM MON., 12 Jan., 1987
*-------------------------------------------------------------------
 
*   ko_dcb(16)      - the KALOBS file dcb buffer (since read
*                   - as type 1 file only the 16 words need
*                   - to be given)
*   values          - the first word in the values block.
*                   - Everything is read through this word.
*                   - Equivalanced to START_NAMES.
 
 
      integer*4 ko_dcb(16), values
 
*     Start the declarations with the variables in the file which we
*     wish to have appear at the start of the file.  In particular the
*     start block numbers, number of blocks, and data base name and
*     description.
 
*   start_names     - Record number for the start of the names
*                   - block
*   start_apr       - Record number for the start of the aprioris
*                   - block
*   start_data      - record number of the start of the data
*                   - records
*   num_values_blocks   - Number of blocks in the VALUES section
*                   - of the file
*   num_names_blocks    - Number of blocks in the NAMES section
*                   - of the file
*   num_apr_blocks  - Number of blocks in the APRIORIS section
*                   - of the file
*   obs_rec_len     - Length of the data records (words)
 
*   num_obs         - Total number of observations in this exp.
*   num_sites       - Number of sites in this experiment
*   num_sources     - Number of sources in this experiment
 
*   version         - Version of the data base being processed
 
*   dum_filv        - filler to get real*8 on correct boundary
 
 
      integer*4 start_names, start_apr, start_data, num_values_blocks,
     .    num_names_blocks, num_apr_blocks, obs_rec_len, num_obs,
     .    num_sites, num_sources, version
 
 
      integer*4 dum_filv
 
*   read_epoch      - Julian epoch when this data was readin
 
      real*8 read_epoch
 
*   values_read     - Indicates that start_@ and num_@_blocks
*                   - are known.
 
 
      logical values_read
 
*   data_base       - Name of the data base from which the data
*                   - was read
 
      character*10 data_base
 
*   expt_title      - concatination of the data base description
*                   - and the names of the sites involved.
*                   - OR input by user. (Length set so that KalObs
*                   - file name and title fit on one line)
 
 
      character*62 expt_title
 
*     End variables to be at the top of the common
*------------------------------------------------------------------------
 
*   avail_baseline  - Indicates that baseline calibration was
*                   - available.  Set during data readin stage
*                   - or when a program updates the contents
*                   - of the data files. See cont_baseline for
*                   - bit masking
*   avail_met(max_sites)    - Indicates which met values are
*                   - available at each site.  The mapping is
*                   - Bit   Meaning
*                   -   1   Pressure
*                   -   2   Temperature
*                   -   3   Relativity humidity
*                   -   4   WVR antenna temperatures
*                   -   5   Data Base WVR delays
*                   -   6   User dry zenith delay (in data records)
*                   -   7   User wet zenith delay (in data records)
 
*   avail_site(max_sites)   - Indicates that the contribution
*                   - was available by site (see above comment)
*                   - see cont_site for bit masking
*   avail_structure(wrd_source) - indicates that structure correction
*                   - is available.
*   av_obs_dur      - Average observation duration to be used if the
*                   - actual duration was not in the data base.
*   cont_baseline   - Contriubtions by baseline (one value for all
*                   - baselines to aviod closure errors)
*                   - Meaning of the bits is:
*                   - Bit   Meaning
*                   -   1   earth tide
*                   -   2   pole tide
*                   -   3   ocean loading
*                   -   4   axis offset contribution
*                   -   5   relativity contribution
*                   -   6   Group delay ion.  See also ION_MASK for
*                   -       bit mapping of editing ion correction
*                   -   7   Phase delay ion.  See also ION_MASK for
*                   -       bit mapping of editing ion correction
*                   -  >8   NOT USED
 
*   cont_site(max_sites)   - saves the bit masked contributions
*                   - for each site
*                   - to be used in this solution.  The meaning of
*                   - the bits is to use each of the following if
*                   - they are available:
*                   - Bit   Meaning
*                   -   1   cable calibration
*                   -   2   feed rotation correction (phase only)
*                   -   3   phase calibration removal
*                   -   4   User defined
*                   -   5   Radial Ocean load avaliable
*                   -   6   Horizontal Ocean load avaiable.
*                   -  >6   NOT USED
 
*   cont_structure(wrd_source) - bit masked for each source
*                   - for which source structure contribution is
*                   - to be used
*   clk_brk_site(max_clk_brk)   - site number for each of the
*                   - clock breaks
 
*   data_mask       - A bit pattern which when and'd with the
*                   - DATA_FLAG generates the KUNW flag.  If the
*                   - result of the and is no zero then the data
*                   - is not used.  Data_Mask = -1 will delete
*                   - all observations with any problem.  Data_mask
*                   - 0 is not allowed, since the 16th bit of DATA_FLAG
*                   - is used to indicate that there is a calibrations
*                   - problem with this observation.
*   data_notes      - Indicates some conditions about the data.  The
*                   - mapping is:
*                   - Bit   Meaning
*                   -   1   Start and stop seconds not found.  Default
*                   -       duration for observations used.
*                   -   2   Data should have phase cal. removed.
*                   -   3   FINAL        values for wobble and UT1 used
*                   -   4   PRELIMINARY  values for wobble and UT1 used
*                   -   5   EXTRAPOLATED values for wobble and UT1 used
*                   -   6   Extended earth tide partials have been
*                   -       added to the KalObs file.
*                   -   7   Maximum temperatures at each site have been
*                   -       added to KalObs file.
*                   -   8   Indicates that the phase_cal_cont has been
*                   -       re-computed in KALUPD (Same sign convention
*                   -       used as GSFC calculation in DBCAL)
*                   -   9   Indicates that the most frequently used
*                   -       number of channels has been set.
*                   -  10   Indicates that clock jumps have been added to
*                           this data and that these should be added to
*                           theoretical model.
*                   -  11   Set if we know the status of UT1 tables (ie if
*                           tides have been added or not)
 
*   data_type       - Type of data to be used in this solution.
*                   - This word is bit mapped.  With the following
*                   - meanings:
*                   - Bit   Meaning
*                   -   1   Use Group-delay data.
*                   -   2   Use Phase-delay data.
*                   -   3   Use SB-delay data.
*                   -   4   Use Phase-delay-rate data.
*                   - There are currently restrictions which allow
*                   - only one delay type to be used at one time.
*   delete_count(max_edit_types)    - Counts the primary reasons
*                   - for not using data in a solution.  The
*                   - ordering determines the primary reasons.
*                   - The meaning of each entry is:
*                   - Word  Meaning
*                   -   1   Database lists as bad FRNGE quality
*                   -   2   ION delay is bad (fails ION_MASK)
*                   -   3   WVR temperatures bad
*                   -   4   WVR Dbase delay bad
*                   -   5   Propagation medium problem
*                   -   6   Cable Cal bad
*                   -   7   SOLVK interactive edit
*                   -   8   SOLVK automatic edit
*                   -   9   Data base grp/rate interactive edit (unw=1)
*                   -  10   phase delay connection problem
*                   -  11   User defined editing flag
*                   -  12   User source delete
*                   -  13   Down-site delete due to time interval
*                   -  14   Below elevation cutoff
*                   -  15   Data do not close within tolerance.
*                   -  16   RESERVED for PLTSL program use.
*                   - These values are reset before each solution
 
*   down_sites(2,max_time_edits)    - site numbers for time range
*                   - editting.  If one of the values is <1 then all
*                   - data to the other site is downweighted
*   down_num        - Number of values in down_sites.
 
*   down_source(wrd_source) - Bit mapped words which indicate which
*                   - sources should not be used in this solution
 
*   ion_mask        - Sets the edit conditions for the ionospheric
*                   - delay.   This value and'd with the ION_BITS
*                   - for each observation should yield zero, if the
*                   - data point in not to be unweighted
*   num_clk_brk     - Number of clock breaks in this experiment
 
*   run_time(7)     - Date and time of the last solution using this
*                   - data file.
*   res_num(max_obep)   - Number of residuals for each baseline
*                   - with delay and rate numbers saved separately.
*   wvr_code_limit  - Limit to be placed on the WVR code quality to
*                   - be used in the solution.  (Value must be
*                   - greater than or equal 0 and less than or
*                   - equal to this value.
 
 
      integer*4 avail_baseline, avail_met(max_sites),
     .    avail_site(max_sites), avail_structure(wrd_source),
     .    av_obs_dur, cont_baseline, cont_site(max_sites),
     .    cont_structure(wrd_source), clk_brk_site(max_clk_brk),
     .    data_mask, data_notes, data_type,
     .    delete_count(max_edit_types), down_sites(2,max_time_edits),
     .    down_num, down_source(wrd_source), ion_mask, num_clk_brk,
     .    run_time(7), res_num(max_obep), wvr_code_limit
 
*   cont_medium(max_sites) - the contributions to be for
*                   - calibrating the propagation medium delays
*                   - The meaning of the bits is:
*                   - Bit   Meaning
*                   - Dry Zenith definition
*                   -   1   Saastamonian dry zenith delay
*                   -   2   Constant dry supplied by user
*                   -   3   User defined dry term (computed from
*                   -       the available whether data)
*                   -   4   User defined dry term in the data records
*                   -   5   SPARE
*                   -   6   Constant dry term based on height
*                   - Wet zenith definition
*                   -   7   Saastamoinen wet zenith delay
*                   -   8   wet delay from WVR anntenna temperatures
*                   -   9   wet delay from database WVR values
*                   -  10   User defined wet delay.
*                   -  11   user defined wet term in the data records.
*                   -  12   SPARE
*                   -  13   constant wet delay supplied by user
*                   - Dry mapping function
*                   -  14   Marini dry mapping function
*                   -  15   CfA-2.2 dry mapping function
*                   -  16   MIT lapse rate/ht trop mapping function
*                   -  17   MIT Temperature mapping function 
*                   -  18   SPARE
*                   -  19   Chao dry mapping function
*                   - Wet mapping function
*                   -  20   Marini wet mapping function
*                   -  21   CfA-2.2 wet mapping function
*                   -  22   cosecant wet mapping function
*                   -  23   user defined wet mapping function
*                   -  24   MIT Wet mapping function (constant)
*                   -  25   MIT Temperature wet mapping function.
*                   -  26   Chao wet mapping function
*                   - Modifiers (Applied to all modern mapping functions)
*                   -  27   Use Max temp at site rather than obs temp
*                   -  28   Use seasonal temperature rather than observed
*                           (This bit gets set if their is no weather
*                            data available.)
 
 
      integer*4 cont_medium(max_sites)
 
*   barometer_hgt(max_sites)    - The height difference between the
*                   - barometer and the intersection of axes of the
*                   - radio telescope. (m).  The sense is height of
*                   - telescope minus height of barometer.
*   calc_ver        - Version of Calc used for theorecticals
*   clk_brk_mag(max_clk_brk)    - the magnitude of the clock
*                   - breaks detected (if these values are set
*                   - in a post-processing operation very large
*                   - values are inserted) (ps)
 
*   correlator_noise(2) - Correlator noise of delay and rates
*                   - (ps and fs/s)
*   down_times(2,max_time_edits)    - Epochs of the downweighted
*                   - time boundaries measured from START_EPOCH.
*                   - Values in fractional days.
 
*   elev_cutoff(max_sites)  - the elevation cutoff angle at each
*                   - site (rads)
*   ellip_hgt(max_sites) - Ellipsoid height of each of the
*                   - sites.  Included here because the
*                   - atmospheric calibration routines will need
*                   - these values (m)
*   latitudes(max_sites)    - Latitudes of the sites (rad).  These
*                   - values are used for atmospheric delay
*                   - calculations
*   longitudes(max_sites)   - The East longitude of the sites.
 
*   n_sigma_limit   - Limit on deviation of postfit resiudals
*                   - for data not to edited by automatic editor in
*                   - KALAN.
 
*   pre_chisq       - Prefit residual chi**sq per number of
*                   - observation (not degrees of freedom)
 
*   res_chisq(max_obep) - the Chi**2/n for each baseline and delay
*                   - and rate (if rates are used in the solution)
 
*   res_wrms(max_obep)  - Weighted wrms scatter of the residuals
*                   - for each baseline and delay and rate. (ps and
*                   - fs/s).  If a back solution is not run, then
*                   - all values are zero.
*   tol_close       - Tolerance on the closure of delays and rate
*                   - around triplets.  The value if given as closure
*                   - error divided by sigma of closure quanity.
*   user_dry_con(max_sites) - Users values for constant dry zenith
*                   - delay (ps)
*   user_wet_con(max_sites) - Users values for constant wet zenith
*                   - delay (ps)
 
 
      real*4 barometer_hgt(max_sites), calc_ver,
     .    clk_brk_mag(max_clk_brk), correlator_noise(2),
     .    down_times(2,max_time_edits), elev_cutoff(max_sites),
     .    ellip_hgt(max_sites), latitudes(max_sites),
     .    longitudes(max_sites), n_sigma_limit, pre_chisq,
     .    res_chisq(max_obep), res_wrms(max_obep), tol_close,
     .    user_dry_con(max_sites), user_wet_con(max_sites)
 
* New atmosphere constants
 
*   lapse_con(max_sites)    - Lapse rate at each site (K/km)
*   ht_con(max_sites)       - Height of tropopoause at each site (km)
*   tbias_con(max_sites)    - Bias in surface temperature (K)
*   adry_con(4,max_sites)   - coefficients of the dry mapping function
*                           - by site
*   bdry_con(4,max_sites)   - Coefficinets of bending due to dry air
*   awet_con(4,max_sites)   - Coeffiencts of wet mapping function
*   bwet_con(4,max_sites)   - Coefficients of bending due to wet air.
 
 
      real*4 lapse_con(max_sites), ht_con(max_sites),
     .    tbias_con(max_sites), adry_con(4,max_sites),
     .    bdry_con(4,max_sites), awet_con(4,max_sites),
     .    bwet_con(4,max_sites)
 
*   clk_brk_epoch(max_clk_brk)  - Julian epochs for the clock
*                   - breaks in this data.
*   clock_epoch     - Julian epoch for the clock polynomials
*   clock_poly(0:max_clk_order,max_sites) - Polynomial coeffs
*                   - for removing the "gross" clock behavior
 
*   end_epoch       - Julian epoch (date+fraction) of last obs
*                   - data set.
*   mid_epoch       - Julian epoch of the middle of the data set
*                   - (No observation is necessarily at this time)
*   start_epoch     - Julian epoch of first observation in data set
*   hp_read_epoch   - Epoch for reading HP files (copied at conversion) 
 
      real*8 clk_brk_epoch(max_clk_brk), clock_epoch,
     .    clock_poly(0:max_clk_order,max_sites), end_epoch, mid_epoch,
     .    start_epoch, hp_read_epoch
 
*   min_frnge_quality   - Minimum FRNGE quality code to be
*                   - used in the solution.
*   dummy_char_FRNGE    - dummy character to ensure that we dont end
*                   - up with an odd byte address.
 
      character*1 min_frnge_quality
 
 
      character*3 dummy_char_FRNGE
 
*   clk_known       - True if the clock polynomials for this
*                   - data are known
*   ecc_known       - True if the eccentricities for this data
*                   - are known.
 
 
      logical clk_known, ecc_known
 
*     START of added item.  Remove count of new values from spare_values.
*   atm_max_temp(max_sites) - Maximum temperature at each site for the
*                   - experiment.  This value may be better for use
*                   - in the CfA mapping function.
 
 
      real*4 atm_max_temp(max_sites)
 
*   usual_num_channels  - Usually number of channels used in this
*                       - data set
*   bclk_site           - Bit mapped word which gives the sites 
*                       - which require baseline dependent clocks.
 
 
      integer*4 usual_num_channels, bclk_site

*   Kal_ver             - Version of KalObs file.  This value starts at 
*                         data base version *100 and is increment by 1 each time
*                         the header is written.
*   ut1_tide            - If value is -1 then no short period tides in UT1,
*                         if zero then tables already include tides.  Once
*                         update pmu is run, then this is zet to -1.  Data_notes
*                         bit 11 is used to make sure a valid value is in this
*                         slot.
*   used_ocean_series   - Ocean loading series number (used to check
*                         against parameter curr_ocean_series in 
*                         Kalman_param.h to see if we need to update.)

      integer*4 Kal_ver, ut1_tide, used_ocean_series
 
*   spare_values(125)   - One spare block, set up at the initial vers
*                   - of READIN.  If more variables are required at
*                   - some later time use this space.
*                   - TAH 870730: Removed 16 words for atm_max_temp.
*                   - TAH 880802: Removed  1 word  for
*                   -             usual_num_channels
*                   - TAH 900701: Reset to 128 for UNIX versions of files
*                   - TAH 900808: Used one word for sites which require
*                                 baseline dependent clocks. BCLK_SITE
*                   - TAH 901011: Used one word for KAL_VER
*                   - TAH 901118: USed one word for UT1_TIDE
*                   - TAH 910905: Used one word for USED_OCEAN_SERIES
* TAH 940923: Reset with new version of KalObs file 
 
      integer*4 spare_values(128)
 
*   last_values_word    - Last word used in the values block.
*                   - This word is proceded by 128 dummies for
*                   - when the file is read.
*   dummy_values(127)   - dummy values at the end of the block
*                   - to ensure that when the block in read that
*                   - no values for later blocks get overwritten
 
 
      integer*4 last_values_word, dummy_values(127)
 
      equivalence ( values, start_names )
 
*--------------------------------------------------------------------
* COMMON DECALRATION
 
      common / values_block / ko_dcb, start_names, start_apr,
     .    start_data, num_values_blocks, num_names_blocks,
     .    num_apr_blocks, obs_rec_len, num_obs, num_sites, num_sources,
     .    version, dum_filv,
     .    read_epoch, values_read, data_base, expt_title,
     .    avail_baseline, avail_met, avail_site, avail_structure,
     .    av_obs_dur, cont_baseline, cont_site, cont_structure,
     .    clk_brk_site, data_mask, data_notes, data_type, delete_count,
     .    down_sites, down_num, down_source, ion_mask, num_clk_brk,
     .    run_time, res_num, wvr_code_limit, cont_medium,
     .    barometer_hgt, calc_ver, clk_brk_mag, correlator_noise,
     .    down_times, elev_cutoff, ellip_hgt, latitudes, longitudes,
     .    n_sigma_limit, pre_chisq, res_chisq, res_wrms, tol_close,
     .    user_dry_con, user_wet_con,
     .    lapse_con, ht_con, tbias_con, adry_con,
     .    bdry_con, awet_con, bwet_con,
     .    clk_brk_epoch, clock_epoch,
     .    clock_poly, end_epoch, mid_epoch, start_epoch, hp_read_epoch,
     .    min_frnge_quality, dummy_char_FRNGE, clk_known, ecc_known,
     .    atm_max_temp, usual_num_channels, bclk_site, Kal_Ver,
     .    ut1_tide, used_ocean_series, 
     .    spare_values, last_values_word, dummy_values
 
*
