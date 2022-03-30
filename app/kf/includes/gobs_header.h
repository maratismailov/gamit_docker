 
*     Gobs Header definition file.  This file contains all of the
*     header information required by the Gobs files.
 
*     The Gobs Files have the following structure which should be
*     transparent to the user when the standard reading/writing/create
*     routines are used (Ghandlers).
 
*     | ---- Header block ------ || --epoch numbers--||-- Epoch 1 -- |
*     | -- Epoch 2 -- | ....... | -- Epoch N-- |
*
*     where
*      | denotes GO_RECL I*4 word boundaries.  The last_header_word
*         will appear within GO_RECL words of boundary.  There is a fixed
*         number of words in the Header Block for a given version of
*         Gobs File System (sys_version is current version).
 
*     The epoch numbers block points to the data numbers of the start
*         of each Epoch Block.  There are a variable number of these
*         records depending on the number of epochs in the specific
*         file.
 
*     The Epoch blocks are further divided into two types of record:
*     | -- Epoch Header-- || -- Obs 1 -- || -- Obs 2 -- |
*     | -- Obs 3 -- | .... | -- Obs M -- |
 
*     The Epoch header contains epoch specific information (similar to
*     the type 4 c-file record).  It is fixed length for a given
*     version of the Gobs File system.
 
*     Each of the Obs records contains one-way observation specific
*     information.  The records here are of a fixed length, but there
*     are variable numbers of Obs records at each epoch.  The internal
*     memory structure of the Obs record is longer than that written
*     to disk.  In order to keep the disk files as small as possible
*     the the memory allocation of the Obs record is longer than the
*     disk version.  In this additional memory, information computed
*     during the SOLVG run is stored here e.g., short_period ut1
*     partials, atmospheric partials, ocean loading corrections are
*     all stored in this extended area and computed at run-time
*     (depending on options specified.)
 
 
*---------------------------------------------------------------------
 
*     File control block.
 
*   go_dcb(16)  - This is the data control block used by
*               - direct access file reading routines.
*   header      - First word of the header block.  The
*               - disk image of the file will start hear.
 
 
      integer*4 go_dcb(16), header
 
*     The next group are the pointers to where various parts of
*     the gobs_file start and there size.
 
*   num_header_recs - Number of GO_RECL I*4 word records in the header
*               - records.
*   num_epoch_recs  - Number of GO_RECL I*4 word records in the epochs
*               - block
*   num_ephead_recs - Number of GO_RECL I*4 word records in each epoch
*               - header block. (This count includes the buffer
*                 entries of one per site).
*   num_eph_only_recs - Number of GO_RECL I*4 word records in just
*                 fixed part of the ephead block.
*   num_eph_buff_words - Number of I*4 words in each 
*                 buffer entry for each site. NOTE: This one is the
*                 number of words not the number of records.
*   num_data_recs    - Number of GO_RECL i*4 words in one-way
*               - observation record.  NOTE: This is the same as
*                 num_gobs (i.e., number of data points).
 
*   sys_version - Version of the Gobs_file definition (multiplied
*               - by 100.  Initial version is 1.00 i.e, 100.
*               - Kalman_param.h contains the current version in
*               - use).
*   Gobs_version    - Version of this file.  This starts with the
*               - ASCII code of the c-file series multiplied by
*               - 100 and the numbers less than 100 are the
*               - incremental version changes.
 
 
      integer*4 num_header_recs, num_epoch_recs, num_ephead_recs, 
     .    num_eph_only_recs, num_eph_buff_words,
     .    num_data_recs, sys_version, Gobs_version
 
*   num_gsites  - Number of sites in this Gobs File
*   num_grcv    - Number of recievers.
*   num_gkd_sites  - number of Kinematic/dynamic sites
*   num_gorb_sites  - number of orbiting sites
*   num_gsvs    - Number of satellites in this file
*   num_gepochs - Number of epochs in this file.
*   num_gobs    - Number of one-way data records in this file
*   num_gchans  - Maximum number of channels used by any
*               - reciever in this Gobs_file.
*   curr_epoch  - Number of the current epoch being read
*               - or written.  The ephead block is automatically
*               - when ever the observation number changes to a
*               - new epoch. Value is set to zero when file opened.
*   curr_obs    - Number of the current observation.  Value is set to
*               - zero when the file is opened.
 
      integer*4 num_gsites, num_grcv, num_gkd_sites, num_gorb_sites, 
     .    num_gsvs, num_gepochs, num_gobs, num_gchans, 
     .    curr_epoch, curr_obs
 
*   header_read - True if the header records of the this
*               - file have already been read (set false when
*               - file opened)
*   ephead_current - Used when files are written and is set false
*               - when ever a one-way data record is written.  When
*               - data outside of the current epoch is read, this
*               - logical will cause the updated ephead block to be
*               - written before the next ephead block is read.
 
 
      logical header_read, ephead_current
 
*   expt_title  - Title of this experiment.  Given by user
*               - when the Gobs_file is created with ctogobs.
 
 
      character*64 expt_title
 
* End of initial header list
*-----------------------------------------------------------------------
 
 
*   avail_site(max_gsites) - Indicates which types of calibration
*               - information are available by site.  It starts
*               - with meteorogical data and progress to
*               - calibrations
*               - Bit   Meaning
*               -   1   Pressure Available
*               -   2   Temperature Available
*               -   3   Relativity Humidity available
*               -   4   WVR calibration available
*               -   5   User defined dry delay (by epoch)
*               -   6   User defined wet delay (by epoch)
*               -   7   Multipath model available
*               -   8   Phase center variation model
*               -   9   Radial Ocean Loading
*               -  10   Horizontal Ocean Loading
*               -  11   Radial Atmospheric pressure loading
*               -  12   Horizontal Atmospheric pressure loading.
*   avail_svs(max_gsvs) - Indicates the types of calibration
*               - information available for the satellites.
*               - Bit   Meaning
*               -   1   Phase center variation (by epoch)
*               -   2   Stochastic orbit variation (by epoch)
 
*   cont_site(max_gsites)   - Bit set to indicate which
*               - to the theoretical model should be applied
*               - for this analysis.
*   cont_svs(max_gsvs)  - Indicates which contributions should
*               - be applied for each satellite.
 
*   edit_mask   - Mask to applied to determine the editing on
*               - the data when it is run.  (Same lay out as
*               - the data_flag variable read from the one-way
*               - data record.  (-1 will delete all data with any
*               - problem.)
*   data_proc(max_gsites)   - Bit mapped word which indicates
*               - the  types of data to be processed in the
*               - solution.
*               - NOTE: this selection is site dependent.
*               - Bit  Meaning
*               -   1  L1 carrier phase
*               -   2  L2 carrier phase
*               -   3  P1 P-code pseudo range L1
*               -   4  P2 P-code pseudo range L2
*               -   5  C1 C/A-code pseudo range L1
*               -   6  D1 L1 Doppler frequency
*               -   7  D2 L2 Doppler frequency
*               -   8  LC carrier phase
*               -   9  LG carrier phase
*               -  10  LC range
*               -  11  LG range
*               - For the range data types, the most accurate
*               - range type will be used (if there is any
*               - debate)
*   data_types(max_gsites) - Bit mapped word giving types of
*               - observations available by site.
*               - Currently coded types are:
*               - Type  Meaning
*               -   1  L1 carrier phase
*               -   2  L2 carrier phase
*               -   3  P1 P-code pseudo range L1
*               -   4  P2 P-code pseudo range L2
*               -   5  C1 C/A-code pseudo range L1
*               -   6  D1 L1 Doppler frequency
*               -   7  D2 L2 Doppler frequency
*   lamba_facts(max_gsvs,max_gdata_types, max_grcv) - Wave-
*               - length factors by satellite, data_type and
*               - site:  Factors are
*               -   0  No Data
*               -   1  unambiguous, undoubled
*               -  -1  ambiguous,   undoubled
*               -   2  unambiguous, doubled
*               -  -2  ambiguous,   doubled
 
*   file_notes  - These are notes on the file which indicate
*               - the status of the file.  Current meanings:
*               - BIT  Meaning
*               -   1  Kinematic data has been processed as
*               -      separate sites.
 
*   delete_count(max_gobs_proc, max_gedt_types) - Records the
*               - numbers of data edited in the solutions by the
*               - type of data.  The categories are (alligned
*               - the same as the data_flag
*               - Element Meaning
*               -
*   down_sites_svs(4,max_down_sites_svs)  - List of the sites and
*                 satellites to be edited
*               - by time.  The order here is
*                 (1) Site number (0 means all)
*                 (2) SVS number (0 means all)
*                 (3) Start epoch
*                 (4) End epoch
*   down_sites_svs_num  - Number of entries to use in the down_sites_svs
*               - array.
 
*   run_time(7) - Runtime of solution (ymdhms and s/100)
 
*   num_res(max_gobs_proc,max_gsites)   - Number of residuals
*               - for each type of data at each site.
*   num_epochs  - Number of epochs in this file (set based on all
*               - c-files merged)
*   obs_interval    - Obsersvation interval by site
*               - (secs)
*   rcv_interval(max_grcv)    - Receivers sampling interval by
*               - site (secs)
 
*   orb_to_rcv(max_gorb_sites)  - Gives the reciever number of the
*               - orbiting receivers i.e., if a receiver is
*               - in this list then it is orbiting.
*   kin_to_rcv(max_gkd_sites)  - Gives the reciever numbers of the
*               - kinematic reveivers.
*   site_to_rcv(max_gsites)    - Gives the site to reciever relation.
 
*   wvr_code_limit  - Limit on the WVR data quality code to
*               - used in the processing.
 
*   cont_medium(max_gsites) - Propagation medium contributions
*               - for calibrating the atmospheric delay
*               - Bit   Meaning
*               - Hydrostatic Zenith Delay
*               -   1   Saastamonian hydrostatic delay
*               -   2   Command specified constant
*               -   3   User defined term in data records
*               -   4   SPARE
*               -   5   SPARE
*               -   6   Constant based on height
*               - Wet Zenith delay
*               -   7   Saastamoinian wet delay
*               -   8   Wet delay from WVR
*               -   9   MTT Wet zenith delay
*               -  10   Seasonal model (climatology)
*               -  11   User defined term in data records
*               -  12   SPARE
*               -  13   SPARE
*               -  14   Constant specfied by user.
*               - Hydrostatic Mapping function
*               -  15   CfA-2.2 mapping function
*               -  16   MIT fully parametrerized mapping fn.
*               -  17   MTT mapping function
*               -  18   User defined function in data records
*               -  19   SPARE
*               -  20   SPARE
*               -  21   Chao dry mapping function
*               - Wet Mapping function
*               -  22   MTT mapping function
*               -  23   User defined function in data records
*               -  24   SPARE
*               -  25   SPARE
*               -  26   Chao Wet mapping
*               - Modifiers
*               -  27   Use maximum temperature at site
*               -  28   Use seasonal temperature
*               -  29   Use seasonal relative humidity
 
 
 
      integer*4 avail_site(max_gsites), avail_svs(max_gsvs),
     .    cont_site(max_gsites), cont_svs(max_gsvs), edit_mask,
     .    data_proc(max_gsites), data_types(max_gsites),
     .    lamba_facts(max_gsvs,max_gdata_types, max_grcv),
     .    file_notes, delete_count(max_gobs_proc, max_gedt_types),
     .    down_sites_svs(max_down_sites_svs), down_sites_svs_num, 
     .    run_time(7), num_res(max_gobs_proc,max_gsites),
     .    num_epochs, obs_interval(max_gsites),
     .    rcv_interval(max_grcv), orb_to_rcv(max_gorb_sites),
     .    kin_to_rcv(max_gkd_sites), site_to_rcv(max_gsites),
     .    wvr_code_limit, cont_medium(max_gsites)
 
*   base_noise(max_gobs_proc)   - Base level noise for each
*               - the data types to processed
*               - (cycles @ appropriate frequency)
*   elev_cutoff(max_gsites) - Elevation angle cutoff for each
*               - site (rads)
*   min_clean_elev  - Minimum elevation used in ctogobs for cleaning
*                 data.  No data below this elevation can be used 
*                 in solvg runs unless furrther cleaning is done.
*   ellip_hgt(max_gsites)   - Ellipsoid height of each site
*               - (used by the atmospheric routines) (m)
*   latitudes(max_gsites)   - Latitudes (used by atmospheric
*               - and ionosheric routines) (rad)
*   longitudes(max_gsites)  - Longtitudes (could be used
*               - by atmosphere and ionosphere routine (rad)
*   n_sigma_limit   - Editing conditions of deviations of
*               - residuals (post-processor input)
*   res_chisq(max_gobs_proc,max_gsites) - Chi**2/n for each
*               - data type at each site
*   res_wrms(max_gobs_proc,max_gsites)  - Weighted RMS scatter
*               - for each data type at each station.
 
 
      real*4 base_noise(max_gobs_proc), elev_cutoff(max_gsites),
     .    min_clean_elev, ellip_hgt(max_gsites), latitudes(max_gsites),
     .    longitudes(max_gsites), n_sigma_limit,
     .    res_chisq(max_gobs_proc,max_gsites),
     .    res_wrms(max_gobs_proc,max_gsites)
 
*     Meteorogical quantities
 
*   user_hyd_con(max_gsites)    - User defined hydrostatic
*               - zenith delay (m)
*   user_wet_con(max_gsites)    - User defined wet zenith delay (m)
*   lapse_con(max_gsites)   - Lapse rates at each site (K/km)
*   ht_con(max_gsites)      - Height of tropopause (km)
*   tbias_con(max_gsites)   - Bias in the surface temp. (K)
*   ahyd_con(3,max_gsites)  - Coefficients in the hydrostatic
*                           - mapping function (constant)
*   awet_con(3,max_gsites)  - Coefficients in the wet
*                           - mapping function (constant)
*   atm_max_temp(max_gsites)    - Maximum atmospheric temperature
*                           - at site during day.
*   atm_mar_sav(max_gsites) - Saved value of the atm. RW process
*                           - noise by site.
*   atm_ataz_sav(2, max_gsites) - Saved values of azimuthal
*                           - assymetry process noise (NS and EW
*                           - directions)
*   atm_map_sav(max_gsites) - Saved values of the mapping
*                           - function process noise.
*   clk_mar_sav(3, max_gsites)  - Saved value of clk. Process
*                           - noise.
*   sit_mar_sav(3, max_gsites)  - Saved values of the site
*                           - stochastic noise in NEU
*                           - coordinates (mainly for dynamic
*                           - receivers).
*   rsw_version(max_gsites) - Receiver software version.
 
 
      real*4 user_hyd_con(max_gsites), user_wet_con(max_gsites),
     .    lapse_con(max_gsites), ht_con(max_gsites),
     .    tbias_con(max_gsites), ahyd_con(3,max_gsites),
     .    awet_con(3,max_gsites), atm_max_temp(max_gsites),
     .    atm_mar_sav(max_gsites), atm_ataz_sav(2, max_gsites),
     .    atm_map_sav(max_gsites), clk_mar_sav(3, max_gsites),
     .    sit_mar_sav(3, max_gsites), rsw_version(max_gsites),
     .    dummy_head_R4
 
*   end_mjd      - Julian date at end of data (plus fraction of
*               - day)
*   mid_mjd      - Julian date at midpoint of data
*   start_mjd    - Julian date at start epoch of measurements.
*   pmu_mjd      - Julian date plus fractional day for polar
*               - position apriori's
*   ephem_mjd    - Julian date plus fractional day for ephemeris
*               - elements.
*   nom_rcv_clk(3,max_gsites)   - Nominal receiver clock offset,
*               - rate, and acceleration (sec, sec/sec,
*   nom_svs_clk(3,max_gsites)   - Nominal satellite clock
*               - offset, rate, and acceleration (sec, sec/sec,
*               - (sec/sec)/sec)
 
 
      real*8 end_mjd, mid_mjd, start_mjd, pmu_mjd, ephem_mjd,
     .    nom_rcv_clk(3,max_gsites), nom_svs_clk(3,max_gsites)
 
*     Apriori geodetic model information
*   ut1_apr(2)  - Value and rate of UT1-AT at pmu_epoch (mas and
*               - mas/day.
*   wob_apr(2,2)    - Pole position values and rate at pmu_epoch
*               - (mas and mas/day).  X and Y change with first
*               - index).
*   nut_ang_apr(2,2)    - Apriori nutation angles and rates
*               - Long, and Obliquity (mas and mas/day)
*   offset_L1(3,max_gsites) - L1 phase center position relative
*               - to ground mark (N,E,U) (m) (Note change in
*               - convension from GAMIT)
*   offset_L2(3,max_gsites) - L2 phase center position relative
*               - to ground mark (N,E,U) (m) (Note change in
*               - convension from GAMIT)
*   site_pos(3,max_gsites)  - Apriori positions of the sites
*               - in geocentric cartesian coordinates (m)
*   svs_elements(9,max_gsvs)    - Apriori orbital elements and
*               - solar radiation values (m, mm/sec, unitless)
*   orb_elements(9,max_gorb_sites)  - Apriori orbital elements and
*               - solar radiation values for orbiting
*               - receivers
 
 
      real*8 ut1_apr(2), wob_apr(2,2), nut_ang_apr(2,2),
     .    offset_L1(3,max_gsites), offset_L2(3,max_gsites),
     .    site_pos(3,max_gsites), svs_elements(9,max_gsvs),
     .    orb_elements(9,max_gorb_sites)
 
*   rcv_sw(max_grcv)  - Receiver software at each site.
*               - e.g. ROG, TRM, etc
*   rcv_motion(max_grcv)  - Defines the type of motion of
*               - the receiver:  Current choice are
*               - STAT  -- Static, fixed receiver
*               - DYNM  -- Dynamic moving receiver
*               - KINE  -- Kinematic
*               - ORBT  -- Receiver is orbiting
 
 
 
 
 
 
      character*4 rcv_sw(max_grcv), rcv_motion(max_grcv)
 
*   site_names(max_gsites)  - Short names of the sites.  These
*               - are the four character codes from the cfile
*               - name with _GPS appended.
*   svs_names(max_gsvs)     - Names of the satellites (PRN
*               - numbers in form PRN_02).
 
 
      character*8 site_names(max_gsites), svs_names(max_gsvs)
 
*   site_long_names(max_gsites)       - Long names of the sites.
*   xfile_name(max_grcv)  - Names of the X-files originally
*                           - used. (Raw Data files)
*   tfile_name              - Name of T-file (Ephemeris file)
*   jfile_name              - Name of J-file (Satellite clock
*                           - file)
 
      character*16 site_long_names(max_gsites), xfile_name(max_grcv),
     .    tfile_name, jfile_name
 
*     Ending declarations of record type:
 
*   spare_header(GO_RECL)   - GO_RECL words left spare for furture use
*               - while maintain backward compatatibility
 
 
      integer*4 spare_header(GO_RECL)
 
*   last_header_word    - Last word of the header.  This is used
*               - to get how large the file is.
 
      integer*4 last_header_word
 
*   dummy_header(GO_RECL)   - Padding a end of record to ensure
*               - that when we read the header we don't overwrite
*               - anything else.
 
 
      integer*4 dummy_header(GO_RECL)
 
*.................................................................
 
 
      common / gobs_head / go_dcb, header, num_header_recs,
     .    num_epoch_recs, num_ephead_recs, num_eph_only_recs,
     .    num_eph_buff_words, num_data_recs, sys_version,
     .    Gobs_version, num_gsites, num_grcv, num_gkd_sites, 
     .    num_gorb_sites, num_gsvs, num_gepochs, num_gobs, num_gchans,
     .    curr_epoch, curr_obs, header_read, ephead_current,
     .    expt_title, avail_site, avail_svs, cont_site, cont_svs,
     .    edit_mask, data_proc, data_types, lamba_facts, file_notes,
     .    delete_count, down_sites_svs, down_sites_svs_num,
     .    run_time, num_res, num_epochs, obs_interval, rcv_interval, 
     .    orb_to_rcv, kin_to_rcv, site_to_rcv, 
     .    wvr_code_limit, cont_medium, base_noise,
     .    elev_cutoff, min_clean_elev, ellip_hgt, latitudes, 
     .    longitudes, n_sigma_limit,
     .    res_chisq, res_wrms, user_hyd_con, user_wet_con, lapse_con,
     .    ht_con, tbias_con, ahyd_con, awet_con, atm_max_temp,
     .    atm_mar_sav, atm_ataz_sav, atm_map_sav, clk_mar_sav,
     .    sit_mar_sav, rsw_version, dummy_head_R4,
     .    end_mjd, mid_mjd, start_mjd, pmu_mjd,
     .    ephem_mjd, nom_rcv_clk, nom_svs_clk, ut1_apr, wob_apr,
     .    nut_ang_apr, offset_L1, offset_L2, site_pos, svs_elements,
     .    orb_elements, rcv_sw, rcv_motion, site_names, svs_names,
     .    site_long_names, xfile_name, tfile_name, jfile_name,
     .    spare_header, last_header_word, dummy_header
 
