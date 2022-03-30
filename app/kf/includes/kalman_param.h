
c The parameter file for SOLVK -- kalman filter software
c
c mod history                 4:16 PM  THU.,  3  JAN., 1985
c.....................................................................
c who  when   what
c TAH 850103 :Input this file
c JLD 860610 :Changed parameters to add SOUC_DEL option
c JLD 870413 :Added MAX_BASELINES
c JLD 870415 :Changed MAX_KAL_RECS from ? to 500
c TAH 870602 :Changed Mal_kal_recs from 500 to 400
c TAH 870720 :Added max_glb_parn, max_glb_sites, max_glb_sources
c TAH 870729 :Added max_glb_trans, max_glb_mar, max_nut_coeff,
c             max_etd_coeff
c TAH 880418 :Added structure_ext, extent and directory for source
c             structure files.
c TAH 910411 :Added versions numbers for globk, htoglb
c TAH 910905 :Added curr_ocean_series and ocean_file to allow
c             computation of ocean loading contributions.
c TAH 920616 :Added version for htoh
c TAH 920817 :Added parameters required by the Kalman filter
c             GPS analysis system (ctogobs,solvg).
c TAH 940205 :Changed max_glb_svs and max_glb_sources to 
c             larger values
c TAH 940608 :Added features for earthquake processing in an
c             automatic fashion; added variance factors for 
c             individual experiments; added true dynamic memory
c             mapping. VERSION 3.2 introduced.
c TAH 940923 :Changed the number of stations in the KalObs files
c             from 8 to 16.  New program convert_kalobs written
c             to convert old kalobs file to the new format.  
c             Introduced max_pts for the maximum number of points
c             to save the polar motion/ut1 fut1_pts and fwob_pts
c             arrays.
c TAH 950124 :Changed max_glb_parm to 4094 and max_glb_sites to 
c             1024.
c TAH 960813 :Increased the max_glb_sites to 2048.
c TAH 960826 :Increased dimensioning of transition rows and colummns
c             and max number of parameters (Vers 4.10 of globk)
c TAH 960906 :Added new command to globk (Ver 4.11; see glb_cntl.h)
c TAH 981112 :Added new feature to rename command. Ver 4.15
c TAH 981020 :Added multiple polar motion/UT1 parameters and allowed
c             time tag block in binary file.
c             Added new satellite radiation parameters and antenna
c             offsets.
c TAH 991110 :Introduced non-secular components to station position
c             evolution and site downweighting.
c TAH 010916 :Increased number of sites to 4096. Increased comfile
c             version
c TAH 061204: Increased to max_gchannels to 16 channels (Also can be set with 
c             max_chan command in autcln)      
c MAF 210629: Increased max_gsvs to 45 to accommodate BeiDou constellation
c-----------------------------------------------------------------------
c
 
      integer*4 max_sites, wrd_site, max_sources, wrd_source, max_bl,
     .    max_obep, max_eprec,max_nut,  max_parn, max_tran, max_clk_brk,
     .    max_cl,   max_av,   max_save, max_headr_size, max_baselines,
     .    kal_cart, sol_cart, site_start, site_end,
     .    souc_start, souc_end, mis_start, mis_end,
     .    krec_len, max_kal_recs, max_deriv,
     .    scr_wrds, scr_forsl, scr_outsl, scr_baksl, max_bclk,
     .    max_pts

c   max_pts         - Maximum number of points in the fut1_inf and
c                     fwob_inf arrays in KalObs files.
 
c
*   max_base_conts  - The number of different baseline dependent
*                   - contributions which can be applied.
*   max_channels    - Maximum number of Mark III channels
*   max_chi_types   - Number of types of parameters in the global
*                   - solution which can have their postfit
*                   - rms computed.
*   max_commands    - Number of commands in the SOLVK program
*   max_clk_order   - Maximum clock order allowed for apriori
*                   - clock model (starts at 0)
*   max_edit_types  - Maximum number of edit types allowed to
*                   - be recorded. Used for DELETE_COUNT.
*   max_etd_coeff   - Maxiumim number extended earth coefficients
*                   - which can be estimated.
*   max_ephem_epoch    - Maxiumum number of epochs for ephemeris
*                   - positions and rates.
*   max_frq_wvr     - Maximum number of channels in the WVR
 
*   max_glb_commands- Number of commands for GLOBK.
 
*   max_glb_mar     - Maximum number of markov parameters allowed
*                   - in a global solution.
*   max_glb_parn    - Maximum number of global parameters which
*                   - can be estimated
*   max_glb_sites   - Maxium number of global sites.  This should be
*                   - large enough to allow for all VLBI sites
*   max_glb_site_wrds   - Maxiumum number of words need for bit mapped
*                   - sites
*   max_etd_sites   - Maximum number of sites allowed in a global solution
*                     if site dependent extended Earth tide coefficients are
*                     to be estimated.
*   max_glb_sol     - Maximum number of experiments which can be
*                   - combined in one solution.
*   max_glb_sources - Maximim number of global sources.  Again should
*                   - be large enough to handle all sources.
*                   - (Even if all these positions can not be estimated
*                   - we need the space for the names)
*   max_glb_sou_wrds    - Maximum number of words for source bit mapping
*   max_glb_svs     - Maximum number of global SVS
*   max_svs_elem    - Maximum number of orbital elements for a satellite
*                     (95/7/29: Value is 17).
*   max_glbapr_sz   - Maxiumum number of real*8 values to be allocated
*                   - the apriori values in Global software
*   max_lcodes      - Number of lcodes to be read from the data base.
*                   - THIS VALUE MUST ACCUARTELY REFLECT THE ACTUAL
*                   - NUMBER OF LCODES USED BECAUSE IT IS USED IN
*                   - ALOCATING SPACE IN THE READIN_USER COMMON.
 
*   max_medium_conts- Number of site dependent propagation medium
*                   - contributions which can be applied
*   max_nut_coeff   - Maxium number of nutation coefficeints which
*                   - which can be estimated in a global solution
*   max_coeff_periods - Maximum number of periods allowed for
*                     coefficiennts
*   max_etd_names   - Max. names fort extended Earth parameters.
*   max_ut1_names   - Max  names for Ut1 diur/semi coefficients
*   max_xy_names    - MAx  names for xy diur/semi coefficients
*   max_ut1_coeff   - Max UT1 diur/semi coefficients
*   max_xy_coeff    - Max xy diur/semi coefficients
 
*   max_site_conts  - Number of site dependent contributions which can
*                   - be applied.
*   max_time_edits  - Maximum number of time ranges which can be
*                   - specified for data not being used.
 
*   max_trans_col   - Maximum number of parameters which can be
*                   - be involved in state transissions
*   max_trans_row   - Maxiumum number of rows in global transission
*                   - matrix.  (Number of parameters which can be
*                   - affected by transissions)
 
*   version_year    - Date at which current version of KalObs files
*   version_month   - were introduced.  This parameter must be updated
*   version_day     - when ever the KalObs files are changed.
 
*   readin_scr_size - Number of scratch common words in readin.
*                   - Must be greater than 4 + 6*max_sources
 
*   substitute      - Substitute value if a quanity is not in
*                   - the database.
 
*   sort_recl       - Record length in bytes for the sort file used
*                     by the global solution software.  Should not
*                     be less than file name length
 
*   curr_ocean_series - Number of the current ocean loading series.
 
 
      integer*4 max_base_conts, max_channels, max_chi_types,
     .    max_commands, max_clk_order, max_edit_types,
     .    max_ephem_epoch, max_frq_wvr, max_glb_commands, max_glb_mar,
     .    max_glb_parn, max_glb_sites, max_glb_site_wrds, max_glb_sol,
     .    max_glb_sources, max_glb_sou_wrds, max_glb_svs, max_svs_elem,
     .    max_glbapr_sz,max_lcodes, max_etd_sites,
     .    max_medium_conts, max_nut_coeff, max_coeff_periods,
     .    max_etd_names, max_ut1_names, max_xy_names,
     .    max_etd_coeff, max_ut1_coeff, max_xy_coeff, max_site_conts,
     .    max_time_edits, max_trans_col, max_trans_row, version_year,
     .    version_month, version_day, version_hour, readin_scr_size,
     .    substitute, sort_recl, curr_ocean_series

*   max_eq    - Maximum number of earthquakes allowed in one solution
*   max_eqfiles - Maximum number of earthquke files allowed 
*                (Added 5.18 090930).
*   max_rn    - Maxumum number of site renames and displacements 
*               (optional allowed).
*   max_eq_wrd - Maxiumum number of 32 bit words for eq entries to mark
*                as used.
*   max_rn_wrd - Maxiumum number of 32 bit words for rename ntries to mark
*                as used.
*   max_eq_cmd - NUmber of earthquake commands

*   max_apr_files - Number of apriori files that can be specified.
*                   Each is read in the sequence given.

      integer*4 max_eq, max_eqfiles, max_rn, max_eq_cmd, max_apr_files,
     .          max_eq_wrd, max_rn_wrd 
 
*   alias_file_default  - default name for alias_file (READIN)
*   default_com_file    - Name of the default name for the SOLVK common
*                       - file.
*   char_sub        - Character subsitute value if quantity is
*                   - not in the database.
*   list_mask_default   - Default name for the list file mask to
*                   - be used in GLRED
*   simobs_help     - Name of the help file for SIMOBS program
*   structure_ext   - Extent and directory for source structure files
*                   - (Set to .mod::sources)
*   ocean_file      - Ocean loading data file.
 
 
      character*(*) alias_file_default, char_sub, list_mask_default,
     .    simobs_help, structure_ext, default_com_file, ocean_file
 
*    globk_version  - Version number for globk (stored as string)
*    globc_version  - Version number for globc (stored as string)
*    glbak_version  - Version number for glbak (stored as string)
*    glist_version  - Version number for glist
*    glorg_version  - Version number for glorg
*    solvk_version  - Version number for solvk.
*    htoh_version   - Version number for htoh
*    solvg_version  - Version numner of solvg
 
 
      character*(*) globk_version, globc_version, 
     .    glist_version, glbak_version, solvk_version,
     .    glorg_version, htoh_version,  solvg_version
 
*   max_vma_space   - Maximum amount of VMA to be used
 
 
      integer*4 max_ema_space, max_vma_space
 
 
      real*8 max_sigma
 
c
*                                 !  maximum number of sites allowed in
c                                   a single experiment
c                                    MOD TAH 940923 Changed to 16
      parameter ( max_sites = 16 )
      parameter ( max_pts   = 10 ) 
c
      parameter (max_baselines = max_sites * (max_sites - 1) / 2)
c                             ! The maximum number of baselines possible
c
*                                                 ! number of 16 bit words
      parameter ( wrd_site = (max_sites-1)/16+1 )
c                                   needed for one bit per site
c
*                                     ! maximum number of sources allowed in
      parameter ( max_sources = 128 )
c                                   a single experiment
c
*                                                     ! number of 16 bit words
      parameter ( wrd_source = (max_sources-1)/16+1 )
c                                   needed for one bit per source
c
*                                                       ! Maximum number of
c                                   baseines which can be processed
      parameter ( max_bl = (max_sites-1)*max_sites/2 )
 
      parameter ( max_bclk = max_bl )
c
c
*                                        ! maximum number of observations
      parameter ( max_obep = 2*max_bl )
c                                   per epoch in an experiment.  If delay and
c                                   rate data are be to used, this number
c                                   will be 2*max_bl
c
*                                         ! maximum number of records from
      parameter ( max_eprec = max_obep+2)
c                                   kalfile file plus 2 which may be needed
c                                   for a single epoch
c
*                                  ! maximum number of nutation periods
      parameter ( max_nut   = 10 )
c                                   which can be estimated (4 for each
c                                   period).
c
*                                  ! maximum number of parameters which
      parameter ( max_parn  = 512)
c                                    can be estimated in one arc of the
c                                    of the Kalman filter.
c
*                                          ! maximum number of state
      parameter ( max_tran = 2*max_sites )
c                                    transission terms. Assummed one for
c                                    for each site clock and atmosphere.
c
*                                       ! the maximum group delay sigma
      parameter ( max_sigma = 1000.d0 )
c                                    that an observation can have for it
c                                    to be included in the solution (ps)
c
*                                    ! maximum number of clock breaks which
      parameter ( max_clk_brk = 10 )
c                                    can be modeled in the solution
c
*                                    ! maxmimum number of station residuals
      parameter ( max_av = 10 )
c                                    which can be averaged while looking for
c                                    clock breaks
*                                       ! the size of the arras used to save
      parameter ( max_save = 2*max_av )
c                                    the station_means.
c
      parameter ( max_headr_size = 4*(max_sites+max_sources)
     .                           + max_parn )
c                                    ! the maximum length of the header
c                                    records in BCK_FILE
c
*                                  ! cartridge for the KALFIL to be written
      parameter ( kal_cart  = 53 )
c                                    to.
c
*                                  ! Cartridge for the SOLVE common files
      parameter ( sol_cart  = 53 )
c
*                                                 ! bounds for the control
      parameter ( site_start = 1 , site_end = 26)
c                                    words which are applicable to stations.
c
*                                                   ! bounds for the control
      parameter ( souc_start = 27 , souc_end = 30 )
c                                    words which are applicable to sources.
c
*                                                 ! bounds for the control
      parameter ( mis_start = 31 , mis_end = 68 )
c                                    words which are miscellaneous.
c
*                                    ! Number of commands in the SOLVK
      parameter ( max_commands = 68)
*                                    ! program.
c
c
*                                    ! maximum number of parameters on
      parameter ( max_deriv = 53 )
c                                    which a single observation can
c                                    depend
*                                               ! record length
*     parameter ( krec_len = max_deriv*5 + 25 )
      parameter ( krec_len = max_deriv*3 + 15 )
c                                    for the KALFIL
c                                    which contains the O-C and derivatives
c                                    to be used by the solution programs
*     NOTE:  Condition below removed in the UNIX version which is using
*     all integer*4 variables.
c     NOTE : krec_len must be an even number **** Hence if max_deriv
c     is changed the constant in krec_len must be change back to 21 (the
c     value actually required).  This procedure is necessary becuase VIS
c     is used to move the records of KALFIL back and forth between main
c     memory and EMA.
c
*                                       ! the maximin number of kalfil
      parameter ( max_kal_recs =  2000 )
c                                    records which can be saved in ema
c                                    before the kalfil is written to
c                                    disk.
c
*                                    ! the size of the 'scratch' common
      parameter ( scr_wrds = 17000 )
c                                    area to be used in the program.  This
c                                    area is used for large dcb buffers
c                                    brious subroutines.  *** Must ***
c                                    be greater than 1296 so that parfile
c                                    can be read. (Changed from 2576 TAH
c                                    870108 for FmpRunProgram to fit)
c
*                                           ! the size of the scratch comm
      parameter ( scr_forsl =16*max_parn )
c                                    to be used in program FORSL.  Used
c                                    in the same manner as scr_wrds but
c                                    do not have the restriction of being
c                                    greater than 1296, but must be greater
c                                    than max_parm*4, because it is used
c                                    as a scratch area.
*                                         ! the size of the scratch common area
      parameter ( scr_outsl =16*max_parn)
c                                    in the program OUTSL.  Should not be
c                                    any smaller than 'scr_forsl'
c
*                                          ! the size of the scratch common
      parameter ( scr_baksl =16*max_parn )
c                                    in the program BAKSL.  Should not be
c                                    any smaller than the 'scr_forsl'
c
*                                         ! maximum number of words available
      parameter ( max_ema_space = 3276800)
c                                    in the smallest partition in which
c                                    the kalman filter will run
c
*                                         ! Maximim number of words of VMA
      parameter ( max_vma_space = 15000000 )
*                                     ! space allowed in BAKSL.  This space
*                                     ! is only used if MAX_EMA_SPACE is
*                                     ! not large enough to run solution
 
*                                      ! Number of baseline independent
      parameter ( max_base_conts = 7 )
*                                    ! contributions.
 
*                                       ! 14 channels used (for saving
      parameter ( max_channels   = 14 )
*                                       ! phase cal and data phases)
 
*                                       ! Number of parameter types which
      parameter ( max_chi_types  = 24 )
*                                       ! can have rms calculated during
*                                       ! global solution.  Add an extra
*                                       ! three to allow for NEU accumulation
 
*                                       ! Clocks can be quadratic
      parameter ( max_clk_order  =  2 )
 
      parameter ( max_cl = (max_clk_order+1)*max_sites)
*                                        ! maximum number of clock parameters
*                                  ! which can be estimated for calculation
*                                  ! of residuals.
 
*                                       ! Upto 16 types of edits.
      parameter ( max_edit_types = 16 )
 
*                                       ! Ephemeris saved at three values
      parameter ( max_ephem_epoch=  3 )
 
*                                       ! Allow upto three channel WVR
      parameter ( max_frq_wvr    =  3 )
 
*                                       ! 80 commands are available
      parameter ( max_glb_commands = 100)
 
*                                       ! 40 global parameters can be markov
* MOD TAH 060516: Increased dimesions for globk
* max_glb_mar    = 4096 
* max_glb_parn   = 8192
* max_glb_sites  = 4096
* max_glb_sol    = 5000
* max_trans_col  = 4096
* max_trans_row  = 4096
* 
* MOD TAH 140809: Increased max_glb_mar from 8192 to 16384.
* MOD TAH 190510: Increased max_glb_mar from 16384 to 40000. 
      parameter ( max_glb_mar    = 40000 )
 
*                                       ! 256 parameters in global solution
* MOD TAH 140809: Increased max_glb_parn 16384 to 32768.
* MOD TAH 190510: Increased max_glb_parn 32768 to 36000
      parameter ( max_glb_parn   = 40000 )

* MOD TAH 190521: Increased number of sites from 8192 to 16384 (needs 14 bits
*                 to represent (more than 16 bits will cause problems).
      parameter ( max_glb_sites  =  16384 )
 
      parameter ( max_glb_site_wrds = (max_glb_sites-1)/32+1 )
 
      parameter ( max_etd_sites  = 32 )
 
*                                       ! Experiments in one solution
      parameter ( max_glb_sol   = 10000)
 
*                                       ! transission elements
* MOD TAH 140809: Increased max_trans_col/row  8192 to 16384.
* MOD TAH 190510: Increased max_trans_col/row 16384 to 18000 (half max_glb_parn)
      parameter ( max_trans_col  = 20000 )
*                                       ! transissions affecting up to 60
      parameter ( max_trans_row  = 20000 )
*                                       ! parameters
*                                       ! Max 512 VLBI sources
      parameter ( max_glb_sources= 512)
 
      parameter ( max_glb_sou_wrds = (max_glb_sources-1)/32+1)
 
      parameter ( max_glb_svs    = 128 ) ! MOD TAH 091002: Allowed upto 128
                                        ! satellites to handle glonass+galielo
*                                       ! Value on 950729.
*                                       ! Updated to 23 981020:
* MOD TAH 190509: In new model ECOMC max 22 could be used at one time.
*                 6 IC, 3 Constant, 10 cos/sin on periodic terms, 3 antenna offsets.
      parameter ( max_svs_elem   = 23 )   

*                                              ! add 64 to allow for
*                                              ! end of record. 
      parameter ( max_glbapr_sz = 6*max_glb_parn+64 ) 
*     parameter ( max_glbapr_sz = 20*max_glb_sites +
*    .                            10*max_glb_sources )
 
*                                       ! Number of LCODES to be read from
      parameter ( max_lcodes     = 106 )
*                                       ! data base.  MUST GIVE ACTUAL NUMBER
*                                       ! BECUASE IT IS USED TO ALLOCATE
*                                       ! SPACE.
 
*                                       ! Number of propagation medium
      parameter ( max_medium_conts=28 )
*                                       ! contributions
 
*                                       ! Maximum number of nutation series
      parameter ( max_nut_coeff  = 30 )
*                                       ! coefficients (Each a pair in and
*                                       ! out of phase).
 
      parameter ( max_coeff_periods = 15 )
      parameter ( max_etd_names  = 12 )
      parameter ( max_ut1_names  =  4 )
      parameter ( max_xy_names   =  6 )
      parameter ( max_ut1_coeff  = max_coeff_periods*max_ut1_names)
      parameter ( max_xy_coeff   = max_coeff_periods*max_xy_names)
      parameter ( max_etd_coeff  = max_coeff_periods*max_etd_names)
*                                       ! Number of site dependent contrib-
      parameter ( max_site_conts =  6 )
*                                       ! utions.
*                                       ! Maxium number of time ranges that
      parameter ( max_time_edits =  5 )
*                                       ! can specified for editing
 
 
*                                       ! Current version date of
*                                       ! the KalObs files.  MUST
*                                       ! BE UPDATED IF FILE STRUCTURE
*                                       ! IS CHANGED.
* MOD TAH 940923: Changed from 90 2 1 to 94 9 23
      parameter ( Version_year  = 1994,
     .            version_month =    9,
     .            version_day   =   23,
     .            version_hour  =    0 )
 
*                                           ! scratch common area.  Must
      parameter ( readin_scr_size  = 1000)
*                                       ! be at least 4+6*max_sources long
 
*                                       ! indicates quanity not in database
      parameter ( substitute  = 000   )

* MOD TAH 151203: Increased file name lengrh from 128 to 256 characters 
      parameter ( sort_recl   = 256   )
 
      parameter ( alias_file_default = 'kalman_alias.dat' )
      parameter ( default_com_file   = 'solvk_com.bin' )
 
*                                       ! indicates character quanity not
      parameter ( char_sub = 'XXXXXXXX')
*                                       ! in database.
 
*                                          ! Default name for the list mask
      parameter ( list_mask_default = 'L@')
*                                       ! in GLRED (file will be purged after
*                                       ! run)
 
      parameter ( simobs_help = 'simobs.help' )
 
      parameter ( structure_ext = '.MOD')
 
* Version number parameters.  Should be updated when ever the stdrelease
* directory is changed.
 
* Mods made from 3.0 to 3.1:  911017:
* -----------------------------------
* SOLVK: Introduced new mapping functions (MTT), seasonal model available
*        for when no met data; ocean loading constributions computed when
*        ever series number is different bewteen KalObs file and
*        curr_ocean_series; mapping function partials fully supported now;
*        Azimuthal assymmetry partial modified.
* GLOBK: Fixed GLSAVE to handle GPS data; introduced file for time-dependent
*        markov statitics on satellite orbits; implemented markov NEU on sites;
*        modded glorg so that parameters can forced to be equal adjustment.
* GLOBK: Incremented to 3.2 for earthquake processing.
* GLOBK: Increased to 4.12 11/27/96: added new commands 
* GLOBK: Increased to 4.13 07/06/97: Added extension to max_chi command for large
*        rotations, increased internal checking of solution compatability. 
* GLOBK: Increased to 4.14 09/09/97: Added new crt/prt/org options to erase files
*        and to stop outputs.
* GLORG: Increased to 4.01 11/27/96: Added new features
* GLORG: Increased to 4.02 05/14/97: Modified site selection in pos_/rate_org
* GLOBK: Increased to 4.15 11/12/97: Added including name substr for global files
*        to used in renames.
* GLOBK: Increased to 4.16 02/16/98: Fixed long standing bug in XYZ_to_GEOD; dNdLat
*                                    derivative in error by factor 2.
* GLOBK: Increased to 4.17 05/14/98: Added rn_used, eq_used to keep track of which
*        renames and earthquakes are actually used.
* GLORG: Introduced assign_number to allow sites to assigned to plates but not
*        used in their estimation when plate poles are computed. Version 4.05 980611.
* GLOBK: Mulitiple pmu estimates, new radiation parameters and sv antenna offset
*        981020 (mods started). Version 5.0
* GLORG: Version increased to make compatible with globk version 
* GLOBK: Version 5.02 (990226) Fixed bugs in multiday PMU and IRW models.
* GLOBK: Version 5.03 (991110) Non-secular site position evolution and site
*        downweighted.  Fixed leap-seconds code, if leap.sec file not found.
* GLORG: Version 5.01 (990226) Replaced use_site with stab_sit command to aviod
*        confusion in commands
* GLORG: Version 5.02 (99110) Propagation of non-secular position evolution.
* GLORG: Version 5.03 (000302) Added constraint sigma to force command.
* GLOBK: Version 5.04 (000404) Cleaned up interaction between globk and glorg
*        when out_glb and out_sol are used.  Now two output binary files
*        are written with the constrained version with a .CON appended to the
*        name.
* GLORG: Version 5.04 (000404) Cleaned up code associated with writing 
*        constrained sol-files in glorg.
* GLOBK: Version 5.05 (000903) Added the sig_neu command which allows components
*        of sites to be downweighted by adding noise to cov_obs
*        Added the _XCL name extent which will stop a site being used.
*        Made changes to eq_name_change but this should not effect anything
*        (needed for new glist features)
* GLORG: Version 5.05 (000903) Added _XHI name extent which will automatically
*        stop the heights being equated for such sites.
* GLIST: Version 5.01 (009093) Added checking and reporting of apriori coordinate
*        agreement 
* GLOBK: Version 5.06: (010525) Fixed problems in the EXTENDED model implemenation 
* GLORG: Version 5.06: (010525) Fixed problems in the EXTENDED model implemenation 
c GLOBK: Version 5.07: (030116) include list parameter before run (PLST)
c GLORG: Version 5.07: (030116):Report on match of parameters from renamed sites (RNRP)
c                               Force arprioris to match on equates (FIXA)
c                               Extended stab_site command to allow time range
c                               and hfile string..
c GLORG: Version 5.08: (030513): Added translation estimation when plate poles
c                               are estimated and a command to turn this feature
c                               off.
c GLOBK: Version 5.08: (030607): Add earthquake log estimates, wild card use_site
c                               decimation, glsave options.
c GlOBK: Version 5.09: (040220): Fixed problem in wild card mar_neu, added SMAR
c                                fixed treatment of satellite antenna offsets.
c                                Changed 'all' return in get_command from 999 to 999999
c                                so there in no overlap with site numbers.
c GLORG; Version 5.10: (040329): Added wild cards and ALL to plate and assign_p
c                                commands, added arguments to stab_min, added
c                                output to frame definition parameters.
c GLOBK: Version 5.10: (040703): Added atmospheric delay estimates (carried foward from
c                                solve vers 10.11) 
c GLOBK: Version 5.11: (040801): Changed memory allocation to be 64-bit compatable
c GLORG: Version 5.11: (040801): Changed memory allocation to be 64-bit compatable
* GLOBK: Version 5.12: (060518): Increased dimensions of numbers of stations and parameters
* GLOBK: Version 5.13: (070823): Introduces MIDP option for outut so that solutions referr
*                                to midpoint of data. Also applied to GLSAVE and
*                                introduced  glsave -M option when run stand alone.
* GLOBK: Version 5.14  (080103): Updated to pass nutation and gravity field information
* GLOBK: Verions 5.15  (080129): Updated to allow orbits to be rotated when different
*                                nutation models have been used in the gamit runs.
* GLOBK: Version 5.16  (080416): Added feature to allow all file names to have wild cards. 
* GLORG: Version 5.16  (080416): Added feature to allow all file names to have wild cards. 
* GLOBK: Version 5.17  (080724): Added reading +REFERENCE_FRAME lines from apr files and
*                                outputing these.
* GLOBK: Version 5.18  (090930): Added BREAK to eq-files, allowed upto 16 eqfiles to be used
*                                Added RENM to output options to output complete rename file.
* GLOBK: Version 5.19  (101015): Allowed full site antex file name and increased receiver
*                                name length to C*20.  curr_glb_com_ver=5.11  
* GLOBK: Version 5.20  (110512): Added IERS2010 mean pole tide; Also changed htoglb to map
*                                etide bits to globk version and initialize sinex as 3+23.
* GLOBK: Version 5.21  (110815): Added multiple options with + in <OPTIONS> string;
* GLORG: Version 5.17  (110815): Added reference frame hiarchy in stab_site command with / 
*                                separator in names
* GLOBK: Version 5.22  (130419): Added APPLOAD_MOD command which allows loading corrections
*                                to be removed during globk run (corrections must be in
*                                binary hfiles).
* GLOBK: Version 5.23  (130716): Added white noise mar_scale term so process noise is indepentent
*                                of time step (allows mixed scale solution combintations).
* COM_FILE: Version 5.13 (131106): Increased max_rn to 65536 (new com_file version and new globk)
* GLOBK: Version 5.24  (131106): MAX_RN increased. 
* GLOBK: Version 5.25  (140327): Added passed of earth radiation and antenna thrust model.
* GLOBK: Version 5.26  (140406): Added atmospheric delay model records and ocean pole tide
* GLBAK: Version 5.02  (140327): Same addition: Needed for correct g-files.
* GLIST: Version 5.03  (140406): Added output of atmoshere and model parameters
* GLOBK: Version 5.27  (140809): Increased # of parameters and markov limits, added apriori sigmas
* GLORG: Version 5.18  (150130): Allow short site names in equate/unequate/force commands so that
*                                all sites that match will be used.
* GLOBK: Version 5.28  (150526): Added I*8 value for total memory needed (allows  >20000 parameters
*                                in input and output solutions).  Fixed VEL MEAN bug. 
* GLORG: Version 5.19  (150526): Same mods as GLOBK 5.28
* GLOBK: Version 5.29  (150825): Modified satellite names to include satellite vechicle number as 
*                                well as PRN i.e., PRN_01 will now be PRN_0149 or PRN_0163 depending
*                                on when data is collected.  Upto 100 satellite are possible in a
*                                single run (max expected currently is 72).
* GLOBK: Version 5.30  (161128): Incread max_rn from 65536 to 262144, max_ss from 4096 to 8192. 
*                                GLORG and curr_glb_com_ver also increases. 
* GLORG: Version 5.20  (161128): Same mods as GLOBK 5.30.
* GLOBK: Version 5.31  (180402): Updated to support GNSS with new satallite naming scheme.  For GPS
*                                the default is to retain the PRN_<pn><sv> form where pn is 2-digit
*                                PRN, <sv> is 2-digit satellite vehicle number.  The new scheme is
*                                <sys><sv3>_<pn> where sys is G--GPS, R-Glonass, E-Galileo, C-Beidou
*                                and the sv number has 3-digits.  The "use_prnn N" command will invoke
*                                using the new scheme for all satellites. 
* GLOBK: Version 5.32  (190510): Updated maximum number of parameters to handle 6000 stations.  
*                                curr_glb_com_ver updated to 5.16 (globk_cntl.h).  Implemented
*                                ECOMC model.
* GLBAK: Version 5.03  (190526): Added frame realization to back solution (allowing loose solutions)
* GLOBK: Version 5.33  (200222): Mean pole model for pole tide, antenna azimuth meta data.
* GLOBK: Version 5.34  (200602): Added LAB<str> option to GLOBK OPTION features.
* GLORG: Version 5.21  (210517): Added BALN option to org_opt to modify station coordinate and 
*                                satellite PCO covariances to account for the number of times a 
*                                station is used. PCO change based on number of ACs.  The operation 
*                                is only done in the output stage of glorg.   If the out_sol command 
*                                is used in glorg to outout the combined binary h-file (.CON extent), 
*                                the balancing operation has been performed.  (Standard out_glb binary 
*                                h-files from globk is not effected).  When using the .CON h-files, 
*                                apr_rot should be used to allow rotation of the frame resolved H-file.  
*                                These .CON files also have the translations resolved.  If converted 
*                                to SINEX, the -d=TR can be used to add back rotation and translation 
*                                uncertainty.  SCALE should not be estimated when frame resolved H-files 
*                                are produced because the scale will be fixed at the reference frame 
*                                values in these files.

      parameter ( globk_version  = '5.34' )
      parameter ( globc_version  = '5.34' )
      parameter ( glbak_version  = '5.03' )
      parameter ( glist_version  = '5.04' )
      parameter ( glorg_version  = '5.21' )
      parameter ( solvk_version  = '4.0' )
      parameter ( htoh_version   = '4.0' )
      parameter ( solvg_version  = '1.0' )
 
* curr_ocean_series - Current ocean loading series to be used (checed
*                     with used_ocean_series to see if needs updating)
*          MOD TAH 911105: Updated series number so that the Q1 argument
*                     would be corrected (ANGFAC(3,8) changed from 0 to 1)
* ocean_file        - Name of the ocean loading data file (same structure
*                     as BLOKQ.)
 
      parameter ( curr_ocean_series = 2 )
      parameter ( ocean_file = '/data3/tah/tables/ocean_910905.dat' )
 
* SOLVG/GPS data processing parameters
 
*   max_gsites  - Maximum number of sites in SOLVG (should be
*               - => max_grcv
*   max_grcv    - Maximum number of recievers in SOLVG.  This
*               - value can differ from the max_gsites because
*               - of kinematic receivers which generate multiple
*               - sites.
*   max_gorb_sites  - Maximum number of orbiting receivers allowed
*   max_gkd_sites  - Maximum number of kinematic receiveres allowed
 
*   max_gsvs    - Maximum number of satellite transmitter in
*               - SOLVG.
*   max_gepoch  - Maximum number of epoch in SOLVG file
*   max_gobs_proc   - Maximum number of data types which can
*               - can processed sumulataneously
*   max_gedt_types  - Maximum number of edit types for SOLVG
*   max_down_sites_svs  - Maximum number of entries allowed for
*               - specifying which epochs at which sites and svs
*               - should not be used
*   max_gchannels   - Maximum number of channels any reveiver
*               - is allowed to have.
*   max_gdata_types - Maximum number of data types allowed in
*               - a cfile.
*   curr_sys_version    - Current version of Gobs file system.
*               - (value *100 i.e. version 1.00 is 100)

*   go_recl     - Record length of the Gobs file.  This value is
*                 set to be a little greater or equal to the one-
*                 way data record (keeps file at smallest possible
*                 size. (With version 1.00 files this value is 64
*                 words.  Actual length is 62 words, I*4)

*   eph_buff_len - Length of the buffer in I*4 words of the record
*                  to save rcv depenedent data from the ephead block
*                  (must be a multiple of go_recl)
*   max_ephem    - Maximum number of days that can be used in orbital
*                  arc. (used when 'Make' option is given on svs_file)
 
      integer*4 max_gsites, max_grcv, max_gorb_sites, max_gkd_sites,
     .    max_gsvs, max_gepoch, max_gobs_proc, max_gedt_types,
     .    max_down_sites_svs, max_gchannels, max_ephem,
     .    max_gdata_types, curr_sys_version, go_recl, eph_buff_len
 
      parameter ( max_gsites = 200 )
      parameter ( max_grcv   = 100 )
      parameter ( max_gorb_sites = 2 )
      parameter ( max_gkd_sites =  6 )
* MOD TAH 200210: Increased from 32 to 35
* MOD TAH 210629: Increased from 35 to 45
      parameter ( max_gsvs = 45 )
      parameter ( max_gepoch = 36000 )
      parameter ( max_gobs_proc =  4 )
      parameter ( max_gedt_types = 32 )
      parameter ( max_down_sites_svs = 50 )
* MOD TAH 011229: Increased to 14 channels (see max_chan in autcln)
* MOD TAH 061204: Increased to 16 channels (Also can be set with max_chan
*                command)      
      parameter ( max_gchannels = 16)
      parameter ( max_ephem  = 100 )
      parameter ( max_gdata_types = 4 )
      parameter ( curr_sys_version = 100 )

      parameter ( go_recl   = 64 )
      parameter ( eph_buff_len = go_recl )
 
* Earthquake and renaming
* MOD TAH 080502: Increased number earthquakes to 256
* MOD TAH 090930: Increaded to 1024
      parameter ( max_eq      =  1024 )
* MOD TAH 090930: Added max_eqfiles
      parameter ( max_eqfiles  =   16 )

* MOD TAH 010526: Increased max renames from 1024 to 8192 to handle
*     tsview outputs
* MOD TAH 080503: Increased again from 8192 to 32768.
* MOD TAH 131106: Increased from 32768 to 65536
* MOD TAH 161128:  Increased from 65536 to 262144
      parameter ( max_rn     = 262144 ) 
      parameter ( max_eq_wrd = (max_eq-1)/32 + 1 )
      parameter ( max_rn_wrd = (max_rn-1)/32 + 1 )
     
      parameter ( max_eq_cmd =   10 )
 
      parameter ( max_apr_files = 10 )

*------------------------------------------------------------
* NEW values SINEX definitions and globk binary file extension

* glbf_version  - Current version of binary global files
*     NOTE: Saved versions of this are in sglb_vers and cglb_vers
*           (check grep on these).

      integer*4 glbf_version

* default_ institute    - Name to be for sinex file creator and owner.
*     SHOULD BE CHANGED WHEN GLOBK INSTALLED

      character*4 default_institute
      
* MOD TAH 981020: Version 102 associated with parameter epoch block.
* MOD TAH 990901: Version 103 associated with extended aprriori model
*     and inclusion of sig_neu command.
* MOD TAH 010526: Version 104: Increased max number of renames to 8192.
* MOD TAH 050622: Version 105: Introduction of satellite information
*                 records.
* MOD TAH 101015: Version 106: Added sti_antmod (C*16) and ssrecv_ty_end
*                 (C*4) to sinf_def to allow full name of antmod file and
*                 receivers to be saved.
* MOD TAH 190606: Version 107: Changed the satellite radiation parameter
*                 mapping to support ECOMC model (see svel_to_code.f and
*                 glb_hdr_def.h)
* NOTE: When version changed, glsave/gw_descript.f should be checked to
*     see that it will recognise the new version.  Looks like mod'd
*     to take any version above 0 (TAH 050622).
      parameter ( glbf_version = 107 )
      parameter ( default_institute   = 'UNK ' )

*------------------------------------------------------------

* NEW VALUES for mulitple PMU estimates 981020

* max_mul_pmu  - Maximum multiple pmu values

      integer*4 max_mul_pmu
      parameter ( max_mul_pmu = 366 )

* max_ptide_hfs -- Maximium number of names for hfiles strings
*     that can be specified in the app_ptid command 
      integer*4 max_ptide_hfs
      parameter ( max_ptide_hfs = 32 )

*-------------------------------------------------------------
* NEW VALUES FOR non-secular position evolution and site down weighting
*
* max_nonsec  -- Maximum number of non-secular times allowed
*     (one non-secular term per entry.  Entry constists of generically
*     site, type, date (JD), parameter (time constant or period),
*     2 coefficients each for NE and U.
* max_ss  -- Maximum number of site sigmas allowed.

      integer*4 max_nonsec, max_ss
* MOD TAH 060518: Increased to 4096 from 2056
* MOD TAH 140225: Increased max_nonsec to 8192 from 4096
      parameter ( max_nonsec = 8192 )
* MOD TAH 161128: Increased max_ss from 4096 to 8192
      parameter ( max_ss     = 8192 )     



 
 
