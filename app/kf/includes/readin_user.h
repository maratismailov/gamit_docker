c
c@USER
c
c     User information
c
c MOD TAH 870130 Put all of the LCODE data and logicals into a common
c     so that at the end of each database we can list the LCODES
c     which were not found.
c
      integer*4 icrt, log_lu, ipass(5), today(7), cart
 
      integer*4 delay_flag
 
c
*   created     - Indicates output obs file created sucessfully
*   default_title   - Indicates defualt experiment title is being
*               - generated (ie., none passed through runstring)
*   kill        - Indicates that a database error has been trapped
*               - by ERRUS and that processing on this data base
*               - should end.
      logical*4 created, default_title, kill
 
c
      character* 8 source_alias(max_sources), site_alias(max_sites)
 
*            KalObs_name    - Name of the output data file
      character*128 KalObs_name
 
      character*40 database_file
 
      character*80 line
 
c
      equivalence (log_lu, ipass(1))
c
      common / user / ipass, icrt, database_file, today, cart,
     .    site_alias, source_alias, created, default_title,
     .    kill, delay_flag, line, KalObs_name
 
c
c
c LCODE common declarations
 
*   lcodes(7,max_lcodes)- equivalence which will have all of the LCODES
*                   - accessible.  Used when LCODES are listed.
*                   - Max_lcodes MUST reflect number of Lcodes because
*                   - it is used in the spacing of the common block
*                   - entries.
      character*14 lcodes(max_lcodes)
 
*   Fcodes(max_lcodes)  - Equivalence to the logicals which indicate
*                   - which LCODES were found (see comments above)
 
 
      logical Fcodes(max_lcodes)
 
* TYPE 1 Record LCODES:
* MOD TAH 891117:  All Lodes become character*8 instead of integer 
*                  arrays.
*   Lapr_etd     - LCODE for earth tide aprioris
*   Lapr_gamma   - LCODE for gamma apriori
*   Laxis_offsets- LCODE for axis offset values
*   Laxis_types  - LCODE for axis offset types
*   LCALC_ver    - LCODE for CALC version
*   Ldelay_flag  - LCODE for delay flag
*   Lellip_hgt   - LCODE for ellipsiod heights of the sites
*   Lfut1_inf    - LCODE for UT1 interpolation information
*   Lput1_inf    - LCODE for UT1 interpolation information
*   Leut1_inf    - LCODE for UT1 interpolation information
*   Lfut1_pts    - LCODE for UT1 values at the table epochs
*   Lput1_pts    - LCODE for UT1 values at the table epochs
*   Leut1_pts    - LCODE for UT1 values at the table epochs
*   Ltai_utc     - LCODE for TAI minus UTC
*   Lfwob_inf    - LCODE for wobble interpolation information
*   Lpwob_inf    - LCODE for wobble interpolation information
*   Lewob_inf    - LCODE for wobble interpolation information
*   Lfwob_pts    - LCODE for wobble values at the table epochs
*   Lpwob_pts    - LCODE for wobble values at the table epochs
*   Lewob_pts    - LCODE for wobble values at the table epochs
*   Linterval    - LCODE for start and stop times in db
*   Lnum_obs     - LCODE for number of observations
*   Lnum_sites   - LCODE for number of sites
*   Lnum_sources - LCODE for number of sources
*   Lsite_names  - LCODE for the names of the sites.
*   Lsite_pos    - LCODE for site positions
*   Lsource_names- LCODE for the names of the sources.
*   Lsource_pos  - LCODE for source positions
*   Lut1_tide    - LCODE for UT1 tidal flag
 
      character*14 Lapr_etd, Lapr_gamma, Laxis_offsets,
     .    Laxis_types, LCALC_ver, Ldelay_flag, Lellip_hgt,
     .    Lfut1_inf, Lput1_inf, Leut1_inf, Lfut1_pts,
     .    Lput1_pts, Leut1_pts, Ltai_utc, Lfwob_inf,
     .    Lpwob_inf, Lewob_inf, Lfwob_pts, Lpwob_pts,
     .    Lewob_pts, Linterval, Lnum_obs, Lnum_sites,
     .    Lnum_sources, Lsite_names, Lsite_pos,
     .    Lsource_names, Lsource_pos, Lut1_tide
 
* TYPE 2 Record LCODES
*   LBtape       - LCODE for BTape number for this observation
*   LFRNGE_qual  - LCODE for FRNGE quality code
*   Lazim        - LCODE for azimuth angles and rates
*   Lbaseline    - LCODE for site names for this observation
*   Lcont_atm    - LCODE for atmospheric (CHAO) delay.  Used only to
*                   - remove contribution from < Ver 6.0 data bases.
*   Lcont_ax_off - LCODE for axis offset contribution
*   Lcont_cable  - LCODE for cable calibration vaules
*   Lcont_etd    - LCODE for earth tide contribution
*   Lcont_gen_rel- LCODE for general relativity contibrution
*   Lcont_ocean  - LCODE for ocean loading contribution
*   Lcont_ptd    - LCODE for pole tide contribution
*   Lcorr_coeff  - LCODE for correlation coefficient
*   Ldata_amp    - LCODE for correlator residual phases and
*                   - amplitudes
*   Lelev        - LCODE for elevation angles and rates
*   Lend         - LCODE for end seconds (past hour) for obs.
*   Lfreq_num    - LCODE for number of frequencies in sequence
*   Lfreq_ref    - LCODE for phase reference frequency
*   Lfreq_seq    - LCODE for frequancy sequence
*   Lgrp_del_amb - LCODE for group delay ambiguity spacing
*   Lion_bits    - LCODE for quality of gion and pion entries
*   Lion_code    - LCODE for quality of gion (OBSELETE)
*   Lion_gcorr   - LCODE for group delay and phase delay rate ion
*                   - correction
*   Lion_pcorr   - LCODE for phase delay ion correction
*   Lion_rmsg    - LCODE for group delay and phase delay rate ion
*                   - sigma.
*   Lion_rmsp    - LCODE for phase delay ion sigma.
*   Lion_sigma   - LCODE for group delay and phase delay rate ion
*                   - sigma to be added to data sigma (OBSELETE)
*   Lmoon        - LCODE for ephemeris of moon
*   Lnum_gamb    - LCODE for number of group delay ambiquities
*   Lnum_pamb    - LCODE for number of phase delay ambiquities (must
*                   - use GETJ -- I*4 get)
*   Lnut_ang     - LCODE for WAHR series nutation angles.
*   Lobs_Gdelay  - LCODE for observed group delay
*   Lobs_NBdelay - LCODE for observed single band delay (OBSELETE)
*   Lobs_SBdelay - LCODE for observed single band delay
*   Lobs_rate    - LCODE for observed phase delay rate
*   Lobs_totphs  - LCODE for observed total phase
*   Lparangle    - LCODE for parralatic angle
*   Lpart_ax_off - LCODE for axis offset partial
*   Lpart_etd    - LCODE for earth tide partials
*   Lpart_gam    - LCODE for gamma partials
*   Lpart_nut    - LCODE for nutation angle partials
*   Lpart_parlx  - LCODE for parralax partials
*   Lpart_ut1    - LCODE for UT1 partials
*   Lpart_wob    - LCODE for wobble partials
*   Lphase_cal   - LCODE for phase cal phases, amplitudes and freq.
*   Lpressure    - LCODE for pressure and rates
*   Lrel_hum     - LCODE for relative humidity and rate
*   Lsec_tag     - LCODE for seconds tag of the epoch
*   Lsigma_SB    - LCODE for single band delay sigma
*   Lsigma_grp   - LCODE for group delay sigma
*   Lsigma_phs   - LCODE for phase delay sigma
*   Lsigma_rate  - LCODE for phase delay rate sigma
*   Lsite_part   - LCODE for site position partial
*   Lsnr         - LCODE for signal to noise ratio
*   Lsource_part - LCODE for source position partial
*   Lstar_ID     - LCODE for source names for this observation
*   Lstart       - LCODE for start seconds (past hour) for obs.
*   Lstop        - LCODE for stop seconds (past hour) for obs.
*   Lsun         - LCODE for ephemeris of sun
*   Ltemp_C      - LCODE for temperature and rate
*   Ltho_delay   - LCODE for theoretical delay
*   Ltho_rate    - LCODE for theoretical rate
*   Lhell_delay  - LCODE for Helling delay
*   Lhell_rate   - LCODE for Helling rate
*   Lshap_delay  - LCODE for Shapiro delay
*   Lshap_rate   - LCODE for Shapiro rate
*   Luncaldly    - LCODE for phase calibration delay to be removed
*                   - from group delay
*   Lunphsdel    - LCODE for new phase-decalibration entry.
*   Lunphsrat    - LCODE for phase cal rate contribution to
*                   - rate observable.
*   Lunw_gflag   - Group delay unweight flag from DB
*   Lunw_pflag   - phase delay unweight flag from DB
*   Lunw_rflag   - phase delay rate unweight flag from DB
*   Lutc_tag     - LCODE for UTC tag of the epoch
*   Lwvr_code    - LCODE for WVR quality code stored in data base
*   Lwvr_delay   - LCODE for WVR delay stored in data base
*   Lwvr_dely    - LCODE for WVR delay stored in data base (OBSELETE)
*   Lut1_val     - LCODE for actual value of UT1 used in Calc.
*   Lpol_val     - LCODE for actual value of XY pole pos. used in Calc.
 
      character*14 LBtape, LFRNGE_qual, Lazim, Lbaseline,
     .    Lcont_atm, Lcont_ax_off, Lcont_cable, Lcont_etd,
     .    Lcont_gen_rel, Lcont_ocean, Lcont_ptd,
     .    Lcorr_coeff, Ldata_amp, Lelev, Lend,
     .    Lfreq_num, Lfreq_ref, Lfreq_seq, Lgrp_del_amb,
     .    Lion_bits, Lion_code, Lion_gcorr, Lion_pcorr,
     .    Lion_rmsg, Lion_rmsp, Lion_sigma, Lmoon,
     .    Lnum_gamb, Lnum_pamb, Lnut_ang, Lobs_Gdelay,
     .    Lobs_NBdelay, Lobs_SBdelay, Lobs_rate,
     .    Lobs_totphs, Lparangle, Lpart_ax_off, Lpart_etd,
     .    Lpart_gam, Lpart_nut, Lpart_parlx, Lpart_ut1,
     .    Lpart_wob, Lphase_cal, Lpressure, Lrel_hum,
     .    Lsec_tag, Lsigma_SB, Lsigma_grp, Lsigma_phs,
     .    Lsigma_rate, Lsite_part, Lsnr, Lsource_part,
     .    Lstar_ID, Lstart, Lstop, Lsun, Ltemp_C,
     .    Ltho_delay, Ltho_rate, 
     .    Lhell_delay, Lhell_rate, Lshap_delay, Lshap_rate,
     .    Luncaldly, Lunphsdel,
     .    Lunphsrat, Lunw_gflag, Lunw_pflag, Lunw_rflag,
     .    Lutc_tag, Lwvr_code, Lwvr_delay, Lwvr_dely,
     .    Lut1_val, Lpol_val
 
**************************************************************************
 
*   Fapr_etd        - FCODE for earth tide aprioris
*   Fapr_gamma      - FCODE for gamma apriori
*   Faxis_offsets   - FCODE for axis offset values
*   Faxis_types     - FCODE for axis offset types
*   FCALC_ver       - FCODE for CALC version
*   Fdelay_flag     - FCODE for delay flag
*   Fellip_hgt      - FCODE for ellipsiod heights of the sites
*   Ffut1_inf       - FCODE for UT1 interpolation information
*   Fput1_inf       - FCODE for UT1 interpolation information
*   Feut1_inf       - FCODE for UT1 interpolation information
*   Ffut1_pts       - FCODE for UT1 values at the table epochs
*   Fput1_pts       - FCODE for UT1 values at the table epochs
*   Feut1_pts       - FCODE for UT1 values at the table epochs
*   Ftai_utc        - FCODE for TAI munus UTC
*   Ffwob_inf       - FCODE for wobble interpolation information
*   Fpwob_inf       - FCODE for wobble interpolation information
*   Fewob_inf       - FCODE for wobble interpolation information
*   Ffwob_pts       - FCODE for wobble values at the table epochs
*   Fpwob_pts       - FCODE for wobble values at the table epochs
*   Fewob_pts       - FCODE for wobble values at the table epochs
*   Finterval       - FCODE for start and stop times
*   Fnum_obs        - FCODE for number of observations
*   Fnum_sites      - FCODE for number of sites
*   Fnum_sources    - FCODE for number of sources
*   Fsite_names     - FCODE for the names of the sites.
*   Fsite_pos       - FCODE for site positions
*   Fsource_names   - FCODE for the names of the sources.
*   Fsource_pos     - FCODE for source positions
*   Fut1_tide       - FCODE for UT1 tidal value
 
      logical Fapr_etd, Fapr_gamma, Faxis_offsets, Faxis_types,
     .    FCALC_ver, Fdelay_flag, Fellip_hgt, Ffut1_inf, Fput1_inf,
     .    Feut1_inf, Ffut1_pts, Fput1_pts, Feut1_pts, Ftai_utc,
     .    Ffwob_inf, Fpwob_inf, Fewob_inf, Ffwob_pts, Fpwob_pts,
     .    Fewob_pts, Finterval, Fnum_obs, Fnum_sites, Fnum_sources,
     .    Fsite_names, Fsite_pos, Fsource_names, Fsource_pos, Fut1_tide
 
* TYPE 2 Record FCODES
*   FBtape          - FCODE for BTape number for this observation
*   FFRNGE_qual     - FCODE for FRNGE quality code
*   Fazim           - FCODE for azimuth angles and rates
*   Fbaseline       - FCODE for site names for this observation
*   Fcont_atm       - FCODE for atmospheric (CHAO) delay.  Used only to
*                   - remove contribution from < Ver 6.0 data bases.
*   Fcont_ax_off    - FCODE for axis offset contribution
*   Fcont_cable     - FCODE for cable calibration vaules
*   Fcont_etd       - FCODE for earth tide contribution
*   Fcont_gen_rel   - FCODE for general relativity contibrution
*   Fcont_ocean     - FCODE for ocean loading contribution
*   Fcont_ptd       - FCODE for pole tide contribution
*   Fcorr_coeff     - FCODE for correlation coefficient
*   Fdata_amp       - FCODE for correlator residual phases and
*                   - amplitudes
*   Felev           - FCODE for elevation angles and rates
*   Fend            - FCODE for end seconds (past hour) for obs.
*   Ffreq_num       - FCODE for number of frequencies in sequence
*   Ffreq_ref       - FCODE for phase reference frequency
*   Ffreq_seq       - FCODE for frequancy sequence
*   Fgrp_del_amb    - FCODE for group delay ambiguity spacing
*   Fion_bits       - FCODE for quality of gion and pion entries
*   Fion_code       - FCODE for quality of gion (OBSELETE)
*   Fion_gcorr      - FCODE for group delay and phase delay rate ion
*                   - correction
*   Fion_pcorr      - FCODE for phase delay ion correction
*   Fion_rmsg       - FCODE for group delay and phase delay rate ion
*                   - sigma.
*   Fion_rmsp       - FCODE for phase delay ion sigma.
*   Fion_sigma      - FCODE for group delay and phase delay rate ion
*                   - sigma to be added to data sigma (OBSELETE)
*   Fmoon           - FCODE for ephemeris of moon
*   Fnum_gamb       - FCODE for number of group delay ambiquities
*   Fnum_pamb       - FCODE for number of phase delay ambiquities (must
*                   - use GETJ -- I*4 get)
*   Fnut_ang        - FCODE for WAHR series nutation angles.
*   Fobs_Gdelay     - FCODE for observed group delay
*   Fobs_NBdelay    - FCODE for observed single band delay (OBSELETE)
*   Fobs_SBdelay    - FCODE for observed single band delay
*   Fobs_rate       - FCODE for observed phase delay rate
*   Fobs_totphs     - FCODE for observed total phase
*   Fparangle       - FCODE for parralatic angle
*   Fpart_ax_off    - FCODE for axis offset partial
*   Fpart_etd       - FCODE for earth tide partials
*   Fpart_gam       - FCODE for gamma partials
*   Fpart_nut       - FCODE for nutation angle partials
*   Fpart_parlx     - FCODE for parralax partials
*   Fpart_ut1       - FCODE for UT1 partials
*   Fpart_wob       - FCODE for wobble partials
*   Fphase_cal      - FCODE for phase cal phases, amplitudes and freq.
*   Fpressure       - FCODE for pressure and rates
*   Frel_hum        - FCODE for relative humidity and rate
*   Fsec_tag        - FCODE for seconds tag of the epoch
*   Fsigma_SB       - FCODE for single band delay sigma
*   Fsigma_grp      - FCODE for group delay sigma
*   Fsigma_phs      - FCODE for phase delay sigma
*   Fsigma_rate     - FCODE for phase delay rate sigma
*   Fsite_part      - FCODE for site position partial
*   Fsnr            - FCODE for signal to noise ratio
*   Fsource_part    - FCODE for source position partial
*   Fstar_ID        - FCODE for source names for this observation
*   Fstart          - FCODE for start seconds (past hour) for obs.
*   Fstop           - FCODE for stop  seconds (past hour) for obs.
*   Fsun            - FCODE for ephemeris of sun
*   Ftemp_C         - FCODE for temperature and rate
*   Ftho_delay      - FCODE for theoretical delay
*   Ftho_rate       - FCODE for theoretical rate
*   Fhell_delay     - FCODE for Helling delay
*   Fhell_rate      - FCODE for Helling rate
*   Fshap_delay     - FCODE for Shapiro delay
*   Fshap_rate      - FCODE for Shapiro rate
*   Funcaldly       - FCODE for phase calibration delay to be removed
*                   - from group delay
*   Funphsdel       - FCODE for new phase-decalibration entry.
*   Funphsrat       - FCODE for phase cal rate contribution to
*                   - rate observable.
*   Funw_gflag      - Group delay unweight flag from DB
*   Funw_pflag      - phase delay unweight flag from DB
*   Funw_rflag      - phase delay rate unweight flag from DB
*   Futc_tag        - FCODE for UTC tag of the epoch
*   Fwvr_code       - FCODE for WVR quality code stored in data base
*   Fwvr_delay      - FCODE for WVR delay stored in data base
*   Fwvr_dely       - FCODE for WVR delay stored in data base (OBSELETE)
*   Fut1_val        - FCODE for actual value of UT1 used in Calc.
*   Fpol_val        - FCODE for actual value of XY pole pos. used in Calc.
 
      logical FBtape, FFRNGE_qual, Fazim, Fbaseline, Fcont_atm,
     .    Fcont_ax_off, Fcont_cable, Fcont_etd, Fcont_gen_rel,
     .    Fcont_ocean, Fcont_ptd, Fcorr_coeff, Fdata_amp, Felev, Fend,
     .    Ffreq_num, Ffreq_ref, Ffreq_seq, Fgrp_del_amb, Fion_bits,
     .    Fion_code, Fion_gcorr, Fion_pcorr, Fion_rmsg, Fion_rmsp,
     .    Fion_sigma, Fmoon, Fnum_gamb, Fnum_pamb, Fnut_ang,
     .    Fobs_Gdelay, Fobs_NBdelay, Fobs_SBdelay, Fobs_rate,
     .    Fobs_totphs, Fparangle, Fpart_ax_off, Fpart_etd, Fpart_gam,
     .    Fpart_nut, Fpart_parlx, Fpart_ut1, Fpart_wob, Fphase_cal,
     .    Fpressure, Frel_hum, Fsec_tag, Fsigma_SB, Fsigma_grp,
     .    Fsigma_phs, Fsigma_rate, Fsite_part, Fsnr, Fsource_part,
     .    Fstar_ID, Fstart, Fstop, Fsun, Ftemp_C, Ftho_delay,
     .    Ftho_rate, Fhell_delay, Fhell_rate, Fshap_delay, Fshap_rate,
     .    Funcaldly, Funphsdel, Funphsrat, Funw_gflag,
     .    Funw_pflag, Funw_rflag, Futc_tag, Fwvr_code, Fwvr_delay,
     .    Fwvr_dely , Fut1_val  , Fpol_val
 
      equivalence ( Lcodes,Lapr_etd )
      equivalence ( Fcodes,Fapr_etd )
 
***** START DECLARATION OF COMMON BLOCK
 
      common /DB_LCODES/ Lapr_etd, Lapr_gamma, Laxis_offsets,
     .    Laxis_types, LCALC_ver, Ldelay_flag, Lellip_hgt, Lfut1_inf,
     .    Lput1_inf, Leut1_inf, Lfut1_pts, Lput1_pts, Leut1_pts,
     .    Ltai_utc, Lfwob_inf, Lpwob_inf, Lewob_inf, Lfwob_pts,
     .    Lpwob_pts, Lewob_pts, Linterval, Lnum_obs, Lnum_sites,
     .    Lnum_sources, Lsite_names, Lsite_pos, Lsource_names,
     .    Lsource_pos, Lut1_tide,  LBtape, LFRNGE_qual, Lazim,
     .    Lbaseline, Lcont_atm, Lcont_ax_off, Lcont_cable, Lcont_etd,
     .    Lcont_gen_rel, Lcont_ocean, Lcont_ptd, Lcorr_coeff,
     .    Ldata_amp, Lelev, Lend, Lfreq_num, Lfreq_ref, Lfreq_seq,
     .    Lgrp_del_amb, Lion_bits, Lion_code, Lion_gcorr, Lion_pcorr,
     .    Lion_rmsg, Lion_rmsp, Lion_sigma, Lmoon, Lnum_gamb,
     .    Lnum_pamb, Lnut_ang, Lobs_Gdelay, Lobs_NBdelay, Lobs_SBdelay,
     .    Lobs_rate, Lobs_totphs, Lparangle, Lpart_ax_off, Lpart_etd,
     .    Lpart_gam, Lpart_nut, Lpart_parlx, Lpart_ut1, Lpart_wob,
     .    Lphase_cal, Lpressure, Lrel_hum, Lsec_tag, Lsigma_SB,
     .    Lsigma_grp, Lsigma_phs, Lsigma_rate, Lsite_part, Lsnr,
     .    Lsource_part, Lstar_ID, Lstart, Lstop, Lsun, Ltemp_C,
     .    Ltho_delay, Ltho_rate,  
     .    Lhell_delay, Lhell_rate, Lshap_delay, Lshap_rate,
     .    Luncaldly, Lunphsdel, Lunphsrat,
     .    Lunw_gflag, Lunw_pflag, Lunw_rflag, Lutc_tag, Lwvr_code,
     .    Lwvr_delay, Lwvr_dely, Lut1_val, Lpol_val
 
      common /DB_FCODES/ Fapr_etd, Fapr_gamma, Faxis_offsets, 
     .    Faxis_types,
     .    FCALC_ver, Fdelay_flag, Fellip_hgt, Ffut1_inf, Fput1_inf,
     .    Feut1_inf, Ffut1_pts, Fput1_pts, Feut1_pts, Ftai_utc,
     .    Ffwob_inf, Fpwob_inf, Fewob_inf, Ffwob_pts, Fpwob_pts,
     .    Fewob_pts, Finterval, Fnum_obs, Fnum_sites, Fnum_sources,
     .    Fsite_names, Fsite_pos, Fsource_names, Fsource_pos, 
     .    Fut1_tide,   FBtape,
     .    FFRNGE_qual, Fazim, Fbaseline, Fcont_atm, Fcont_ax_off,
     .    Fcont_cable, Fcont_etd, Fcont_gen_rel, Fcont_ocean,
     .    Fcont_ptd, Fcorr_coeff, Fdata_amp, Felev, Fend, Ffreq_num,
     .    Ffreq_ref, Ffreq_seq, Fgrp_del_amb, Fion_bits, Fion_code,
     .    Fion_gcorr, Fion_pcorr, Fion_rmsg, Fion_rmsp, Fion_sigma,
     .    Fmoon, Fnum_gamb, Fnum_pamb, Fnut_ang, Fobs_Gdelay,
     .    Fobs_NBdelay, Fobs_SBdelay, Fobs_rate, Fobs_totphs,
     .    Fparangle, Fpart_ax_off, Fpart_etd, Fpart_gam, Fpart_nut,
     .    Fpart_parlx, Fpart_ut1, Fpart_wob, Fphase_cal, Fpressure,
     .    Frel_hum, Fsec_tag, Fsigma_SB, Fsigma_grp, Fsigma_phs,
     .    Fsigma_rate, Fsite_part, Fsnr, Fsource_part, Fstar_ID,
     .    Fstart, Fstop, Fsun, Ftemp_C, Ftho_delay, Ftho_rate,
     .    Fhell_delay,  Fhell_rate, Fshap_delay, Fshap_rate,
     .    Funcaldly, Funphsdel, Funphsrat, Funw_gflag, Funw_pflag,
     .    Funw_rflag, Futc_tag, Fwvr_code, Fwvr_delay, Fwvr_dely,
     .    Fut1_val, Fpol_val
 
