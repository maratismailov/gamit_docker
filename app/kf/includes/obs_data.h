 
*------------------------------------------------------------------
*                                                     OBS_DATA.FTNI
*     The data part definition of the new Kalman filter files.
*     There will be one of these groups of records for each type 2
*     record in the data base.
*
*     T.Herring                   10:08 PM MON., 12 Jan., 1987
*------------------------------------------------------------------
 
*   Btape       - Btape number for this observation
*   data_flag   - Bit masked which indicates any problem
*               - with the data at this epoch.  Not all of these
*               - bits are turned on in the data stored in the
*               - data file, but are included here for
*               - consistency with DELETE_COUNT priority scheme.
*               - Those actually not set in data file are shown
*               - with a star (*) next to bit.  The bit masking
*               - is:
*               - Bit   Meaning
*               -   1   Database lists as bad FRNGE quality
*               - * 2   ION delay is bad (fails ION_MASK)
*               - * 3   WVR temperatures bad
*               - * 4   WVR Dbase delay bad
*               - * 5   Propagation medium problem
*               - * 6   Cable Cal bad
*               -   7   SOLVK interactive edit
*               -   8   SOLVK automatic edit
*               -   9   Data base grp/rate interactive edit (unw=1)
*               -  10   phase delay connection problem
*               -  11   User defined editing flag
*               - *12   User source delete
*               - *13   Down-site delete due to time interval
*               - *14   Below elevation cutoff
*               - *15   Data do not close within tolerance.
*               -  16   RESERVED for PLTSL program use.
*   data_amp(2,max_channels)  - Amplitude and Phase of the
*               - signal from frnge (d-6, 0.001 deg)
*   ion_flag    - Bit mappped word giving the status of the group
*               - and phase ionospheric corrections. (Corresponds
*               - to ICORR and LCODE ION_BITS in GSFC software).
*               - The meaning of the bits is:
*               - Bit   Meaning
*               -   1   matching observation has quality code of
*               -       8 or 9, but unweighted in SOLVE
*               -   2   No matching group delay data
*               -   3   Matching observation has quality code of
*               -       1 to 7.
*               -   4   group ION correction available
*               -   5   Down weight flag for group ion correction
*               -       (Can only be reset if bit 2 not turned on)
*               -   6   Matching observation has quality code 0
*               -   PHASE ION VALUES
*               -   7   matching observation has quality code of
*               -       8 or 9, but unweighted in SOLVE
*               -   8   No matching group delay data
*               -   9   Matching observation has quality code of
*               -       1 to 7.
*               -  10   phase ION correction available
*               -  11   Down weight flag for group ion correction
*               -       (Can only be reset if bit 2 not turned on)
*               -  12   Matching observation has quality code 0
*               -       (no Fringes)
*   num_channels    - Number of channels in the data (VLBI
*               - channels for BWS)
*   num_grp_amb - number of group delay ambiquities
*               - (see also I*4 num_phs_amb)
*   num_wvr_frq(2)  - number of WVR frequencies at each site
*   medium_flag(2)  - Records the quality of the atmopsheric
*               - data at the two sites.  Bit masked
*               - words mean:
*               - Bit   Meaning
*               -   1   Pressure bad
*               -   2   Temperature bad
*               -   3   relative humidity bad
*               -   4   WVR antenna temps bad
*               -   5   Data Base WVR delay bad
*               -   6   Cable Cal bad.
*   phase_cal(2,2,max_channels) - phase cal amplitude and phase
*               - in each channel at each site (d-6, .001 deg)
 
*   site(2)     - site numbers for this observation
*   source      - source number for this observation
 
*   wvr_code(2) - Code for WVR values read from data base
*               - The meaning of these values is:
*               - First 2 digits:
*               - -1  No WVR data available
*               -  0  WVR data with in 2 deg of obs position, and
*               -     taken during the observation
*               -  1  WVR data taken during observation, but > 2 deg
*               -     from obs position
*               -  2  WVR data mapped from zenith point taken before
*               -     observation
*               -  3  WVR data interpolated from 2 closest available
*               -     time points, one on each side of observation.
*               -     These points are with in 30 minutes of mid-epoch
*               -     of observation.
*               -  4  Same as 3 but data points are > 30 minutes from
*               -     observation
*               -  5  WVR data extrapolated (data only before or after
*               -     observation)
*               -  6  WVR data taken from zenith delays measured during
*               -     observation period.
*               -
*               - In addition if more than one WVR data point then
*               - number is added to code by (CODE+N*100).
 
 
 
      integer*4 Btape, data_flag, data_amp(2,max_channels), ion_flag,
     .    num_channels, num_grp_amb, num_wvr_frq(2), medium_flag(2),
     .    phase_cal(2,2,max_channels), site(2), source, wvr_code(2)
 
*   num_phs_amb - Number of phase delay ambiquities from
*               - the phase delay generated at read in time
 
 
 
      integer*4 num_phs_amb
 
*   FRNGE_qual  - FRNGE quality code from data base.
*   dum_fil_char - Charcater filler to get words on correct boundaries
 
 
      character*2 FRNGE_qual, dum_fil_char
 
*   azimuth(2,2)    - Azimuth and rate of change (rad, rad/s)
*   atm_part(2,2)   - Atmospheric delay partial for site 1 and 2
*               - and for delay and rate.  This value is not
*               - stored in the file.  It is generated during the
*               - KALGN run to match the current calibration being
*               - used.
*   AXO_part(2,2) - Axis offset partials sites 1 and 2,
*               - and delay and rate
*   AXO_cont(2) - Axis offset contribution to delay and rate (ps.fs/s)
 
*   cable_cont(2,2) - Cable calibration delay and rate
 
*   corr_coeff  - The correlation coefficient (1d-4 units)
 
*   db_sigma(4) - Sigmas for the grp, phs and SB delays and rates as
*               - read from the data base (ps and fs/s).  These values
*               - should never be changed.  Sigma is used for
*               - solutions.)
*   duration    - duration of the observation (sec)
*   elev(2,2)   - elevation angles and rate (rad, rad/sec)
 
*   etd_cont(2) - Earth tide contribution (delay,rate)
*   etd_part(3,2,2)  - Earth tide partials for h,l, and lag,
*               - at sites 1&2 for delay and rate.
*   etd_ext_part(12,2,2) - Full earth tide partials.  These
*               - values are for h and ln and le, in and out-of-phase in
*               - diurnal (K1) and semidiurnal (M2) bands for
*               - each site and delay and rate.
*   feed_rot_cont(2,2)  - Feed rotation contribution at site 1 and
*               - and 2 for phase delay and rate.
 
*   freq_seq(max_channels)  - difference from phs_ref_freq
*               - of the frequencies of the channels
 
*   gamma_part(2) - Gamma partial for delay and rate
 
*   gen_rel_cont(2) - General relativity contribution
*   ion_corr(3) - ionospheric correction (grp, phs, rate)
*   ion_sigma(3) - ionospheric correction sigma (grp, phs,
*               - rate)
*   nut_part(2,2)   - Nutation partials for dspi and deps for
*               - delay and rate
 
*   ocean_cont(2)   - Ocean loading contribution
*   parlx_part(2)   - Parallax partials for delay and rate
*   ptd_cont(2) - Pole tide contribution (delay,rate)
 
*   phase_cal_cont(2,2) - Contribution of phase cal to for sites 1
*                   - &2 for delay and rate. The signs are set so
*                   - these values are the contribution to the
*                   - model. (ps, fs/s)
 
*   pmu_part(3,2)   - Polar motion, UT1 partials for delay
*               - and rate. (ps/mas for delay, all components;
*               - (fs/s)/mas for rate, all components)
 
*   pressure(2,2)   - Pressure value and rate at the two
*               - sites (mbar and mbar/sec)
*   rel_hum(2,2)    - relative humidity and rate (0-1 and
*               - 0-1/s)
 
*   sigma(4)    - Data sigmas for grp,phs, SB delay and rate
*               - (ps and fs/s)
*   snr         - Data signal-to-noise ratio
 
*   site_part(3,2,2)    - site partial derivatives for XYZ for
*               - sites 1 and 2, and delay and rate
*   source_part(2,2)    - source position partials RA and Dec,
*               - delay and rate.
*   source_struc_cont(3)    - Contribution from source structure
*               - for the grp, phs and rate observable
 
*   temp_C(2,2) - Temperature and rate at the two sites
*               - (C and C/sec)
 
*   user_cont(2)    - User contribution to delay and rate. You
*               - can put anything you like here (ps and fs/s)
*   user_dry_zen(2) - User dry zenith delay.  Again any thing you
*               - like (ps and fs/sec)
*   user_wet_zen(2) - User wet zenith delay. (ps and ps/s)
 
 
*   WVR_azim(2) - Azimuth of the WVR when measurement made (rad)
*   WVR_dbase(2,2)  - WVR values from the data base (zenith site
*               - 1 and 2 and delay and rate (ps fs/s)
*   WVR_elev(2) - Elevation of the WVR when measurement made (rad)
*   WVR_frq(max_frq_wvr,2)  - WVR frequencies
*   WVR_tant(max_frq_wvr,2,2)   - WVR antenna temperatures and
*               - rates at the two sites
*   WVR_tant_sig(max_frq_wvr,2,2)   - sigmas for the WVR antenna
*               - temperatures and rates
 
 
 
      real*4 azimuth(2,2), atm_part(2,2), AXO_part(2,2), AXO_cont(2),
     .    cable_cont(2,2), corr_coeff, db_sigma(4), duration,
     .    elev(2,2), etd_cont(2), etd_part(3,2,2),
     .    etd_ext_part(12,2,2),
     .    feed_rot_cont(2,2), freq_seq(max_channels), gamma_part(2),
     .    gen_rel_cont(2), ion_corr(3), ion_sigma(3), nut_part(2,2),
     .    ocean_cont(2), parlx_part(2), ptd_cont(2),
     .    phase_cal_cont(2,2), pmu_part(3,2), pressure(2,2),
     .    rel_hum(2,2), sigma(4), snr, site_part(3,2,2),
     .    source_part(2,2), source_struc_cont(3), temp_C(2,2),
     .    user_cont(2), user_dry_zen(2), user_wet_zen(2), WVR_azim(2),
     .    WVR_dbase(2,2), WVR_elev(2), WVR_frq(max_frq_wvr,2),
     .    WVR_tant(max_frq_wvr,2,2), WVR_tant_sig(max_frq_wvr,2,2)
 
* New variables for atmospheric conditions and general Earth rotation
 
*   lapse_usr(2)    - Lapse rate at each site (K/km)
*   ht_usr(2)       - Height of tropopause bt site (km)
*   tbias_usr(2)    - Bias in surface temperature by site (K)
*   adry_usr(4,2)   - Coefficients for dry mapping function by site
*                   - (dimensionless)
*   bdry_usr(4,2)   - Coefficienst for bending due to dry air
*   awet_usr(4,2)   - Coefficinets for wet mapping function
*   bwet_usr(4,2)   - Coefficients for bending due to wet air
 
*   eor_ut1_part(4,2)   - Partials for diurnal cos and sin and
*                   - Semidiurnal cos and sin variations in UT1
*   eor_xy_part(6,2)    - Partials for diurnal retrograde PM,
*                   - and semidiurnal Prograde and retrograde PM
*                   - (paired in quadrature components)
 
 
      real*4 lapse_usr(2), ht_usr(2), tbias_usr(2), adry_usr(4,2),
     .    bdry_usr(4,2), awet_usr(4,2), bwet_usr(4,2),
     .    eor_ut1_part(4,2), eor_xy_part(6,2)
 
*   db_theoretical(4)   - The theorectical grp, phase, ans SB delays
*               - and rates.  (Value read from database and never
*               - changed.) (All calibrations are applied to a copy
*               - of db_theoretical)
*   epoch       - Julian epoch with fractional part of a day of
*               - epoch of this measurement
*   grp_amb     - Value of the group delay ambiquity
*   observable(4)  - observed value for grp delay, phase delay
*               - SB delay and phase delay rate (ps, ps, fs/sec)
*   phs_amb     - Value of the phase delay ambiquity
*   phs_ref_freq    - Phase delay reference frequency (Hz)
*   theoretical(4)    - theoretical values for grp delay,
*               - phase delay, SB delay and phase delay rate.
 
 
 
      real*8 db_theoretical(4), epoch, grp_amb, observable(4), phs_amb,
     .    phs_ref_freq, theoretical(4)

*   Additional variables added for Calc 7.0 Data bases
*   pmu_calc(3) - Actual values of Pole x,y and UT1-AT used in the
*                 calc solution for this observation (mas for all)
*   clock_jump(2,2) - Clock jump and rate at the sites in the baselines.  Sign is
*                 NOT APPLIED.  Therefore site 1 should be substacted and
*                 and site 2 added. (ps and fs/s)

      real*8 pmu_calc(3), clock_jump(2,2)

*   ocean_rad(2,2) - Site dependent ocean radial loading by sites
*                    (index 1) and delay and rate (2) (ps and fs/s)
*   ocean_horz(2,2) - Site dependent ocean horizontal loading 
*                     contribution by site and delay/rate.

      real*4 ocean_rad(2,2), ocean_horz(2,2)
 
*   last_data_word  - last word of data record
*   dummy_data(127) - dummy space to stop overwritting records
*               - when file is read.
* MOD TAH 900711: Removed 6 words for pmu_calc variable (allows
*                 compatability with old Calc 6.0 KalObs files)
*                 Decrement dummy_data to 121
* MOD TAH 901001: Removed 4 words for the clock_jump
* MOD TAH 910103: Removed 4 words for the clock_jump rate terms
* MOD TAH 910905: Removed 8 words for the ocean_rad and ocean_horz
*                 variables.(Was 113 for dummy_data)

* WARNING:  When dummy data is changed, make sure that HPK_to_UNIX is also
*     changed in routine convert_data.
 
      integer*4 last_data_word, dummy_data(105)
 
********************
* COMMON DECLARATION
********************
 
 
      common / data_block / site, source, Btape, data_flag, data_amp,
     .    ion_flag, num_channels, num_grp_amb, num_wvr_frq,
     .    medium_flag, phase_cal, wvr_code, num_phs_amb, FRNGE_qual,
     .    dum_fil_char,
     .    azimuth, atm_part, AXO_part, AXO_cont, cable_cont,
     .    corr_coeff, db_sigma, duration, elev, etd_cont, etd_part,
     .    etd_ext_part, feed_rot_cont, freq_seq, gamma_part,
     .    gen_rel_cont, ion_corr, ion_sigma, nut_part, ocean_cont,
     .    parlx_part, ptd_cont, phase_cal_cont, pmu_part, pressure,
     .    rel_hum, sigma, snr, site_part, source_part,
     .    source_struc_cont, temp_C, user_cont, user_dry_zen,
     .    user_wet_zen, WVR_azim, WVR_dbase, WVR_elev, WVR_frq,
     .    WVR_tant, WVR_tant_sig,
     .    lapse_usr, ht_usr, tbias_usr, adry_usr,
     .    bdry_usr, awet_usr, bwet_usr, eor_ut1_part, eor_xy_part,
     .    db_theoretical, epoch, grp_amb,
     .    observable, phs_amb, phs_ref_freq, theoretical,
     .    pmu_calc, clock_jump, ocean_rad, ocean_horz,
     .    last_data_word, dummy_data
 
*
