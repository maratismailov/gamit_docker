 
*     This is the individual one-way record definitions for the Gobs
*     files.
 
*   site        - Ground site number.  Relates to position and
*                 atmospheres.  (Clocks and Biases are related to
*                 the rcv number.
*   rcv         - Receiver number.  (May be differerent from
*                 site number if a kinematic site has been 
*                 processed as mulitple sites.
*   kd_site     - Kinematic/dynamic site number (0 if this site is not
*                 kinematic/dynamic).  This site number indicates
*                 that there is a time dependent position available
*                 for this site.
*   orb_site    - Orbital site number if this is an orbital
*                 receiver.
*   svs         - Satellite number from the list of satelllites
*   data_flag   - Bit mapped word indicating status of this
*                 observation.  If this arrangements of flags is
*                 changed then cerrfl_to_df should be updated.
*                 The meaning of the BITs is
*                 Bit Meaning
*                   1  Low SNR (recorded in SOLVG)
*                   2  Low SNR from Cfile (GAMIT Error flag 3)
*                   3  Rclock deviates more than limit (ctogobs)
*                   4  Bad data in GAMIT (GAMIT Error flag 2).  Also
*                      used if L1 and L2 range deviate by more than 100 m.
*                   5  Bias too close to end of arc (autcln)
*                   6  Not enough data between biases (autcln)
*                   7  Insufficient data in normal point (ctogobs)
*                   8  Too much scatter in NPoint (ctogobs)
*                   9  Too much multipath (solvg)
*                  10  WVR data bad
*                  11  Met. data bad (when available)
*                  12  Prefit phase residual edit (ctogobs)
*                  13  Prefit range residual edit (ctogobs)
*                  14  Flagged bad in GAMIT (error flag -1)
*                  15  Flagged low-elev in GAMIT (error flag 4)
*                  16  Below min. elev cutoff in SOLVG
*                  17  Edit with DOWN_SITE command
*                  18  Edit with DOWN_SVS command
*                  19  Autcln max_scan_dd edit (to many DDScan flags).
*                  20  AUTCLN Postfit residual edit (autcln)
*                  21  BAD Model Deleted Data and Unknown GAMIT edit flag
*                  22  Bad L2 data, but L1 is possibly good.
*                  23  Good data but no double differences for this
*                      one way measurent. 
*                  24  Data below min_ctog_elev (not cleaned)
*                  25  Pre-edited in ctogobs using pre-edit command
*                  26  MINIMac range data (ignore these ranges)
*                  27  Data below final elevation cutoff (mask out except
*                      when normal points are formed).
*                  28  DCB applied to C1 range
*                  29  DCB applied to both C1 and P2 ranges.
*                  ......
*                  32  GAMIT Bias flag
*                  31  CTOGOBS Bias flag.
*                  30  No data (used for internal book keeping)
*   kine_flag   - Kinematic data flag (indicates status of kinematic
*                 data.
*                 BIT  Meaning
*                   1  Receiver is moving.  Do not use data.
 
 
 
      integer*4 site, rcv, kd_site, orb_site, svs, data_flag, kine_flag

 
*   elev(2)     - Elevation angle and rate of change (rad and
*                 rad/s)
*   azimuth(2)  - Azimuth and rate of change (rad and rad/s)
*   snr(max_gdata_types) - SNR for each observable type
*   multipath_cont(max_gdata_types) - multipath constribution of
*                 each obervable.
*   wvr_los     - Line of site WVR calibration (meters)
 
 
      real*4 elev(2), azimuth(2), snr(max_gdata_types),
     .    multipath_cont(max_gdata_types), wvr_los
 
*     Primary partial deviatives
 
*   clk_partial - L1 carrier phase partial with respect to the
*                 ground receiver clock (L1 cycles per nsec)
*   site_part(3)    - L1 carrier phase partial with respect to
*                 cartesina XYZ coordinates of the receiver/
*   svs_partials(9) - L1 carrier phase partial with respect to
*                 6 orbit elements and 3 radiation parameters
*                 (L1 cycles per m, per mm/s and per unitless)
 
 
      real*8 clk_partial
      real*4 site_part(3), svs_partials(9)
 
*     Observables and theoreticals
 
*   observable(max_gdata_types) - Observables dependent on the
*                 data_types for this receiver (all units cycles
*                 at the appropriate frequency.
 
*   cycle_biases(max_gdata_types)   - Apriori values of the number
*                 of cycles needed for this data type.  Although
*                 this is real*8 it should contain an integer for
*                 full wavelength receivers.  (Real*8 format is
*                 capable of represprenting integers over the full
*                 range of I*4 (as far as I can tell) and therefore
*                 there should be loss of precision using real*8 and
*                 half wavelengths can be easily included.  For
*                 range measurements, millsecond clock jumps are
*                 also stored here.  (These are removed from the 
*                 clock model.)
*   cf_theoreticals(max_gdata_types)    - Theoretical values
*                 dependendent on data type (all cycles at
*                 appropriate frequency)
*   cf_sigmas(max_gdata_types)   - Sigmas of each data type.
*                 (Computed in ctogobs) (Cycles at appropriate
*                 frequency)
 
 
      real*8 observable(max_gdata_types),
     .    cycle_biases(max_gdata_types),
     .    cf_theoreticals(max_gdata_types) 

      real*4 cf_sigmas(max_gdata_types)
 
*   last_data_word   - Last word (to be saved to disk) of the
*                 obs records.
 
 
      integer*4 last_data_word
 
***** Extended part of observation record (Computed at run-time).
 
*   ocean_rad_cont(max_gdata_types) - Radial ocean tide
*                 contribution to observables (cycles)
*   ocean_hor_cont(max_gdata_types) - Horizontal ocean tide
*                 contribution to observables (cycles)
*   dum_dat1    - Dummy word for alignment
 
 
      real*4 ocean_rad_cont(max_gdata_types),
     .    ocean_hor_cont(max_gdata_types), dum_dat1
 
*   atm_part(max_gdata_types)  - Atmospheric delay partial
*                 (Computed based on the atmospheric model used).
*   atm_az_part(max_gdata_types)  - Atmospheric azimuthal
*                 assymetry partial
*   atm_map_part(max_gdata_types)  - Mapping function partials
*   ion_los_part(max_gdata_types)   - Ionsoheric delay line
*                 of site partial.
*   iob_mod_part(max_gdata_types)   - Ionspheric delay
*                 model partial (paremtrized as zenith
*                 ionopshere)
*   pmu_part(3) - X and Y wobble and UT1 partials (all cycles
*                 per milliarcsec, including UT1)
*   eor_ut1_part(4,max_gdata_types) - Diurnal and
*                 semidiurnal quaduature (cos and sin) partials
*   eor_xy_part(6,max_gdata_types)  - Diurnal and
*                 semidiurmal quadrature partials for polar
*                 motion.
*   nut_part(2,max_gdata_types) - Nutation partials
*                 (implemented as diurnal polar motion)
*   align_r4     - Real*4 for allignmnet reasons
*   theoreticals(max_gdata_types) - Theorectical delay with all
*                 calibrations and model changes applied.
*   sigmas(max_gdata_types) - Sigmas with all of the modifications
*                 applied.
 
 
      real*4 atm_part(max_gdata_types), atm_az_part(max_gdata_types),
     .    atm_map_part(max_gdata_types), ion_los_part(max_gdata_types),
     .    iob_mod_part(max_gdata_types), pmu_part(3),
     .    eor_ut1_part(4,max_gdata_types),
     .    eor_xy_part(6,max_gdata_types), nut_part(2,max_gdata_types),
     .    sigmas(max_gdata_types), align_r4

      real*8 theoreticals(max_gdata_types)
 
*..................................................................
 
 
 
      common / gobs_data / site, rcv, kd_site, orb_site, svs,
     .    data_flag, kine_flag, elev, azimuth, snr,
     .    multipath_cont, wvr_los,
     .    clk_partial, site_part, svs_partials,
     .    observable,  cycle_biases, cf_theoreticals, cf_sigmas,
     .    last_data_word, ocean_rad_cont, ocean_hor_cont, dum_dat1,
     .    atm_part, atm_az_part,
     .    atm_map_part, ion_los_part, iob_mod_part, pmu_part,
     .    eor_ut1_part, eor_xy_part, nut_part, align_r4,
     .    theoreticals, sigmas
 
