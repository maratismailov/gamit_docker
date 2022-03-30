*------------------------------------------------------------------
*                                                      OBS_APR.FTNI
*     The apriori's part of the header records of the new Kalman
*     filter files.
*
*     T.Herring                   09:11 PM MON., 12 Jan., 1987
*------------------------------------------------------------------
 
*   aprioris    - First word of the aprioris block.  This header
*               - word can be used to read the file and is
*               - eqivalanced to AXIS_TYPES.
 
      integer*4 aprioris
 
*   axis_types(max_sites)    - types of axes configuration
*               - at each site (CALC definition)
*   tai_utc         - Difference between tai and utc (seconds) so
*               - then ut1-utc can be computed.
 
 
      integer*4 axis_types(max_sites), tai_utc
 
*   atm_mar_apr(3,max_sites)    - site/experiment dependent values
*               - the atmosphere markov values
*   atm_az_mar_apr(2,max_sites) - site/exp dependent values
*               - for NS-EW azimutal assymentric atmosphere
*               - markov values
*   clk_mar_apr(3,max_sites)    - site/exp dependent values
*               - for the marckov values of the site clocks
 
*   ecc_change(3, max_sites)  - Ecentricity corrections at
*               - each site in global XYZ coordinates(m) .

*   dum_fila    - filler to get real*8 in correct place
 
 
      real*4 atm_mar_apr(3,max_sites), atm_az_mar_apr(2,max_sites),
     .    clk_mar_apr(3,max_sites), ecc_change(3, max_sites)

      integer*4 dum_fila
 
*   axis_offsets(max_sites) - axis offset values at each of
*               - the sites
*   gamma_apr   - apriori value of gamma
*   nut_ang_apr(2,2)- apriori values of the nutation angles
*               - at about the center of the data set (mas)
 
*   ephem_epoch(max_ephem_epoch)  - epochs for the
*               - ephemeris entries (beginning, approx mid,
*               - end)
 
*   ephem_moon(3,2,max_ephem_epoch) - the ephemeris for the
*               - the moon (X,Y,Z, Xdot,Ydot,Zdot) at each
*               - ephem epoch
*   ephem_sun(3,2,max_ephem_epoch) - the ephemeris for the
*               - the sun (X,Y,Z, Xdot,Ydot,Zdot) at each
*               - ephem epoch
 
*   fut1_inf(4)     - The intepolation information for UT1.
*                   - These values are:
*                   - 1 -- Julian date of first point in table
*                   - 2 -- Spacing between data in days
*                   - 3 -- Number of data points
*                   - 4 -- units in terms of seconds of time.
*   fut1_pts(max_pts)     - the UT1 values at the epoch given in
*                   - Ut1_inf.  To obtain seconds multiply the
*                   - values by UT1_inf(4).
*   fwob_inf(3)     - The intepolation information for wobble.
*                   - These values are:
*                   - 1 -- Julian date of first point in table
*                   - 2 -- Spacing between data in days
*                   - 3 -- Number of data points
*   fwob_pts(2,max_pts)   - the wobble values at the epoch given in
*                   - Wob_inf, in mas.
 
*   pmu_epoch   - Epoch for the Polar motion, UT1 values.
 
*   site_pos(3,max_sites)   - apriori site position used in
*               - in the theoretical delay and rate (m)
*   source_pos(2,max_sources)   - apriori source positions
*               - used in the theorectical delay and rate.
*               - (mas for both RA (at equator) and Dec)
*   etd_apr(3) - Apriori values for h,l, and lag angle
*               - (dimensionless, deg)
 
*   ut1_apr(2)  - Values of UT1-AT and rate at the pmu
*               - epoch (mas and mas/s)
 
*   wob_apr(2,2)- Values of x and y position and rates at
*               - the pmu epoch (mas) (X and Y pole change with
*               - first index)

*   eor_ut1_val(4,2) - Apriori values for UT1 diurnal and semi
*               - diurnal corrections (value at midepoch and rate)
*   eor_xy_val(6,2)  - Apriori values for prograde and retrograde
*               - PM corrections (value at midepch and rate)
*   ext_etd_val(12,max_sites,2) - Apriori values for extended ETD
*               - corrections (value at midepoch and rate of change)
 
 
      real*8 axis_offsets(max_sites), gamma_apr, nut_ang_apr(2,2),
     .    ephem_epoch(max_ephem_epoch),
     .    ephem_moon(3,2,max_ephem_epoch),
     .    ephem_sun(3,2,max_ephem_epoch), fut1_inf(4), 
     .    fut1_pts(max_pts), fwob_inf(3), fwob_pts(2,max_pts),
     .    pmu_epoch, site_pos(3,max_sites),
     .    source_pos(2,max_sources), etd_apr(3), ut1_apr(2),
     .    wob_apr(2,2),
     .    eor_ut1_val(4,2), eor_xy_val(6,2), 
     .    ext_etd_val(12,max_sites,2)
 
*   user_atm_dsc    - Description of the user defined atmospheric
*               - delay conributions
*   user_cont_dsc   - Description of the user defined contributions
*               - to the theoretical delays.
*   user_edit_dsc   - Description of the user defined edit flag
 
*   user_pmu_dsc    - Description of user inserted Polar motion/UT1
*                   - tables
 
 
      character*80 user_atm_dsc, user_cont_dsc, user_edit_dsc,
     .    user_pmu_dsc

* MOD TAH 910905: Added saving of the atmospheric mapping function
*     markov process values.

*   atm_map_mar_apr(max_sites) - Saved values of the mapping function
*     markov process statistics. Removed max_sites words from spare

      real*4 atm_map_mar_apr(max_sites)
 
*   spare_apr(128)  - One spare block to used if needed later
*				 910905: Removed max_sites words
*   Reset with 93/09/23 version of KalObs file
 
 
      integer*4 spare_apr(128)
 
*   last_apr_word   - Last word of the aprioris block
*   dummy_apr(127)  - dummy 127 words to stop overwritting
*               - of files during reading.
 
 
      integer*4 last_apr_word, dummy_apr(127)
 
      equivalence ( aprioris, axis_types )
 
*-------------------------------------------------------------------
* Common declaration 
 
 
      common / apr_block / axis_types, tai_utc, atm_mar_apr,
     .    atm_az_mar_apr, clk_mar_apr, ecc_change, dum_fila,
     .    axis_offsets,
     .    gamma_apr, nut_ang_apr, ephem_epoch, ephem_moon, ephem_sun,
     .    fut1_inf, fut1_pts, fwob_inf, fwob_pts, pmu_epoch, site_pos,
     .    source_pos, etd_apr, ut1_apr, wob_apr, 
     .    eor_ut1_val, eor_xy_val, ext_etd_val, user_atm_dsc,
     .    user_cont_dsc, user_edit_dsc, user_pmu_dsc, 
     .    atm_map_mar_apr, spare_apr,
     .    last_apr_word, dummy_apr
 
*
