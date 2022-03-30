 
*     Common block for modear progrma
 
*  max_sat       - Maximum number of satellites
*  max_eph_eps   - Maximum number of ephemeris epochs in the sp3 files.
*  max_epoch     - Maximum number of epochs in cfile
 
      integer*4 max_sat, max_eph_eps, max_epoch
 
      parameter ( max_sat = 100 )
      parameter ( max_eph_eps = 1000 )
      parameter ( max_epoch   = 3000 )

* SP3 File information:     
*  sp3_xyz(3,max_sat, max_eph_eps) - Earth fixed coordinates of satellites
*                          (XYZ, m)
*  sp3_time(max_ephs_eps) - MJD's for the ephemeris entries
*  sp3_clk(max_sat, max_eph_eps) - The clock corrections for the satellites

* Nav file information (if need)
*  nav_clk(3,max_sat, max_eph_eps) - Clock poly nominal information from
*                                    the navigation file
*  nav_time(max_eph_eps)           - MJD's for the times given in the nav
*                                    file

* Satellite position information
*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_clk(max_sat)      - Satellite clock errors (secs)
*  svs_ang(3,max_sat)    - Zenith angle, Azimuth and range to sat.

*  site_xyz(3)           - Site XYZ
*  site_loc(3)           - lat, long of site

*  loc_rot(3,3)          - Rotation matrix from XYZ to NEU
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  site_clock(max_epoch) - Estimated site clock for each epoch (second)

*  data_epoch            - Epoch of the data (mjd)

*  omc(max_sat)          - Oobs minus C (m) (single diff)
*  omc_data(max_sat)     - Oobs minus C (m) for data file
*  dry_zen               - Dry zenith delay (m) computed from station
*                        - station height
*  dXYZ_L1(3), dXYZ_L2(3) - Change in position from mark to phase
*                         - center at L1 and L2
*  dXYZ_tide(3)           - Change in position due to tides.


      real*8  sp3_xyz(3,max_sat, max_eph_eps), 
     .    sp3_time(max_eph_eps), sp3_clk(max_sat, max_eph_eps),
     .    nav_clk(4,max_eph_eps, max_sat),
     .    svs_xyz(3,max_sat),
     .    svs_clk(max_sat), svs_ang(3,max_sat),
     .    site_xyz(3), site_loc(3), loc_rot(3,3),
     .    P1o(max_sat), P2o(max_sat),
     .    L1o(max_sat), L2o(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, omc_data(max_sat), site_clock(max_epoch),
     .    dry_zen, dXYZ_L1(3), dXYZ_L2(3), dXYZ_tide(3)
     
 
*  num_sat       - Number of satelites found
*  num_sp3       - Number of epoch entries in the sp3 file
*  num_nav(max_sat)  - Number of entries in the NAV file for each PRN
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
*  num_chan      - NUmber of channels observed at this epoch
*  num_data      - Number of channels in data file
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)

      integer*4 num_sat, num_sp3, num_nav(max_sat)

*  omc_OK(max_sat) - true if omc is OK else set false.
*  svclk_OK        - Set true when we have values for the clocks
*                    either from the SP3 file or from the Nav file

      logical omc_OK(max_sat), svclk_OK
 
*  in_cf         - Name of input cfile
*  out_cf        - Name of output cfile (same as input name, after the
*                  input has been renamed to name ending in .o)
*  sp3_file      - Name of SP3 file
*  nav_file      - NAme of Rinex navigation file


      character*256 in_cf, out_cf, sp3_file, nav_file


*---------------------------------------------------------------

      common / modear_com / sp3_xyz, 
     .    sp3_time, sp3_clk, nav_clk, svs_xyz, svs_clk, svs_ang,
     .    site_xyz, site_loc, loc_rot,
     .    P1o, P2o, L1o, L2o, P1c, P2c, L1c, L2c,
     .    data_epoch, omc_data, site_clock,
     .    dry_zen, dXYZ_L1, dXYZ_L2, dXYZ_tide, 
     .    num_sat, num_sp3, num_nav, 
     .    omc_OK, svclk_OK,
     .    in_cf, out_cf, sp3_file, nav_file
