 
*     Common block for svinert progrma
 
*  max_sat       - Maximum number of satellites
 
      integer*4 max_sat
 
      parameter ( max_sat = 100 )
 
*  toe_jd(max_sat)  - Ephemeris times for each satelite (as jd)
*   toe(max_sat)        - Ephemeris time in seconds from start
*                       - GPS week.
*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_loc(3,max_sat)    - Phi, Lambda, and height of sattellite
*  svs_ang(3,max_sat)    - Zenith angle, Azimuth and range to sat.
*  start, stop, step     - Start, stop and step for position
*  aode(max_sat)         - Age of ephemeris
*  af0(max_sat), af1(max_sat), af2(max_sat)     - Clock corrections
*  crs(max_sat), crc(max_sat)    - Radial corrections
*  dn(max_sat)   - Change to mean motion
*  m0(max_sat)   -Mean anomaly
*  cuc(max_sat), cus(max_sat)    - Correction for latitude
*  ecc(max_sat)          - Eccentricity
*  art(max_sat)          - root semimajor axis
*  cic(max_sat), cis(max_sat)    - Inclination corrections
*  om0(max_sat)          - longitude of node
*  i0(max_sat)           - Inclination
*  w(max_sat)            - argument of perigee
*  omd(max_sat)          - Omega dot
*  idt(max_sat)          - Inclination rate
*  cflg12(max_sat), pflg12(max_sat)      - Flags
*  weekno                - Week number
*  svacc, svhealth, tgd - Miscellaneous
*  aodc(max_sat)         - Age of clokc
*  site_xyz(3)           - Site XYZ (Initial value)
*  site_apr(3)           - Apriori position to be used as we move
*                          along the trajectory.
*  site_loc(3)           - lat, long of site
*  diff_xyx(3)           - difference site XYZ coordinates
*  loc_rot(3,3)          - Rotation matrix from XYZ to NEU
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  data_epoch            - Epoch of the data (jd)
*  diff_epoch            - Epoch of difference data (jd)
*  omc(max_sat)          - Oobs minus C (m) (single diff)
*  omc_data(max_sat)     - Oobs minus C (m) for data file
*  omc_diff(max_sat)     - Oberved minus C for difference data set
*  apart(4,max_sat)      - Partial derivatives
*  apart_data(4,max_sat) - direct data omc 
*  kgain(22,max_sat),temp_gain(22,max_sat) - Kalman Gain and
*                         a temp intermeridate comp.
*  Extended Kalman Gain for inertial input.
*  The state vector is now:
*  X    Y    Z
*  Site clock.
*  Xd   Yd   Zd
*  Xdd  Ydd  Zdd
*  Xddb Yddb Zddb
*  T1   T2   T3
*  T1d  T2d  T3d
*  T1b  T2b  T3b
* where d is time deriative, and T is orientation (to be applied
* accelerations) and b denotes bias in rate and accelerometers.

*  cov_parm(22,22), sol_vec(22) - Covariance matrix and solution
*                         vector for the parameters
*  wn(22)                  - Process noise model per step
*  acat_mat(max_sat,max_sat) - Covariance matrix of data plus
*                          predicted parameter noise
*  dx(22)                 - Adjustmenst to the parameters
*  out_spacing  - Spacing for output values (secs)
 
      real*8 toe_jd(max_sat), toe(max_sat), svs_xyz(3,max_sat),
     .    svs_loc(3,max_sat), svs_ang(3,max_sat),
     .    start, stop, step, aode(max_sat),
     .    af0(max_sat), af1(max_sat), af2(max_sat),
     .    crs(max_sat), crc(max_sat), dn(max_sat), m0(max_sat),
     .    cuc(max_sat), cus(max_sat), ecc(max_sat), art(max_sat),
     .    cic(max_sat), cis(max_sat), om0(max_sat), i0(max_sat),
     .    w(max_sat), omd(max_sat), idt(max_sat),
     .    cflg12(max_sat), pflg12(max_sat), weekno,
     .    svacc, svhealth, tgd, aodc(max_sat), site_xyz(3),
     .    site_apr(3), site_vel(3), site_acc(3), diff_xyz(3), 
     .    site_loc(3), loc_rot(3,3), lat_deg, long_deg,
     .    P1o_data(max_sat), P2o_data(max_sat),
     .    L1o_data(max_sat), L2o_data(max_sat),
     .    P1o_diff(max_sat), P2o_diff(max_sat), 
     .    L1o_diff(max_sat), L2o_diff(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, diff_epoch, omc_data(max_sat),omc_diff(max_sat),
     .    apart_data(4,max_sat),
     .    omc(max_sat), apart(4,max_sat), kgain(22,max_sat),
     .    temp_gain(22,max_sat), cov_parm(22,22), sol_vec(22),
     .    wn(22), acat_mat(max_sat,max_sat), dx(22), data_noise,
     .    chi, out_spacing

*  Additional variables needed
*  Td_data(3)    - Theta rate data
*  Td_sig(3)     - Theta rate sigma
*  Xdd_data(3)   - X acceleration data 
*  Xdd_sig(3)    - Sigma for acceleration data

*  X_apr_sig(3) - Apriori sigma for position
*  Xd_apr_sig(3) - Apriori sigma for velocity
*  Xdd_apr_sig(3)  - Apriori sigma for accelration
*  Xddb_apr_sig(3)  - Apriori sigma for accelration bias
*  T_apr_sig (3)    - Apriori sigma for orientation
*  Td_apr_sig(3)    - Apriori sigma for orienation rate
*  Tdb_apr_sig(3)   - APriori sigma for orientation rate bias.

*  ep_inert         - Epoch for the inertial measurement.

*  Xdd_obs(3)       - Three observed accelarations
*  Td_obs(3)        - Thres observed rotation rates
 
      real*8 Td_data(3), Td_sig(3) , Xdd_data(3), Xdd_sig(3) ,
     .       X_apr_sig(3), Xd_apr_sig(3) , Xdd_apr_sig(3) ,
     .       Xddb_apr_sig(3) , T_apr_sig (3), Td_apr_sig(3),
     .       Tdb_apr_sig(3), ep_inert,
     .       Xdd_obs(3), Td_obs(3) 

*  num_sat       - Number of satelites found
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
*  num_chan      - NUmber of channels observed at this epoch
*  num_data      - Number of channels in data file
*  num_diff      - Number of channels in difference file
*  num_used      - Actual number of channels used.
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)
*  chan_data(max_sat) = channels numbers form data file
*  chan_diff(max_sat) = channels in difference file (PRN numbers)
*  num_data_types - NUmber of data types 
*  num_diff_types - number of data types in difference file
 
      integer*4 num_sat, prn(max_sat), num_chan,  chan(max_sat),
     .          nchi, num_used, num_data_types, num_diff, num_data,
     .          chan_data(max_sat), chan_diff(max_sat), num_diff_types

*  omc_OK(max_sat) - true if omc is OK else set false.
*  omc_OK_data(max_sat), omc_OK_diff(max_sat) - Omc OK for data and diff
*  Inert_OK   - Set true if inert data file OK.
      logical omc_OK(max_sat), omc_OK_data(max_Sat), 
     .        omc_OK_diff(max_sat), Inert_OK
 
*  nav_file      - NAme of Rinex navigation file
*  data_file     - Name of Rinex data file
*  diff_file     - Names of rinex data file to be differenced from
*                  first one
*  command_file  - Name of command file
*  inert_file    - Name of inertial data file

      character*256 nav_file, data_file, diff_file, command_file ,
     .              inert_file 

* site_name      - Name of the site
* diff_name      - Names of site in difference names
* data_types(10)  - Type of data in the channels
* diff_types(10)  - TYpe of data in difference file
* out_type       - Sets output coordinate type (XYZ or NEU)

      character*8 site_name, data_types(10), diff_name, diff_types(10)
      character*4 out_type

*---------------------------------------------------------------

      common / svinert_com /  toe_jd, toe, svs_xyz,
     .    svs_loc, svs_ang, start, stop, step, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12, weekno,
     .    svacc, svhealth, tgd, aodc, site_xyz, site_apr, 
     .    site_vel, site_acc, diff_xyz,
     .    site_loc,  loc_rot, lat_deg, long_deg,
     .    P1o_data, P2o_data, L1o_data, L2o_data,
     .    P1o_diff, P2o_diff, L1o_diff, L2o_diff, P1c, P2c, L1c, L2c,
     .    data_epoch,  omc, apart, kgain, temp_gain, cov_parm, sol_vec,
     .    diff_epoch, omc_diff, omc_data, apart_data,
     .    wn, acat_mat, dx, data_noise, chi, out_spacing,
     .    Td_data, Td_sig , Xdd_data, Xdd_sig ,
     .    X_apr_sig, Xd_apr_sig , Xdd_apr_sig ,
     .    Xddb_apr_sig , T_apr_sig , Td_apr_sig,
     .    Tdb_apr_sig, ep_inert,
     .    Xdd_obs, Td_obs,
     .    num_sat, prn,  num_chan,  chan, chan_data, nchi, num_used,
     .    num_data,
     .    num_data_types, num_diff, chan_diff, num_diff_types, omc_ok,
     .    omc_OK_data, omc_OK_diff, inert_OK,
     .    nav_file, data_file, site_name, data_types, out_type,
     .    diff_file, diff_types, diff_name, command_file, inert_file
 
