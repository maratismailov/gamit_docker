 
*     Common block for svsp3 program

* MOD TAH 180312: Removed SP3 file information and put in includes/sp3_def.h

      include '../includes/sp3_def.h' 

* MOD TAH 151021:Updated code to handle large number of data types
      integer*4 max_data_types
      parameter ( max_data_types = 36 ) ! FIXED: If > 18 code needs mod

*  sp3_xyz(3,max_sat, max_eph_eps) - Earth fixed coordinates of satellites
*                          (XYZ, m)
*  sp3_time(maxmax_ephs_eps) - JD's for the ephemeris entries
*  svs_clk(max_sat,max_eph_eps)       - Clock error by time (us)
*  svs_dt(max_sat)       - Clock error (us) 
*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_loc(3,max_sat)    - Phi, Lambda, and height of sattellite
*  svs_ang(3,max_sat)    - Zenith angle, Azimuth and range to sat.
*  start, stop, step     - Start, stop and step for position
*  site_xyz(3)           - Site XYZ
*  site_geod(3)          - geodetic co-lat, long of site (rad), and heighy
*  diff_xyz(3)           - difference site XYZ coordinates
*  diff_geod(3)          - geodetic co-lat, long of ref site (rad), and heighy
*  loc_rot(3,3)          - Rotation matrix from XYZ to NEU
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  data_epoch            - Epoch of the data (jd)
*  diff_epoch            - Epoch of difference data (jd)
*  omc(max_sat)          - Oobs minus C (m) (single diff)
*  omc_data(max_sat)     - Oobs minus C (m) for data file
*  omc_diff(max_sat)     - Oberved minus C for difference data set
*  apart(6,max_sat)      - Partial derivatives
*  apart_data(6,max_sat) - direct data omc 
*  kgain(6,max_sat),temp_gain(6,max_sat) - Kalman Gain and
*                         a temp intermeridate comp.
*  cov_parm(6,6), sol_vec(6) - Covariance matrix and solution
*                         vector for the parameters
*  wn(6)                  - Process noise model per step
*  acat_mat(max_sat,max_sat) - Covariance matrix of data plus
*                          predicted parameter noise
*  dx(6)                 - Adjustmenst to the parameters
*  out_spacing  - Spacing for output values (secs)
 
      real*8 svs_xyz(3,max_sat), svs_dt(max_sat), 
     .    svs_loc(3,max_sat), svs_ang(3,max_sat),
     .    start, stop, step, weekno,
     .    site_xyz(3), diff_xyz(3), site_geod(3),
     .    diff_geod(3), loc_rot(3,3), lat_deg, long_deg,
     .    P1o_data(max_sat), P2o_data(max_sat),
     .    L1o_data(max_sat), L2o_data(max_sat),
     .    P1o_diff(max_sat), P2o_diff(max_sat), 
     .    L1o_diff(max_sat), L2o_diff(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, diff_epoch, omc_data(max_sat),omc_diff(max_sat),
     .    apart_data(6,max_sat),
     .    omc(max_sat), apart(6,max_sat), kgain(6,max_sat),
     .    temp_gain(6,max_sat), cov_parm(6,6), sol_vec(6),
     .    wn(6), acat_mat(max_sat,max_sat), dx(6), data_noise,
     .    chi, out_spacing

      real*8 atm_apr_sig(2)  ! Aprori sigmas for the atmospheric
                ! delays at the two sites.  Set with -atm option.


      real*8 svres_sum(max_sat)    ! Mean residual by satellite (also
                  ! used as summation array)
     .,      svres_var(max_sat)   ! RMS scatter of satellite residuals
                  ! (also used for summation).

      real*4 data_rxver  ! Version of data rinex file
      real*4 diff_rxver  ! Version of diff rinex file


      integer*4 svres_num(max_sat)  ! Number of values for each aatellite
 
*  num_chan      - NUmber of channels observed at this epoch
*  num_data      - Number of channels in data file
*  num_diff      - Number of channels in difference file
*  num_used      - Actual number of channels used.
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)
*  chan_data(max_sat) = channels numbers form data file
*  chan_diff(max_sat) = channels in difference file (PRN numbers)
*  num_data_types - NUmber of data types 
*  num_diff_types - number of data types in difference file
 
      integer*4 num_chan,  chan(max_sat),
     .          nchi, num_used, num_data_types, num_diff, num_data,
     .          chan_data(max_sat), chan_diff(max_sat), num_diff_types

*  debug_start    = Start epoch for debug
*  debug_end      - End epoch for debug
*  proc_start, proc_end - Processing start and end

      integer*4 debug_start, debug_end, proc_start, proc_end

*  omc_OK(max_sat) - true if omc is OK else set false.
*  omc_OK_data(max_sat), omc_OK_diff(max_sat) - Omc OK for data and diff
      logical omc_OK(max_sat), omc_OK_data(max_Sat), 
     .        omc_OK_diff(max_sat)
* Nodiff -- Set true if no difference data file is given and
*     point positioning will be used.
      logical Nodiff
 
      logical rep_dels   ! Set true to report deletss (set with 
                         ! -rep_dels options),

*  nav_file      - NAme of Rinex navigation file
*  data_file     - Name of Rinex data file
*  diff_file     - Names of rinex data file to be differenced from
*                  first one

      character*256 nav_file, data_file, diff_file

* site_name      - Name of the site
* diff_name      - Names of site in difference names
* data_types(max_data_types)  - Type of data in the channels
* diff_types(max_data_types)  - TYpe of data in difference file
* out_type       - Sets output coordinate type (XYZ or NEU)

      character*8 site_name, diff_name
      character*2 data_types(max_data_types),diff_types(max_data_types)
      character*4 out_type
      character*8 anal_typ    ! String with analysis type (LC,L1,L2)

      character*10 gnsstouse  ! String with codes for GNSS to use
                            ! GRECJSI (set with file name)

*---------------------------------------------------------------

      common / svtrack_com /  svs_xyz,
     .    svs_loc, svs_ang, start, stop, step,  weekno,
     .    site_xyz,  diff_xyz, svs_dt, 
     .    site_geod, diff_geod,   loc_rot, lat_deg, long_deg,
     .    P1o_data, P2o_data, L1o_data, L2o_data,
     .    P1o_diff, P2o_diff, L1o_diff, L2o_diff, P1c, P2c, L1c, L2c,
     .    data_epoch,  omc, apart, kgain, temp_gain, cov_parm, sol_vec,
     .    diff_epoch, omc_diff, omc_data, apart_data,
     .    wn, acat_mat, dx, data_noise, chi, out_spacing, atm_apr_sig,
     .    svres_sum, svres_var,
     .    data_rxver, diff_rxver,  
     .    num_chan,  chan, chan_data, 
     .    nchi, num_used, num_data, debug_start, debug_end,
     .    proc_start, proc_end, svres_num, 
     .    num_data_types, num_diff, chan_diff, num_diff_types, omc_ok,
     .    omc_OK_data, omc_OK_diff, Nodiff, rep_dels, 
     .    nav_file, data_file, site_name, data_types, out_type,
     .    diff_file, diff_types, diff_name, anal_typ, gnsstouse
 
