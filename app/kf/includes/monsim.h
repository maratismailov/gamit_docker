*
*     Include file for MONSIM program.
*
*   com_len         - Length of the commands in this program
*   max_sv          - Maximum number of SV
*   max_par         - Maximum number of parameters
*   num_mon_comms   - Number of commands
 
      integer*4 com_len, max_sv, max_par, num_mon_comms
 
      parameter ( com_len      = 16 )
      parameter ( max_sv       = 24 )
      parameter ( max_par      =  5 )
      parameter ( num_mon_comms= 12 )
 
 
*   start_time      - Time to start (seconds past the epoch of the
*                   - orbital elements)
*   end_time            - Time of last obs (sec)
*   sampling_interval   - Sampling interval in seconds (should not
*                   - be less than about 60 seconds)
*   obs_interval        - Spacing between obs in batch (sec)
*   tref            - reference time (start of current batch of
*                   - data)
*   ctime           - current time (sec)
*   prev_time       - Epoch of previous observation
*   sv_inclination  - Inclination of satellite orbits (rad)
*   sv_a            - Semimajor axis (m)
*   mean_motion     - Mean motion of SV (rads/sec)
*   arg_perigee(max_sv) - argument of perigee for each SV
*   long_node(max_sv)       - Long of node for each SV
*   az_up(max_sv)       - Azmith of visible SV (rad)
*   el_up(max_sv)       - Eleve of visible SV (rad)
*   min_elevation   - Minimum elevation to use (rad)
*   lowest_elev     - Lowest elevation used in batch (rad)
 
*   site_loc(3)     - Site coordinates as lat,long, rad
*   site_xyz(3)     - XYZ coordinates of site
 
*   mar_atm         - Atmosphere RW (mm^2/s)
*   apr_atm         - Apriori variance on atm (mm^2)
*   mar_clk         - CLock RW (mm^2/s)
*   apr_clk         - Apriori variance on clk (mm^2)
 
*   apr_site(3)     - Apriori sigma on site (mm^2) in
*                   - N,E, and UP coordinates
*   var_obs         - Data variance (mm^2)
 
*   cov_par(max_par,max_par)    - Parmeter covariance
*   cov_obs(max_sv,max_sv)  - data covariance
*   k_gain(max_sv, max_par) - Kalman gain
*   temp_gain(max_sv, max_par)  - Temporary calc
*   a_part(max_par,max_sv)  - Partial derivatives
 
*   sig_par(max_par)        - Sigmas on parameter estimates
*   condition(max_par)  - Condition numbers
 
      real*8 start_time, end_time, sampling_interval, obs_interval,
     .    tref, ctime, prev_time, sv_inclination, sv_a, mean_motion,
     .    arg_perigee(max_sv), long_node(max_sv), az_up(max_sv),
     .    el_up(max_sv), min_elevation, lowest_elev, site_loc(3),
     .    site_xyz(3), mar_atm, apr_atm, mar_clk, apr_clk, apr_site(3),
     .    var_obs, cov_par(max_par,max_par), cov_obs(max_sv,max_sv),
     .    k_gain(max_sv, max_par), temp_gain(max_sv, max_par),
     .    a_part(max_par,max_sv), sig_par(max_par), condition(max_par)
 
*   num_par     - Number of parameters estimated
*   num_sv      - NUmber of satellites
*   num_up      - number of visible SV
*   num_in_batch    - Number of data in each batch
*   num_orb_planes  - Number of orbital planes
*   svnum_up(max_sv) - Numbers of the satellites currently up
 
      integer*4 num_par, num_sv, num_up, num_in_batch, num_orb_planes,
     .    svnum_up(max_sv)

*   azel_file_open  - Indicates that azel_file is open

      logical azel_file_open
 
*   mon_commands(num_mon_comms) - Commands in this
*                   - program.
*   azel_file  - Name of the Azimith elevation file
 
 
      character*(com_len) mon_commands(num_mon_comms)

      character*128 azel_file
 
*----------------------------------------------------------------
* COMMON DECLARATION
*
 
      common / monsim_com / start_time, end_time, sampling_interval,
     .    obs_interval, tref, ctime, prev_time, sv_inclination, sv_a,
     .    mean_motion, arg_perigee, long_node, az_up, el_up,
     .    min_elevation, lowest_elev, site_loc, site_xyz, mar_atm,
     .    apr_atm, mar_clk, apr_clk, apr_site, var_obs, cov_par,
     .    cov_obs, k_gain, temp_gain, a_part, sig_par, condition,
     .    num_par, num_sv, num_up, num_in_batch, num_orb_planes,
     .    svnum_up, azel_file_open, mon_commands, azel_file
 
*---------------------------------------------------------------------
 
