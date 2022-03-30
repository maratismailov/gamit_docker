*
*     Include file for MONSIM program.
*
*   com_len         - Length of the commands in this program
*   max_sv          - Maximum number of SV
*   max_par         - Maximum number of parameters
*   num_mon_comms   - Number of commands
*   max_files       - Maximum number of input files
*   max_clk_poly    - Dimension for the latest clock polynomial
*                   - coefficients
*   max_epochs      - Maxiumum number of epochs allowed
 
      integer*4 com_len, max_sv, max_par, num_mon_comms, max_files,
     .          max_clk_poly, max_epochs
 
      parameter ( com_len      = 16 )
      parameter ( max_sv       = 24 )
      parameter ( max_par      =  5 )
      parameter ( num_mon_comms= 24 )
      parameter ( max_files    = 100 )
      parameter ( max_clk_poly = 10 )
      parameter ( max_epochs   = 5000 )
 
 
*   ctime           - current time (sec)
*   next_time		- Time of next observation
*   prev_time       - Epoch of previous observation
*   file_times(max_files) - Current time for each of the data files
*                       - being read
*   cycles_to_mm        - Converstion from cycles to mm (one values at the
*                       - moment but maybe one per file)
*   epoch_to_sec        - Conversion from epcoh number to seconds
*   az_up(max_sv)       - Azmith of visible SV (rad)
*   el_up(max_sv)       - Eleve of visible SV (rad)
*   min_elevation   - Minimum elevation to use (rad)
*   lowest_elev     - Lowest elevation used in batch (rad)
 
*   site_loc(3)     - Site coordinates as lat,long, rad
 
*   mar_atm         - Atmosphere RW (mm^2/s)
*   apr_atm         - Apriori variance on atm (mm^2)
*   mar_clk         - CLock RW (mm^2/s)
*   apr_clk         - Apriori variance on clk (mm^2)
 
*   apr_site(3)     - Apriori sigma on site (mm^2) in
*                   - N,E, and UP coordinates
*   var_obs         - Data variance (mm^2)
 
*   cov_par(max_par,max_par)    - Parmeter covariance
*   sol_par(max_par)        - Parameter estimates
*   cov_obs(max_sv,max_sv)  - data covariance
*   prefit(max_sv)          - Prefit data residuals
*   k_gain(max_sv, max_par) - Kalman gain
*   temp_gain(max_sv, max_par)  - Temporary calc
*   a_part(max_par,max_sv)  - Partial derivatives
 
*   sig_par(max_par)        - Sigmas on parameter estimates
*   condition(max_par)  - Condition numbers

*   bias(max_files)     - Constant value which must be added to all
*                    - data from a given file
*   add_ambig(max_files) - Additive ambiguity as a function of file number
*   clock_poly(0:max_clk_poly-1)  - Estimates of the site clock polynomial
*   lambda   - Wavelength factor (2 for doubling recievers)
*   samping_interval - sampling interval of data
 
      real*8 ctime, next_time, prev_time, file_times(max_files),
     .    cycles_to_mm,  epoch_to_sec, az_up(max_sv),
     .    el_up(max_sv), min_elevation, lowest_elev, site_loc(3),
     .    mar_atm, apr_atm, mar_clk, apr_clk, apr_site(3),
     .    var_obs, cov_par(max_par,max_par), sol_par(max_par), 
     .    cov_obs(max_sv,max_sv), prefit(max_sv),
     .    k_gain(max_sv, max_par), temp_gain(max_sv, max_par),
     .    a_part(max_par,max_sv), sig_par(max_par), condition(max_par),
     .    bias(max_files), clock_poly(0:max_clk_poly-1), lambda,
     .    add_ambig(max_files), 
     .    sampling_interval, epoch_range(2), atm_apr , bias_minel
 
*   num_par     - Number of parameters estimated
*   num_sv      - NUmber of satellites
*   num_up      - number of visible SV
*   num_in_batch    - Number of data in each batch
*   clk_order       - Order of the clock polynomial to fit to the
*                     data before processing.
*   svnum_up(max_sv) - Numbers of the satellites currently up
*   num_files   - Number of data files to be read
*   unitd(max_files) - Unit numbers for open files (if value negative then file
*               - is not open
 
      integer*4 num_par, num_sv, num_up, num_in_batch, clk_order,
     .    svnum_up(max_sv), num_files, unitd(max_files), col_data

*   azel_file_open  - Indicates that azel_file is open

      logical azel_file_open, no_bias
 
*   mon_commands(num_mon_comms) - Commands in this
*                   - program.
*   azel_file  - Name of the Azimith elevation file
*   data_files(max_files)	- Names of the input files
*   buffers(max_files)		- Lines read from current files
 
      character*(com_len) mon_commands(num_mon_comms)

      character*80  buffers(max_files)

      character*128 azel_file, data_files(max_files)
 
*----------------------------------------------------------------
* COMMON DECLARATION
*
 
      common / monsim_com / ctime, next_time, prev_time, file_times, 
     .    cycles_to_mm, epoch_to_sec, az_up, el_up,
     .    min_elevation, lowest_elev, site_loc, mar_atm,
     .    apr_atm, mar_clk, apr_clk, apr_site, var_obs, cov_par,
     .    sol_par, cov_obs, prefit, k_gain, temp_gain, 
     .    a_part, sig_par, condition, bias, clock_poly, lambda,
     .    add_ambig, sampling_interval, epoch_range, atm_apr,
     .    bias_minel,
     .    num_par, num_sv, num_up, num_in_batch, clk_order,
     .    svnum_up, num_files, unitd, col_data, 
     .    azel_file_open, no_bias, mon_commands, 
     .    buffers, azel_file, data_files
 
*---------------------------------------------------------------------
 
