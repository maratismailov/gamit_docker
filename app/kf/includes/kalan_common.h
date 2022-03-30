 
*     Main common block for the KALAN program.  This contains
*     most of the working information we need to communicate
*     within the KALAN program
*
*                             10:51 PM SUN., 31 May., 1987
 
*   max_stack_entry - Maximum number of command files which
*               - can be transferred to.
*   num_kalan_commands  - Number of kalan commands
*   num_mar_types       - Number of markov types allowed
*   max_est             - Maxiumnum of estimated parameters
*                       - for atm and clk statistics
 
      integer*4 max_stack_entry, num_kalan_commands, num_mar_types,
     .    max_est
 
      parameter ( max_stack_entry = 5 )
 
      parameter ( num_kalan_commands = 13 )
 
      parameter ( num_mar_types   = 46 )
 
      parameter ( max_est = 2*max_sites )
 
*------------------------------------------------------------------
 
*   bnum_est        - Number of parameters being estimated
*                   - for atm and clk statistics from the delay
*                   - rate data.
*   command_unit    - Command file unit number currently being
*                   - used
*   ed_count(8,2)   - Editting count for 8 catelgories and for
*               - for delay and rate (if both used).  The
*               - categories are:
*               - #   res>limit  data used  bit 8 status
*               - 1     yes           no           1
*               - 2     yes          yes           1
*               - 3      no           no           1
*               - 4      no          yes           1
 
*               - 5     yes           no           0
*               - 6     yes          yes           0
*               - 7      no           no           0
*               - 8      no          yes           0
*               - where bit 8 status is SOLVK automatic edit
*               - bit from data_flag in the KalObs File.
 
*   ed_rec(max_obep)    - Record numbers which should be
*               - read from the KalObs file for editing
*   ed_stat(max_obep)   - Value of the ed_count category for
*               - each observation at an epoch
 
*   max_atm_iter    - Maximum number of iterations on atm
*               - statistics estimation
*   nedits      - Number of records of KalObs file which
*               - should be read at each epoch
 
*   nep_obs     - Number of observations (rates) at the current
*               - epoch
*   next_char   - Next character in command buffer to be read
*   num_allan   - Number of values used to compute the allan
*               - variance at each interval
*   num_est(max_est)  - Number of estimates for chisq for each
*               - site during atmospheric delay sigma estimation.
*   pel         - Program control command number (used when going
*               - between segments)
 
*   stack_entry - Current position on the command stack
*   stack_unit(max_stack_entry ) - stack containing the unit numbers
*               - of the command files
 
*   unit_out    - Unit number for the output file
*   unit_log    - Unit number fot the log file
 
      integer*4 bnum_est, command_unit, ed_count(8,2),
     .    ed_rec(max_obep), ed_stat(max_obep), max_atm_iter, nedits,
     .    nep_obs, next_char, num_allan, num_est(max_est), pel,
     .    stack_entry, stack_unit(max_stack_entry ), unit_out,
     .    unit_log
 
*   iallan_epochs   - start index in ema_data for allan_epochs
*   iallan_mar      - start index in ema_data for allan_mar
*   iallan_sig      - start index in ema_data for allan_sig
 
*   ibak_array  - start index in ema_data for the bak_array
*   irate_data  - start index in ema_data for the rate_data
*   irate_sites - start index in ema_data for the rate_sites
*   ia_part     - start index in ema_data for the a_part
*   icov_obs    - start index in ema_data for the cov_obs
*   icov_atm    - start index in ema_data for the cov_atm
*   ikgain      - start index in ema_data for the kgain
*   irate_res   - start index in ema_data for the rate_residual
 
      integer*4 iallan_epochs, iallan_mar, iallan_sig, ibak_array,
     .    irate_data, irate_sites, ia_part, icov_obs, icov_atm, ikgain,
     .    irate_res
 
*   atm_chisq(max_est)    - Chisq values for the estimates
*               - of the atmospheric zenith delay statistics at
*               - each site
*   atm_start   - Start guess for the atmospheric random walk
*               - statistics (ps**2/2) Usual value is about 0.75
*   atm_tol     - Tolerance on match of atmspheric delay observed
*               - scatter to apriori estimate of scatter (all in
*               - Chisq space).  Typical value would 1.1.  Covers
*               - all values from 0.9 to 1.1.
*   atm_mar(max_est)  - The markov statistics for the atmosphere
*               - zenith delay random walk.
 
*   clk_start   - Start variance for the clocks (ps**2/s)
 
*   sum_chisq(max_est)    - Sum of the chisq values for the
*               - atmospheric statistics at each site
 
*   sum_weight(max_est)    - Weights used in computing the
*               - Atmosphere statistics
*   sum_wgh_all - Summ of all of the weights
*   kaldum      - Dummy real*4 to allign real*8
 
      real*4 atm_chisq(max_est), atm_start, atm_tol, atm_mar(max_est),
     .    clk_start, sum_chisq(max_est), sum_weight(max_est),
     .    sum_wgh_all, kaldum
 
*   al_interval - data spacing to be used for curent allan variance
*               - calculation
*   al_start    - Start spacing for the intevals to be used
*               - for computing the allan variances (sec)
*   al_step     - Multiplying factor to be used for generating
*               - allan varaiance intervals.  Must be greater
*               - than one.
*   al_tol      - tolerance for matching nomial interval to
*               - actual data interval (0-1)
*   allan_dt    - Averge time interval actually used for the
*               - calculation of the allan variances
*   allan_err   - Estimated sigma for the allan standard derivation
*               - Assumes markov estimates are independent, therefore
*               - tends to be an upper bound.
 
 
*   allan_sums(3)   - Summation values for the allan varaince
*               - calculations.  The first element contains the
*               - sum for the allan variance itself, the
*               - second, the sum for the varaiance of the estimate,
*               - the third, the average interval actually used
*               - for the calculation.
*   allan_std   - sqrt of allan variance for the current interval
*               - (units of markov element/sec)
 
*   atm_est(max_est)  - Estimates of the zenith delay rate
*               - for each site (computed for each epoch)
*   total_span  - Total span of the data set in seconds (used for
*               - computing allan variances)
 
 
      real*8 al_interval, al_start, al_step, al_tol, allan_dt,
     .    allan_err, allan_sums(3), allan_std, atm_est(max_est),
     .    total_span
 
*   allan_variances - Indicates that the Allan variances should
*               - be computed from the specified list of
*               - markov elements.
*   BakFile_open    - Indicates that the bak file is open
 
*   atmos_stats - Indicates that we are to compute the atmosphere
*               - statistics from the rate residuals
 
 
*   KalObs_open - Indicates that the KalObs file is open
 
*   write_data  - Indicate KalObs file should be written
*   write_header- Indicates KalObs header should be written
 
 
      logical allan_variances, BakFile_open, atmos_stats, KalObs_open,
     .    write_data, write_header
 
*   mar_types( num_mar_types)    - the names of the
*               - markov elements which can be acessed
*   mar_units( num_mar_types)    - Units of each of the markov
*               - elements
 
      character*8 mar_types( num_mar_types),
     .            mar_units( num_mar_types) 
 
*   kalan_commands(num_kalan_commands)  - The kalan commands
*               - for KALAN
 
      character*12 kalan_commands(num_kalan_commands)
 
*   BakFile_name    - Name of the current BakFile
 
*   command_file    - Name of the Command file
*   log_file        - Name of the log file.
*   output_file     - Name of the output file (may be an LU)
 
*   stack_file(max_stack_entry)   - Names of the command files on the
*                   - stack
 
      character*64 BakFile_name, command_file, log_file, output_file,
     .    stack_file(max_stack_entry)
 
*   buffer          - Command line read from command file
 
 
      character*80 buffer
 
*     Start common declartaion
 
*   bnum_est        - Number of atm/clk parmaters
*   command_unit    - Current command file unit number
*   ed_count   - Editting count for 8 catelgories and for
*               - for delay and rate .  The
*               - categories are:
*               - #   res>limit  data used  bit 8 status
*               - 1     yes           no           1
*               - 2     yes          yes           1
*               - 3      no           no           1
*               - 4      no          yes           1
 
*               - 5     yes           no           0
*               - 6     yes          yes           0
*               - 7      no           no           0
*               - 8      no          yes           0
*               - where bit 8 status is SOLVK automatic edit
*               - bit from data_flag in the KalObs File.
 
*   ed_rec    - Record numbers which should be
*               - read from the KalObs file for editing
*   ed_stat   - Value of the ed_count category for
*               - each observation at an epoch
 
*   max_atm_iter    - Maximum number of iterations on atm
*               - statistics estimation
*   nedits      - Number of records of KalObs file which
*               - should be read at each epoch
 
*   nep_obs     - Number of observations at current epoch
*   next_char   - next character number in command buffer
*   num_allan   - Number of values in allan variacalc.
*   num_est     - Number of estimates for chisq for each
*               - site during atmospheric delay sigma estimation.
*   pel         - Program control command number
 
 
*   stack_entry - Current position on the command stack
*   stack_unit   - stack containing the unit numbers
*               - of the command files
 
*   unit_out    - Unit number for the output file
*   unit_log    - Unit number fot the log file
 
*   iallan_epochs   - start index for allan_epochs
*   iallan_mar      - start index for allan_mar
*   iallan_sig      - start index fot allan_sig
 
*   ibak_array  - start index in ema_data for the bak_array
*   irate_data  - start index in ema_data for the rate_data
*   irate_sites - start index in ema_data for the rate_sites
*   ia_part     - start index in ema_data for the a_part
*   icov_obs    - start index in ema_data for the cov_obs
*   icov_atm    - start index in ema_data for the cov_atm
*   ikgain      - start index in ema_data for the kgain
*   irate_res   - start index in ema_data for the rate_residual
 
*   atm_chisq   - Chisq values for the estimates
*               - of the atmospheric zenith delay statistics at
*               - each site
*   atm_start   - Start guess for the atmospheric random walk
*               - statistics (ps**2/2) Usual value is about 0.75
*   atm_tol     - Tolerance on match of atmspheric delay observed
*               - scatter to apriori estimate of scatter (all in
*               - Chisq space).  Typical value would 1.1.  Covers
*               - all values from 0.9 to 1.1.
*   atm_mar     - The markov statistics for the atmosphere
*               - zenith delay random walk.
*   clk_start   - Start value for clocks
*   sum_chisq   - Sum of the chisq values for the
*               - atmospheric statistics at each site
*   sum_weight  - Weights used in computing atm_stats
*   sum_wgh_all - Sum of all the weights in the atm_stats calculation
 
*   al_interval - Current interval to be for allan var calc (sec)
 
*   al_start    - Start spacing for the intevals to be used
*               - for computing the allan variances (sec)
*   al_step     - Multiplying factor to be used for generating
*               - allan varaiance intervals.  Must be greater
*               - than one.
*   al_tol      - tolerance for matching nomial interval to
*               - actual data interval (0-1)
*   allan_dt    - Averge time interval actually used for the
*               - calculation of the allan variances
*   allan_err   - Estimated sigma for the allan standard derivation
*               - Assumes markov estimates are independent, therefore
*               - tends to be an upper bound.
 
*   allan_sums  - Summation values for the allan varaince
*               - calculations.  The first element contains the
*               - sum for the allan variance itself, the
*               - second, the sum for the varaiance of the estimate,
*               - the third, the average interval actually used
*               - for the calculation.
*   allan_std   - sqrt of allan variance for the current interval
*               - (units of markov element/sec)
 
*   atm_est     - Estimates of the zenith delay rate
*               - for each site (computed for each epoch)
*   total_span  - Total span of the data set in seconds (used for
*               - computing allan variances)
 
*   allan_variances - Indicates that the Allan variances should
*               - be computed from the specified list of
*               - markov elements.
*   BakFile_open    - Indicates that the bak file is open
 
*   atmos_stats - Indicates that we are to compute the atmosphere
*               - statistics from the rate residuals
 
 
*   KalObs_open - Indicates that the KalObs file is open
 
*   write_data  - Indicate KalObs file should be written
*   write_header- Indicates KalObs header should be written
 
*   mar_types    - Types of markov elements which can be
*                   - processed.
 
*   kalan_commands  - The kalan commands
*               - for KALAN
 
*   BakFile_name    - Name of the current BakFile
*   Command_file - NAme of the current command file
 
*   log_file        - Name of the log file.
*   output_file     - Name of the output file
 
*   stack_file   - Names of the command files on the
*                   - stack
*   buffer          - Command line
 
      common / Kalan_common / bnum_est, command_unit, ed_count, ed_rec,
     .    ed_stat, max_atm_iter, nedits, nep_obs, next_char, num_allan,
     .    num_est, pel, stack_entry, stack_unit, unit_out, unit_log,
     .    iallan_epochs, iallan_mar, iallan_sig, ibak_array,
     .    irate_data, irate_sites, ia_part, icov_obs, icov_atm, ikgain,
     .    irate_res, atm_chisq, atm_start, atm_tol, atm_mar, clk_start,
     .    sum_chisq, sum_weight, sum_wgh_all, kaldum,
     .    al_interval, al_start,
     .    al_step, al_tol, allan_dt, allan_err, allan_sums, allan_std,
     .    atm_est, total_span, allan_variances, BakFile_open,
     .    atmos_stats, KalObs_open, write_data, write_header,
     .    mar_types, mar_units, 
     .    kalan_commands, BakFile_name, Command_file,
     .    log_file, output_file, stack_file, buffer
 
*-------------------------------------------------------------------------
 
 
