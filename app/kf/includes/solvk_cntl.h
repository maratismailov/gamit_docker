 
*---------------------------------------------------------------------
*                                                     SOLVK_CONTROL.FTNI
*     This is the main part of the common block for the Kalman filter
*     programs.  This file will contain information about how the
*     current filter solution is to run.  All of the programs in the
*     solution suite will use this common.
*     The common is saved in a files ('comfile') during the solution,
*     but this file should be treated as a scratch file and other
*     software should not access this file.
*
*     A second part of the common 'MARKOV' block contains all of the
*     information about the Markov and apriori statistics.
*
*     T.Herring                   10:59 PM SAT., 17 Jan., 1987

* MOD TAH 930219: Added cov_mar_neu to allow for markov NEU site
*     positions.  Replaces the XYZ markov process.
*--------------------------------------------------------------------
 
 
*   isoldc(16)  - the 16 word DCB used to read and write this
*               - common as a type 1 file
*   control     - first word in the main common block.  Used to read
*               - and write the file.
 
      integer*4 isoldc(16), control
 
*     This set of declarations must go at the top of the file for decoding
*     purposes
 
*   start_control  - First block number of the control part of SOLVK
*               - common
*   start_markov- First block number of the markov part of the common
 
*   num_control_blocks - Number of blocks in the control part of the common
*   num_markov_blocks   - Number of blocks in the markov part of the
*                   - common
 
      integer*4 start_control, start_markov, num_control_blocks,
     .    num_markov_blocks
 
*   control_read   - Indicates that the control part of file has been read
*   markov_read - Indicates that the markov part of the file has been
*               - read.
 
      logical control_read, markov_read
 
*     End top variable declarations
*-----------------------------------------------------------------------
 
*   clk_parn(2,max_sites)   - the parameter numbers of the
*               - clock parameters (offset and rate) for
*               - each site.  If parameter is not being
*               - estimated, the entry in clk_parn is set
*               - to zero.
*   icrt        - User's session LU number.
*   iprnt       - The print device LU number for this solution.
*   ilog        - The log device LU number for this solution.  A
*               - single line per solution is normally output to
*               - this device.
*   ind_mar(max_parn)   - the parameter numbers of the markov
*               - parameters in this solution.
*   nblocks     - size of the covariance matrix and solution
*               - vector in 128 word blocks (used when writing
*               - and reading data files)
*   num_av      - Number of prefit residuals to be average while
*               - looking for clock breaks.
*   num_krec    - Actual number of records written to KALFILE.
*               - This value in increment for each data type
*               - written out, and is not incremented if data
*               - is out of time order.
*   num_mar     - number of markov parameters in the current
*               - solution.
*   num_tran    - number of state transition parameters in the
*               - current solution
 
*   par_num     - Total number of parameters to be estimated in
*               - current solution.  Computed in PARAM_LIST
 
*   pre_num     - Number of prefit residuals in prefit chi**2
*   solution_number -   Number of the current solution.  Used when
*               - multiple solutions are run on the one data set.
*   state_tran(max_tran)    - array which points to the parameter
*               - number of the rate component of the state
*               - transition parameters.  There is an implicit
*               - assumption in the software that the parameter
*               - imediately prceeding the rate term is the offset
*               - term.  A warning message is printed if this is
*               - not the case.
*   total_epochs    - Total number of epochs of data to be used in
*               - this solution.
*   used_epochs - Number of epochs of data which contain at least
*               - some weighted data.
 
      integer*4 clk_parn(2,max_sites), icrt, iprnt, ilog,
     .    ind_mar(max_parn), nblocks, num_av, num_krec, num_mar,
     .    num_tran, par_num, pre_num, solution_number,
     .    state_tran(max_tran), total_epochs, used_epochs
 
*   atm_mar(3,max_sites)    - random, random walk and integrated
*               - random walk statistics for the atmospheric
*               - zenith delay parameters parameters.
*               - Units: ps**2, ps**2/sec, ps**2/sec**2
*   atm_apr(2,max_sites)    - Apriori variances for the atmosphere
*               - zenith delay and rate.
*               - Units: ps**2 and (fs/sec)**2
*   clk_mar(3,max_sites)    - random, random walk and integrated
*               - random walk statistics for the clock markov
*               - parameters.
*               - Units: ps**2, ps**2/sec, ps**2/sec**2
*   clk_apr(2,max_sites)    - apriori variances for the clock
*               - offset and rate parameters.
*               - Units: ps**2 and (fs/sec)**2
*   cov_mar(max_parn)   - the step variances for each of the
*               - Markov elements.
*   cov_apr(max_parn)   - the apriori variances for each of the
*               - parameters being estimated.
 
      real*4 atm_mar(3,max_sites), atm_apr(2,max_sites),
     .    clk_mar(3,max_sites), clk_apr(2,max_sites),
     .    cov_mar(max_parn), cov_apr(max_parn)

* MOD TAH 930219: Added markov step covariance matrix for NEU
*   cov_mar_neu(3,3,max_sites) -Step covariance matrix for a NEU
*         stochastic process.

      real*4 cov_mar_neu(3,3,max_sites)
 
*   clk_brk_criteria    - Jump in the prefit residuals averaged
*               - over NUM_AV observations which will be treated
*               - as a clock break (ps).
      real*4 dum_allign
 
      real*8 clk_brk_criteria
 
*   bck_soln    - Indicate back solution is to be run.  Set by
*               - giving BCK_FILE name (unless name of file NONE)
 
      logical bck_soln
 
*   apr_file    - Name of the file containing the aproiri site and
*               - source positons to be used in the current
*               - solution
*   bck_file    - Name of the file for the back solution results
*               - to be put in.  (Computing the back solution will
*               - slow the Kalman filter soltions considerably)
*   com_file    - Name of the Kalman filter control common file.
*               - May be passed through the runstring or defaults
*               - to value given in Kalman_blockdata.
*   data _file  - Name of the file containing the VLBI data.  This
*               - file is greated using READIN
*   ecc_file    - Name of the file containing the eccentricity
*               - data for the VLBI data set.
*   glb_file    - Name of file in which the covariance matrix
*               - solution will be stored for later processing
*               - in GLOBK.
*   kal_file    - Name of the file to be passed to the solution
*               - program.  This file contains only pre-fit
*               - residuals and partial derivatives.
*   mar_file    - Name of the SOLVK control markov file (in CI
*               - this will be the root)
*   sol_file    - Name of the file used to store the covariance
*               - matrices for the solution.
*   sd_file     - File containing the semidiurnal and diurnal UT1,
*               - pole position and tide corrections.
*   log_file    - Name of log file for solution
*   out_file    - Name of output file (may be 6)
 
      character*64 apr_file, bck_file, com_file, data _file, ecc_file,
     .             glb_file, kal_file, mar_file, sol_file, 
     .             sd_file, log_file, out_file
 
*------------------------------------------------------------------------
*     Ending declarations
 
*   last_control_word  - Last word in the control part of the common block
*   dummy_control(127) - 127 words to make sure that we donot overwrite
*                   - anything when we read the block
 
      integer*4 last_control_word, dummy_control(127)
 
      equivalence ( control, start_control )
 
*-----------------------------------------------------------------------
*     Start of common declaration
 
 
      common / solvk_control / isoldc, start_control, start_markov,
     .    num_control_blocks, num_markov_blocks, control_read,
     .    markov_read, clk_parn, icrt, iprnt, ilog, ind_mar, nblocks,
     .    num_av, num_krec, num_mar, num_tran, par_num, pre_num,
     .    solution_number, state_tran, total_epochs, used_epochs,
     .    atm_mar, atm_apr, clk_mar, clk_apr, cov_mar, cov_apr,
     .    cov_mar_neu, 
     .    dum_allign,
     .    clk_brk_criteria, bck_soln, apr_file, bck_file, com_file,
     .    data _file, ecc_file, glb_file, kal_file, mar_file, sol_file,
     .    sd_file, log_file, out_file, 
     .    last_control_word, dummy_control
 
