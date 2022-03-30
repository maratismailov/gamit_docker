CTITLE GEN_GOBS
 
      subroutine gen_gobs(L1_cyc_cse, L2_cyc_cse, data_flag_cse,
     .    L1r_rng_cse, L2r_rng_cse,  params_cse, par_flag_cse  )

      implicit none
 
*     Routine to generate the Gobs file for processing in SOLVG.
*     Here all of the cfiles are re-opened and read as the gobs
*     file is written.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for
*                 each observation
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep)
 
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement  (may be
*                     fractional
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*                   - by epoch.
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)    - L1 rnage residuals
*                 passed becuase millisecond jumps are needed.
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)    - L1 rnage residuals
*                 passed becuase millisecond jumps are needed.
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT errors and Gobs file errors
*   trimlen     - Length of string routine.
 
 
c      integer*4 ierr, trimlen
 
****  Dummy stub:

      return
      end


