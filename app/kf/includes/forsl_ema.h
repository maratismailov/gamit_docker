c
c The ema declaration for program FORSL
c
c .....................................................................
c Mod history                12:27 PM  FRI., 29  MAR., 1985
c WHO  WHEN   WHAT
c .....................................................................
c TAH 850329 : Typed in declaration
c
c .....................................................................
c
      integer*4 for_ema(max_ema_space)
 
      common / forsl_ema /  for_ema
 
c
c Variables
c ---------
c The array 'for_ema' is an integer array which will fill most
c     of the ema space available to the program.  To allow
c     efficent use of the ema space, the variables which are
c     used in ema are 'equivalenced' to this array in software,
c     i.e., the locations of the variables depend on the size
c     of the arrays at run time.  This procedure should ensure
c     the minimum amount of swapping between main memory and EMA.
c
c     The following equivalence between 'for_ema' and the true
c     arrays is used.
c
c     Index in              Effective variable
c     for_ema
c
c     icov_parm            = cov_parm(1) -> cov_parm(par_num*par_num)
c
c     isol_parm            = sol_parm(1) -> sol_parm(par_num)
c
c     irec_parm            = rec_count(1) -> rec_count(nbl+2)
c
c     ia_part              = a_part(1)   -> a_part(max_deriv*nbl)
c
c     itemp_gain           = temp_gain(1) -> temp_gain(par_num*nbl)
c
c     ik_gain              = k_gain(1)   -> k_gain(par_num*nbl)
c
c     icov_obs             = cov_obs(1)  -> cov_obs(nbl*nbl)
c
c     ipre_fit             = prefit(1)   -> prefit(nbl)
c
c     ikal_array           = kal_array(1) -> kal_array(krec_len*max_krecs)
c
c
c     where
c cov_parm  -- the covariance matrix of the estimated parameters.
c     This matrix is updated for each epoch of observations made
c sol_parm  -- the estimated parameters, also updated for each epoch
c     of observation
c rec_count -- the start and stop records in Kalfile for each epoch
c     plus the record numbers of the observations actually used in
c     the kalman filter estimation.  These values are used in BAKSL.
c a_part    -- the partial derivatives of the delay/rate model with
c     respect to the paramters.  This matrix is in compressed form
c     (zeros suppressed), and must be used in conjuction with
c     the main memory array 'a_pointer'.
c temp_gain -- is temporary storage of some matrix products which
c     should spead the program considerable
c k_gain    -- is the Kalman filter gain matrix and gives the
c     effect of each observation on the parameter matrix
c cov_obs   -- the covariance matrix of the observations.  We also
c     include here the effects of the 'white noise' components of
c     of the Markov parameters.
c prefit    -- the prefit residual from kalfil
c
c par_num   -- the number of parameters being estimated from the
c     data set (see &kalcn).
c nbl       -- the number of baselines possible for this data
c     set (true_nsite*(true_nsite+1)/2) (see &kalcn)
c max_deriv -- Maximum number of parameters on which an observation
c     can depend (see &kalpa)
c kal_array -- The array into which kalfil is read from disk
c
