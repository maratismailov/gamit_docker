c
c the SOLVK ema common block
c
c......................................................................
c mod history                 9:45 AM  FRI., 18  JAN., 1985
c
c who  when   what
c......................................................................
c TAH 850118 : Typed common block in
c
c
c......................................................................
c
c Variables
c ---------
c
c .........................
c Normal equation variables
c .........................
c
c norm_eq  -- the normal equations used for estimating the clock
c        offset and rate in SOLVK.
c
c sol_vec -- the left hand side of the normal equations which will
c     be replaced by the solution vector during the inversion of
c     normal equations using 'invert_vis'
c
c cl_deriv -- the partial derivatives of the clocks parameters
c        (Needs to be in ema so that VIS can be used).
c
c
      real*8 norm_eq(max_cl,max_cl), sol_vec(max_cl),
     .       cl_deriv(max_cl)
 
c
c
c ...................................
c Storage array for generating KALFIL
c ...................................
c
c kal_array -- the array in ema which will contain the KALFIL records
c        until the array is full, at which time the records will be
c        written to disk.  This procedure is used to minimize
c        'disk-thrashing'.
c
      integer*4 kal_array(krec_len,max_kal_recs)
 
c
      common /solvk_ema/  norm_eq, sol_vec, cl_deriv, kal_array
c

