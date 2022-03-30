*
* Include file for processing non-secular aproiri coordinates in 
* globk utility packages corcom, merge_apr.f and unify_apr.f
*
      integer*4 num_nonsut_types
      parameter ( num_nonsut_types = 4 )

      integer*4 max_nonsut   ! Maximum of non-secular entries 
                             ! (should >= max_nonec in kalman_param.h)
      parameter ( max_nonsut = 8192 )

* NOTE: Site number is generated from a site list that is contained
*       in someother include (depends on program).
*
* Entries for non-secular position changes
*-----------------------------------------
* num_nonsut -- Number of non-secular terms
* param_nonsut(2,max_nonsut) -- Gives the site number (1) and the
*     type of non-secular term.  The type defined are:
*     1 -- Offset and Rate change (applied after reference epoch)
*     2 -- Periodic (cosine and sine terms)
*     3 -- Exponential (exp(-t/tau))
*     4 -- Logarithmic (ln(1+t/tau))
*    >5    Not defined.

      integer*4 num_nonsut, param_nonsut(2,max_nonsut)

* apr_sig_nonsut(2,max_nonsut) -- Apriori sigmas for the estimated 
*     non-secular terms

      real*4 apr_sig_nonsut(2,max_nonsut) 

* Real*8 values
*   apr_val_nonsut(8,max_nonsut)  - Non-secular parameters for station
*                  positions.  The order is:
*                  JD -- Julian date either for periodic 0 phase
*                        of start of exponential/logarithm/steps
*                  parameter -- Either period or delay times (days)
*                  Paired for X Y and Z (read in as NEU and converted)
*                  coef 1    -- Exponential/logarithm constant;
*                               cosine term or offset (m)
*                  coef 2    -- skipped for exponential and log,
*                               sine term (m) or rate change (m/yr)
*                  (see also I*4 param_nonsut(2,max_nonsut)
*   est_val_nonsut(8,max_nonsut) -- Parameters for the estimated values
*       of the non-secular terms.  (param_est_nons defines the other 
*       parameters of the non-secular estimates).

      real*8 apr_val_nonsut(8,max_nonsut), est_val_nonsut(8,max_nonsut) 

      character*8 nonsut_types(num_nonsut_types)

      common / nonsut_CH / nonsut_types
      common / nonsut_4B / num_nonsut, param_nonsut, apr_sig_nonsut
      common / nonsut_8B / apr_val_nonsut, est_val_nonsut
