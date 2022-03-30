 
*     Include file for compare_pmu
 
*   max_sm      - Maximum number of smoothed values allowed
*   max_pr      - Maximum number of values allowed in the
*               - primary PMU series.
*   max_tab     - Maximum number of values in table for UT1
*               - (Needed to take out discontinuities in
*               - UT1-UTC)
 
      integer*4 max_sm, max_pr, max_tab
 
*                                     ! Over 20 years of daily values
      parameter ( max_sm = 10000 )
*                                     ! 4 times the current VLBI data
      parameter ( max_pr = 10000 )
*                                 ! set.  But allow for smooth data
*                                 ! to be re-read.
      parameter ( max_tab = 10000 )
 
*     Start the declarations
 
*   sm_jd(max_sm)   - Julian dates of the smoothed values
*   sm_xp(max_sm)   - X-pole position for smoothed values (asec)
*   sm_yp(max_sm)   - Y-pole position for smoothed values (asec)
*   sm_UT(max_sm)   - UT1-AT for smoothed values (tsec)
*   sm_xsig(max_sm) - Sigma for smooth X-pole position (asec)
*   sm_ysig(max_sm) - Sigma for smooth Y-pole position (asec)
*   sm_usig(max_sm) - Sigma for smooth UT1-AT (tsec)
*   sm_xfwhm(max_sm) - FWHM For x pole (days)
*   sm_yfwhm(max_sm) - FWHM For y pole (days)
*   sm_ufwhm(max_sm) - FWHM For ut1     (days)
*   sm_start, sm_stop   - Start and stop JD for smooth series
*   sm_step         - Step size in days for smooth series
*   sm_fwhm_scale   - Scaling to be used on the fwhm of filter
*                   - for smoothing primary series.
 
*   pr_jd(max_pr)   - Julian dates of the primary values
*   pr_xp(max_pr)   - X-pole position for primary values (asec)
*   pr_yp(max_pr)   - Y-pole position for primary values (asec)
*   pr_UT(max_pr)   - UT1-AT for primary values (tsec)
*   pr_xsig(max_pr) - Sigma for primary X-pole position (asec)
*   pr_ysig(max_pr) - Sigma for primary Y-pole position (asec)
*   pr_usig(max_pr) - Sigma for primary UT1-AT (tsec)
*   pr_start, pr_stop   - Start and stop JD for primary series
 
*   prx_stats(4)        - Statistics accumulation for primary-smooth
*   pry_stats(4)        - Statistics accumulation for primary-smooth
*   pru_stats(4)        - Statistics accumulation for primary-smooth
*   secx_stats(4)   - Statistics accumulation for secondary-smooth
*   secy_stats(4)   - Statistics accumulation for secondary-smooth
*   secu_stats(4)   - Statistics accumulation for secondary-smooth
 
*   sec_jd          - Julian dates of the primary values
*   sec_xp          - X-pole position for primary values (asec)
*   sec_yp          - Y-pole position for primary values (asec)
*   sec_UT          - UT1-AT for primary values (tsec)
*   sec_xsig            - Sigma for primary X-pole position (asec)
*   sec_ysig            - Sigma for primary Y-pole position (asec)
*   sec_usig            - Sigma for primary UT1-AT (tsec)
*   sec_start, sec_stop - Start and stop JD for primary series
 
      real*8 sm_jd(max_sm), sm_xp(max_sm), sm_yp(max_sm),
     .    sm_UT(max_sm), sm_xsig(max_sm), sm_ysig(max_sm),
     .    sm_usig(max_sm), sm_start, sm_stop, sm_step, sm_fwhm_scale,
     .    pr_jd(max_pr), pr_xp(max_pr), pr_yp(max_pr), pr_UT(max_pr),
     .    pr_xsig(max_pr), pr_ysig(max_pr), pr_usig(max_pr),
     .    pr_start, pr_stop, prx_stats(4), pry_stats(4), pru_stats(4),
     .    secx_stats(4), secy_stats(4), secu_stats(4), sec_jd, sec_xp,
     .    sec_yp, sec_UT, sec_xsig, sec_ysig, sec_usig,
     .    sec_start, sec_stop,
     .    sm_xfwhm(max_sm), sm_yfwhm(max_sm),  sm_ufwhm(max_sm) 
 
*   num_sm          - Number of values in the smooth series
*   num_pr          - Number of values in the primary series
*   in_unit         - Unit number for reading files
*   out_unit            - Unit number for writing files.  Only one
*                   - of each needed.
 
      integer*4 num_sm, num_pr, in_unit, out_unit
 
*   gen_sm          - Tells whether to generate smoothed series
*                   - .false. if out_sm = 'USEPR', .true.
*                   - otherwise.
 
      logical gen_sm
 
*   UT1_sm          - Definition of UT1 for smoothed ('R' will
*                   - make series regularized, tides removed)
*   UT1_pr          - Definition of UT1 for priminary ('R' if
*                   - series is regularized, tides removed)
*   UT1_sec         - Definition of UT1 for secondary ('R' if
*                   - series is regularized, tides removed)
 
      character*8 UT1_sm, UT1_pr, UT1_sec
 
*   out_sm          - Name of output file for smoothed series
*   in_pr           - Name of input file for input primary
*   out_pr          - Name of output file for primary residuals
*   in_sec          - Name of input file for input secondary
*   out_sec         - Name of output file for secondary
*                   - residuals
 
      character*128 out_sm, in_pr, out_pr, in_sec, out_sec
 
*-----------------------------------------------------------------
* COMMON DECLARATION
 
      common / compare_com / sm_jd, sm_xp, sm_yp, sm_UT, sm_xsig,
     .    sm_ysig, sm_usig, sm_start, sm_stop, sm_step, sm_fwhm_scale,
     .    pr_jd, pr_xp, pr_yp, pr_UT, pr_xsig, pr_ysig, pr_usig,
     .    pr_start, pr_stop, prx_stats, pry_stats, pru_stats,
     .    secx_stats, secy_stats, secu_stats, sec_jd, sec_xp, sec_yp,
     .    sec_UT, sec_xsig, sec_ysig, sec_usig, sec_start, sec_stop,
     .    sm_xfwhm, sm_yfwhm,  sm_ufwhm,
     .    num_sm, num_pr, in_unit, out_unit, gen_sm, UT1_sm, UT1_pr,
     .    UT1_sec, out_sm, in_pr, out_pr, in_sec, out_sec
 
*--------------------------------------------------------------------
 
 
 
 
