
******************************************************
*     Declarations for the ppb program
 
*         max_sites     - Maximum number of sites allowed
*       max_bl          - Max number of baselines
*       max_soln            - Max number of solutions per baseline
*       max_ent         - Max number of entries allowed
*                       - (max_bl*max_soln)
 
      integer*4 max_sites, max_bl, max_soln, max_ent
 
*   Re                  - Radius of the Earth (m)
 
 
      real*8 Re
 
      parameter ( max_sites = 256 )
      parameter ( max_bl = max_sites*(max_sites-1)/2 )
      parameter ( max_soln = 20 )
      parameter ( max_ent  = max_bl*max_soln)
      parameter ( Re = 6378145.d0)
 
***** Global type declarations
 
*   num_sites       - Number of sites found
*   min_num         - Min. number of determinations
*   num_ent         - Number of entries found
*   num_est(max_ent)        - Number of values in determination
*   num_soln(max_ent)       - Solution number for each entry
*   peak_soln           - Largest solution number found.
*   site(2,max_ent)     - Baseline (site numbers)
 
      integer*4 num_sites, min_num, num_est(max_ent),
     .    num_soln(max_ent), peak_soln, site(2,max_ent),
     .    num_ent, num_used(4)
 
*   length(max_ent)     - Length of each baseline (m)
*   wrms(max_ent)           - WRMS scatter of baseline (mm)
*   chi(max_ent)            - Chi (nrms) of baseline
*   max_len             - Maximum length to be considered. (km)
*   min_len             - Minimum length (km)
 
      real*8 length(max_ent), wrms(max_ent), chi(max_ent), max_len,
     .       min_len 
 
****  Estimation normal equations
 
*   bvec_av(2), norm_av(2,2)    - Normal equations and solution
*                       - vector for average results
*                           - Normal equations and solution
*                           - vector for site results
*   hs_av(max_soln), us_av(max_soln)    - Average horizontal and
*                   - Up noise (mm)
*   hs_st(max_soln,max_sites), us_st(max_soln,max_sites) - Site
*                   - dependent  horizontal and Up noise (mm)
*   hs_apr_av(max_soln), us_apr_av(max_soln)    - Average horizontal and
*                   - Up noise (mm) Apriori values.
*   hs_apr_st(max_soln,max_sites), us_apr_st(max_soln,max_sites) - Site
*                   - dependent  horizontal and Up noise (mm)
*                   - Apriori values.
 
      real*8 bvec_av(2), norm_av(2,2),
     .    bvec_ch(2), norm_ch(2,2),
     .    bvec_st(2*max_sites), norm_st(2*max_sites,2*max_sites),
     .    hs_av(max_soln), us_av(max_soln), hs_con(max_soln),
     .    hc_av(max_soln), uc_av(max_soln), hc_con(max_soln),
     .    hs_st(max_soln,max_sites), us_st(max_soln,max_sites),
     .    hs_apr_av(max_soln), us_apr_av(max_soln), 
     .    hs_apr_con(max_soln),
     .    hc_apr_av(max_soln), uc_apr_av(max_soln), 
     .    hc_apr_con(max_soln),
     .    hs_apr_st(max_soln,max_sites), us_apr_st(max_soln,max_sites)
 
*   site_names(max_sites)   - Names of the sites
 
 
      character*8 site_names(max_sites)
 
*             sum_file  - Name of the summary file.
 
 
      character*256 sum_file
 
*---------------------------------------------------------------
* COMMON DECLARATION.
 
      common / ppb_com / length, wrms, chi, max_len, min_len,
     .    bvec_av, norm_av, bvec_ch, norm_ch,
     .    bvec_st, norm_st, hs_av, us_av, hs_con, hs_st, us_st, 
     .    hs_apr_av, us_apr_av, hs_apr_con, hs_apr_st, us_apr_st, 
     .    hc_av, uc_av, hc_con, hc_apr_av, uc_apr_av, hc_apr_con,
     .    num_sites, num_used, 
     .    min_num, num_est, num_soln, peak_soln, site, num_ent,
     .    site_names, sum_file
 
