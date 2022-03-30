 
*     VELROT.H Include file for velrot.f

* MOD TAH 001010: Vers 1.01 -- Added AV_DIST command which will 
*     average the velocity over the distance specified.
 
*   max_sites   - Maximum number of sites allowed (also applies
*               - to ties).
*   max_links   - Maxiumum of linked fundamental sites
*   velrot_ver  - Version number
 
      integer*4 max_sites, max_links, max_swrds
      character*(*) velrot_ver
 
      parameter ( max_sites = 4096 )
      parameter ( max_swrds = max_sites/32+1 )
      parameter ( max_links = 4*max_sites)
      parameter ( velrot_ver = '1.01' )
 
***** Common block declarations
 
*   sys1_coord(3,2,max_sites)   - Positions and velocities of
*                           - each site in system 1 (m, m/yr)
*   sys2_coord(3,2,max_sites)     - Positions and velocities of
*                           - each site in system 2 (m, m/yr)
*   sys1_cov(3,3,max_sites) - Covariance matrix for NEU velocities
*   sys2_cov(3,3,max_sites) - Covariance matrix for NEU velocities

*   height_weight           - Weight to be given to heights in
*                           - transformation solution.
*   eq_dist - Distance over which velocities should be linked 
*             (meters)
*   cp_dist - Comparison distance (will cause sites to be marked
*             on output (m)
*   av_dist - Distance over which the velocities will be averaged
 
      real*8 sys1_coord(3,2,max_sites), sys2_coord(3,2,max_sites),
     .       sys1_cov(3,3,max_sites), sys2_cov(3,3,max_sites),
     .       height_weight, eq_dist, cp_dist, av_dist
 
*     Transformation quantities
*   norm_eq(7,7), bvec(7)       - Transformation normal equations
*                           - and solution vector.
*   trans_parm(7)           - Actual transformation parmeters
*   trans_parm_out(7)           - Estimated transformation
*                           - parameters in units of the
*                           - output (X-pole,Y-pole, UT1,
*                           - [all mas], XYZ translation (m),
*                           - Scale (ppb).
*   trans_sigma(7)          - Uncertainties in the parameters
*                           - (scaled by RMS fit)
*   sum_prefit              - Sum of prefit residual**2.
*   sum_weight              - Sum of the weights 
*   sum_postfit             - Estimates of sum of postfit
*                           - residuals**2
*   rms_fit                 - Overall RMS of fit (m)
*   chi_fit                 - Chi of fit
*   stats_fund(3,3)         - Summation statics for NE and U
*                           - (number, sum res, sum res**2)
*   stats_all(3,3)          - Summation statistics for all
*                           - sites.
 
      real*8 norm_eq(7,7), bvec(7), trans_parm(7), trans_parm_out(7),
     .    trans_sigma(7), sum_prefit, sum_postfit, sum_weight, rms_fit,
     .    stats_fund(3,3), stats_all(3,3), frame_epoch, chi_fit
 
*   num_sys1, num_sys2      - Number of sites in systems 1 & 2
*                           - (number in sys 1 will typically
*                           - increase after the ties have
*                           - been applied)
*   num_fund                    - Number of "fundamental" stations
*                           - to be used in transforming sys 1
*                           - to sys2
*   fund_links(2,max_links) - Pairs of sites from sys 1 and
*                           - sys 2 to used for transformation
*   num_parn                    - Number of transformation
*                           - parameters to be estimated (6 if
*                           - no scale, 7 if scale.)
*   out_unit                    - Unit number for output (200 if
*                           - file, 6 if screen output)
*   num_data                    - Number of components used in
*                           - estimating the transformation.
*   av_sys1(max_swrds), av_sys2(max_swrds) - Bit mapped variable
*                             which shows we have averaged site
*   av_u1(max_swrds), av_u2(max_swrds) - Records the site used in
*                             specific averaging case.
 
      integer*4 num_sys1, num_sys2, num_fund,
     .    fund_links(2,max_links), num_parn, out_unit, num_data,
     .    av_sys1(max_swrds), av_sys2(max_swrds), 
     .    av_u1(max_swrds), av_u2(max_swrds) 
 
*   sys1_frame              - Frame for system 1 (either NAFD
*                           - AM0-2)
*   sys2_frame              - Frams for system 2 (as above)
*   out_frame               - Frame for output (as above)
*                           - except entried with deciminal
*                           - year for the reference of the
*                           - frame).
 
      character*10 sys1_frame, sys2_frame, out_frame
 
*   sys1_names(max_sites)       - Names of sites in system 1
*   sys2_names(max_sites)       - Names of sites in system 2
*   fund_names(max_sites)       - Names of the sites to be used in
*                           - transforming sys1 to sys2
 
 
      character*8 sys1_names(max_sites), sys2_names(max_sites),
     .     fund_names(max_sites)
 
*   sys1_file, sys2_file        - Names of the system 1 and system
*                           - 2 files.
*   out_file                    - Name of the output file.
*   fund_file               - Name of the file containing the
*                           - list of fundamental sites to be
*                           - used in transforming sys1 to
*                           - sys2
 
 
      character*256 sys1_file, sys2_file, out_file, 
     .    fund_file 

*   param_opt  -- Options for parameters to be estimated
*                 T -- Translation
*                 R -- Rotation
*                 S -- Scale
*                 L -- Local, two parameters only

      character*8 param_opt
 
*----------------------------------------------------------------------
*     COMMON BLOCK DECLARATION
 
*   sum_prefit              - Sum of prefit residual**2.
*   sum_postfit             - Estimates of sum of postfit
*   rms_fit                 - Overall RMS of fit (m)
 
      common / corcom_com / sys1_coord, sys2_coord,
     .    sys1_cov, sys2_cov, frame_epoch,
     .    height_weight, norm_eq, bvec, trans_parm, trans_parm_out,
     .    trans_sigma, sum_prefit, sum_postfit, 
     .    sum_weight, rms_fit, chi_fit, eq_dist, cp_dist,
     .    av_dist, stats_fund,
     .    stats_all, num_sys1, num_sys2,  num_fund,
     .    fund_links, num_parn, out_unit, num_data, 
     .    av_sys1, av_sys2, av_u1, av_u2, sys1_frame,
     .    sys2_frame, out_frame, sys1_names, sys2_names, 
     .    fund_names, sys1_file, sys2_file, out_file,
     .    fund_file, param_opt
 
*---------------------------------------------------------------------
 



