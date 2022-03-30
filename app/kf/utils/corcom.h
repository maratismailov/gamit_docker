 
*     CORCOM.H Include file for corcom.f
 
*   max_sites   - Maximum number of sites allowed (also applies
*               - to ties).
 
      integer*4 max_sites
 
      parameter ( max_sites = 9192 )
 
 
***** Common block declarations
 
*   sys1_coord(3,2,max_sites)   - Positions and velocities of
*                           - each site in system 1 (m, m/yr)
*   sys2_coord(3,2,max_sites)     - Positions and velocities of
*                           - each site in system 2 (m, m/yr)
*   sys1_epoch(max_sites), sys2_epoch(max_sites)    - Epochs for the
*                           - positions in systems 1 and 2
*                           - (Deciminal years)
*   ties_coord(3,max_sites) - Coordinate differences for the
*                           - ties.
*   frame_epoch             - Epoch for output frame reference
*                           - (defines epoch at orientation of
*                           - sys 1 transformed and sys 2 will
*                           - be the same.
*   height_weight           - Weight to be given to heights in
*                           - transformation solution.
 
      real*8 sys1_coord(3,2,max_sites), sys2_coord(3,2,max_sites),
     .    sys1_epoch(max_sites), sys2_epoch(max_sites),
     .    ties_coord(3,max_sites), frame_epoch, height_weight
 
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
*   sum_postfit             - Estimates of sum of postfit
*                           - residuals**2
*   rms_fit                 - Overall RMS of fit (m)
*   stats_fund(3,3)         - Summation statics for NE and U
*                           - (number, sum res, sum res**2)
*   stats_all(3,3)          - Summation statistics for all
*                           - sites.
 
      real*8 norm_eq(7,7), bvec(7), trans_parm(7), trans_parm_out(7),
     .    trans_sigma(7), sum_prefit, sum_postfit, rms_fit,
     .    stats_fund(3,3), stats_all(3,3)
 
*   num_sys1, num_sys2      - Number of sites in systems 1 & 2
*                           - (number in sys 1 will typically
*                           - increase after the ties have
*                           - been applied)
*   num_ties                    - Number of ties
*   num_fund                    - Number of "fundamental" stations
*                           - to be used in transforming sys 1
*                           - to sys2
*   fund_links(2,max_sites) - Pairs of sites from sys 1 and
*                           - sys 2 to used for transformation
*   num_parn                    - Number of transformation
*                           - parameters to be estimated (6 if
*                           - no scale, 7 if scale.)
*   out_unit                    - Unit number for output (200 if
*                           - file, 6 if screen output)
*   num_data                    - Number of components used in
*                           - estimating the transformation.
 
      integer*4 num_sys1, num_sys2, num_ties, num_fund,
     .    fund_links(2,max_sites), num_parn, out_unit, num_data
 
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
*   ties_names(2,max_sites) - Names of the ties sites (order
*                           - read from input.  Second site is
*                           - mapped to first site name).
*   ties_type(max_sites)    - The type of tie given.  If not NEU
*                             then assummed to be XYZ.
*   fund_names(max_sites)       - Names of the sites to be used in
*                           - transforming sys1 to sys2
 
 
      character*8 sys1_names(max_sites), sys2_names(max_sites),
     .    ties_names(2,max_sites), fund_names(max_sites),
     .    ties_type(max_sites)
 
*   sys1_file, sys2_file        - Names of the system 1 and system
*                           - 2 files.
*   out_file                    - Name of the output file.
*   ties_file               - Name of the ties files
*   fund_file               - Name of the file containing the
*                           - list of fundamental sites to be
*                           - used in transforming sys1 to
*                           - sys2
 
 
      character*256 sys1_file, sys2_file, out_file, ties_file,
     .    fund_file
 
*----------------------------------------------------------------------
*     COMMON BLOCK DECLARATION
 
*   sum_prefit              - Sum of prefit residual**2.
*   sum_postfit             - Estimates of sum of postfit
*   rms_fit                 - Overall RMS of fit (m)
 
      common / corcom_com / sys1_coord, sys2_coord,
     .    sys1_epoch, sys2_epoch, ties_coord, frame_epoch,
     .    height_weight, norm_eq, bvec, trans_parm, trans_parm_out,
     .    trans_sigma, sum_prefit, sum_postfit, rms_fit, stats_fund,
     .    stats_all, num_sys1, num_sys2, num_ties, num_fund,
     .    fund_links, num_parn, out_unit, num_data, sys1_frame,
     .    sys2_frame, out_frame, sys1_names, sys2_names, ties_names,
     .    fund_names, sys1_file, sys2_file, out_file, ties_file,
     .    fund_file, ties_type
 
*---------------------------------------------------------------------
 
