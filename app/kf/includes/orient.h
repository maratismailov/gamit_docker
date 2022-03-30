 
*     This include file containes the variables used in the
*     orient_system group of subroutines.  These routines are used to
*     "fixed" the rotational and translation orientation of the
*     coordinate systems in GPS and VLBI.  For these routines to work
*     properly all site corrdinates (and rates) should be estimated with
*     weak constraints.  Baseline length results and sigmas should
*     not be affected by the applaication of these conditions.
*
*   nes     - Number of sites whose positions have been
*           - estimated
*   nos     - number of sites in the origin definition
 
*   irotran - Start location of the rotran (rotation and
*           - translation) transformation matrix
*   ixtol   - Start location of the XYZ global to local NEU
*           - transformation matrix
*   icov_con    - Start of the covariance matrix of the local NE
*           - coordinates for the sites used in the conditiion
*   isol_con    - Start of NE solution vector
*   itr     - start of the combined rotran and xtol
*           - transformation matrix
*   iscratch    - Array space for scratch use in saving intermediate
*           - prodocts
 
      integer*4 nes, nos, irotran, ixtol, icov_con, isol_con, itr,
     .    iscratch
 
*   norm_theta(6,6) - Normal equations for the rot/tran
*           - parameter estimates
*   b_theta(6)  - Solution vector for rot/tran
*   theta_pos(6)    - Six estimates of the rotation/translation for
*               - position
*   theta_rate(6)   - Six estimates of the rotation/translation
*               - for rate estimates
*   chi_pos, chi_rate   - Chi**2/f for the position and rate
*               - estimates.
 
      real*8 norm_theta(6,6), b_theta(6), theta_pos(6), theta_rate(6),
     .    chi_pos, chi_rate
 
 
      common / orient_com / norm_theta, b_theta, theta_pos, theta_rate,
     .    chi_pos, chi_rate, nes, nos, irotran, ixtol, icov_con,
     .    isol_con, itr, iscratch
 
