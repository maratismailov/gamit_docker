 
*     Include file for fitsnr which will compute phase variations from SNR
*     oscillations.
 
* Parmeters
 
*   max_epochs  - Maximum number of epochs allowed
*   max_sat     - Maximum number of satellites
*   max_poly        - Maximum number of sin(el) poynomial fit
*   max_rf      - Maxium number of reflectors
 
      integer*4 max_epochs, max_sat, max_poly, max_rf
 
      parameter ( max_epochs = 6000 )
* MOD TAH 200210: Increased from 32 to 35
      parameter ( max_sat    =   35 )
      parameter ( max_poly   =    5 )
      parameter ( max_rf     =   10 )
 
 
*----------------------------------------------------------------------
 
*     Real*8 declarations
 
*   norm_eq(max_poly, max_poly) - Normal equations for fiting
*                               - gain curves
*   bvec(max_poly)              - Solution vector.
 
      real*8 norm_eq(max_poly, max_poly), bvec(max_poly)
 
*   epoch(max_epochs)               - JD of measurements plus fraction
*                               - of day
*   snr_L1o(max_epochs, max_sat)    - L1 oberved SNR values
*   snr_L2o(max_epochs, max_sat)    - L2 oberved SNR values
*   gain_L1o(max_epochs, max_sat)   - L1 oberved gain values (SNR with
*                               - elevation angle dependent
*                               - polynomial removed).
*   gain_L2o(max_epochs, max_sat)   - L2 oberved gain values
*   gain_L1s(max_epochs, max_sat)   - L1 smoothed gain (smoothed to
*                               - different levels during processing
*   gain_L2s(max_epochs, max_sat)   - L2 smoothed gain (smoothed to
*                               - different levels during processing
*   gain_L1r(max_epochs, max_sat)   - L1 residual gain (difference
*                               - between observed and smoothed)
*   gain_L2r(max_epochs, max_sat)   - L2 residual gain.
 
*   phs_L1a(max_epochs, max_sat)    - Phase adjustments at L1 (cycles)
*   phs_L2a(max_epochs, max_sat)    - Phase adjustments at L1 (cycles)
 
*   gain_poly_L1(max_poly)      - Polynomial coefficients for L1 fit
*   gain_poly_L2(max_poly)      - Polynomial coefficienst for L2 fit
*                               - (In X-correlation mode, L1 gain
*                               - removed first).
*   zen_L1, zen_L2              - Zenith values for the expected SNR
*                               - at L1 and L2
*   L2_code_scale               - Scaling to be applied to L2 code
*                               - tracking SNR values (after L2 SNR
*                               - removed).
 
      real*8 epoch(max_epochs), snr_L1o(max_epochs, max_sat),
     .    snr_L2o(max_epochs, max_sat), gain_L1o(max_epochs, max_sat),
     .    gain_L2o(max_epochs, max_sat), gain_L1s(max_epochs, max_sat),
     .    gain_L2s(max_epochs, max_sat), gain_L1r(max_epochs, max_sat),
     .    gain_L2r(max_epochs, max_sat), phs_L1a(max_epochs, max_sat),
     .    phs_L2a(max_epochs, max_sat), gain_poly_L1(max_poly),
     .    gain_poly_L2(max_poly), zen_L1, zen_L2, L2_code_scale
     
*   h(max_rf)  - Heights of the refectors in L1 wavelengths
*   r(max_rf)  - Reflection coefficients
*   an(max_rf) - Angle reflectors make to ground (entered in degrees
*                saved in radians)

      real*8 h(max_rf),  r(max_rf),  an(max_rf) 
 
*   az(max_epochs, max_sat)     - Azimuth (rads)
*   el(max_epochs, max_sat)     - Elevation angle (rads)
 
      real*4 az(max_epochs, max_sat), el(max_epochs, max_sat)
 
*   flags(max_epochs, max_sat)      - Data flag for SNR: Bits are:
*                               - Bit Meaning
*                               -   1 SNR is measured
*                               -   2 L1 SNR is code measurement
*                               -   2 L2 SNR is code measurement
 
*   num_epochs      - Number of epochs
*   num_sat         - Number of satellites
*   prn_list(max_sat)  - PRN numbers that go with channels
*   num_poly_L1, num_poly_L2    - Polynomial orders for L1 and L2
*                   - sin(el) dependence fit
*   num_rf          - Number of reflectors
 
 
      integer*4 flags(max_epochs, max_sat), num_epochs, num_sat,
     .    prn_list(max_sat), num_poly_L1, num_poly_L2, num_rf
 
*   infile      - Name of input file
 
      character*128 infile
 
 
 
*----------------------------------------------------------------------
 
*     Common declaration

      common / fitsnr_comm /  norm_eq, bvec,
     .    epoch, snr_L1o,
     .    snr_L2o, gain_L1o,
     .    gain_L2o, gain_L1s,
     .    gain_L2s, gain_L1r,
     .    gain_L2r, phs_L1a,
     .    phs_L2a, gain_poly_L1,
     .    gain_poly_L2, zen_L1, zen_L2, L2_code_scale,
     .    h, r, an,
     
     .    az, el,
 
     .    flags, num_epochs, num_sat,
     .    prn_list, num_poly_L1, num_poly_L2, num_rf,
 
     .    infile
 
 
 
