 
*     Include file for fitsnr which will compute phase variations from SNR
*     oscillations.
 
* Parmeters
 
*   max_epochs  - Maximum number of epochs allowed
*   max_sat     - Maximum number of satellites
*   max_poly        - Maximum number of sin(el) poynomial fit
*   max_seq     - Maximum length of sequence allowed for spectral
*                 fitting
*   max_cmp     - Maximum number of spectral components allowed
 
      integer*4 max_epochs, max_sat, max_poly, max_seq, max_cmp
 
      parameter ( max_epochs = 20000 )
* MOD TAH 200210: Increased from 32 to 35
      parameter ( max_sat    =   35 )
*     Increased Polynomial order from 5 to 6 when offset added.
      parameter ( max_poly   =    5 )
      parameter ( max_seq    = 5000 )
      parameter ( max_cmp    = 5000 )
 
 
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
*   min_elev                    - Minimum elevation angle to use.
*   min_snr                     - Minumim SNR to consider in the processing.
 
      real*8 epoch(max_epochs), snr_L1o(max_epochs, max_sat),
     .    snr_L2o(max_epochs, max_sat), gain_L1o(max_epochs, max_sat),
     .    gain_L2o(max_epochs, max_sat), gain_L1s(max_epochs, max_sat),
     .    gain_L2s(max_epochs, max_sat), gain_L1r(max_epochs, max_sat),
     .    gain_L2r(max_epochs, max_sat), phs_L1a(max_epochs, max_sat),
     .    phs_L2a(max_epochs, max_sat), gain_poly_L1(max_poly),
     .    gain_poly_L2(max_poly), zen_L1, zen_L2, L2_code_scale,
     .    min_elev, min_snr

* gain_L1q, gain_L2q -- Gain sequence for one-sat from rise to set.
* per_cmp(max_cmp)   - Period of component variation
* cos_cmp(max_cmp)   - Cosine amplitude in SNR
* sin_cmp(max_cmp)   - Sine amplitude in SNR
     
      real*8 gain_L1q(max_seq), gain_L2q(max_seq), per_cmp(max_cmp),
     .       cos_cmp(max_cmp) , sin_cmp(max_cmp)
 
*   az(max_epochs, max_sat)     - Azimuth (rads)
*   el(max_epochs, max_sat)     - Elevation angle (rads)
 
      real*4 az(max_epochs, max_sat), el(max_epochs, max_sat)
      
      real*4 seq(2,max_seq), work(2,max_seq)
  
*   flags(max_epochs, max_sat)      - Data flag for SNR: Bits are:
*                               - Bit Meaning
*                               -   1 SNR is measured
*                               -   2 L1 SNR is code measurement
*                               -   3 L2 SNR is code measurement
*                               -   4 L1 SNR OK, non-zero
*                               -   5 L2 SNR OK, non-zero
 
*   num_epochs      - Number of epochs
*   num_sat         - Number of satellites
*   prn_list(max_sat)  - PRN numbers that go with channels
*   num_poly_L1, num_poly_L2    - Polynomial orders for L1 and L2
*                   - sin(el) dependence fit
*   num_L2_code(2)  - Number of L1 measurements in code tracking
*                     mode (1) and number in codeless mode (2)
*   out_opts        - Output options.  Bit mapped to indicate
*                     which things should be output.
*                     Mapping is:
*                     BIT  Name   Meaning
*                       1  GAIN   elevation angle gain values
*                       2  OBSG   Observed gain by satellite and time
*                                 after removing the gain curves
*                       3  FFTC   Components of the FFT fits to the 
*                                 SNR values
*                       4  DPHS   Corrections to phase values (main result
*                                 of program)
*                       5  SMTH   Gain values determined from the FFT
*                                 coefficients used in getting the DPHS values
*                       6  RMSE   RMS of LC phase corrections as funnction of
*                                 elevation angle.

 
      integer*4 flags(max_epochs, max_sat), num_epochs, num_sat,
     .    prn_list(max_sat), num_poly_L1, num_poly_L2,
     .    num_L2_code(2), out_opts

*   start_seq    - Start epoch number of sequence of data
*   end_seq      - End  epoch number of sequence of data
*   num_seq      - Number of data in the sequence
*   mf           - Current max frequency.
*   num_cmp      - NUmber of components needed to represent SNR.
*   bnd_cmp(max_cmp)  - Band (i.e., index from i=1,mf) for coefficients

      integer*4 start_seq, end_seq, num_seq , num_cmp,
     .          bnd_cmp(max_cmp), mf
      
*   infile      - Name of input file
 
      character*128 infile, in_cf, out_cf

*   out_cmds    - List of output options (see out_opts)

      character*4 out_cmds(32)
 
 
 
*----------------------------------------------------------------------
 
*     Common declaration

      common / fitsnr_comm /  norm_eq, bvec,
     .    epoch, snr_L1o,
     .    snr_L2o, gain_L1o,
     .    gain_L2o, gain_L1s,
     .    gain_L2s, gain_L1r,
     .    gain_L2r, phs_L1a,
     .    phs_L2a, gain_poly_L1,
     .    gain_poly_L2, zen_L1, zen_L2, L2_code_scale, min_elev,
     .    min_snr, 
     .    gain_L1q, gain_L2q, per_cmp, cos_cmp, sin_cmp,
     
     .    az, el,
     .    seq, work,
    
     .    flags, num_epochs, num_sat,
     .    prn_list, num_poly_L1, num_poly_L2, num_L2_code,
     .    out_opts,
     .    start_seq, end_seq, num_seq, num_cmp, bnd_cmp, mf,
     
     .    infile, in_cf, out_cf

      common / cmds_comm / out_cmds
 
 
 
