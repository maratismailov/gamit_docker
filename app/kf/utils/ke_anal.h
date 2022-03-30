*
*     Common include block for KE_ANAL
*
*   max_bins        - Maximum number of bins for distributions
*   max_elev        - Maxumium number of elevation bins
*   max_exper   - Maximum number of experiments
 
      integer*4 max_bins, max_elev, max_exper
 
      parameter ( max_bins  = 100 )
      parameter ( max_elev  =  12 )
      parameter ( max_exper = 100 )
 
*   sig_sums(5, max_bins)   - Summation of statistics as function
*           - of sigmas.  Values are:
*           - (1) -- Number of values
*           - (2) -- Sums of weights
*           - (3) -- Sum of difference by weight
*           - (4) -- Sum of diff**2 by weight
*           - (5) -- Sum of index (eg elevation angle) by wght
*   elev_sums(5,max_bins)   - Summation of statistics as function
*           - of Elevation angle
*   rate_sums(5,max_bins)   - Summation of statistics as function
*           - of rate bins
*   frng_sums(5,max_bins)   - Fringe quality sums (numerical values
*           - only)
*   time_sums(5,max_bins)   - Statistics by UTC (1hr bins)
*   overal_sums(5)      - Overall statistics
 
      real*8 sig_sums(5, max_bins), elev_sums(5,max_bins),
     .    rate_sums(5,max_bins), frng_sums(5,max_bins),
     .    time_sums(5,max_bins), overal_sums(5)
 
*     The following are the final stats arrays
*   sig_fin(6, max_bins)    - final  statistics as function
*           - of sigmas.  Values are:
*           - (1) -- Weighted mean
*           - (2) -- sigma of mean (scaled)
*           - (3) -- WRMS about 0 offset
*           - (4) -- Chi**2 about 0 offset
*           - (5) -- Number of values
*           - (6) -- Average value in bin
*   elev_fin(6,max_bins)    - final statistics as function
*           - of Elevation angle
*   rate_fin(6,max_bins)    - final statistics as function
*           - of rate bins
*   frng_fin(6,max_bins)    - final statistics for Fringe
*           - (numerical values only)
*   time_fin(6,max_bins)    - final Statistics by UTC (1hr bins)
*   overal_fin(6)       - Overall final statistics
 
      real*8 sig_fin(6, max_bins), elev_fin(6,max_bins),
     .    rate_fin(6,max_bins), frng_fin(6,max_bins),
     .    time_fin(6,max_bins), overal_fin(6)
 
*     The following are local summations before rmoveing mean
*     difference
*   sig_lsum(5, max_bins)   - Summation of statistics as function
*           - of sigmas.  Values are:
*           - (1) -- Number of values
*           - (2) -- Sums of weights
*           - (3) -- Sum of difference by weight
*           - (4) -- Sum of diff**2 by weight
*           - (5) -- Sum of index (eg elevation angle) by wght
*   elev_lsum(5,max_bins)   - Summation of statistics as function
*           - of Elevation angle
*   rate_lsum(5,max_bins)   - Summation of statistics as function
*           - of rate bins
*   frng_lsum(5,max_bins)   - Fringe quality sums (numerical values
*           - only)
*   time_lsum(5,max_bins)   - Statistics by UTC (1hr bins)
*   overal_lsum(5)      - Overall statistics
 
*   elev_bins(max_elev+1)   - Lower side of bins for elevation
*           - angle bins
 
*   exper_mean(max_exper)   - Weighted means offset between group
*           - and phase for each experiment
*   exper_sig(max_exper)    - Sigma of the mean for each exp
*                       - (scaled by chi**2)
*   exper_chi(max_exper)    - Chi**2 for each experiment
*   exper_epoch(max_exper)  - Mid_epoch of each experiment
 
      real*8 sig_lsum(5, max_bins), elev_lsum(5,max_bins),
     .    rate_lsum(5,max_bins), frng_lsum(5,max_bins),
     .    time_lsum(5,max_bins), overal_lsum(5), elev_bins(max_elev+1),
     .    exper_mean(max_exper), exper_sig(max_exper),
     .    exper_chi(max_exper), exper_epoch(max_exper)
 
*   exper_num(max_exper)    - Number of good data in experiment
*   ne      - Number of experiments processed
 
 
      integer*4 exper_num(max_exper), ne
 
*   ko_names(max_exper) - Names of the KalObs files in this
*           - analysis
 
      character*128 ko_names(max_exper)

*   bvec(4), norm_eq(4,4) - Normal equations and solution  
*   --added by rwk 960730

      real*8 bvec(4), norm_eq(4,4), exper_soln(4,100), exper_snsg(4,100)
 
*-----------------------------------------------------------------
* COMMON DECARATION
*
      Common / ke_anal_comm / sig_sums, elev_sums, rate_sums,
     .    frng_sums, time_sums, overal_sums, sig_fin, elev_fin,
     .    rate_fin, frng_fin, time_fin, overal_fin, sig_lsum,
     .    elev_lsum, rate_lsum, frng_lsum, time_lsum, overal_lsum,
     .    elev_bins, exper_mean, exper_sig, exper_chi, exper_epoch,
     .    bvec, norm_eq, exper_soln, exper_snsg,
     .    exper_num, ne, ko_names
 
