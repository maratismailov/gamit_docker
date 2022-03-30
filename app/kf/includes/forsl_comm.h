c
c Control common block for FORSL
c
c ..................................................................
c mod history                10:42 AM  MON.,  1  APR., 1985
c WHO  WHEN   WHAT
c ..................................................................
c TAH 850401  Typed in the common block
c TAH 870211  Updated for new file structures and the addition of
c             higher order polynomials for the apriori clocks.
c TAH 870415  Modified to allow space for the KalObd file record
c             number, and the duration of the observation.
c...................................................................
c
c
c
*   save_kalobs_rec(max_obep)   - The record numbers from the
*           - KalObs file for each observation at this epoch
*   save_dur(max_obep)      - The duration of each of the
*           - observations at this epoch (sec)
 
*   save_unw(max_obep)  - The unweight flags for the observations
*           - at this epcoh
 
      integer*4 save_kalobs_rec(max_obep), save_dur(max_obep),
     .    save_unw(max_obep)
 
c .......................
c epoch control variables
c .......................
c
c current_time -- the epoch of the current set of observations (Julian
c     day + fraction of day) days.
c previous_time -- the epoch of the previous observations (days)
c delta_time  -- the difference between current_time and previous_time
c     (days).  This value must always be positive i.e., the filter
c     must always go forward
c
      real*8 current_time, previous_time, delta_time
 
c
c ...................................
c Partials matrix pointer information
c ...................................
c part_pnt -- the pointers to the partial matrix 'a_part'.  The numbers
c     in this array are arranged in pairs as (start partial number,
c     number of contigous partial numbers at this partial number) i.e.,
c     the partial pointers (3,4,5, 8,9, 11) would be represented as
c     (3,3 8,2, 11,1)
c index_pnt -- an array which indicates the number of pairs of pointers
c     in part_pnt for each observation.
c max_krec  -- the maximun number of kalfile records which can be read
c     into ema for this run
c nbl  -- the maximum number of observations per epoch which could be
c     expected during this run
c nepoch_obs -- the number of (weighted) data points in the epoch
c     being processed
c nepoch_tot -- the total number of data points in the epoch being
c     processed.
c
      integer*4 part_pnt(max_deriv,max_obep), index_pnt(max_obep),
     .    max_krec, nbl, nepoch_obs, nepoch_tot
 
c
c ....................................
c Information about data in this epoch
c ....................................
c save_site -- the site involved for each observation in the current
c     epoch.
c save_souc -- the source observed in the current epoch
c elev_ang  -- the elevation angles for the observations at the
c     current epoch
c save_rec  -- the record numbers from kalfile (saved for writing in
c     bak_file
c
      integer*4 save_site(2,max_obep), save_souc, save_rec(max_obep)
 
c
      real*4 elev_ang(max_sites)
 
c
c .........................
c Back solution information
c .........................
c bak_recl -- the record length of the bak_file (used in BAKSL)
c
      integer*4 bak_recl
 
c
c ..........................
c Prefit residual statistics
c ..........................
c sum_chi_prefit -- the chi**2 of the prefit residuals
c sum_data_prefit -- the number of data used in the calculation of
c     sum_chi_prefit
c
      real*8 sum_chi_prefit, sum_data_prefit
c
      common / forsl_control / current_time, previous_time, delta_time,
     .    part_pnt, index_pnt, max_krec, nbl, nepoch_obs, nepoch_tot,
     .    save_site, save_souc, elev_ang, bak_recl,
     .    save_rec, save_unw, save_dur, save_kalobs_rec,
     .    sum_chi_prefit, sum_data_prefit
 
