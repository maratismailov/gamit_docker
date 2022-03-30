c
c  The KALFIL definition common block
c
c.......................................................................
c Mod history                12:07 PM  WED., 27  MAR., 1985
c WHO  WHEN   WHAT
c ......................................................................
c TAH 850327: Typed in definitions
c TAH 850401: Put the elevation angles (real*4) into kalfil
c             so that the atmospheric noise contribution could be
c             computed in FORSL.
c TAH 870415: Added observation duration, KalObs file record
c             record number.
c ......................................................................
c
c
*   kalobs_rec  - The record number from the KalObs file for this
*               - record.
*   kobs_dur    - the duration of this observation in secs (obtained
*               - at READIN time).
 
*   kunw_flag   - the unweight flag obtained by anding the data_mask
*               - and data_flag.
 
      integer*4 kalobs_rec, kobs_dur, kunw_flag
 
c Variables
c ---------
c
c ksite -- the site number for this observation
c kstar -- the star number for this observation
c kunw_flag -- the unweight flag.  Now Integer*4.
c kal_parn -- the number of partial derivatives for this observation
c kelev   -- the elevation angles at each site in radians (real*4)
c kepoch  -- the Julian day and fraction of a day for this obs.
c prefit_res -- the observed-theoretical delay or rate for this
c     observation.  The observable type of a record in kalfil
c     can be determined by 'delay_type' which indicates which
c     observables will be used in the solution.
c sigma_res -- the SNR sigma of the prefit_res
c pointer -- an array which tells which parameter partial derivatives
c     are stored in the deriv array
c deriv -- the partial derivative array (containing only non-zero
c     partial derivatives)
c padding -- one word to ensure the file record length is even
c
      integer*4 ksite(2), kstar, kal_parn, padding
 
c
      real*4 kelev(2)
 
c
      real*8 kepoch, prefit_res, sigma_res
 
c
      integer*4 pointer(max_deriv)
 
c
      real*8 deriv(max_deriv)
 
c
c
      common /kal_file_des/ ksite, kstar, kunw_flag, kal_parn,
     .    kelev, pointer, kepoch, prefit_res, sigma_res, deriv,
     .    kobs_dur, kalobs_rec, padding
c

