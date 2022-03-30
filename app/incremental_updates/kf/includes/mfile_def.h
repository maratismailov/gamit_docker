*
*     Definition of the mfile contents for use in the
*     autcln and ctogobs.
* MOD TAH 980915: Updated m-file declaration for the addition
*     of multiple gradients and new satellite orbital elements.

* VERSIONS: Original  931
*           980915    940 -- Added multiple gradients.
*           180601   1061 -- Changed slot numbers for ECOMC model
*                            (Satellite PCO moved from 2000 to 2400)
*           200126   1071 -- Slot mumber changes for ECOMC
*           210710   1072 -- Slot number changes for 45 Beidou satellites.
*                            (changes limit maximum number of satellites
*                            to 45.  Can't be increased more without more
*                            changes to the slot numbers.
*
* PARMATERS needed in the definition of mfile
 
*   mf_maxprm   - Maximum number of parameters allowed.
*   mf_maxprm_cf  - Maximum number of parameters allowed
*                 in cfile partials.
*   mf_maxcfl	- Maximum number of c-files allowed.
*   mf_maxtfl   - Maximum number of t-files (and hence
*               - sessions).
*   mf_maxatm   - Maximum number of zenith delays allowed.
*   mf_maxsat   - Maximum number of satellites.
*   mf_maxorb   - Maximum orbital parameters
*                 MOD TAH 980915: Changed value to 6IC+3Rad+6one-per-rev+
*                                 3 offsets = 18 (was 15).
*   mf_maxgrad  - Maximum number of gradient parameters per station.
 
      integer*4 mf_maxprm, mf_maxcfl, mf_maxtfl, mf_maxatm,
     .          mf_maxsat, mf_maxorb, mf_maxprm_cf, mf_maxgrad

* MOD SIMON 050907: Changed from 50 to 100 max cfiles. 
      parameter ( mf_maxcfl    =  100 ) 
      parameter ( mf_maxtfl    =   5 )
* MOD TAH 990915: Changed maxatm and maxgrad from 24/24 to 49/25.
* MOD RWK 170324: Changed maxatm from 49 to 113
      parameter ( mf_maxatm    = 113 )
      parameter ( mf_maxgrad   =  25 )
* MOD TAH 200210: Increased from 32 to 35
* MOD TAH 201031: Increased from 35 to 45 for Beidou
      parameter ( mf_maxsat    =  45 )
* MOD TAH 980915: changed 15 to 18
* MOD TAH 190520: Changed from 18 to 22 to accomomdate ECOMC model      
      parameter ( mf_maxorb    =  22 )

c     parameter ( mf_maxprm_cf =  13 + mf_maxorb*mf_maxsat + 30 )
      
c     parameter ( mf_maxprm = (6+mf_maxatm+mf_maxgrad)*mf_maxcfl +
c    .                        (mf_maxorb+3)*mf_maxsat +
c    .                         mf_maxcfl*mf_maxsat*2  + 6 )
* MOD TAH 000221: Changed algorithm for computing number of slots and
*      fixed small bug in mf_max_prm)
*      (+2 on orbit is to account for bias flags, +6 is for EOP)
      parameter( mf_maxprm_cf =  3 + mf_maxatm + 2*mf_maxgrad +
     .                          (mf_maxorb+2)*mf_maxsat + 6 ) 

      parameter ( mf_maxprm = (6+mf_maxatm+2*mf_maxgrad)*mf_maxcfl +
     .                        (mf_maxorb+3)*mf_maxsat +
     .                         mf_maxcfl*mf_maxsat*2  + 6 )

* The six at the end are for the eop parameters.

*---------------------------------------------------
* BLOCK 1: Parameter estimates and aprioris

*   mf_nversn   - M-file version number 
*   mf_ndy      - Number of sessions.  Initially
*               - autcln will handle only single session
*   mf_nepch    - number of epochs in this analyis
*   mf_mtpart   - Total number of entries in the solution
*               - vectors define defined below.  Max size
*               - is mf_maxprm
 
*   mf_idms(mf_maxprm)  - Sets which parameters could be degrees
*               - minitues and seconds. (not used).
*   mf_islot_all(mf_maxprm) - Parameter codes for solution list.
*               - (see GAMIT fills1 for definition)
*   mf_nsat     - Number of satellites
*   mf_isat(mf_maxsat)   - List of PRN numbers for the sattellites
*               - used.
*   mf_nsite    - Number of sites in mfile
*   mf_nrfile   - Number of cfiles (total for this mfile).
*   mf_ntfile   - Number of tfiles
*   mf_norb     - Number of orbital parameters in cfile.
*               - (15 with Berne model, 9 in baseline mode).
 
      integer*4 mf_nversn, mf_ndy, mf_nepch, mf_mtpart,
     .     mf_idms(mf_maxprm), mf_islot_all(mf_maxprm), mf_nsat, 
     .     mf_isat(mf_maxsat),mf_nsite, mf_nrfile, mf_ntfile, mf_norb
 
*   mf_aprval(mf_maxprm)    - Apriori values of parameter
*   mf_adjust(mf_maxprm)    - Adjustments to the aprioris.  (Units
*               - match partials so that partial*adjust gives L1
*               - cycles.
      real*8 mf_aprval(mf_maxprm), mf_adjust(mf_maxprm)
 
*   mf_alabel(mf_maxprm)    - Labels for the parameters
 
      character*20 mf_alabel(mf_maxprm)
 
*   mf_sitet(mf_maxcfl)    - Monument names for the sites in the
*               - cfiles
      character*12 mf_sitet(mf_maxcfl)
 
*   mf_rfname(mf_maxcfl)   - Names of cfiles used in this mfile
*   mf_tfname(mf_maxtfl)    - Names of the tfiles (multisession
*               - option (not coded currently).
 
      character*16 mf_rfname(mf_maxcfl), mf_tfname(mf_maxtfl)
 
      common / mf_type1 / mf_aprval, mf_adjust, mf_nversn, mf_ndy, 
     .    mf_nepch, mf_mtpart, mf_idms, mf_islot_all, mf_nsat, mf_isat, 
     .    mf_nsite, mf_nrfile, mf_ntfile, mf_norb, mf_alabel, mf_sitet,
     .    mf_rfname, mf_tfname
 
*.......................................................................
* BLOCK 2: Session information
 
*   mf_idy      - Session number
*   mf_mdy(3)   - Month, day, year of start of session
*  mf_nepch_sess    - Number of epochs in this session
*   mf_inter    - Interval between epochs (seconds)
*   mf_skip     - Solve decimation factor
*   mf_ncfile   - Number of cfiles in this session
*   mf_nsat_sess    - Number of satellites in this session.
 
 
      integer*4 mf_idy, mf_mdy(3), mf_nepch_sess, mf_inter, mf_skip,
     .    mf_ncfile, mf_nsat_sess
 
*   mf_hms(3)   - Hour,minute, second for start of session
*   mf_stawgh(mf_maxcfl)   - Station weights for this session
*   mf_satwgh(mf_maxsat)    - Satellite weights for this session.
 
      real*8 mf_hms(3), mf_stawgh(mf_maxcfl), mf_satwgh(mf_maxsat)
 
*   mf_cfname(mf_maxcfl)   - Names of the cfiles in this session.
 
      character*16 mf_cfname(mf_maxcfl)
 
      common / mf_type2 / mf_hms, mf_stawgh, mf_satwgh, mf_idy, mf_mdy,
     .    mf_nepch_sess, mf_inter, mf_skip, mf_ncfile, mf_nsat_sess,
     .    mf_cfname
 
*....................................................................
* BLOCK 3: Cfile information for this session (1-record per cfile)
*
*   mf_mdy_cf(3)    - First epoch for this cfile (same for all)
*   mf_kpart        - Number of partials in cfile (same for all)
*   mf_islot_cf(mf_maxprm)  - Slots in the cfile header (same for
*                   - all)
*   mf_numzen       - Number of zenith delays (assumed same for
*                   - all sites).
*   mf_idtzen(mf_maxatm)    - Tabular points for the zenith delays
*                   - (By epoch number)
* MOD TAH 980915: Updated gradient entries
*   mf_numgrad        - Number of gradient parameters
*   mf_idtgrad(mf_maxgrad)  -  Tabular points for the gradient 
*                     parameters (by epoch number).
 
      integer*4 mf_mdy_cf(3), mf_kpart, mf_islot_cf(mf_maxprm),
     .    mf_numzen, mf_idtzen(mf_maxatm), mf_numgrad,
     .    mf_idtgrad(mf_maxgrad)
 
*   mf_hms_cf(3)    - Start hrs, min and sec.
*   mf_elvcut(mf_maxcfl)   - Elevation cutoff by site. (degrees)
 
      real*8 mf_hms_cf(3), mf_elvcut(mf_maxcfl)
 
*   mf_zenmod       - Zenith delay mode (CON or PWL).
*   mf_gradmod      - Gradient delay mode (CON or PWL) 
 
      character*3 mf_zenmod, mf_gradmod
 
      common / mf_type3 / mf_hms_cf, mf_elvcut, mf_mdy_cf, mf_kpart,
     .    mf_islot_cf, mf_numzen, mf_idtzen, mf_numgrad,
     .    mf_idtgrad, mf_zenmod, mf_gradmod
 
*----------------------------------------------------------------------
