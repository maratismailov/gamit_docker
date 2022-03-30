 
*     This file contains the definition of the record which describe
*     nature of each of the solutions in the global solution.

*MOD TAH 941018: Changed the number stations for the secc_change
*     back to 8 for backwards compatability in the global files.
*     No software actually uses these entries and they have never
*     been set GPS global files. 
* MOD TAH 950826: Added information for SINEX compatability

      integer*4 max_sln_sites

      parameter ( max_sln_sites = 4 )
 
*   sdelete_count(max_edit_types)  - Delete counts for this
*               - experiment
*   snum_obs    - number of observations in this solution
*   stai_utc    - Value of tai minus UTC for this eperiment
*               - (Used to compute UT1-UTC)
*   sversion    - Version of data base for this solution
*   snum_parn   - Number of paramters in solutions (>= ver 1.0)          
 
      integer*4 sdelete_count(max_edit_types), snum_obs, stai_utc,
     .    sversion, snum_parn
 
*   sepoch      - Mid epoch for this experiment
*   ssvs_epoch  - Epoch for the initial conditions of the satellite
*                 orbits in this solution.
*   sut1_apr(2) - apriori values of UT1-AT for this experiment
*   swob_apr(2,2)   - Aproiri value for x and y pole position
*   snut_ang_apr(2,2)   - Apriori values for nutation angles
 
      real*8 sepoch, ssvs_epoch, sut1_apr(2), swob_apr(2,2), 
     .       snut_ang_apr(2,2)

*   secc_change(3,max_sln_sites) - Eccentricity change for the sites
*               - in this experiment.

      real*4 secc_change(3,max_sln_sites)
      
* sgpst  - Time offset between GPS time and UTC

      integer*4 sgpst_utc

* sgframe  - Frame used in gamit
* sgprec   - Precession model
* sgsrpmod - Solar Radiation model
* sgtime   - Time system used

      character*8 sgframe, sgprec, sgsrpmod, sgtime

* sprog_gen -- Name of file type from which this solution came.  Normally
*              will be SNX or GLK.  The first character is + or = if
*              added or combined global.
      character*4 sprog_gen

* srun_time -- Time solution was generated

      real*8 srun_time
 
*   sdata_base  - Name of the data base for this experiment.  With
*                 version 1.00+ files this entry is the SNX program
*                 name +version
 
      character*10 sdata_base
 
*   sKalObs_file    - Name of the KalObs file for this experiment
 
      character*128 sKalObs_file
           
*   sanal_type      - Analysis type (X E etc.) Last character is
*                     enttry (+ or =) from Sinex File record
*                     (For version 1.0 and greater)

      character*6   sanal_type
      
* Currently 104 I*4 are used, leaving 24 I*4 words (96 bytes free).

*----------------------------------------------------------------
* New entries added to support SINEX.

*   sepoch_end  -  End epoch of data in this solution
*   sepoch_start - Start epoch of data in the solution

      real*8 sepoch_end, sepoch_start

*   sglb_vers  -  Version of global file from which this record
*                 came.  (If less than 100, then no station information
*                 records are avaiable)
*   srec_inf   -  Starting record for station information record
*   snum_inf   -  Number of records in the station information.
*   scons_type -  Indicates strength of constraints in solution
*                 (0--tight, 1--significant, 2--loose)
*   ssys_type  -  Bit set for each type of system used in the solution
*                 Bit 1 -- VLBI (R)
*                 Bit 2 -- GPS  (P)
*                 Bit 3 -- SLR  (L)
*                 Bit 4 -- LLR  (M)

*   snum_sites -  Number of sites in this solution.

      integer*4 sglb_vers, srec_sinf, snum_sinf, scons_type, 
     .          ssys_type, snum_sites

*   sowner    -  Owner of the file from which this solution came
*   screator  -  Creator of the file from which this solution came

      character*4 sowner, screator
*
*   sexpt_title -  A Description of this experiment

      character*32 sexpt_title 

*   sgamit_mod  - Bit mapped word with GAMIT tide and E-rot models.

      integer*4 sgamit_mod

* SINEX extension: Used 20 I*4 words, leaving us with 4 i*4 words 
*     free 
* MOD TAH 980514 : Used  1 I*4 words, leaving us with 3 i*4 words
*                  (sgamit_mod variable introduced).
*----------------------------------------------------------------
*   sln_dummy(127)  - Padding at end of record
 
      integer*4 sln_dummy(127)
 
*-----------------------------------------------------------------
*     Common Declaration
 
      common / sln_block / sdelete_count, snum_obs, stai_utc, sversion,
     .    snum_parn, secc_change, sgpst_utc,
     .    sgframe, sgprec, sgsrpmod, sgtime,
     .    sprog_gen, srun_time, 
     .    sepoch, ssvs_epoch, 
     .    sut1_apr, swob_apr, snut_ang_apr, sdata_base,
     .    sKalObs_file, sanal_type, 
     .    sepoch_end, sepoch_start,
     .    sglb_vers, srec_sinf, snum_sinf, scons_type, ssys_type,
     .    snum_sites, sowner, screator,
     .    sexpt_title, sgamit_mod,
     .    sln_dummy

*---------------------------------------------------------------- 
