 
*     Include file containing the definition of the header records
*     for the Kalman filter global solution files.  The full layout
*     of the files is:
*
*     Global Header  Names list  solution records.. 
*     ^-----------^  ^---------^ ^--- ........ ---^ 
*
*     Station inf records    Satellite Inf Records
*     ^-----------------^    ^-------------------^
*
*     parameter types   parameter epoch range
*     ^--------------^  ^--------------------^
*
*     Aproiri types   Apriori epoch range^
*     ^-----------^  ^-------------------^  
*     Apriori values  
*     ^------------^
*
*     Covariance matrix and solution.
*     ^----------- .....  ----------^
*
*     where:
*     Global header contains the information given in this common
*         block.
*     Names list contains the site and source names compressed together
*         This block can take a variable number of records.
*     Solution records contains the descriptive information about each
*         of the KalObs files in this solution.  There is one of these
*         records for each file. (See sln_def.ftni). ** NOTE: The number of
*         records for the solutions records is assumed to equal the number of
*         solutions in the global **.  
*         In verion 1.0 there is an additional solution record
*         for each of the input combined global files.
*     Station inf records are station information for each of the solutions
*         in the combined binary. 
*     Parameter type record(s) contains the code numbers for each of the
*         paramters in this global solution.  The definitions of these
*         codes in given for varibable GLB_TYPES.
*     Parameter epoch range records constains times of validity of 
*         parameters (1.02 version)
*     Apriori types record(s) contains the code numbers for the site and
*         source positions which were held fixed in the SOLVK or GLOBK
*         solutions.  The code definitions are the same as those for the
*         parameters.
*     Apriori epoch range records constains times of validity of 
*         parameters (1.02 version)
*     Apriori values record(s) contains the aproiri values of the
*         fixed parameters in the same order as the codes are given.
*     Covariance matrix and solution records contain the solution and
*         covariance matrix from the SOLVK or GLOBK solutions.  The
*         real*8 values are saved as par_num*par_num covarinve matrix
*         values plus par_num solution vector values.  There is no
*         record gap between the matrix and the solution.
*
*     The codes for the parameter types are: (The values through
*     to 29 match those used in PLTSL)
*     Code   Meaning
*        7   Site X coordinate (site number in high 16 bits)
*        8   Site Y coordinate (site number in high 16 bits)
*        9   Site Z coordinate (site number in high 16 bits)
*       10   Axis offset (site number in high 16 bits)
*       11   Source right ascension (source number in high 16 bits)
*       12   Source declination (source number in high 16 bits)
*       13   X/Y pole position (high 16 bits contains:
*               1 -- X Pole value (mas)
*               2 -- Y Pole value (mas)
*               3 -- X random walk component (mas/year)
*               4 -- Y random walk component (mas/year)
*               5 -- X seasonal offset (mas)
*               6 -- X seasonal rate (mas/year)
*               7 -- Y seasonal offset (mas/day)
*               8 -- Y seasonal rate (mas/year)
*       14   UT1-AT (high 16 bits contain:
*               1 -- UT1-AT (mas)
*               2 -- Rate (mas/year)
*               3 -- annual Seasonal offset (mas)
*               4 -- annual Seasonal rate   (mas/year)
*               5 -- semiannual seasonal offset (mas)
*               6 -- semiannual seasonal rate   (mas/year)
*       15   Nutation angles (High 16 bits are the same as for wobble
*              with DEPS for X, and DPSI for Y)
*       16   Tide love number l
*       17   Tide love number h
*       18   Tide love lag angle
*       19   Extended earth tide parameter 1
*       20   Extended earth tide parameter 2
*       21   Extended earth tide parameter 3
*       22   Extended earth tide parameter 4
*       23   Extended earth tide parameter 5
*       24   Extended earth tide parameter 6
*       25   Extended earth tide parameter 7
*       26   Extended earth tide parameter 8
*       27   Extended earth tide parameter 9
*       28   Extended earth tide parameter 10
*       29   Extended earth tide parameter 11
*       30   Extended earth tide parameter 12
*       31   UT1 Diurnal Cosine
*       32   UT1 Diurnal Sine
*       33   UT1 semidiurnal Cosine
*       34   UT1 semidiurnal Sine
*       35   XY retrograde diurnal Cosine
*       36   XY Retrograde diurnal Sine
*       37   XY Prograde semidiurnal cosine
*       38   XY Prograde semidiurnal sine
*       39   XY retrograde semidiurnal cosine
*       40   XY retrograde semidiurnal sine
*       41   Gamma
*       42   Site position X coordinate rate
*       43   Site position Y coordinate rate
*       44   Site position Z coordinate rate
*       45   Source position right ascension rate
*       46   Source position declination rate
*       47   Nutation coefficent 1 (high 7 bits gives the coefficient
*               number. Bit 16 gives in (0) and out of phase (1)
*       48   Extended Earth tide coefficients: High 7 bits gives the
*               coefficient number.  Bit 16 gives in (0) and out of
*               phase (1).  High 16 bit word gives the station number
*               if site dependent, and 0 here if global tide.)
*       49   UT1 coeffients (high 7 bits give frequency, Bit 16 is
*               0 for cos, and 1 for sin)
*       50   XY coefficients (high 7 bits give frequency, Bit 16 is
*               0 for cos, and 1 for sin)
*       51   SV orbits (high 7 bits of low word gives the element, 
*                and high 16 bits gives the SV number.)
* CODES Changed 190606 to support ECOMC model which is combined 
*       ECOM and ECOM2 model.
*                Element codes ares:
*                     NEW (ECOMC)         (OLD)                    Map
*                1  - X POS             - X POS                    1
*                2  - Y POS             - Y POS                    2
*                3  - Z POS             - Z POS                    3
*                4  - X VEL             - X VEL                    4
*                5  - Y VEL             - Y VEL                    5
*                6  - Z VEL             - Z VEL                    6
*                7  - RAD PRES DIRECT   - RAD PRES DIRECT          7
*                8  - Y AXIS BIAS       - Y AXIS BIAS              8
*                9  - B AXIS BIAS       - Z AXIS BIAS (never used) 0
*               10  - COS DIRECT        - B AXIS BIAS              9
*               11  - SIN DIRECT        - X AXIS BIAS (never used) 0
*               12  - COS Y BIAS        - COS DIRECT               10
*               13  - SIN Y BIAS        - SIN DIRECT               11
*               14  - COS B BIAS        - COS Y BIAS               12
*               15  - SIN B BIAS        - SIN Y BIAS               13
*               16  - COS 2U DIRECT     - COS B BIAS               14
*               17  - SIN 2U DIRECT     - SIN B BIAS               15
* Added New radiation parmeter types (981020)                      0
*               18  - COS 4U DIRECT     - SIN X1 BIAS (never used) 0
*               19  - SIN 4U DIRECT     - SIN X3 BIAS (never used) 0
*               20  - <OPEN>            - SIN Z1 BIAS (never used) 21
*               21  - SVANT X OFF       - SVANT X OFF              22
*               22  - SVANT Y OFF       - SVANT Y OFF              23
*               23  - SVANT Z OFF       - SVANT Z OFF
*       52   GPS Bias parameter.  High 16 bit word gives the double
*               difference number, computed from Station and SV numbers
*       52   Translation parameters.  High 16 bits gives X Y or Z (as
*               1,2 or 3.
*       53   Translation rate parameters.  High 16 bits gives X Y Z 
*       54   Scale parameter
*       55   Scale rate parameter.
* MOD TAH 981020: Introduced new parameter number types for multi-epoch
*            polar motion and UT1.
*       56   Multi-epoch X-polar motion, epoch number in top I*2 part,
*            Bits 7-16 give 1-offset, 2 = rate
*       57   Multi-epoch Y-polar motion, epoch number in top I*2 part,
*            Bits 7-16 give 1-offset, 2 = rate
*       58   Multi-epoch UT1, epoch number in top I*2 part,
*            Bits 7-16 give 1-offset, 2 = rate
*       59   Generic rotation terms about X, Y, Z axes
*       60   Generic rotation rate terms about X,Y,Z axes
*       61   Atmospheric Zenith days.
* MOD TAH 030615: Add log parameters
*       62   Log N component (m)
*       63   Log E component (m)
*       64   Log U component (m)
   
 
*   cglb_dcb(16)    - DCB buffer for reading and writing global
*                   - file.
 
      integer*4 cglb_dcb(16)
 
*   glb_header      - First word in the glb header file
 
*   cnum_header     - Num records for header block
*   cnum_names      - Num records for names block
*   cnum_full       - Num records for the full names block 
*   cnum_soln_recs  - Num records for solution records
*   cnum_par_types  - Num records for parameter types
*   cnum_apr_types  - Num records for apriori types
*   cnum_apr_vals   - Num records for apriori values
*   cnum_par_vals   - Num records for parameter
*                   - covariance matrix and solution vector
 
*   crec_names      - Record number for start of names block
*   crec_full       - Record number for full names block
*   crec_solutions  - Record number for start of solution records
*   crec_par_types  - Record number for start of parameter types
*   crec_apr_types  - Record number for start of apriori types
*   crec_apr_vals   - Record number for start of apriori values
*   crec_par_vals   - Record number for start of parameter
*                   - covariance matrix and solution vector
 
      integer*4 glb_header, cnum_header, cnum_names, cnum_full,
     .    cnum_soln_recs,
     .    cnum_par_types, cnum_apr_types, cnum_apr_vals, cnum_par_vals,
     .    crec_names, crec_full, 
     .    crec_solutions, crec_par_types, crec_apr_types,
     .    crec_apr_vals, crec_par_vals
 
*   cfile_type      - Contains the string 'GLOBAL' so that we can
*                   - tell this is a global file.
 
      character*6 cfile_type
 
*   cexpt_title     - Title of the global file
 
      character*62 cexpt_title
 
*   cnum_apr_codes  - Number of apiori values saved in this
*                   - solution
*   cnum_parn       - Number of parameters in this global solution
*                   - file.
*   cnum_sites      - Number of sites in this global solution
*   cnum_sources    - Number of sources in this global solution
*   cnum_svs        - Number of satellites in this gobal solution
*   crun_time(7)    - Runtime of this global solution.
*   ctai_utc        - TAI - UTC for last experiment in solution
*   cchisq          - Chi squared for solution
 
      integer*4 cnum_apr_codes, cnum_parn, cnum_sites, cnum_sources,
     .    crun_time(7), ctai_utc, cnum_svs
     
      real*4 cchisq
 
*   cdelete_count(max_edit_types)  - Accumulated sum of the
*                   - delete counts for this global soltion
*   cnum_obs        - Total number of observations in this global
 
      integer*4 cdelete_count(max_edit_types), cnum_obs
 
*   cepoch_end      - Epoch of last experiment in solution (JD)
*   cepoch_expt     - Epoch to which the solution parmeter values
*                   - are referred (Mainly for when rates are
*                   - estimated (JD))
*   cepoch_start    - Epoch of first experiment in solution (JD)
 
*   cetd_apr(3)     - Aprior values for lag, h and l.
*   cnut_ang_apr(8) - apriori value for nutation angles and
*                   - random walk and seasonal model (mas and mas/year)
*   cut1_apr(6)     - Apriori value of UT1 and its rate, and its
*                   - seasonal model (mas and mas/year)
*   cwob_apr(8)     - Apriori value of wobble and random walk and
*                   - seasonal model. (mas and mas/year)
*   csvs_epoch      - Ephemeris epoch for this experiment.  (Last
*                     on for multiple day solutions)
 
      real*8 cepoch_end, cepoch_expt, cepoch_start, cetd_apr(3),
     .    cnut_ang_apr(8), cut1_apr(6), cwob_apr(8), csvs_epoch 

* MOD TAH 950824: Added new variables to support writing SINEX format
*     files.

*   cglb_vers -- Version number of binary h-file by 100.
*   cnum_sinf -- Total number of records in the station info records
*                (two stations per record, and each global uses an
*                 integer number of records)
*   crec_sinf -- Record number of the first station information record.
*                (The individual solution records contain the record
*                 numbers for there station information).
*   cnum_acvc      -- Number of apriori constraints saved at the end of
*                 file.
*   cnum_acvc_recs -- Number of records in the apriori covariances for 
*                sites.
*   crec_acvc -- Record number for the start of the apriori covariances
*   cgpst_utc -- Difference between gpst minus utc for this solution
*   crec_comb_soln - Number of soln records that are for combined solutions
*                   (rather than single day solutions)
*   cnum_comb      - Number of combined solution records in this file.
*   ccons_type     - Indicates the strengths of the constraints (0--tight,
*                    1--signficant, 2--loose).
*   csys_type      - Type of system type used: Bit mapped
*                    1 -- VLBI
*                    2 -- GPS
*                    3 -- SLR

      integer*4 cglb_vers, cnum_sinf, crec_sinf, cnum_acvc, 
     .          cnum_acvc_recs, 
     .          crec_acvc, cgpst_utc, crec_comb_soln, cnum_comb, 
     .          ccons_type, csys_type
      
*   cowner    -- Owner of this file (4 characters)
*   ccreator  -- Creater of this binary file (4 character)
*   cprog_gen -- Generating format type (=GLK for globk binary files)
*   canal_type -- Type of analsys ('X V E')

      character*4 cowner, ccreator, cprog_gen
      character*8 canal_type
*
*   Additional GAMIT variables
*   cgframe  - Frame used in gamit
*   cgprec   - Precession model
*   cgsrpmod - Solar Radiation model
*   cgtime   - Time system used
*   cgnut    - Nutation model used
*   cggrav   - Gravity field
*   ceradmod  - Earth-radiation model (Added 140327) 
*   cantradmod  - Antenna-radiation mode (Added 140327) 
*   cionsrc   - Source of ionospheric delay; if not NONE then 2nd order
*               ioncorrection has been applied (140403)
*   cmagfield - Magnet field model used with 2nd order correction (140403)


*   cdryzen   - Dry Zenith delay model (140403)
*   cwetzen   - Wet Zenith delay model (140403)
*   cdrymap   - Dry Mapping function (140403)
*   cwetmap   - Wet Mapptng function (140403)

      character*8 cgframe, cgprec, cgsrpmod, cgtime, cgnut, cggrav,
     .            ceradmod, cantradmod, cionsrc, cmagfield

      character*4 cdryzen, cwetzen, cdrymap, cwetmap
      
* MOD TAH 981020: New record pointers for the epoch range variables

*   cnum_par_ep  -- Number of records needed to save the parameter
*               epoch range
*   cent_par_ep  -- Number of entries in the parameter epoch array
*   crec_par_ep  -- Start of record numbers parameter epoch numbers

*   cnum_apr_ep  -- Number of records needed to save the apriori
*               epoch range
*   cent_apr_ep  -- Number of entries in the apriori epoch array
*   crec_apr_ep  -- Start of record numbers apriori epoch numbers
*   cgamit_mod   -- Bit mapped word for gamit tide and E-rot models
*                   The e-tide model is mapped into the upper 16-bit
*                   words.  E-rot model is in lower 16 bits.  Value
*                   taken from sln_def.h, if not present.
*                   See BIT definition in ggamit_mod globk_cntl.h.
* MOD TAH 110512: Added MPT2010 (bit 23) for IERS2010 pole tide model


      integer*4 cnum_par_ep, cent_par_ep, crec_par_ep, 
     .          cnum_apr_ep, cent_apr_ep, crec_apr_ep,
     .          cgamit_mod

* MOD TAH 050622: Added new record type (Version 1.05 binary files)
*   cload_mod  -- Bit mapped loading models applied in gamit.  Models split 
*                 into 8 bit blocks: 
*                  1- 8 - otidemod;
*                        
*                  9-16 - atmloadmod;
*                      9 -- 1  Atm load applied
*                 17-24 - atmtidemod;
*                     17 -- 1  Atm S1/S2 tide applied
*                 25-32 - hydro load mod;
*                     25 -- 1  Hydrographic load applied
*   cnum_svinf -- Total number of records satellite information records.
*                 Four records paper 128*word record are saved
*   crec_svinf -- Record number of the first satellite information record.
*                 (If this value is zero then there are no satellite records)

      integer*4 cload_mod, cnum_svinf, crec_svinf

*   cspeopmod   -- Name of short period eop model
*   cetidemod   -- Name of earth tide model
*   cotidemod   -- Name of ocean tide model
*   coatmlmod   -- Name of atmospheric load model
*   catmtdmod   -- Name of atmospheric tidal loading model
*   chydromod   -- Name of hydroloading model.

      character*8 cspeopmod ,cetidemod, cotidemod, coatmlmod, 
     .            catmtdmod, chydromod

*   last_glb_header_wrd - Last word in the global file header
*                   - block
*   Glb_header_dummy(127)   - Padding for end of header section.
*       MOD TAH 900912: Used two words of header dummy for csvs_epoch
*       MOD TAH 950824: Used thirteen words for version extension information
*       MOD TAH 950919: Used additional 3 words for cprog_gen and canal_type
*       MOD TAH 950920: Used 4 character*8 variables for cgframe, cgprec, 
*                       cgsrpmod, cgtime
*   MOD TAH 960109: Because of change in the number of blocks used to
*     write the header records; the current required size of Glb_header_dummy
*     is 109.  The variable has keep it dimension at 127 to ensure that
*     reading the header does not overwrite any memory:
*   MOD TAH 981020: Used 6 more words and so size of block is now 103 words.
*   MOD TAH 981112: Used 1 word for cgamit_mod.  Size now 102 words.
*   MOD TAH 050622: Used 3 words for cload_mod, cnum_svinf, crec_svinf. 
*                   12 words for cspeopmod ,cetidemod, cotidemod, coatmlmod, 
*                   catmtdmod, chydromod
*                   For a total of 15 words: Size now 87 Words
*   Checked TAH 071227: Actually only 72 words remaining.
*   MOD TAH 080102: Added cgnut and cggrav (leaving 68 words)
*   MOD TAH 140327: Added ceradmod and cantradmod (leaving 64 words)
*   MOD TAH 140403: Added cionsrc and cmagfield (leaving 60 words)
*   MOD TAH 140403: Added cdryzen, cwetzen, cdrymap, cwetmap all C*4
*                   (leaving 56 words)

*     Glb_header_dummy(64) -- Current required size.
 
      integer*4 last_glb_header_wrd, Glb_header_dummy(127)
 
*----------------------------------------------------------
*     Common block declation
 
      common / GLB_header_block / cglb_dcb, glb_header, cnum_header,
     .    cnum_names, cnum_full, 
     .    cnum_soln_recs, cnum_par_types, cnum_apr_types,
     .    cnum_apr_vals, cnum_par_vals, crec_names, crec_full, 
     .    crec_solutions,
     .    crec_par_types, crec_apr_types, crec_apr_vals, crec_par_vals,
     .    cfile_type, cexpt_title, cnum_apr_codes, cnum_parn,
     .    cnum_sites, cnum_sources, cchisq, cnum_svs, 
     .    crun_time, ctai_utc, cdelete_count,
     .    cnum_obs, cepoch_end, cepoch_expt, cepoch_start, cetd_apr,
     .    cnut_ang_apr, cut1_apr, cwob_apr, csvs_epoch,
     .    cglb_vers, cnum_sinf, crec_sinf, cnum_acvc, cnum_acvc_recs, 
     .    crec_acvc, cgpst_utc, crec_comb_soln, cnum_comb,
     .    ccons_type, csys_type,
     .    cowner, ccreator, cprog_gen, canal_type,
     .    cgframe, cgprec, cgsrpmod, cgtime, 
     .    cnum_par_ep, cent_par_ep, crec_par_ep, 
     .    cnum_apr_ep, cent_apr_ep, crec_apr_ep, cgamit_mod,
     .    cload_mod, cnum_svinf, crec_svinf, cspeopmod, cetidemod,
     .    cotidemod, coatmlmod, catmtdmod, chydromod, cgnut, cggrav, 
     .    ceradmod, cantradmod, cionsrc, cmagfield, 
     .    cdryzen, cwetzen, cdrymap, cwetmap,
     .    last_glb_header_wrd,
     .    Glb_header_dummy
     
