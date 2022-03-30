 
 
*     The htoglb_comm.h include file.  This file is needed for running
*     the H file to global file conversion.

*-----------------------------------------------------------------
* MOD for Version 3.21: 941005
*     Three main changes were made for this version:
*     (1) get_svs_inf routine added because of complex treatment needed
*         to fix a bug in Solve 9.28 in which the satellite information
*         (PRN and number of measurements) was split over multiple lines
*         if more than 24 satellites were used.
*     (2) Added code to allow JPL stacov matrices to be read and converted
*         to globk comapatible binary-hfiles
*     (3) Moved the version number to this file (from kalman_param.h) and
*         moved this file to the local directory.  (Allows stand-alone 
*         updating of the file.
* MOD for Version 3.22: 941212
*     (1) Add capability to read SLR(GSFC) covariance matrix files.
* MOD for Version 3.23: 941230 
*     (1) Fixed bug in skip_htoh so that correction position in file
*         is returned
*     (2) Updated for new Gamit hfile format.
* MOD for Version 3.24: 950112
*     (1) Added reading GSFC VLBI covariance files.
* MOD for version 3.25: 950301
*     (1) Added reading CTOGOBS Kalman filter hfiles.
* MOD for version 3.26: 950729
*     (1) Added reading extra radiation parameters from GAMIT h-files
*     (2) Changed the name generation of the binary hfiles to include
*         4-char code from the hfile.
* MOD for version 3.27: 950829
*     (1) Added additional variables for SINEX compatability using the
*         version 1.00 binary files
*     (2) Allowed more robust treatment of parameters in the input files
*         which htoglb does not know how to interprett.
* MOD for version 4.00:950901 
*     (1) Full comparablity with version 1.0 global files
*     (2) Upgraded reading sinex files and new h-file formats
* MOD for version 4.01: 951014
*     (1) Added features and variables to support glbtosnx.  Fixed
*         some small problems with reading sinex and hfiles.
* MOD for version 4.02: 960225
*     (1) Add check in dec_site_id to check for duplicate names.
*         When a duplicate is found the occ number is incremented.
*     (2) Add report_stat calls.
* MOD for version 4.03: 950430
*     (1) Made compatible with version 1.0 of SINEX.  Mainly changes
*         to the size of the decscription arrays for receivers and
*         antennas and the parameter name (now CH*6)
* MOD for version 4.04: 960618
*     (1) Modified reading the correlation matrix for Sinex 1.0 which
*         puts the sigma in the diagonal elements.
* MOD for version 4.05: 960808
*     (1) Changed to decoding of the correlation matrix entries to fix
*         problems introduced in 4.04
*     (2) Allowed for a second reading of the sinex file to get the
*         aproiri values if these are given before the estimates. 
*         (Parameter list based on estimates is needed to decode apriori).
*  4.06 - 960926: Changed the reading of the antenna phase model
*                 names so that the IGS_01 name is encoded to the
*                 4 character name.
*  4.07 - 961129: Changed reading of translation and translatation rate
*                 plus added scale and scale rate.
*  4.08 - 970312: Changed the handling of the antenna and reciver names to
*                 getmore characters. 
*  4.09 - 970716: Checked if pt is blank, and replaced with ' A' if it is. 
*  4.10 - 971122: Added checking eigenvalues in make_snx when constraints
*                 are removed from the matrix.  Added routine jacobi.f from
*                 numerical recipes.
*  4.11 - 980123: Added user specifiable memory size; and extents based on
*                 type of system (D=DOR, S=SLR, R=VLB, G/C=GPS)
*  4.12 - 981020: Added muliple polar motion parameters

*  5.00 - 990101: Increased version for new version of globk
*  5.01 - 990315: Added extra trimlead's to string to make sure names are
*                 left justified in fields.  Also checked the deconstraint
*                 sigma so that rotation loosening is only applied if the
*                 SINEX file is not being deconstrained.
*  5.02 - 991131: Added optional naming feature (GPSW, Y2DOY, Y4DOY and Y4),
*                 fixed problem with length of sn for receiver and antenna.
*  5.03 - 991206: Added option not to check eigenvalues (comp_eigen)
*  5.04 - 000308: Fixed problem with end time being interpretted as not given
*                 because year was set to zero. (Added additional check of
*                 components of the date being zero).
*  5.05 - 020624: Better error checking when number od solutions records is
*                 large
*  5.06 - 040510: Added -D=TR option to runstring to apply 10 mas and/or 1 m
*                 rotation/translation deconstraint.
*  5.07 - 040703: Added reading of atmospheric delays from the hfiles.
*  5.08 - 040801: Added support for 64-bit memory machines with addressing using
*                 integer*8.
*  5.09 - 050214: Modified to read hfiles generated by SOLVE 10.17 and greater.
*  5.10 - 050516: Added -F=<string> feature to allow control of expt code
*                 included in hfile name
*  5.11 - 050616: Increased size of apriori arrays to avoid over flows

*  6.00 - 050622: Added satellite feautures to support version 2.01 of SINEX.  
*                 Added svinf_def.h declaration of satellite information.  Changed
*                 version of binary Hfiles.
*  6.01 - 070120: Added -t=AC options to remove constraints on satellite phase 
*                 center estimates and center of mass estimates in sinex files
*                 Updared phase center model code (SMod_to_full)
*  6.02 - 070504: Fixed code to allow multiple SINEX files to be read.  Problem
*                 was with keepinf track of which SOLUTION blocks had been read
*  6.03 - 080103: Updated for nutation and gravity field models being passed
*  6.04 - 080416: Fixed bug due to atmospheric delay parameters in reading
*                 satellite apriori covariance matrix.
*  6.05 - 090504: Minor change to way CODE + OCC + PT is converted to 8-character
*                 gsite_name for output
*  6.06 - 100511: Introduced -s option which translate PT code directly into 
*                 8-character name.  PT 001 _1PS, 002 _2PS, 010 _APS etc
*                 If the OCC code is not A, then PS is changed to BS, etc; 
*                 if name is non-unique, it is incremented until unique
*                 Also introduced -v to compute eigenvalues (-e removed)
*  6.07 - 101015: Increased qrecv_ty to C*20 and allowed long antex file 
*                 name. 
*  6.08 - 140327: Added reading of eradmod and antradmod model names from h-file.
*  6.09 - 150501: Updated to handle svnav.dat version 2.0 and h-file version 3.0 
*                 (gamit_hf_ver/100)
*  6.10 - 171003: Updated to handle CATREF IGS SINEX files with no Soln entries in some
*                 blocks with values in other block (caused duplicate sites in older
*                 versions).
*  6.11 - 180202: Added default feature to apply -d=TR to GAMIT BASELINE mode h-files
*                 Added option (-autooff/-a*) to stop this application and run like
*                 pre-6.11 versions of HTOGLB.  Add F to -D=option to force even in
*                 non-BASELINE model solutions..
*  6.12 - 180401: Updated to handle GNSS processing.  Satellite name format changed
*                 {GREC}<SVM>_(PRN> where GREC system, <SVN> is 3-digit satellite 
*                 vehicle number and <PRN> is two digit PRN (internally stored with
*                 factors of 100 added for sys type (see conoff and offcon functions)
*  6.13 - 190627: Added svi_lnum and ss_lnum information so see if site/satellite 
*                 actually used in solution.  Need for writing SINEX from BACK GLX
*                 files generated by glbak.
*  6.14 - 191004: Same change to csec value when station estimates read to make them
*                 nearest 10-minites of hour (To handle GYPSY-X 150 second time tag
*                 offset).
*  6.15 - 200222: Updated to set gamit_mod bits for IERS20 secular pole and antenna
*                 azimuth information.
*  8.16 - 200727: Modified to remove additional days from COD 3-day orbit SINEX files.
*  6.16 - 210112: Added -L=<Var> option to limit the variance of site written to the 
*                 SINEX file. (Poorly determined sites with small eigenvalues in 
*                 normal equations can become negative).

* VERSION INFORMATION:

*   htoglb_version  - Version number for htoglb
      character*(*) htoglb_version
      parameter ( htoglb_version = '6.16' )

*--------------------------------------------------
*   glbtosnx_version  - Version number for glbtosnx
*  1.02 - Changed method of printing file names so that
*         end of string would show if it is longer than
*         Sinex field.
*  1.03 - corrected =SNX line output for the last global 
*         so that verions 0.05 is printed.
*  1.04 - Changed format in the solution block so that Sun 
*         version will not extend past columns.
*  1.05 - 960809: Updated to write Sinex Version 1.0.
*  1.06 - 960926: Changed the reading of the antenna phase model
*                 names so that the IGS_01 name is decoded from
*                 4 char variables.
*  1.07 - 961129: Similar mods to 4.07 of htoglb to account for
*                 scale and translation
*  1.08 - 970312: Changed the handling of the antenna and reciver names to
*                 getmore characters. 
*  1.09 - 980209: Fixed bug associated with small files with no
*                 orbit parameters estimated (amount of memory to read
*                 the aprioris needed to be increased).
*  1.10 - 980218: Improved handling of null and non-printable characters
*                 in binary files (for use with old binary hfiles).
*  1.11 - 981020: Added multiple polar motion parameters
*  1.12 - 990315: Improved some of the outputing to get as much information
*                 in sinex file as possible.  Re-structured libaries for
*                 use with program hfupd.
*  1.13 - 000122: Y2K updates to fix 00 instead of 100 be written into
*                 files.
*  1.14 - 020624: Better error checking when number od solutions records is
*                 large
*  1.15 - 030201: Changed automatic naming to match day-of-week convention
*                 of IGS.  For multiday files, automatic name is generated 
*                 with 7.
*  1.16 - 040801: Modified memory allocation to be 64-bit compatable.
*  1.17 - 050616: Imcreased size of apriori arrays to aviod overflows

*  2.00 - 050622: Added satellite feautures to support version 2.01 of SINEX.  
*                 Added svinf_def.h declaration of satellite information.  Changed
*                 version of binary Hfiles.
*  2.01 - 060501: Modified name_to_id to better output pt and occ codes for renamed
*                 sites.  Specifically XPS sites get an X for PT. 
*  2.02 - 090504: More modifcation to all unique generatation of 8-character site
*                 name to code+occ+pt.  (Lists of occ and pt now saved).
*  2.03 - 100511: Introduced -s option which translates 8-character name back to
*                 PT code directly.  EQ-end characters are ignored.
*  2.04 - 101015: Added long site antex file name and receiver type to C*20
*  2.05 - 110912: Added miscellaneous parameter class so that normal equations can
*                 decode when they contain unknown parameter types
*  2.06 - 130321: Added secondary comment file with + name in comment file name
*                 i.e., head.snx+head.pbo will read +SNX/-SNX comments from head.pbo
*  2.07 - 180401: Added to GNSS sateliite information output
*  2.08 - 190627: Added svi_lnum and ss_lnum testing to see if site/SVS actually used
*                 solution.
*  2.09 - 210112: Added -V=<Var> option to limit the variance of site written to the 
*                 SINEX file. (Poorly determined sites with small eigenvalues in 
*                 normal equations can become negative).

      character*(*) glbtosnx_version
      parameter ( glbtosnx_version = '2.09' )

* Values based on glsave_comm.h
*   max_apr_saves   - Maxium number of apriori values which
*                   - may need to be saved
      integer*4 max_apr_saves

* MOD TAH 090622: Inccreased size by max_svs_elem*max_glb_svs to account for 
*     satellites.

      parameter ( max_apr_saves = 8*max_glb_sites + 4*max_glb_sources +
     .                            6*max_mul_pmu + 
     .                            max_svs_elem*max_glb_svs)


*-----------------------------------------------------------------
 
*   llr(3,max_glb_sites)    - apriori Lat, long, radius of the
*                       - sites
*   site_pos(3,max_glb_sites )  - XYZ coordinates of sites.(m)
*   site_vel(3,max_glb_sites) - XYZ velocities of the sites (m/yr)
*   site_axo(max_glb_sites) - Apriori value for axis offset (m)
*   site_atm(max_glb_sites) - Apriori values for atmospheric delays (m)
*   source_pos(2,max_glb_sources) - Radio Source positions (mas)
*   pmu_pos(2,3)        - Apriro values of polar motion and UT1
*                         offset and rate read from the parameter
*                         list
*   svs_pos(max_svs_elem,max_glb_svs)  - 6 orbital elements XYZ and
*                       - XYZ dot (m and m/s) plus three
*                       - solar radiation dependent bias
*                       - terms. (Ratio of effect to direct
*                       - solar radiation).  There are 11 radiation
*                         parameter types.
*   qstart_epoch, qend_epoch    - Start and stop times of data
*   qmid_epoch          - Mid-point epoch time.
*   qscale(max_glb_parn)    - Scaling factors from GPS coordinates
*                       - to GLOBK units.  (For sites also
*                       - includes conversion from radians to
*                       - NEU displacements at surface)
*   rot_mat(3,3, max_glb_sites) - Rotation matrices from NEU
*                       - XYZ for each site
 
*   qapr_list(max_glb_parn) - Apriori values for those
*                       - parameters which were not estimated.
*   qsol(max_glb_parn)  - Solution vector read from input file.
*                       - These need to be scaled and, for
*                       - sites rotated.
*   qsig(max_glb_parn)  - Sigmas for solutions.  Read from input
*   qref_ep(max_glb_parn) - Reference epochs for the parameters
*                         (These need to be moved to the solution
*                          reference epoch if possible)
*   qapr(max_glb_parn)  - Apriori values for solution
*   qapr_sig(max_glb_parn) - Apriori sigmas


*   qmul_par_ep(max_glb_parn) - Epoch list for multiday parameter 
*                           estimates
*   qmul_apr_ep(max_glb_parn) - Epoch list for multiday apriori 
*                           estimates

*   snx_constraint      - Constraint level to be applied to SNX files
*                         if more heavily constrained
 
 
      real*8 llr(3,max_glb_sites), site_pos(3,max_glb_sites ),
     .    site_vel(3,max_glb_sites), site_axo(max_glb_sites), 
     .    site_atm(max_glb_sites),
     .    source_pos(2,max_glb_sources), pmu_pos(2,3),
     .    svs_pos(max_svs_elem,max_glb_svs), 
     .    qstart_epoch, qend_epoch, qmid_epoch, 
     .    qscale(max_glb_parn), rot_mat(3,3, max_glb_sites),
     .    qapr_list(max_glb_parn), qsol(max_glb_parn),
     .    qsig(max_glb_parn), qref_ep(max_glb_parn),
     .    qapr(max_glb_parn), qapr_sig(max_glb_parn),
     .    snx_constraint, 
     .    qmul_par_ep(max_glb_parn), qmul_apr_ep(max_glb_parn) 

*   qnut_ang_apr(2,2) - Apriori values for the nutation angles
*   qut1_apr(2)       - Apriori vlaues for UT1
*   qwob_apr(2,2)     - Apriori values for wobble parameters

      real*8 qnut_ang_apr(2,2), qut1_apr(2), qwob_apr(2,2)

*   qtai_utc

      integer*4 qtai_utc
 
*   qnum_sites  - Number of sites in the Q file
*   qnum_sources - Number of radio sources
*   qnum_svs        - Number of satellites in the Q file.
*   qnum_parn   - number of estimated parameters in the solution
*   qnum_misc   - number miscellaneous parameters
*   qparn_sites(3,max_glb_sites)    - Parameter numbers for each
*               - site.
*   qparn_vel(3,max_glb_sites)  - Paramater numbers for the 
*                 velocities at each site.
*   qparn_axo(max_glb_sites)    - Parameter numbers for Axis offset.
*   qparn_atm(max_glb_sites)    - Parameter numbers for atm delays
*   qparn_sources(2,max_glb_sources) - Parameter numbers for radio sources
*   qparn_svs(max_svs_elem,max_glb_svs)    - parmaeter numbers for
*               - satelllites
*   qparn_pmu(2,3)              - parameter numbers for polar motion
*                 and UT1 offset and rate
*   qparn_tran(3,2)  - parameter numbers for translations and rate
*   qparn_scale(2)   - parameter numbers for scale and rate
*   qparn_misc(max_glb_parn)    - Miscellaneous parameters that we
*                      know nothing about
*   qglb_codes(max_glb_parn)    - codes for the parameters in the
*               - solution
*   qapr_codes(max_apr_saves)    - codes for apriori values which
*               - were not estimated. Use the max_apr_saves values
*                 from glsave.
*   qnum_apr_codes  - Number apriori codes
*   qnum_par_codes  - Number of parameter codes
*   qnum_soln_recs  - Number of solution records including the
*                     entries for combined solutions
*   qnum_comb       - Number of combined solution records (set to
*                     1 if there are more than data set in the solution)
*   qrun_time(7)    - Run time for the GPS solution
*   rad_num         - Number of radiation parameters in current
*                     hfile.
*   svant_num       - Number of antenna offset parameters in hfile.
*   atm_num         - Number of atmospheric delays in hfile
*   itoo(max_glb_parn)  - Mapping from the input h-file parameter
*                     number to its output number (0 if parameter is
*                     not to be output)
*   atoo(max_glb_parn)  - Mapping of the apriori parameter number
*                     to output.  For GAMIT hfiles atoo is same as itoo,
*                     for sinex they can be different.
*   atos(max_glb_parn)  - Mapping of apriori parameter number to 
*                     site number.  These are the only constraints we
*                     try to remove.
*   num_itoo, num_atoo  - Number of values in the above lists.
* MOD TAH 190627: Added gtol_sites and gtol_svs arrays to save mapping
*    from global to local sites/satellites
      integer*4 gtol_sites(max_glb_sites) ! Mapping from global to local
*                     site number (-1 means site not used)
      integer*4 gtol_svs(max_glb_svs)     ! Mapping from global to local
*                     satellite number (-1 means satellite not used).
      integer*4 qnum_sites, qnum_sources, qnum_svs, qnum_parn,
     .    qnum_misc, qparn_sites(3,max_glb_sites), 
     .    qparn_vel(3,max_glb_sites), qparn_axo(max_glb_sites), 
     .    qparn_atm(max_glb_sites), 
     .    qparn_sources(2,max_glb_sources), 
     .    qparn_svs(max_svs_elem,max_glb_svs),
     .    qparn_pmu(2,3), qparn_tran(3,2), qparn_scale(2),
     .    qparn_misc(max_glb_parn), 
     .    qglb_codes(max_glb_parn), qapr_codes(max_apr_saves),
     .    itoo(max_glb_parn),  
     .    atoo(max_apr_saves), atos(max_apr_saves),
     .    num_itoo, num_atoo,
     .    qrun_time(7), qnum_apr_codes, qnum_par_codes, 
     .    qnum_soln_recs, qnum_comb, rad_num, svant_num, atm_num

* MOD TAH 200727: Added logical array rm_param to denote parameters to
*     be removed from normal equations.
      logical rm_param(max_glb_parn) ! Set true to remove parameters from
              ! day different from center day of SINEX
     .,       rm_site(max_glb_sites) ! Set true to remove site name for 
              ! sites on extra days.
      logical sngday  ! Set true for single day sinex only. (default true,
              ! Set false with -mday option).

* MEMORY Management

* usr_vma_space  -- User assigned memory space (given with -m <size> 
*      option in runstring.
* istart_vma     -- Address in vma_data to start putting data

      integer*4 usr_vma_space 
      integer*8 istart_vma
     

* ---------------------------------------------------------
* New  variables for reference frames:

*   qgpst_utc  -  Time difference between GPST and UTC 
*   qcons_type -  Type of constraint used in this globk run. 
*                 Based on station constraints.  0--tight constraints;
*                 1--signficant constraint; 2--loose
*   qsys_type      - Type of system type used: Bit mapped
*                    1 -- VLBI
*                    2 -- GPS
*                    3 -- SLR
*   solve_ver     Solve version of ascii hfile (*100 to get integer)
*   gamit_hf_ver  Version of gamit hfile (introduced 2.0 101015) 

      integer*4 qgpst_utc, qcons_type, qsys_type, solve_ver,
     .          gamit_hf_ver 
      

*   rad_known   - Logical that indicates that we know the radiations
*                 in the current hfile.  Set to False for each new
*                 hfile.
*   svant_known - Logical to indicate if we know number of svant offset
*                 parameters.

*   apr_missed  - Logical to indicated that apr_missed while reading
*                 the sinex file because it appears before the solution
*                 estimate in the file.
*   comp_eigen  - Logical set true to check the eigenvalues
*   igs_ptname  - Set true when IGS style occupancy codes will be generated
*   auto        - automatic applying -d=TR deconstraints to GAMIT BASELINE
*                 hfiles (verson 6.11). Default is true.  Set false with -a 
*                 option.  Also set false with when orbits estimated.  Can
*                 be overriden with F in the -D= option string. 
*   orb_est     - Set true if orbit estimated.


      logical rad_known, svant_known, apr_missed, comp_eigen, 
     .        igs_ptname, auto, orb_est

* estread, aprread, epread -- Set true once the matrix est and apriori and 
*     estimates themselvves have been read.   We then can remove the constraints.
* neq_sys, neq_read, vec_read, neq_inv  -- Set true when normal equations 
*     system found, the components are read, and the system inverted

      logical estread, aprread, epread, neq_sys, neq_read, 
     .        vec_read, neq_inv

 
*   glb_dir     - Directory for the outpur global files
*   hfile       - Name of the Q file
*   glb_file        - Name of the output global file
*   ephem_file      - Name of output ephemeris file

*   hfile_type  - Character string denoting the type of
*                 file.  This string is the one which will
*                 added to the end of the site names when
*                 new names are generated.  There are two
*                 types: GPS -- GAMIT GPS Solution,
*                        TER -- SOLVEM Terrestrial data
*   datum_type  - Type of datum for this hfile.  Currently
*                 two are accepted:
*                 LLR - Lat, Long, Radius Geocentric
*                 XYZ - Cartesian XYZ

      character*3 hfile_type, datum_type
 
      character*128 glb_dir, hfile, glb_file, ephem_file
 
*   qsite_names(max_glb_sites)      - Names of the sites
*   qsource_names(max_glb_sources)  _ Radio source names
*   qsvs_names(max_glb_svs)     - NAmes of the satellites
 
      character*8 qsite_names(max_glb_sites), 
     .            qsource_names(max_glb_sources),
     .            qsvs_names(max_glb_svs)

*   qgframe    -  Frame used in gamit.
*   qgprec     -  Precesion constant system used in gamit
*   qgsrpmod   -  Solar radiation pressure model used in gamit
*   qgtime     -  Time system used in gamit

       character*8 qgframe, qgprec, qgsrpmod, qgtime
 
*   qfull_names(max_glb_sites)      - Full names of the sites.
 
      character*32 qfull_names(max_glb_sites)

*   hf_ext  - Ooptional extent for the hfile name which is passed
*             on the keys line of GPS hfiles.

      character*4 hf_ext

*   decon_str - Deconstraint string (R for rotation, T for translation)
*   gamit_str - Deconstraint string for GAMIT BASELINE h-files
      character*8 decon_str, gamit_str
      
*---------------------------------------------------------------------
* MOD TAH 981020: New variables to support multiple polar motion
*     parameters

*  qparn_mul_pmu(2,3,max_mul_pmu) -- Parameter numbers for muliple
*                                    pmu values
*  qpmu_mul_apr(2,3,max_mul_pmu)  -- Apriori values for multiple pmu
*                                    values
*  qpmu_mul_asig(2,3,max_mul_pmu) -- Apriori sigmas for multi-pmu.  At
*                                    some point need to upgrade to 
*                                    covariance matrix.
*  qnum_mul_pmu(2,3)              -- Number of mutiple pmu values
*  qent_par_codes  -- Number of multi-epoch parameter codes
*  qent_apr_codes  -- Number of multi-epoch apriori codes
       
      integer*4 qparn_mul_pmu(2,3,max_mul_pmu), qnum_mul_pmu(2,3),
     .          qent_par_codes, qent_apr_codes
      real*8 qpmu_mul_apr(2,3,max_mul_pmu),
     .       qpmu_mul_asig(2,3,max_mul_pmu)

*  name_format  -- String containing naming format to be used
*  expt_code    -- String to includes on binary hfile name rather
*                  the default extracted from the input file name

      character*8 name_format, expt_code

*---------------------------------------------------------------------
* Save of code, occ and pt designations
      character*4 code_all(max_glb_sites)  ! all 4-character codes
      character*1 pt_all(max_glb_sites)    ! all 1-characeter pt codes
      integer*4   occ_all(max_glb_sites)   ! all interger occ codes/

      integer*4 pad01
                              

*---------------------------------------------------------------------
*
*     COMMON DECLARATION
*
 
      common / qtoglb_com / llr, site_pos, site_vel, site_axo, 
     .    site_atm, source_pos, pmu_pos, svs_pos,
     .    qstart_epoch, qend_epoch, qmid_epoch, qscale, 
     .    rot_mat, qapr_list, qsol,
     .    qsig, qref_ep, qapr, qapr_sig, snx_constraint,
     .    qmul_par_ep, qmul_apr_ep, 
     .    qnut_ang_apr, qut1_apr, qwob_apr, qpmu_mul_apr, 
     .    qpmu_mul_asig,
     .    qtai_utc, 
     .    qnum_sites, qnum_sources, qnum_svs, qnum_parn, 
     .    qnum_misc, qparn_sites, 
     .    qparn_vel, qparn_axo, qparn_atm, qparn_sources, 
     .    qparn_svs, qparn_pmu, qparn_tran,  qparn_scale, qparn_misc,
     .    qglb_codes, qapr_codes, itoo, atoo, atos,
     .    num_itoo, num_atoo, qrun_time, 
     .    qnum_apr_codes, qnum_par_codes, qnum_soln_recs,
     .    qnum_comb, rad_num, svant_num, atm_num, gtol_sites, gtol_svs, 
     .    usr_vma_space, pad01, istart_vma,
     .    qparn_mul_pmu, qnum_mul_pmu, qent_par_codes, qent_apr_codes,
     .    rad_known, svant_known, apr_missed, comp_eigen, igs_ptname, 
     .    auto, orb_est, estread, aprread, epread, 
     .    neq_sys, neq_read, vec_read, neq_inv,
     .    qgpst_utc, qcons_type, qsys_type, solve_ver, gamit_hf_ver, 
     .    glb_dir, hfile, glb_file,
     .    ephem_file, qsite_names, qsource_names, qsvs_names, 
     .    qgframe, qgprec, qgsrpmod, qgtime, qfull_names,
     .    hfile_type, datum_type, hf_ext, name_format, expt_code,
     .    decon_str, gamit_str
      common / descrpt / occ_all, code_all, pt_all

      common / redparam / sngday, rm_param, rm_site
 
*----------------------------------------------------------------------
* Varibles for SINEX comapatibilty.
*
*     qdata_st  - Start JD for the data in this solution for this
*                  site.  This nee
*     qdata_en  - End JD for the data in this solution.

*     qrecv_st  - Start JD for the receiver characteristics (not
*                  the data spanned in this solution)
*     qrecv_en  - End JD for  receiver characteristics
*     qante_st  - Start JD for the antenna characteristics (not
*                  the data spanned in this solution) 
*     qante_en  - End JD for antenna characteristics

      real*8 qdata_st(max_glb_sites), qdata_en(max_glb_sites),
     .       qrecv_st(max_glb_sites), qrecv_en(max_glb_sites), 
     .       qante_st(max_glb_sites), qante_en(max_glb_sites)
 
*      qarp_ecc(3) - Single station eccentricity of the
*           - antenna with respect ground mark
*           - (NEU m)
*     qL1a_ecc(3)  - L1 phase center to ARP (NEU m).
*     qL2a_ecc(3)  - L2 phase center to ARP (NEU m).
*     qelev_cut    - Elevation cutoff angle (rads)
* MOD TAH 200205: Added qantdaz for deviation of antenns from
*     true north (3.10 version of h-file).
 
      real*4 qarp_ecc(3,max_glb_sites), qL1a_ecc(3,max_glb_sites), 
     .       qL2a_ecc(3,max_glb_sites), qelev_cut(max_glb_sites)
      real*4 qantdaz(max_glb_sites)  ! Antenna azimuth error (deg)
 
*         qnum_zen - Number of zenith delay parameters

      integer*4 qnum_zen(max_glb_sites)

*   qatmload(3,max_glb_sites) -- Average daily load applied at this site (NEU mm)
*   qhydload(3max_glb_sites) -- Average hydrographic load (NEU mm)

       real*4 qatmload(3,max_glb_sites), qhydload(3,max_glb_sites)
 
 
*   qant_mod   - Antenna model used for this station
* MOD TAH 101015: Increased length of C*16
 
      character*16 qant_mod(max_glb_sites)

*   qcode      - Code name for site from sinex file. 
*   qante_sn   - Antenna serial number (----- if not known)
*   qrecv_sn   - Serial number for the receiver (----- if
*           - not known)
 
      character*8 qcode, qante_sn(max_glb_sites), 
     .            qrecv_sn(max_glb_sites)
 
*   qrecv_ty  - Type of receiver
*   qante_ty   - Antenna type
*   qrecv_fw   - Firm ware version for the receiver.
*   qexpt_title - Title for the experiment
*   qradome_ty  - Radome type on antenna
      character*20 qrecv_ty(max_glb_sites)  
      character*16 qante_ty(max_glb_sites), 
     .             qrecv_fw(max_glb_sites)
      character*4 qradome_ty(max_glb_sites)
     
      character*32 qexpt_title
 
*-----------------------------------------------------------------------
 
      common / qsinf_com / qdata_st, qdata_en, qrecv_st, qrecv_en, 
     .       qante_st, qante_en, qarp_ecc, qL1a_ecc, qL2a_ecc, qantdaz,
     .       qelev_cut, qnum_zen , qatmload, qhydload, qant_mod,
     .       qcode, qante_sn, qrecv_sn,
     .       qrecv_ty, qante_ty, qradome_ty, qrecv_fw, qexpt_title

*-----------------------------------------------------------------------
* MOD TAH 050622: Added block for satellite information
*
      integer*4 qsvi_prn(max_glb_svs)    ! PRNs of saved entries
     .,         qsvi_svn(max_glb_svs)    ! SV number
     .,         qsvi_block(max_glb_svs)  ! Block number

      real*8 qsvi_antpos(3,2,max_glb_svs)  ! Antenna positions read from input files
* MOD TAH 200130: Added qsvi_antapr(3,max_glb_svs) which is apriori LC value used
*     in estimates of PCO.
      real*8 qsvi_antapr(3,max_glb_svs)    ! Apriori PCO for estimates 
      real*8 qsvi_launch(max_glb_svs)      ! Launch date for satellite
      character*16 qsvi_antmod(max_glb_svs)  ! Antenna model for each satellite (SINEX
                                         ! only need 10 charcaters in name)
      character*4 qsvi_ocode(max_glb_svs)  ! Codes for satellites. 
                                         ! 1/2 Frequencues; A/R and E/F

      common / qsvinf_com / qsvi_prn, qsvi_svn, qsvi_block
     .,      qsvi_antpos, qsvi_antapr, qsvi_launch 
     ,,      qsvi_antmod, qsvi_ocode


*-----------------------------------------------------------------------

* Information to be saved for the solution information records

* max_sln_save - Maximum number of solutions from which information
*                can be saved
* MOD TAH 190603: Increased to handle GAGE NAM14/IGS14 solutions.
      integer*4 max_sln_save

      parameter ( max_sln_save = 4000000 )

* qsepoch_start, qsepoch_end, qsrun_time - Times for data start, end and
*      the runtime.
* 

      real*8 qsepoch_start(max_sln_save), qsepoch_end(max_sln_save),
     .       qsrun_time(max_sln_save)

* qapr_cov(3,3,max_glb_sites)  - Apriori constraints applied to the station
*            positions
* qvel_cov(3,3,max_glb_sites)  - Apriori constraints applied to the station
*            velocities
* qapr_vel(3,3,max_glb_sites)  - Apriori constraints on velocities
* qapr_tran(3,2)               - Apriori constraints on transaltions and rate
* qapr_scale(2)                - Apriori constraints on scale and rate
* qapr_misc(max_glb_parn)      - Apriori constraints on miscellaneous parameter
* qapr_pmu(2,3)                - Aprioti constriants on Polar motion/UT1
* qapr_svs(6,6,max_glb_svs)    - Apriori covariance matrix for satellite
*                                position and velocity
* qapr_svs_ant(3,max_glb_svs)  - Apriori sigma on antenna offsets.
* qapr_diag(max_glb_parn)      - Diagonal element apriori sigmas

      real*8 qapr_cov(3,3,max_glb_sites), qvel_cov(3,3,max_glb_sites),
     .       qapr_svs(6,6,max_glb_svs),qapr_svs_ant(3,max_glb_svs),
     .       qapr_diag(max_glb_parn), qapr_vel(3,3,max_glb_sites),
     .       qapr_pmu(2,3), qapr_tran(3,2), qapr_scale(2), 
     .       qapr_misc(max_glb_parn)

* qnum_apr_cov, qnum_apr_svs, qnum_apr_diag - Number of entries in
*      each of the above categories.

      integer*4  qnum_apr_cov, qnum_apr_svs, qnum_apr_diag,
     .           qnum_vel_cov

* qsnum_parn - number of parameters
* qscons_type  - Types of constrainted applied
* qssys_type   - Bit mapped data system type
* qsglb_ver    - Version number of input files
* qsnum_save   - Number of solution inf values saved


      integer*4 qsnum_parn(max_sln_save), qscons_type(max_sln_save),
     .          qssys_type(max_sln_save), qsglb_ver(max_sln_save), 
     .          qsnum_save

* qsowner, qscreator - Owner and creator
* qsprog_gen   - Generator format
* qowner       - Owner from %=SNX line
* qcreator     - Creator from %=SNX line

      character*4 qsowner(max_sln_save), qscreator(max_sln_save),
     .            qsprog_gen(max_sln_save), qowner, qcreator

* qsanal_type - Type of solution

      character*6 qsanal_type(max_sln_save)

* qskalobs_file - Name of the input file
* qsexpt_title  - Title of the experiment

      character*128 qskalobs_file(max_sln_save)
      character*32  qsexpt_title (max_sln_save)

      common / qsln_save / qsepoch_start, qsepoch_end, qsrun_time,
     .         qapr_cov, qvel_cov, qapr_svs, qapr_svs_ant, 
     .         qapr_diag, qapr_vel, 
     .         qapr_pmu, qapr_tran, qapr_scale, qapr_misc, 
     .         qsnum_parn, qscons_type, qssys_type, qsglb_ver,
     .         qsnum_save, qnum_apr_cov, qnum_vel_cov,
     .         qnum_apr_svs, qnum_apr_diag,
     .         qsowner, qscreator, qowner, qcreator, qsprog_gen,
     .         qsanal_type, qskalobs_file, qsexpt_title 

*--------------------------------------------------------------------

*  Variables for glbtosnx

*   glbtosnx_opt - Options for glbtosnx.  Bit mapped current to options
*      Bit  Option
*        1  Output on header information.  No covariance matrix.
*   otoi(max_glb_parn) - Mapping of the output parameter number to
*      the input parameter number

      integer*4 glbtosnx_opt, otoi(max_glb_parn)

      real*8 Var_Limit ! Variance limit for outout of site to SINEX 
                       ! files (TAH 210112: Ver 2.09)


*   snx_comfile  - Name of the file woth comments for sinex files

      character*128 snx_comfile

      common  / glbtosnx_comm / Var_limit, glbtosnx_opt, otoi, 
     .          snx_comfile

*--------------------------------------------------------------------
