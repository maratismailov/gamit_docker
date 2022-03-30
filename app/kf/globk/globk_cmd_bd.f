      block data globk_cmd_bd

      implicit none
 
 
*     Block data containing the commands for GLOBK.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cmds.h'
      include '../includes/globk_common.h'
 
*                                     !  1--Common block file name
      data glb_commands / 'COM_FILE'
*                                     !  2--Name of sorted solution list
     .,                   'SRT_FILE'
*                                     !     file
*                                     !  3--Name of ascii file which
*                                     !     containes the experiment list
     .,                   'LST_FILE'
*                                     ! 3b--Contains list of earthquakes
     .,                   'EQ_FILE'  
*                                     ! 3c--Tells glinit to make the svs_file
     .,                   'MAKE_SVS'  
*                                     !  4--Direction for sorting (+/-1)
     .,                   'SRT_DIR '
* MOD TAH 180402: Added new pre-command to specify that the old PRN names 
*      should be used for GPS (Default backwards comparabiloty)
     .,                   'USE_PRNN'
 
*     Any of the followwing commands will cause GLSORT to be run
*     if it has not yet been run
 
*                                     !  5--Apriori file with site and
     .,                   'APR_FILE'
*                                     !     source positions
*                                     !  6--Name of intermediate file for
     .,                   'SOL_FILE'
*                                     !     covariance matrices
*                                     !  7--Back solution file name
     .,                   'BAK_FILE'
*                                     !     (forces bak solution to run)
*                                     !  8--Sigma for apriori value of axis
     .,                   'APR_AXO '
*                                     !     offset (m)
*                                     !  9--Apriro sigmas for extended earth
     .,                   'BAK_PRTS'
*                                     !     tide coefficients
*                                     ! 10--apriori sigma for gamma
     .,                   'APR_GAM '
*                                     ! 11--apriori sigmas for nutation
     .,                   'APR_NANG'
*                                     !     angles (mas)
*                                     ! 12--Apriori sigmas for nutation
     .,                   'APR_NCOE'
*                                     !     series coefficients (mas)
*                                     ! 13--Apriori sigma for right
     .,                   'APR_RAO '
*                                     !     ascension orgin (mas)
*                                     ! 14--apririo sigmas for site
     .,                   'APR_SITE'
*                                     !     positions and velocities [XYZ (m)
*                                     !     and XYZ dot (m/y)]
*                                     ! 15--apriori sigmas for source
     .,                   'APR_SOUR'
*                                     !     position and velocity [RA,dec (mas)
*                                     !     and RA,Dec dot (mas/y)]
*                                     ! 16--apriori sigmas for h,l, lag
     .,                   'APR_TID '
*                                     ! 17--apriori sigmas for network
     .,                   'AMR_TRAN'
*                                     !     translation [XYZ (m)] and
*                                     !     markov process (m**2/yr)
*                                     ! 18--apriori sigmas for X an Y wobble
     .,                   'APR_WOB '
*                                     !     rates, ans seasonal terms.
*                                     ! 19--apriori sigmas fo UT1, rates and
     .,                   'APR_UT1 '
*                                     !     seasonal terms
*                                     ! 20--Markov variances for axis
     .,                   'MAR_AXO '
*                                     !     offset (m)
*                                     ! 21--markov variances for nutation
     .,                   'MAR_NANG'
*                                     !     angles (mas)
*                                     ! 22--markov variances for site
     .,                   'MAR_SITE'
*                                     !     positions and velocities [XYZ (m)
*                                     !     and XYZ dot (m/y)]
*                                     ! 23--markov variances for source
     .,                   'MAR_SOUR'
*                                     !     position and velocity [RA,dec (mas)
*                                     !     and RA,Dec dot (mas/y)]
*                                     ! 24--markov variances for X an Y wobble
     .,                   'MAR_WOB '
*                                     !     rates, ans seasonal terms.
*                                     ! 25--markov variances fo UT1, rates and
     .,                   'MAR_UT1 '
*                                     !     seasonal terms
*                                     ! 26--Value for axis offset (m)
     .,                   'VAL_AXO '
*                                     ! 27--Values for etd coefficients
     .,                   'VAL_ETD '
*                                     ! 28--values for nutation series
     .,                   'VAL_NCOE'
*                                     !     coefficients
*                                     ! 29--Values for X an Y wobble
     .,                   'VAL_WOB '
*                                     !     rates, ans seasonal terms.
*                                     ! 30--Values for UT1, rates and
     .,                   'VAL_UT1 '
*                                     !     seasonal terms
*                                     ! 31--Indicates that we should compute
     .,                   'COMP_RES'
*                                     !     should compute postfit residuals
*                                     !     (Default is not to)
*                                     ! 32--Name of output file containing
     .,                   'OUT_PMU '
*                                     !     Polar motion/UT1 values (only if
*                                     !     back solution)
*                                     ! 33--Name of output file for nutation
     .,                   'OUT_NUT '
*                                     !     angles (only if back solution)
*                                     ! 34--Name of output file for binary
     .,                   'OUT_GLB '
*                                     !     global solution.  May be later
*                                     !     used in GLOBK run
*                                     ! 35--Name of output ASCII file which
     .,                   'OUT_APR '
*                                     !     can be used as input for apriori
*                                     !     site and source positions.
*                                     ! 36--Options for output of back
     .,                   'BAK_OPTS'
*                                     !     solution (not implemented)
*                                     ! 37--Period of FCN in sidereal days
     .,                   'FCN_PER '
*                                     ! 38--Description line for solution
     .,                   'DESCRIPT' 
     .,                   'CRT_OPT '
     .,                   'PRT_OPT ' 
     .,                   'USE_SITE' 
     .,                   'USE_SOUR' 
     .,                   'APR_UANG'
     .,                   'APR_XYAN' 
     .,                   'APR_ETD '
     .,                   'APR_SVS '
     .,                   'MAR_UANG'
     .,                   'MAR_XYAN'
     .,                   'MAR_ETD '
     .,                   'MAR_SVS '
     .,                   'APR_UCOE'
     .,                   'APR_XYCO'
     .,                   'APR_ECOE'
     .,                   'SVS_FILE'
     .,                   'GLB_TIDE' 
*                                      ! 56--Nutation coefficient file
     .,                   'IN_NUT  '
     .,                   'IN_PLAN '
     .,                   'IN_PMU  '
     .,                   'IN_SD   '
     .,                   'APR_NEU '
     .,                   'MAR_NEU '
     .,                   'SVS_MARF'
*                                      ! New commands for entering radiation
*                                      ! parameters
     .,                   'APR_RAD '
     .,                   'MAR_RAD '
*                                      ! New parameter for max chi**2 inc.
*                                      ! and the max_prefit_diff
     .,                   'MAX_CHII'
     .,                   'USE_POS '
     .,                   'USE_NUM '
     .,                   'APR_TRAN'
     .,                   'MAR_TRAN'
     .,                   'APR_SCAL'
     .,                   'MAR_SCAL'
     .,                   'ORG_OPT '
     .,                   'ORG_CMD '
     .,                   'ORG_OUT '
     .,                   'NO_DIRCP'
     .,                   'SOURCE  '
     .,                   'RAD_RESE'
     .,                   'APR_SVAN'
     .,                   'MAR_SVAN'
     .,                   'VAL_SVAN'
     .,                   'MUL_PMU '
     .,                   'APP_PTID'
     .,                   'APR_ROT '
     .,                   'MAR_ROT '
     .,                   'EP_TOL  '
     .,                   'IRW_MOD '
     .,                   'SIG_NEU'
     .,                   'UNI_WGHT'
     .,                   'DECIMATE'
     .,                   'GLB_OPT '
     .,                   'APR_ATM '
     .,                   'MAR_ATM '
     .,                   'DEL_SCRA'
     .,                   'FREE_LOG'
     .,                   'APP_MODL'
     .,                   'DUMMY099'
     .,                   'DUMMY100' /

      data etd_names /    'RADIAL_DIUR_A+ '
     .,                   'RADIAL_DIUR_A- '
     .,                   'EAST_DIUR_A+   '
     .,                   'EAST_DIUR_A-   '
     .,                   'NORTH_DIUR_A+  '
     .,                   'NORTH_DIUR_A-  '
     .,                   'RADIAL_SEMI_A+ '
     .,                   'RADIAL_SEMI_A- '
     .,                   'EAST_SEMI_A+   '
     .,                   'EAST_SEMI_A-   '
     .,                   'NORTH_SEMI_A+  '
     .,                   'NORTH_SEMI_A-  ' /
 
      data nut_names /    'SLOPE_PSI '
     .,                   'SLOPE_EPS '
     .,                   'FCN_PSI   '
     .,                   'FCN_EPS   '
     .,                   '18Y_PSI   '
     .,                   '18Y_EPS   '
     .,                   '182.6D_PSI'
     .,                   '182.6D_EPS'
     .,                   '13.7D_PSI '
     .,                   '13.7D_EPS '
     .,                   '9YR_PSI   '
     .,                   '9YR_EPS   '
     .,                   '365.3D_PSI'
     .,                   '365.3D_EPS'
     .,                   '27.6D_PSI '
     .,                   '27.6D_EPS '
     .,                   '121.7D_PSI'
     .,                   '121.7D_EPS'
     .,                   '13.6D_PSI '
     .,                   '13.6D_EPS '
     .,                   '9.1D_PSI  '
     .,                   '9.1D_EPS  '
     .,                   '365.2D_PSI'
     .,                   '365.2D_EPS'
     .,                   '31.8D_PSI '
     .,                   '31.8D_EPS '
     .,                   '177.8D_PSI'
     .,                   '177.8D_EPS'
     .,                   '27.1D_PSI '
     .,                   '27.1D_EPS ' /

      data coeff_periods/ 'CONSTANT  '
     .,                   'FCN       '
     .,                   '18.6_YR   '
     .,                   '182.6_DAY '
     .,                   '13.66_DAY '
     .,                   '9.3_YR    '
     .,                   '365.26_DAY'
     .,                   '27.6_DAY  '
     .,                   '121.7_DAY '
     .,                   '13.62_DAY '
     .,                   '9.1_DAY   '
     .,                   '365.2_DAY '
     .,                   '31.8_DAY  '
     .,                   '177.8_DAY '
     .,                   '27.1_DAY  ' /

      data ut1_names   /  'UT1_DIUR_A+  '
     .,                   'UT1_DIUR_A-  '
     .,                   'UT1_SEMI_A+  '
     .,                   'UT1_SEMI_A-  ' /

      data xy_names    /  'XY_PRO_DIUR_A+ '
     .,                   'XY_PRO_DIUR_A- ' 
     .,                   'XY_RET_SEMI_A+ '
     .,                   'XY_RET_SEMI_A- '
     .,                   'XY_PRO_SEMI_A+ '
     .,                   'XY_PRO_SEMI_A- '  /

      data eq_commands /  'EQ_DEF  ',
     .                    'EQ_RENAM',
     .                    'EQ_COSEI',
     .                    'EQ_POST ',
     .                    'EQ_PRE  ',
     .                    'RENAME  ',
     .                    'EQ_LOG  ',  
     .                    'BREAK   ',  
     .                    'RESET   ',  
     .                    'DUMMY_10'  /

      data nonsec_types / 'OFFSET  ',
     .                    'PERIODIC',
     .                    'EXP     ',
     .                    'LOG     '  /

      data load_types   / 'USED    ',
     .                    'ATMLOAD ',
     .                    'HYDROLD '  /


      data glb_svs_file  / ' ' /
     .,    nut_inp_file  / ' ' /
     .,    plan_inp_file / ' ' /
     .,    pmu_inp_file  / ' ' /
     .,    svs_mar_file  / ' ' /

 
      end
