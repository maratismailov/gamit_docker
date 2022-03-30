*
* Include data for the MAKEXK commands
*
*......................................................................
* mod history                  
* G. Chen, MIT,  July 20, 1996
* MOD AZ 190305: change from USE_GPTGMF to ATM_MODELC
*......................................................................
*
*
*     Include file which contains the names of the commands used
*     by MAKEXK. Commands contains the basic command set, the @_conts
*     values contain the names of the contributions which can
*     be applied.  
 
c
*   max_xk_cmds    - Number of commands in the MAKEXK program
      integer*4  max_xk_cmds, file_start, file_end,
     .    site_start,site_end,
     .    mis_start,mis_end
  
      parameter ( max_xk_cmds = 60)
*                                                 ! bounds for the control
      parameter ( file_start = 1 , file_end = 2)
c                                   words which are applicable to obs. files.
      parameter ( site_start = 3 , site_end = 9)
c                                   words which are applicable to stations.
c
*                                                 ! bounds for the control
      parameter ( mis_start = 10 , mis_end = max_xk_cmds )
c                                    words which are miscellaneous.
c

*                                    

*   commands(max_xk_cmds)  -  The array containing the commands
*               - for the MAKEXK.  
c

      character*16 commands(max_xk_cmds)


c
*.... initialise the control name array
*                                       !  1  -------
      data commands      / 'OBS_FILE',
     .                     'RDX_FILE',
                                        !  2  --------
     .                     'ATM_BIAS',
     .                     'ATM_STATS',
     .                     'SITE_POS',    
     .                     'TIMEDEP_PROCNS',
     .                     'SITE_STATS',        
     .                     'RCV_TYPE',
     .                     'ANTE_OFF',   
                                        !  3  ------
     .                     'APR_FILE',       ! Coordinates of site
     .                     'NAV_FILE',
     .                     'SVJ_FILE',
     .                     'RES_ROOT',
     .                     'BF_SET',
     .                     'DEBUG',
     .                     'ION_STATS',
     .                     'START_TIME',          ! YY MM DD HR MIN seconds
     .                     'INTERVAL',            ! seconds
     .                     'NUM_EPOCHS',
     .                     'OUT_TYPE',
     .                     'DATA_NOISE',
     .                     'DATA_TYPES',          ! # L1 L2 .. 
     .                     'EDIT_SSV',
     .                     'AMBIN_FILE',
     .                     'CHANGE_WN',         
     .                     'AMB_CYCLE',         
     .                     'GEOD_OUT',
     .                     'CUT_OFF' ,
     .                     'SEARCH_TYPE',             
     .                     'ATM_FILE',
     .                     'POS_ROOT' ,
     .                     'MET_MODEL',
     .                     'FLOAT_TYPE',
     .                     'MODE',
     .                     'BACK_TYPE',
     .                     'OUT_SIG_LIMIT',
     .                     'RMS_EDIT_TOL',
     .                     'EXCLUDE_SVS',
     .                     'WLS_ROOT',
     .                     'RWL_ROOT',
     .                     'SUM_FILE',
     .                     'STOPGO_MODE',
     .                     'MWWL_JUMP',
     .                     'ANTMOD_FILE',
     .                     'USR_ADDBF',
     .                     'USR_DELBF',    
     .                     'REF_NEU',    
     .                     'ATM_MODELC',    
     .                     'TIME_UNIT',    
     .                     'SET_MSEC_BF',
     .                     'DCB_FILE',
     .                     'IONEX_FILE',
     .                     'IONLOS_FILE',
     .                     'MIN_TOLS',
     .                     'RM_CSLIP',
     .                     'TR_GNSS',
     .                     'USE_BLQ',
     .                     'DUMMY58',
     .                     'DUMMY59',
     .                     'DUMMY60'    /

