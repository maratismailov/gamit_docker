 
      block data glorg_bd

      implicit none  
 
*     Block data for glorg commands
 
      include '../includes/kalman_param.h'
      include '../includes/glorg_common.h'
 
 
*                                         ! 1 -- Name of GLOBK command file
      data glorg_commands / 'COM_FILE'
*                                         ! 2 -- Name of the apririo file to
*                                         !      use
     .,                     'APR_FILE'
*                                         ! 3 -- Defination of sites to be
*                                         !      used to get get plate rotation
*                                         !      vectors
     .,                     'PLATE_DE'
*                                         ! 4 -- Name of a saved global file
*                                         !      to be used
     .,                     'GLB_FILE'
*                                         ! 5 -- Gives the names of the sites
*                                         !      to use for fixing origin
*                                         !      all must have been estimated
     .,                     'USE_SITE' 
*                                         ! 6 -- Sites to be used for cov.
*                                         !      matrix computations
     .,                     'COV_SITE'
*                                         ! 7 -- Causes parameters to be
*                                         !      to be forced to be equal
*                                         !      adjustments.
     .,                     'EQUATE  '  
*                                         ! 8 -- Local_eq rotate positions into
*                                         !      NEU before applying equates
     .,                     'LOCAL_EQ'     
*                                         ! 9 -- Force paprameter to value
     .,                     'FORCE   '     
     .,                     'POS_ORG ' 
     .,                     'RATE_ORG' 
     .,                     'CND_HGTV' 
     .,                     'EQ_DIST ' 
     .,                     'UNEQUATE' 
     .,                     'FIRST_EQ' 
     .,                     'CONSTRAI'    
     .,                     'COND_SIG' 
     .,                     'STAB_ITE' 
     .,                     'STAB_MIN' 
     .,                     'SOURCE  ' 
     .,                     'ASSIGN_P' 
     .,                     'OUT_SOL ' 
     .,                     'STAB_SIT' 
     .,                     'NOPLATET'
     .,                     'EQ_4CHAR'
     .,                     'DUMMY26 '
     .,                     'DUMMY27 '
     .,                     'DUMMY28 '
     .,                     'DUMMY29 '
     .,                     'DUMMY30 '
     .,                     'DUMMY31 '
     .,                     'DUMMY32 '  /

      data  org_types    /  'XROT    ' 
     .,                     'YROT    '
     .,                     'ZROT    '
     .,                     'XTRAN   ' 
     .,                     'YTRAN   '
     .,                     'ZTRAN   '
     .,                     'SCALE   '     /
      end
 
