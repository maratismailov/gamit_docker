 
*     This include file contains the commands for GLOBK and the
*     names of some quantities which must be specified by name
*     (Nutation and ETD coefficient names)
 
 
*   glb_commands(max_glb_commands)  - Commands for GLOBK
*   eq_commands(max_eq_cmd)         - Commands for reading the
*                                     earthquake file.
 
      character*8 glb_commands(max_glb_commands),
     .    eq_commands(max_eq_cmd)
 
*   nut_names(max_nut_coeff)        - Names of the nutation series
*                                   - parameters
*   coeff_periods(max_coeff_periods) - Periods which are allowable
*                                      coefficient estimation eg.,
*                                      182.6_D
 
      character*10 nut_names(max_nut_coeff),
     .    coeff_periods(max_coeff_periods)
 
*   etd_names(max_etd_names)        - Names of the earth tide parameters
*                                   - eg, RADIAL_DIUR_COS
*   ut1_names(max_ut1_names)        - Names for the UT1 coefficent
*                                     parameters eg UT1_DIUR_COS
*   xy_names(max_xy_names)          - Names for the x and y polar motion
*                                   - coefficients eg, XY_RET_DIUR_COS

      character*16 etd_names(max_etd_names), ut1_names(max_ut1_names),
     .             xy_names(max_xy_names)

* MOD TAG 120714: Added LOAD_TYPES. set with the APP_MODL command.  Values
*     set in the globk_cmd_bd.f block data



      character*8 load_types(3)  
 
 
      common / glb_com_list / glb_commands, coeff_periods, eq_commands,
     .      nut_names, etd_names, ut1_names, xy_names, load_types
 
