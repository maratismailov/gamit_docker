      block data pltsl_bd

      implicit none 
c
c     Block data which will intialise the commands and other plot
c     default values
c
c Include files
c -------------
*                             ! the plot parameter file
      include 'plot_param.h'
c
*                             ! the plot common block
      include 'pltsl_com.h'
c
c
c EMA area mapping defaults (needed in case data is processed before a
c     file is read)
*                               ! time field in bak_array
      data ibepoch    /  1  /
*                               ! markov field in bak_array
      data imar_val   /  5  /
*                               ! markov field in bak_array
      data ipost_res  /  9  /
*                               ! elevation angles in bak_array
      data ibelev     / 13  /
*                               ! record numbers in bak_array
      data ibak_recs  / 17  /
*                               ! source number in bak_array
      data ibak_sou   / 18  /
*                               ! the unweight flags in bak_array
      data ibak_unw   / 19  /
c
c
c intialise the number of commands and number of plot types
      data num_types    / max_types    /
c
c
      data plot_types / 'CLK_OFFS'
     .,       'CLK_RATE'
     .,       'ATM_OFFS'
     .,       'ATM_RATE'
     .,       'NS_AZATM'
     .,       'EW_AZATM'
     .,       'N_SITE'
     .,       'E_SITE'
     .,       'U_SITE'
     .,       'AXIS_OFF'
     .,       'RA_SOURC'
     .,       'DEC_SOUC'
     .,       'X_POLE'
     .,       'XP_RAT'
     .,       'Y_POLE'
     .,       'YP_RAT'
     .,       'UT1-AT'
     .,       'UT_RAT'
     .,       'DPHI'
     .,       'DEPS'
     .,       'L_TIDE'
     .,       'H_TIDE'
     .,       'TIDE_LAG'
     .,       'RAD_1_CS'
     .,       'RAD_1_SN'
     .,       'EAS_1_CS'
     .,       'EAS_1_SN'
     .,       'SOU_1_CS'
     .,       'SOU_1_SN'
     .,       'RAD_2_CS'
     .,       'RAD_2_SN'
     .,       'EAS_2_CS'
     .,       'EAS_2_SN'
     .,       'SOU_2_CS'
     .,       'SOU_2_SN'
     .,       'UT1_1_CS'
     .,       'UT1_1_SN'
     .,       'UT1_2_CS'
     .,       'UT1_2_SN'
     .,       'XY_-1_CS'
     .,       'XY_-1_SN'
     .,       'XY_+2_CS'
     .,       'XY_+2_SN'
     .,       'XY_-2_CS'
     .,       'XY_-2_SN'
     .,       'GAMMA   '
     .,       'MAP_FUNC'
*                            ! Time is type 48
     .,       'TIME    '
     .,       'ELEV_ANG'
     .,       'DEL_RES '
     .,       'RATE_RES' /



      data unit_label /
     .        '(ns)'
     .,       '(D-14)'
     .,       '(ps)'
     .,       '(ps/s)'
     .,       '(ps)'
     .,       '(ps)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(us/s)'
     .,       '(mas)'
     .,       '(us/s)'
     .,       '(mts)'
     .,       '(us/s)'
     .,       '(mas)'
     .,       '(mas)'
     .,       ' '
     .,       ' '
     .,       '(deg)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(m)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       '(mas)'
     .,       ' '   
     .,       '(C) '
     .,       ' '   
     .,       '(deg)' 
     .,       '(ps) '
     .,       '(fs/s)' /

c
c Units of the for plot labels
c     data unit_label / '(ns)  ',
c    .                  '(D-14)',
c    .                  '(ps)  ',
c    .                  '(m)   ',
c    .                  '(mas) ',
c    .                  '(mts) ',
c    .                  '(deg) ',
c    .                  '(fs/s)' /
c
c Units look-up table. (For each plot type these entries give the units
c     to be used on the plots of this entry)
c     data unit_table / 1, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
c    .   6, 5, 5, 0, 0, 7, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 7, 3, 8 /
c
c Conversion factors from internal units to output units
*                               !  1 - convert ps to ns
      data conv_mar  / 0.001,
*                               !  2 - convert fs/s to D-14
     .                   0.1,
*                               !  3 - leave as ps
     .                   1.0,
*                               !  4 - leave as ps
     .                   1.0,
*                               !  5 - leave as ps
     .                   1.0,
*                               !  6 - leave as ps
     .                   1.0,
*                               !  7 - leave as m
     .                   1.0,
*                               !  8 - leave as m
     .                   1.0,
*                               !  9 - leave as m
     .                   1.0,
*                               ! 10 - leave as m
     .                   1.0,
*                               ! 11 - leave as mas
     .                   1.0,
*                               ! 12 - leave as mas
     .                   1.0,
*                               ! 13 - leave as mas
     .                   1.0,
*                               ! 14 - leave as mas
     .                   1.0,
     .                   1.0,
     .                   1.0,
*                               ! 15 - convert mas to mts
     .             0.0666667,
     .             0.0666667,
*                               ! 16 - leave as mas
     .                   1.0,
*                               ! 17 - leave as mas
     .                   1.0,
*                               ! 18 - leave dimensionless
     .                   1.0,
*                               ! 19 - leave dimensionless
     .                   1.0,
*                               ! 20 - convert rads to deg
     .             57.295780,
*                               ! 21 - leave in m
     .                   1.0,
*                               ! 22 - leave in m
     .                   1.0,
*                               ! 23 - leave in m
     .                   1.0,
*                               ! 24 - leave in m
     .                   1.0,
*                               ! 25 - leave in m
     .                   1.0,
*                               ! 26 - leave in m
     .                   1.0,
*                               ! 27 - leave in m
     .                   1.0,
*                               ! 28 - leave in m
     .                   1.0,
*                               ! 29 - leave dimension less
     .                   1.0,
     .                15*1.0  /
 
c
      end
 
