      block data plot_bd

      implicit none 

c
c     Block data which will intialise the commands and other plot
c     default values
c
c Include files
c -------------
*                        ! the plot parameter file
      include 'plot_param.h'
c
*                        ! the plot common block
      include 'plot_com.h'
c
c Initialise the view_size
      data view_size / 0.125,  0.975,  0.075, 0.975 /
c
c Initalise signs of axes
      data sign_x / 1.0 /,  sign_y / 1.0 /
c
c Initalise the line types, piont type, and error bar type
      data line_type / 0,0 /, point_type / 4 /, errb_type / 1 /,
     .     pen_type  / 1   /, font_type  / 2 /, 
     .     font_name / '6x10' /

      data errb_scale / 1.0, 1.0 /
c
c Initalise character size (mm)
      data charsz_x / 2.0 /, charsz_y / 2.0 /
      data font_size / 6, 8 / 
c
c Initalise field defaults
      data x_field / 0, 1, 5 /, y_field / 1, 6, 7 /,
     .     p_field / 0, 0 /,
     .     v_field / 1, 2, 3, 4, 5, 6, 7, 0 /
c
c Initalise number of data
      data num_data / 0 /, start_data / 1 /, step_data / 1 /
c
c Initalize the header information
      data  input_num_headers  / 3 /,
     .      actual_num_headers / 3 /
 
c
c EMA area mapping defaults (needed in case data is processed before a
c     file is read)
*                              ! start of bak_array
      data ibak_array /  1 /
c
*                              ! start of x_data
      data ix_array   / 20 /
*                              ! start of y_data
      data iy_array   / 24 /
*                              ! start of pt_data
      data ipt_array  / 28 /
c
c LU information
      data num_device / max_device /
c
      data avail_lu / 4,  9, 24, 26 /
      data lu_ids   / 1,  1,  2,  3 /
c
c intialise the number of commands and number of plot types
      data num_commands / max_plt_cmds /
c
c Set the logical variables
c
      data file_read / .false. /
      data view_set  / .false. /
      data scale_set / .false. /
      data pop_scale / .false. /
      data use_def_scale / .true. /
      data reset_scales  / .true. /
      data map_mode  / .false. /
c
      data xdiff / .false. /
      data ydiff / .false. /
 
c     Polynomial fitting information
 
*                                     ! Use units of plot
      data conv_poly  / 1.d0, 1.d0 /
      data poly_units / 'Xunit', 'Yunit' /

      data vscale   / 1.0 /
      data vconfid  / 0.95 /
 
c The commands for PLTSL
*                                    ! 1  - End program command
      data commands / 'END     ',
*                                    ! 2  - draw plot
     .                'DRAW    ',
*                                    ! 3  - clear screen
     .                'ERASE   ',
*                                    ! 4  - get name of data file
     .                'FILE    ',
*                                    ! 5  - read the data file
     .                'READ    ',
*                                    ! 6  - defines X data field in file
     .                'X_FIELD ',
*                                    ! 7  - defines Y data field in file
     .                'Y_FIELD ',
*                                    ! 8  - defines the view area
     .                'VIEW    ',
*                                    ! 9  - defines user scales of view area
     .                'SCALE   ',
*                                    ! 10 - draw X axis at y min
     .                'XMN_AXIS',
*                                    ! 11 - draw X axis at y max
     .                'XMX_AXIS',
*                                    ! 12 - draw Y axisat x min
     .                'YMN_AXIS',
*                                    ! 13 - draw Y axis at x max
     .                'YMX_AXIS',
*                                    ! 14 - line type
     .                'LINE    ',
*                                    ! 15 - point type
     .                'POINT   ',
*                                    ! 16 - error bar type
     .                'ERRBARS ',
*                                    ! 17 - put label on plot
     .                'LABEL   ',
*                                    ! 18 - character size
     .                'CHARSZ  ',
*                                    ! 19 - sign of x axis
     .                'SIGN_X  ',
*                                    ! 20 - sign of y axis
     .                'SIGN_Y  ',
*                                    ! 21 - pen number
     .                'PEN     ',
*                                    ! 22 - font number to be used
     .                'FONT    ',
*                                    ! 23 - psuedo paper size
     .                'PS_SIZE ',
*                                    ! 24 - set X scale for plot
     .                'X_SCALE ',
*                                    ! 25 - set Y scale for plot
     .                'Y_SCALE ',
*                                    ! 26 - define the point field
     .                'P_FIELD ',
*                                    ! 27 - if scales reset for new data
     .                'RESET_SC',
*                                    ! 28 - draw all axes
     .                'AXES    ',
*                                    ! 29 - Help command
     .                'HELP    ',
*                                    ! 30 - Help command
     .                '?????   ',
*                                    ! 31 - Number of headers
     .                'HEADERS ',
*                                    ! 32 - status of plot parameters
     .                'STATUS  ',
*                                    ! 33 - pause for n milliseconds
     .                'PAUSE   ',
*                                    ! 34 - Define window for fitting poly.
     .                'WINDOW  ',
*                                    ! 35 - Tells to fit polynomial
     .                'FIT_POLY',
*                                    ! 36 - Unit labels for polynomial.
     .                'POLY_UNI',
*                                    ! 37 - Draw polynomial
     .                'PDRAW   ',
*                                    ! 38 - Mark the window for poly fitting
     .                'MARKWIND',
*                                    ! 39 - Identify point
     .                'IDENTIFY',
*                                    ! 40 - Difference x data
     .                'XDIFF   ',
*                                    ! 41 - Difference x data
     .                'YDIFF   ',
     .                'MPROJECT',
     .                'MSET    ',
     .                'MDRAW   ',
     .                'MRESET  ',
     .                'V_FIELD ',
     .                'VREAD   ',
     .                'VDRAW   ',
     .                'NETWORK ',
     .                'SFIELD  ',
     .                'SREAD   ',
     .                'SDRAW   ',
     .                'KEYBOARD',
     .                'RETURN  ', 
     .                'NEW_WIND',
     .                'POP_SCAL',
     .                'CPLOTX  ',    ! Used just to get help
     .                'DUMMY58 ',
     .                'DUMMY59 ',
     .                'DUMMY56 '  /

 
c
c Information for time axis plots
c -------------------------------
*                                     ! no turn over on years, start at 1950
      data tax_data /  9999,    50,
*                                     ! 12 months in a year,   start at 1
     .                   13,     1,
*                                     ! 30 days in a month,    start at 1
     .                   32,     1,
*                                     ! 24 hours in a day,     start at 0
     .                   24,     0,
*                                     ! 60 minutes in an hour, start at 0
     .                   60,     0,  
                                      ! 60 seconds in an min, start at 0 
     .                   60,     0 /
c
      end
 
