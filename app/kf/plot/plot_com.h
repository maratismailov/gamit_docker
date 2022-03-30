c
c     The main labelled common block for the PLOT program
c
c     T. Herring             10:22 AM  TUE., 19  NOV., 1985
c
* MOD TAH 910312: added Errbar_scale to rescale error bars
 
      common / plot_com / userlu, termlu, plotlu, 
     .    terminal, pcontrol,
     .    pel, buffer, runstring, input_file, saved_files, 
     .    meta_file,  file_read,
     .    num_data, start_data, step_data,
     .    input_num_headers, actual_num_headers, headers,
     .    typical_record,
     .    view_set, scale_set, use_def_scale, reset_scales,
     .    x_size_mm, y_size_mm, aspect_ratio,
     .    view_size, scale_size, default_scale, xdiff, ydiff,
     .    sign_x, sign_y, ref_valx, ref_valy,
     .    x_field, y_field, p_field,
     .    line_type, point_type, errb_type, errb_scale,
     .    pen_type, font_type,
     .    charsz_x, charsz_y, avail_lu, lu_ids, num_device, device_id,
     .    commands, num_commands, xaxis_label, yaxis_label, tax_data,
     .    num_epochs, ibak_array, ix_array, iy_array, ipt_array,
     .    bak_recl,
     .    num_labels, label_orient, label_just, label_char, label_pos,
     .    label_all, cori, just, pos, lfc, lf, wh,
     .    poly_order, num_poly_data, poly_window, conv_poly,
     .    poly_vec, poly_mat, sum_poly, poly_labels, poly_units,
     .    plt_pad,
     .    plim, ppos, rota, map_ou, map_proj, map_limit, map_cd,
     .    map_grid_space, map_sa, map_la, map_mode,
     .    vpos, vvel, vsig, vrho, vpnt, vname, v_field, vscale,
     .    vconfid, 
c         this is new for strain clovers: kurt 
     .    spos, strain, strainsig, s_field, num_strains, 
c         end new stuff
     .    num_velvec, ignore_col1, pushlu, pop_scale, num_files, 
     .    font_size,  font_name
 
*   label_char(max_labels)  - Charater number of the last character
*               - in each of the labels which have been saved.
*   num_labels  - number of labels to output
*   wh(2)       - The width and height of the labels which have
*               - been put on axes. Saved using SWH.
 
      integer*4 label_char(max_labels), num_labels, wh(2)
 
*   cori(2)     - Orientation of labels currently in effect.  Updated
*               - by call to SCJ.
*   just(2)     - Justification of labels currently in effect.  Updated
*               - by call to SCJ
*   label_just(2,max_labels)    - Justification for each of the labels
*               - to be output.
*   label_orient(2,max_labels)  - Orientation of each of the labels to
*               - be output
*   label_pos(2,max_labels)     - Positions of each of the labels to be
*               - output
*   lf(2)       - Current delta x and delta y which produce a line feed
*               - (used for labels and titles on axes)
*   lfc(2)      - Line feed in characters.  ROTOR will convert this
*               - value to lf in world coordiantes.
*   pos(2)      - Postion of labels currently in effect.  Updated by
*               - call to S2MOV or SR2MV.
 
      real*4 cori(2), just(2), label_just(2,max_labels),
     .    label_orient(2,max_labels), label_pos(2,max_labels), lf(2),
     .    lfc(2), pos(2)
 
*   label_all   - The concatination of all labels to be output.
 
      character*(max_labels_len) label_all
 
c
c Descriptions
c ------------
c
c IO devices
c ----------
c userlu -- the user's logical unit number (may be a file)
c pushlu -- LU number saved when keyboard command is used TAH 131015
c termlu -- the log device logical unit number (may be a file)
c plotlu -- the plot device logical unit number (can not be a file)
c terminal -- tell program if user is using a terminal (rather than
c     file input)
c pcontrol -- controls the execution of segments.  If pcontrol
c        equals 1 -- run segment one.
c        equals 2 -- run segment two. etc.
c pel  -- determines the operation of the program segement when
c     it is run.
c buffer -- buffer used to read user terminal or file
c runstring  -- 4 runstring parameters read with RCPAR
 
      integer*4 userlu, pushlu, termlu, plotlu, pcontrol, pel
 
c
      logical terminal
 
c
      character*256 buffer
 
c
 
      character*256 runstring(5)
 
c Data file information
c ---------------------
c input_file -- the name of the input file containing data
* meta_file  -- name to be given to meta file.
c file_read -- logical to indicate that data hase been read
c num_data  -- number of data read so far
* start_data -- Data point to start plot with (default 1)
* step_data  -- Number of data to step between points on plot (default 1)
c input_num_headers -- the number of header records in the file to be read.
c     (This number is given by the user and defaults to 3)
c actual_num_headers -- the number of header records which can be saved
c     for later use.
c headers -- the chartacter array of headers
c typical record  -- a copy of the first data record incase user doesnot
c     know what is in the file.
c
      character*256 input_file, meta_file
 
c saved_files(max_files) -- Names of saved files passed in runstring 
c     Accesed with #N
      character*256 saved_files(max_files)
      
c
      logical file_read
 
c
      integer*4 num_data, start_data, step_data, 
     .          input_num_headers, actual_num_headers

      integer*4 num_files  ! Number of files passed in runstring
c
      character*(char_recl) headers(max_headers), typical_record
 
c
c Logicals for contolling plots
c -----------------------------
c view_set -- indicate that JVIEW has been called
c scale_set -- indicates that JWIND has been called
c use_def_scale -- tell set_scale that it should use default scales
c reset_scales -- determines whether the scales will be reset when new
c     data are read.
c
      logical view_set, scale_set, use_def_scale, reset_scales
      logical pop_scale   ! Set true to pop to previous scale.
 
c
c xdiff, ydiff -- Tells program to difference values from current read
c     from those already in the data arrays.  NO CHECKS ARE MADE TO SEE
c     DATA CAN BE LEGITIMALLY DIFFERENCED.  The errors bars are summed
c     in an rms sense
 
 
      logical xdiff, ydiff
 
c
c Plot size information
c ---------------------
c x_size_mm -- the physical size of the plotting area in mm (obtained
c     from call to JDPMM) in the x direction (long side)
c y_size_mm -- same as x_size_mm except for y direction
c aspect_ratio -- the aspect ratio of plotting area (set for each device
c     in subroutine setup)
c
c view_size -- four element array containing the minimum and maximum values
c     of the area of the screen to be plotted on. Default values are
c     0.05 0.95 0.05 0.95 (order is for x and y axes). Values may be
c     changed by command VIEW
c scale_size -- the user coordinates of the corners of the view area.
c     Order is x min, x max, ymin, y max.  If the maximum values of
c     either x or y is less than the minimium value then the direction
c     of the corresponding axis will be reversed. Scale size can be
c     changed by SCALE command
c default_scale -- the default values for the user coordinates of the
c     view area.  These values are computed when the data is read in
c sign_x -- gives the sign of the x axis (+1 default, -1 indicates
c     axis direction is reversed)
c sign_y -- same as sign_x but for y_axis
c
      real*4 x_size_mm, y_size_mm, aspect_ratio, view_size(4),
     .    scale_size(4), default_scale(4)
 
c
      integer*4 sign_x, sign_y
 
c
c Reference values
c ----------------
c ref_valx -- the 0th order value of the x plot quanitity. This value
c     is used to stop lost of precision for large mean x values.
c ref_valy -- same as ref_valx for the y values
c
      real*8 ref_valx, ref_valy
 
c
c Field information
c -----------------
c x_field  -- the information about how to read the x data from file
c     x_field(1) -- data type  0 = time, 1 = normal data
c     x_field(2) -- type of data to be read combined with station:
c        (see plot type commands).
c     x_field(3) -- number of columns to be read.  For type 0 data (time)
c        x_field(3) can be any value between 1 and 6 with the order
c        being yr, month, day, hour, min, second. For type 1 data
c        x_field (3) can be one or two with the second argument
c        being taken as the standard deviation.
c
c y_field  -- same as x_field except for y data.
c
c p_field  -- defines the field in the data file to be used to generate
c     the plot symbols.
c     p_field(1) -- column number for plot character (1=A,2=B .... )
c     p_field(2) -- column number for edit flag.  A non zero value
c                   in this field will cause a symbol to be lower case
c
      integer*4 x_field(3), y_field(3), p_field(2)
 
c
c Plotting information
c --------------------
c line_type -- the type of line connecting point (see Driver manual for
c     line types for each device. If line type = 0, no line will connect
c     pionts. The second entry in line type give the point number to
c     to be connected if the point type is set (if its value is non-zero)
c point_type -- the type of point to be plotted. See Driver manual.
c     If point_type is zero, no point will be plotted
c errb_type -- type of error bar.  If errb_type = 0 no error bars will be
c     drawn.  No error bars will be drawn if x/y_field(3) = 1.
c errb_scale -- Rescales the x and y error bars (default 1)
c pen_type  -- type of pen to be used in pen.
c font_type -- the font to be for labels.  Currently only three fonts
c     are supported:
c     font_type 1 -- Simplex Roman font (fast)
c     font_type 2 -- Full Times Roman font
c     font_type 3 -- Graphics 1000 Symbol set
c     The default font is font_type 3
c charsz_x - the size of the characters in the x direction in mm.
c charsz_y - the size of the characters in the y direction in mm.
c
      integer*4 line_type(2), point_type, errb_type, pen_type,
     .    font_type

c font_name  - Name of font (for X version)

      character*256 font_name 
 
c
      real*4 charsz_x, charsz_y, errb_scale(2)
      integer*4 font_size(2)   ! Font size in pixels (returned by jfont)
 
c
c Device tables
c -------------
c avail_lu -- the lu numbers of the plot devices in the CfA system
c lu_id    -- the type of graphics devices the associated with each
c     lu.  There are three types currently:
c     type 1 -- HP2648A graphics terminal
c     type 2 -- HP2635/HP150 grahpics terminal
c     type 3 -- HP7470A flat bed plotter
c     The association of lu's and types is set in PLOBD
c num_device -- number of devices in the system.
c device_id -- the type of device being used (1,2 or 3)
c
      integer*4 avail_lu(max_device), lu_ids(max_device), num_device,
     .    device_id
 
c
c Command information
c -------------------
c commands -- the list of commands available to the system  (see PLKDB)
c num_commands -- the number of commands available
c
c xaxis_label -- the default axis label to be used for the xaxis
c yaxis_label -- the default axis label to be used for thr yaxis
c     These labels are constructed in GET_FIELD
c
      character*8 commands(max_plt_cmds)
 
c
      integer*4 num_commands
 
c
      character*36 xaxis_label, yaxis_label
 
c
c Information for time axis
c -------------------------
c tax_data -- gives the turnover values and start values for various
c     spacings in the time axis labels (see &plkbd)
* Increased to 6 to allow seconds
c
      integer*4 tax_data(2,6)
 
c
c num_epochs -- the number of epochs of data in the plot file
c
      integer*4 num_epochs
 
c
* Polynomial fitting variables
 
*   poly_order  - Order of the polnomial to be fitted (See MAX_POLY
*               - for maximum value)
*   num_poly_data   - Number of data in polynomial fit
 
      integer*4 poly_order, num_poly_data
 
*   poly_window(4)  - Min x, max x, min y and max y ranges to be used
*                   - when polynomials fitted
 
      real*4 poly_window(4)
 
*   conv_poly(2)    - Conversion factors for X and Y data for output
*                   - (ONLY) of the polynomial coefficients.
*   poly_vec(0:max_poly)    - Solution vector for polynomial
*   poly_mat(max_poly_matrix)   - Normal equationd for fitting
*                           - polymonial
*   sum_poly(2) - Summation statistics for polynomial.
*               - Sum_poly(1) = sum {res**2*wgh}
*               - Sum_poly(2) = sum {whg}
 
      real*8 conv_poly(2), poly_vec(0:max_poly),
     .    poly_mat(max_poly_matrix), sum_poly(2)
 
*   poly_labels(max_poly_label)  - Labels for polynomial.  May be
*               - Reference with \p<n>, where n is
*               - 1 -- Reference value and number of data
*               - 2 -- WRMS and Chi**2/f
*               - 3 -- Polynomial 0 th order value and sigma
*               - 4 -- Polynomial 1 st order value and sigma
*               - 5> and do on.
*               - (NOTE: Sigmas are scaled by Chi**2/f)
      character*80 poly_labels(max_poly_label)
 
*   poly_units(2)   - Unit labels for X and Y data types
 
      character*5 poly_units(2)
 
c EMA mapping information
c -----------------------
c ibak_array -- the start address of the bak_array in ema_data
c ix_array   -- the start address of x_data
c iy_array   -- the start address of y_data
c ipt_array  -- the start address of pt_data
c
      integer*4 ibak_array, ix_array, iy_array, ipt_array
 
c
c Bak_file record mapping
c -----------------------
c bak_recl  -- the record length of bak file
 
      integer*4 bak_recl
 
*
* Mapping parameters
* ------------------
*
*  plim(2,4)  -- Limits on map area. Usually only first index needed.
*  ppos(2)    -- lat and long of projection point (degs)
*  rota       -- Rotation angle of map (deg)

*  map_ou     -- Outline set to be used (default 'CO')
*  map_proj   -- Projection to be used.
*  map_limit  -- Type of limit.
*  map_cd     -- Dummy to get words aligned

*  plt_pad    -- Padding need by the sun4/60

      integer*2 plt_pad

      real*4  plim(2,4), ppos(2), rota

      character*2 map_ou, map_proj, map_limit, map_cd

*     Additional control on map
 
*   map_grid_space  - Spacing of the grif lines on the map
*   map_sa          - Altitude of satellite (earth radii)
 
      real*4 map_grid_space, map_sa
 
*   map_mode        - Indicates that we are in mapping mode
*   map_la          - Indicates that equation, and poles should be
*                   - labeled.
 
      logical map_mode, map_la
 
*     Error ellipise and velocity vector variables
 
*   vpos(2, max_velvec) - Positions for starts of the velocity vectors
*   vvel(2, max_velvec) - Values for the velocities
*   vsig(2, max_velvec) - X and Y sigmas on the velocities
*   vrho(max_velvec)    - correlation bewteen sigmas.
*   vscale              - Conversion from units of velocity input to
*                       - mm on the page (laser printer output, 200 mm
*                       - wide)
*   vconfid             - Confidence of error ellipse (0-1)
 
      real*4 vpos(2, max_velvec), vvel(2, max_velvec),
     .    vsig(2, max_velvec), vrho(max_velvec), vscale, vconfid

*   v_field(8)          - Field information which controls reading of the
*                       - input file.  The entries are:
*                       - 1 - column for x position (longitude)
*                       - 2 - column for y position (latitude)
*                       - 3 - column for velocity in x (East)
*                       - 4 - column for velocity in y (north)
*                       - 5 - column for sigma velocity in x (East)
*                       - 6 - column for sigma velocity in y (north)
*                     - 7 - column for correlation between x and y
*                     - 8 - column for name of site.
*   num_velvec          - Number of velocity vectors to be processed
*   vpnt(max_velvec)    - Point types to be output at the base of each
*                       - of the vectors (use the P_field command to enter)
 
      integer*4 v_field(8), num_velvec, vpnt(max_velvec)
 
*     vname(max_velvec)   - Names to be put at the base of vector.
      character*16 vname(max_velvec)
 
*     STRAIN RATE INFORMATION
 
*     Positions for center of cloverleaf (long, lat)
      real*4 spos(2, max_strains)
*     components of strain rates (eps1,eps2,spin,theta) 
      real*4 strain(4, max_strains)
*     uncertainties of strain rates  (eps1,eps2,spin,theta) 
      real*4 strainsig(4, max_strains)
*     Field information which controls reading of the Strain Rate File
*                       - 1 - column for x position (longitude)
*                       - 2 - column for y position (latitude)
*                       - 3 - column for 1st eigenvalue of strain rate
*                       - 4 - column for sigma of 1st eigenvalue of strain rate 
*                       - 5 - column for 2nd eigenvalue of strain rate
*                       - 6 - column for sigma of 2nd eigenvalue of strain rate 
*                       - 7 - column for spin rate (1/yr) 
*                       - 8 - column for sigma of spin rate (1/yr) 
*                       - 9 - column for azimuth of second eigenvector (deg)
*                       -10 - column for uncertainty of azimuth of second eigenvector (deg)
      integer*4 s_field(10)    
*     Number of strain cloverleaves to be processed
      integer num_strains
*    

*   ignore_col1         - Indicates that column1 is not special. (Default
*                         is true, and set false by passing non-zero value
*                         as 5 runstring paremeter.

      logical ignore_col1
 
