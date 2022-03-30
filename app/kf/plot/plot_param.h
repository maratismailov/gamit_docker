c
c     Parameter file for the plotting program  PLOT
c
c     T. Herring             10:18 AM  TUE., 19  NOV., 1985
c
c
      integer*4 max_data, max_plt_cmds, max_device, scr_size, max_plot,
     .    headr_size, max_site, max_souc, max_mar,
     .    char_recl, int_recl, max_headers, max_velvec, max_ellipse,
     .    max_clover, max_strains
 
c
      integer*4 max_plt_space
 
c
      character*(*) help_file, pltsl_help_file
 
c
*   dcb_size    - Size of the DCB to used in reading the bakfile
*               - (PLTSL only)
*   max_labels  - Maximum number of labels which can be saved
*   max_labels_len  - Maximum length of all of the labels to be
*               - saved and output in OUTPUT_LABELS
*   max_mar_types   - Maximum number of markov parameters allowed
*                   - (Used in PLTSL)
 
*   max_poly    - Maximum order of polynomial which can be fitted
*               - to the data.
*   max_poly_label  - Maximum number of labels allowed for polynomials
*               - (number of coefficients plus two for reference epoch
*               - and rms)
*   max_poly_matrix - Size of lower triangluar matrix for normal
*               - equations for polynomial.
*   max_types   - Maxium number of plot types.  Four greater than
*               - max_mar_types to allow for time, elevation, and
*               - delay and rate residuals (PLTSL only)
*   max_units   - Maximum number of units for the quantities in the
*               - the BAKFILE (PLTSL only)
*   max_files   - Maximum number of files passed in runstring
 
      integer*4 dcb_size, max_labels, max_labels_len, max_mar_types,
     .    max_poly, max_poly_label, max_poly_matrix, max_types,
     .    max_units, max_files
 
      parameter ( dcb_size       = 2576 )
      parameter ( max_labels     =  100 )
      parameter ( max_labels_len = 1024 )
      parameter ( max_mar_types  =   47 )
*                                           ! Upto quadradic allowed
      parameter ( max_poly       =   10 )
      parameter ( max_poly_label = max_poly + 3 )
      parameter ( max_poly_matrix= (max_poly+1)*(max_poly+2)/2 )
      parameter ( max_types      = max_mar_types + 4 )
      parameter ( max_units      =    8 )
      parameter ( max_files      =   10 )
 
c Parameters
c ----------
c max_data -- Maximum number of data which can be plotted from one file
c
      parameter ( max_data = 5000 )
c
c max_plt_cmds -- Maximum number of commands in plotting package
! MOD TAH 131015: Added KEYBOARD and RETURN commands (52-54 commands)
! MOD TAH 131015: Added NEW_WIND new window to set new size and position.
c
      parameter ( max_plt_cmds = 60 )
c
c max_device   -- maximum number of devices in system
c
      parameter ( max_device = 4 )
c
c scr_size -- the size of the scratch common area for the plot program
c
*                                          ! 200 words = 400 characters
      parameter ( scr_size   = 200000 )
c
c max_plot -- the maximum number of points which can be plotted
c     at one pass
c
*                                           ! two real*4 values
      parameter ( max_plot   = scr_size/4 )
c
*                                 ! maximum number of sites allowed
      parameter ( max_site =  8 )
c                                   in the plot files
c
*                                 ! maximum number of sources allowed
      parameter ( max_souc = 32 )
c                                   in the plot files
c
*                                 ! maximum number of markov parameters
      parameter ( max_mar  =100 )
c                                   which can be read from plot files
c
*                                                                  ! the
      parameter ( headr_size = 4*max_site + 4*max_souc + max_mar )
c                                   largest size the headr record can be
c
c
*                                    ! the number of bytes in each record
      parameter ( char_recl = 500  )
c                                   of the file.  The actual number may
c                                   be less than than this value
c
*                                          ! the maximum file record length
      parameter ( int_recl  = char_recl/4)
c                                   in worrds.
c
*                                     ! maximum number of header records
      parameter ( max_headers = 10)
c                                   which can be saved.
c
*                                     ! Maximum number of velocity vectors 
*                                     ! that can be saved
      parameter ( max_velvec = 1024 ) 
*                                     ! number of points in error ellipse
      parameter ( max_ellipse = 100 )
c     maximum number of points in strain clover leaf
      parameter ( max_clover = 1000) 
c     maximum number of strain records
      parameter ( max_strains = 1024)
*
*                                         ! the available ema space.
      parameter ( max_plt_space = 20000000)
c                                   Some of this space will be VMA
c
c
*                                            ! the name of the help
      parameter ( help_file = 'cplotx.hlp' )
c                                   file to be used in PLOT
c
*                                                   ! the name of the help
      parameter ( pltsl_help_file = 'pltsl.hlp' )
c                                   file to be used in PLTSL
 
