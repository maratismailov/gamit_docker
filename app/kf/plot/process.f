CTITLE PROCESS
      subroutine process( iel, bak_array, x_data, y_data, pt_data)

      implicit none 
c
c     Routines to process the commands given by the user.  iel
c     is the command number and determines the lines to code
c     to be executed (see block data)
c
c Include files
c -------------
*                        ! the parameter file
      include 'plot_param.h'
c
*                        ! the common block file
      include 'plot_com.h'
c
c Varaibles
c ---------
c iel    -- the command number based on its position in 'commands'
c bak_array -- the array containing the solution
c x_data -- the ema data to be plotted
c y_data -- the ema data to be plotted
c pt_data -- the ema lists of point types
c
      integer*4 iel
 
c
      integer*4 bak_array(bak_recl,1), pt_data(2,1)
 
c
      real*4 x_data(2,1), y_data(2,1)
 
c
*   cdum        - Dummy character for READ_LINE
 
      character*1 cdum, cd
 
 
*   ierr        - READ_LINE error
*   indx        - Pointer in string for readline
*   pause_secs  - Pause time in seconds
 
      integer*4 ierr, jerr, indx, pause_secs, i, id

* MOD TAH 131015: User position and window size
*     wind_pos is: x0, y0,  wdt0, hgt0
      integer*4 wind_pos(4)  ! Initial window position and size
      
 
c
c Functions
c ---------
c trimlen -- HP string length utility
c
      integer*4 trimlen, dum

* xp_err(max_ellipse), yp_err(max_ellipse)    ! points defining the error ellipse
* dxp_err(max_ellipse), dyp_err(max_ellipse)  ! points defining the error ellipse
*                                             ! Differential from end of ellipse.
*
      real*4 xp_err(max_ellipse),  yp_err(max_ellipse), 
     .       dxp_err(max_ellipse), dyp_err(max_ellipse)
c
c.... Goto the required subroutine (see PLOBD for commands)
      goto    ( 100,  200,  300,  400,  500,  600,  700,  800,  900,
     .   1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
     .   2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
     .   3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900,
     .   4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900,
     .   5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700 ) iel
c
c.... The END command: turn off
  100 continue
         pcontrol = 1
         return
c
c.... The DRAW command
  200 continue

*        see if step and start have been passed in the command (if
*        none are passed set values back to 1 and 1 )
         indx = 9
         step_data = 1
         start_data = 1
         call read_line(buffer, indx, 'I4', ierr, step_data, cdum)
         if( ierr.ne.0 .or. step_data.le.0 ) then 
             step_data = 1
         end if
         call read_line(buffer, indx, 'I4', ierr, start_data, cdum)
         if( ierr.ne.0 .or. start_data.le.0 ) then 
             start_data = 1
         end if
  
         call draw(x_data, y_data, pt_data)
         return
c
c.... The ERASE command
  300 continue
*                     ! Graphics 1000 II routine to clear screen
         call jclr
         return
c
c.... The FILE command
  400 continue
         call get_file
c
c....    Now set pcontrol to go to segment one to read the file
*                       ! force file read in segment 1
         iel = -4
         pcontrol = 1
c
         return
c
c.... The READ command
  500 continue
         pcontrol = 9
c        call read_data(bak_array, x_data, y_data, pt_data)
         return
c
c.... The X_FIELD command
  600 continue
         call get_field(x_field, xaxis_label)
         return
c
c.... The Y_FIELD command
  700 continue
         call get_field(y_field, yaxis_label)
         return
c
c.... The VIEW command
  800 continue
         call get_view(buffer)
         return
c
c.... The SCALE command
  900 continue
         call get_scale
         return
c
c.... The XAXIS_MN command
 1000 continue
c
c....    Add the defualt label to the buffer if not time axis
         if( x_field(1).ne.0 ) then
            buffer(trimlen(buffer)+2:) = xaxis_label
         end if
c
         pcontrol = 3
         return
c
c.... The XAXIS_MX command
 1100 continue
         pcontrol = 3
         return
c
c.... The YAXIS_MN command
 1200 continue
c
c....    If not a time axis add the default label
         if( y_field(1).ne.0 ) then
            buffer(trimlen(buffer)+2:) =  yaxis_label
         end if
c
         pcontrol = 3
         return
c
c.... The YAXIS_MX command
 1300 continue
         pcontrol = 3
         pel = iel
         return
c
c.... The LINE command
 1400 continue
C        call get_int(buffer, line_type, 1)
         indx = 9
         call read_line(buffer,indx, 'I4', ierr, line_type,cdum)
         if( ierr.ne.0 ) line_type(1) = 1
         call report_error('IOSTAT',ierr,'decod',buffer,0,'LINE_TYPE')
 
*        See if second argument
         call read_line(buffer,indx, 'I4',ierr, line_type(2),cdum)
*                                          ! Set for no point connection
         if( ierr.ne.0 ) line_type(2) = 0
 
*                                      ! Force line change, if we are not
         if( line_type(1).gt.0 ) then
*                                      ! turning off line drawing.
            call jlstl(line_type)
         end if
         return
c
c.... The POINT command
 1500 continue
         call get_int(buffer, point_type, 1)
         return
c
c.... The ERRBARS command
 1600 continue
         indx = 9
         call read_line(buffer, indx, 'I4', ierr, errb_type, cd)

*        Now get the scaling factors (if they are there)
         call read_line(buffer, indx, 'R4', ierr, errb_scale(1), cd)
         if( ierr.eq.0 ) then
             call read_line(buffer, indx, 'R4', jerr, errb_scale(2), cd)
          else
              errb_scale(1) = 1.0
          end if

*        See if error on second value (if there is then set equal to first)
         if( jerr.ne.0 ) then
             errb_scale(2) = errb_scale(1)
         end if

         return
c
c.... The LABEL command
 1700 continue
         pcontrol = 3
         return
c
c.... The CHARSZ command
 1800 continue
         call get_charsz
         call set_charsz( charsz_x, charsz_y )
         return
c
c.... The SIGN_X command
 1900 continue
         call get_int(buffer, sign_x, 1)
         return
c
c.... The SIGN_Y command
 2000 continue
         call get_int(buffer, sign_y, 1)
         return
c
c.... The PEN command
 2100 continue
         call get_int(buffer, pen_type, 1)
         call jcolr(pen_type)
         return
c
c.... The FONT command
 2200 continue
         indx = 9
         call getword(buffer, font_name, indx)
         read(font_name,'(i3)',iostat=ierr ) font_type
         if( ierr.ne.0 ) font_type = 1
*                    ! The font will be selected when labels are wriiten
         return
c
c.... The PS_SIZE command.  This is psuedo page size to make characters
c     come out correctly for a specific aspect ratio and size device
 2300 continue
          read(buffer(9:),*, iostat=ierr, err=2310 )
     .        x_size_mm, y_size_mm
 2310     continue
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'decod',buffer,0,
     .            'PS_SIZE')
          end if
c
          return
c
c.... The X_SCALE command
 2400 continue
         call get_xy_scale(x_field, scale_size(1), poly_window(1),
     .                     ref_valx )
         return
c
c.... The Y_SCALE command
 2500 continue
         call get_xy_scale(y_field, scale_size(3), poly_window(3),
     .                     ref_valy )
         return
c
c.... The P_FIELD command
 2600 continue
         call get_int(buffer, p_field, 2)
         return
c
c.... The RESET_SC command
 2700 continue
c
c....    See if string contains 'off'
         call get_cmd(buffer(9:), 'OFF     ', 1, iel)
c
*                             ! off found
         if( iel.eq.1 ) then
            reset_scales = .false.
*                             ! turn on reset scale
         else
            reset_scales = .true.
         end if
c
         return
c
c.... The AXES command
 2800 continue
*                       ! execute in segment 3
         pcontrol = 3
         return
c
c.... The HELP command
 2900 continue
c
c.... The ???? command  ! the same as help
 3000 continue
C        call help
*                       ! Schedule segment 9
         pcontrol = 9
c
         return
c
c..... The HEADERS command
 3100  continue
          indx = 9
          call read_line(buffer, indx, 'I4', ierr, input_num_headers,cd)
c
c         See if valid
          actual_num_headers = min(input_num_headers, max_headers)

*         Now see if ignore col one passed
          call read_line(buffer, indx, 'I4',ierr, id, cd)
          if( ierr.eq.0 ) then
*             Set so that column1 is not ignored and tested.
              ignore_col1 = .false.
          else
              ignore_col1 = .true.
          end if
c
          return
c
c.... The STATUS command
 3200 continue
C         call status
          pcontrol = 9
          return
c
c...  The PAUSE command
 3300 continue
*                                              ! Pause in millisec
          call get_int(buffer, pause_secs, 1)
C         call susp(2, pause_secs)
          call sleep( pause_secs )
          return
 
*.... The WINDOW command
 3400 continue
C         call get_poly_window
          pcontrol = 4
          return
 
*.... The FIT_POLY command
 3500 continue
C         call fit_poly
*                             ! A New fourth segment for poly fitting
          pcontrol = 4
          return
 
*.... The POLY_UNI command. Get labels and any conversion of units.
 3600 continue
          indx = 9
          call read_line(buffer, indx, 'CH', ierr,dum,poly_units(1))
          call read_line(buffer, indx, 'CH', ierr,dum,poly_units(2))
 
          call read_line(buffer, indx, 'R8', ierr, conv_poly(1), cdum)
          if( ierr.ne.0 ) conv_poly(1) = 1.d0
          call read_line(buffer, indx, 'R8', ierr, conv_poly(2), cdum)
          if( ierr.ne.0 ) conv_poly(2) = 1.d0
          return
 
*.... The PDRAW command
 3700 continue
          pcontrol = 4
          return
 
*.... The MARKWIND command (Display the window for the polynomial fitting
 3800 continue
          pcontrol = 4
          return
 
*.... The IDENTIFY command (Gives information about data point)
 3900 continue
          pcontrol = 4
          return
 
*.... The XDIFF command (differening invoked)
 4000 continue
          xdiff = .true.
          return
 
*.... The YDIFF command
 4100 continue
          ydiff = .true.
          return
 
*.... The MPROJECT command
 4200 continue
 
*         Command to read the map projection to be used. (Arguments start
*         at character 9
          indx = 9
          call read_line(buffer,indx, 'CH', ierr, id, map_proj)
          call casefold(map_proj)
          call multiread(buffer,indx, 'R4', ierr, ppos, cdum, 2)
          call read_line(buffer,indx, 'R4', ierr, rota, cdum )
 
*         values will be set when we draw map
          call set_map( 1 )
          return
 
*.... The MSET    command
 4300 continue
 
*        This command gets the limits on the the part of the map to
*        to be viewed.
 
         indx = 9
         call read_line(buffer, indx, 'CH', ierr, id, map_limit)
 
*        Reaction here depends on what is passed.
         call casefold( map_limit)
         if( map_limit.eq.'MA' ) return
         if( map_limit.eq.'LI' ) then
*            we need eight values.
             call multiread(buffer,indx,'R4', ierr, plim, cdum,8)
         else
*            Read four values
             do i = 1,4
                call read_line(buffer,indx,'R4',ierr, plim(1,i), cdum)
             end do
         end if
 
*        Agains thats all.  Will get the set the values later
         call set_map( 1 )
         return
 
*.... The MDRAW   command
 4400 continue
 
*        See if boundaries are passed.
         call set_map( 1 )
         indx = 9
         call decode_mdraw( buffer, indx)    
 
*        Now draw the map
         call draw_map
         return

*....  The MRESET command, gets us out of mapping mode
 4500  continue
         call set_map( 0 )
         return

*.... The V_FIELD command, read the velocity field information
 4600 continue 
         indx = 9
         call multiread(buffer, indx, 'I4', ierr, v_field, cdum,8)
*        set the x and y fields to look like normal data
         x_field(1) = 1
         y_field(1) = 1
         return

*     The VREAD command, reads the velocity data
 4700 continue
         call read_vfile ( bak_array )
         return

*     The VDRAW command
 4800 continue
         call vdraw( buffer, xp_err, yp_err, dxp_err, dyp_err,
     .               max_ellipse )
         return

*     The NETWORK command
 4900 continue
         call network
         return

*     The SFIELD (Strain Field) command
 5000 continue
c        think about this
         indx = 9
         call multiread(buffer, indx, 'I4', ierr, s_field, cdum,10)
*        set the x and y fields to look like normal data
         x_field(1) = 1
         y_field(1) = 1
         return
    
*     The SREAD command: read a strain file
 5100 continue
         call read_sfile ( bak_array )
         return

*     The SDRAW command: draw a strain rate clover leaf
 5200 continue
         call sdraw (buffer)
         return

*     KEYBOARD comamnd: Sends control to keyboard for interactive
*     commands
 5300 continue
*        If we are not in terminal mode, save the file lu number 
*        can set to 5
         if( .not.terminal ) then
             pushlu = userlu
             userlu = 5
             terminal = .true.
         else
             write(*,5310)
 5310        format('KEYBOARD command can only be used under file',
     .              ' control. Command ignored')
         end if
         return

*     RETURN command: Sends control back to file with commands
 5400 continue
*        Make sure we are in correct state to return to file control
         if( terminal .and. pushlu.ne.0 ) then
             terminal = .false.
             userlu = pushlu
         else
             write(*,5410)
 5410        format('RETURN can only be used when file control was',
     .              ' transferred with KEYBOARD command')
         endif
         return

*     NEW_WIND new window command;
 5500 continue
*        Read values from buffer
         indx = 9
         call read_line(buffer, indx, 'I4', jerr, wind_pos(1), cdum )
         if( wind_pos(1).ne.0 .and. jerr.eq. 0 ) then  ! Read remaining values
              call multiread(buffer, indx, 'I4', jerr, wind_pos(2), 
     .                       cdum,3)
         end if
         if( jerr.eq.0 ) then          
             call jbegn( meta_file, wind_pos(1), wind_pos(2),
     .                   wind_pos(3), wind_pos(4) )
         end if
         return

*     POP_SCALE
 5600 continue
         scale_set = .false.
         pop_scale = .true.
         call set_scale
         return
         
*     CPLOTX
 5700 continue 
*        Command to get help;
         write(*,5710) 
 5710    format('* CPLOTX command is used just to get help')
         return



c.... Thats all the command
      end
 
