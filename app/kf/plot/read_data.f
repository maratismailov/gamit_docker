CTITLE    .................................................................
 
      subroutine read_data(bak_array, x_data, y_data, pt_data)
c
c     Routine to read data from the file ema into data ema area
c     The file is read according to the field declarers given
c     in the X_FIELD and Y_FIELD commands
c     A reference value is removed from both the x and the
c     y data to aviod rounding error problems with storing
c     the data in single precision.
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c bak_array -- the ema storage of the input file
c x_data, y_data the x and y data and there errors to be plotted
c pt_data -- the symbol to be used for plot
c
      real*4 x_data(2,1), y_data(2,1)
 
c
      integer*4 bak_array(bak_recl,1), pt_data(2,1)
 
c
c
c Local variables
c ---------------
c ierr -- general error variable
c inbuffer -- the buffer into which the data from the file is
c     read.
c values -- values to be read from file
c dx, dy -- range of the default scales (used to increase boundaries by
c     one percent)

c     no longer used
c     integer*4 ierr
 
c     not longer used
c     character*(char_recl) inbuffer
 
c
      real*8 values(2)
 
c
      real*4 dx, dy
 
c
c Variables used in reading file
c ------------------------------
c point -- the point to be used (returned from decode -- only the y_value
c     used), and it edit flag (0 is good data)
c
      integer*4 point(2), i
 
c
*   old_ref_valx    - Old reference value.  Used when differencing
*   old_ref_valy    - Old reference value.  Used when differencing
      real*8 old_ref_valx, old_ref_valy
 
      real*4 t1, t2
 
c
c.... Set file read false
      scale_set = .false.
      if( reset_scales ) then
*                                   ! use the default scales
         use_def_scale = .true.
      else
         use_def_scale = .false.
      end if
c
c.... Initalize the default boundaries
      do i = 1, 3, 2
*                                       ! minimum value
         default_scale(i)   =  1.d20
*                                       ! maximum value
         default_scale(i+1) = -1.d20
      end do
c
c.... see if the input file has been given and that the field information
c     is valid
      if( .not.file_read .or. x_field(1).lt.0 .or.
*                                 ! either no file or no fields
     .    y_field(1).lt.0 ) then
c
*                                    ! no file named or read
         if( .not.file_read  ) then
            write(termlu,'(2a)') 
     .       ' READ_DATA Error: No data file has been read,'
     .      ,' use FILE <file name> command")'
         end if
c
*                                            ! no valid x_field
         if( x_field(1).lt.0 ) then
            write(termlu,'(2a)') 'READ_DATA Error: No valid X_FIELD '
     .        ,'has been given, use X_FIELD <field data> command'
         end if
c
         if( y_field(1).lt.0 ) then
            write(termlu,'(2a)') 'READ_DATA Error: No valid Y_FIELD '
     .           ,'has been given, use Y_FIELD <field data> command'
         end if
c
*                              ! copy from file ema to data ema
      else
c
         num_data = 0
c
         old_ref_valx = ref_valx
         old_ref_valy = ref_valy
c
c....    Now read the data
         do i = 1, num_epochs
c
C           read(200,'(a)') inbuffer
c
c....       increment number of data
            num_data = num_data + 1
c
c....       get x_data firstly
            call decode_data(x_field, bak_array(1,i), values, point)
C           call decode_data(x_field, inbuffer      , values, point)
c
            if( num_data.eq.1 .and. reset_scales ) ref_valx = values(1)
 
c
            IF( point(1).ge.0 ) THEN
*                                 ! Difference
            if( xdiff ) then
               x_data(1,num_data) = x_data(1,num_data) - values(1) +
     .                                               old_ref_valx
               x_data(2,num_data) = sqrt(x_data(2,num_data)**2 +
     .                                   values(2)**2)
               values(1) = x_data(1,num_data)
               values(2) = x_data(2,num_data)
               ref_valx  = 0
            else
               x_data(1,num_data) = values(1) - ref_valx
               x_data(2,num_data) = values(2)
            end if
c
c....       Check the boundaries (if we have to)
            call check_bound(values, default_scale(1), ref_valx)
 
c
c....       get y_data next
            call decode_data(y_field, bak_array(1,i), values, point)
C           call decode_data(y_field, inbuffer      , values, point)
c
            if( num_data.eq.1 .and. reset_scales ) ref_valy = values(1)
 
c
*                                 ! Difference
            if( ydiff ) then
               t1 = y_data(1,num_data)
               t2 = y_data(2,num_data)
 
               y_data(1,num_data) = y_data(1,num_data) - values(1) +
     .                                               old_ref_valy
               y_data(2,num_data) = sqrt(y_data(2,num_data)**2 +
     .                                   values(2)**2)
               values(1) = y_data(1,num_data)
               values(2) = y_data(2,num_data)
               ref_valy  = 0
            else
               y_data(1,num_data) = values(1) - ref_valy
               y_data(2,num_data) = values(2)
            end if
 
            pt_data(1,num_data)  = point(1)
*                                            ! Edit flag
            pt_data(2,num_data)  = point(2)
*                         ! if point(1) => 0
            END IF
c
c....       Check boundaries
*                                            ! Check point really there
            if( point(1).ge.0 ) then
*                                            ! Get scales only for "good"
                if ( point(2).eq.0 ) then
*                                            ! data
                    call check_bound(values, default_scale(3), ref_valy)
                end if
            else
                num_data = num_data - 1
            end if
c
*                              ! looping over data
         end do
c
c....    Now increase the scales by 1%
         dx = default_scale(2) - default_scale(1)
         default_scale(1) = default_scale(1) - 0.01*dx
         default_scale(2) = default_scale(2) + 0.01*dx
c
         dy = default_scale(4) - default_scale(3)
         default_scale(3) = default_scale(3) - 0.01*dy
         default_scale(4) = default_scale(4) + 0.01*dy
c
c....    Assign the default scales to the scale
         if( reset_scales ) then
            do i = 1,4
               scale_size(i)  = default_scale(i)
               poly_window(i) = default_scale(i)
 
            end do
         end if
c
c....    Report the scales
         call report_scales(default_scale)
c
c....    Turn off differences
         xdiff = .false.
         ydiff = .false.
 
*                              ! name has been given
      end if
c
      return
      end
 
