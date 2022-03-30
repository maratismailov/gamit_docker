CTITLE    .................................................................
 
      subroutine draw(x_data, y_data, pt_data)
c
c     Routine to draw the data in ema using the line type given
c     in line_type, point type in point_type
c
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
c x_data, y_data the x and  data and there errors to be plotted
c pt_data -- the type of point to be plotted if point_type<0
c
      real*4 x_data(2,1), y_data(2,1)
 
c
      integer*4 pt_data(2,1)
 
c
c
c Local variables
c ---------------
c xp  -- the array of x values to be plotted
c yp  -- the array of y values to be plotted
c symbol -- the character to be plotted
c ofx, ofy -- the offset in the position of charaters to make them centered
c i,j,k    -- loop counters
* pt_min, pt_max -- Min and max point types when second argument of line
*          type is set negative
* slat, slong  -- Used in converting data when plotting in map mode
c
      real*4 xp(max_plot), yp(max_plot), slat, slong
c
      integer*4 symbol, i, j, k, pt_min, pt_max
c
c Functions
c ---------
c ifbrk  -- logical function to detect break (no longer used)
c incpt  -- True if a point should be included in the line drawing
c line_break -- Indicates a break in a line so stop connecting
 
      logical incpt, line_break
c     logical ifbrk --no longer used
 
c
c Scratch common
c --------------
      common xp, yp
 
c
c.... Make sure we are set up
      if( .not.map_mode ) then
          call set_view
          call set_scale
      end if
      call set_charsz( charsz_x, charsz_y )
c
c.... Check to see if the have any data
*                                 ! No data.
      if( num_data.eq.0 ) then
          call report('DRAW Error: No data read yet, use READ command')
          return
      end if
 
*     Turn on clipping
      if( .not.map_mode ) then
         call jmcur
         call gsclip(1)
      end if
c
c.... Loop over data plotting the symbols
      do i = start_data, num_data, step_data
         if( .not.map_mode ) then
            xp(1) = sign_x*x_data(1,i)
            yp(1) = sign_y*y_data(1,i)
         else
            slong = x_data(1,i) + ref_valx
            slat  = y_data(1,i) + ref_valy
            call maptrn(slat, slong, xp(1), yp(1))
         end if
c
c....    Check which point type
         if( point_type.gt.0 ) then
c
c....       Plot the point with a mark
            call j2mrk(xp,yp,-point_type)
 
         end if
c
c....    see if characters to be plotted
         if( point_type.eq.-1 ) then
c
*           Directly call GKS routine.
            symbol = (pt_data(1,i)+64)
            call j2mrk( xp(1), yp(1), symbol ) 
c           call points(xp(1), yp(1), 1, symbol, 0 )
*                                                ^
*                                                | No Line
         end if
c
c....    see if graphics symbol
         if( point_type.eq.-2 ) then
            symbol = pt_data(1,i)
            call j2mrk(xp,yp, -symbol)
         end if
c
c....    Check for break
C        if( ifbrk() ) then
C           call jmcur
C           call report('User break: while drawing data')
*           Turn off clipping before we return
C           call gsclip(0)
C           return
C        end if
c
      end do
c
c.... Now do error bars
      if( errb_type.gt.0 .and. .not.map_mode) then
         do i = start_data, num_data, step_data
            call draw_errb(x_data(1,i), y_data(1,i))
c
c....       Check for break
C           if( ifbrk() ) then
C              call jmcur
C              call report(' User break: while drawing error bars')
*              Turn off clipping before we return
C              call gsclip(0)
C              return
C           end if
 
         end do
      end if
c
c
c.... Check on drawing line.  line_type has two argumnent.  The first sets
*     line style and the second the type of point to be connected with a
*     line.  If the first argument is negative then, decreasing x coordinate
*     will cause a line break.  If the second argument is negative then 
*     all points of differenct types will be connected with a separate line 
*     for each.
      if( Line_type(1).ne.0 .and. line_type(2).ge.0 ) then
*                                ! set the line type
         call jlstl(abs(line_type(1)))
c
c....    Loop over data in blocks of the largest we can fit in core
*               ! counter over data
         i = start_data-step_data
*               ! counter over data within size of max_plot
         j = 0
         do while ( i.lt. num_data-step_data+1 )
            i = i + step_data
c
c....       Put the data in main memory arrays
*           See if only seclected points are to be connected
*                             ! Default to include
            incpt  = .true.
*           if line_type(2) set, then see if point matches
            if( line_type(2).ne.0 .and.
     .          line_type(2).ne.pt_data(1,i) ) incpt = .false.
 
*           See if we should add
            if( incpt ) then
               j = j+1
               if( .not.map_mode ) then
                  xp(j) = sign_x*x_data(1,i)
                  yp(j) = sign_y*y_data(1,i)
               else
                  slong = x_data(1,i) + ref_valx
                  slat  = y_data(1,i) + ref_valy
                  call maptrn(slat, slong, xp(j), yp(j))
               end if
            end if
 
c****       See if break in data
            line_break = .false.
            if( j.gt.1 ) then
*                                            ! Break in line
                if( xp(j).lt.xp(j-1) .and.
     .              line_type(1).lt.0  ) then
                    line_break = .true.
                    j = j - 1
                else 
                    line_break = .false.
                end if
            end if
 
            if( (j.ge.max_plot .or. i.ge.num_data-step_data+1)
*                                                     ! plot this batch
     .          .or. line_break )                then
c
               call j2ply(j,xp,yp)
               if( j.ge.max_plot .or.line_break ) i = i -step_data
c
*                        ! reset to zero
               j = 0
*                        ! go back one point so that line will be continous
c
c....          Check to see if we have reached end of data
               if ( i.ge.num_data-1) i = num_data
c
            end if
c
c....       Check for break
C           if( ifbrk() ) then
C              call jmcur
C              call report(' User break: while drawing lines')
*              Turn off clipping before we return
C              call gsclip(0)
C              return
C           end if
 
         end do
*                 ! line to be drawn
      end if

***** Check for drawing lines bewteen all points of the same type.

      if( line_type(1).ne.0 .and. line_type(2).lt.0 ) then
*                                ! set the line type
         call jlstl(abs(line_type(1)))

*        Now scan for min and max point types
         pt_min = 99999
         pt_max = -99999
         do i = start_data, num_data, step_data
            pt_min = min(pt_min, pt_data(1,i))
            pt_max = max(pt_max, pt_data(1,i))
         end do

         pt_max = min(pt_max, abs(line_type(2)))                           
         write(*,'(a,i3,a,i5)') ' Drawing lines between point type '
     .        ,pt_min,' to ',pt_max
c        below replaced by above to avoid splitting Hollerith -- rwk 970920
c         write(*,'('' Drawing lines between point type '',i3,
c     .             '' to '',i5)') pt_min, pt_max

*        Now loop over all of the point types.
         do k = pt_min, pt_max
            write(*,'('' Drawing line '',i3)') k
c
c....       Loop over data in blocks of the largest we can fit in core
*               ! counter over data
            i = start_data - step_data
*               ! counter over data within size of max_plot
            j = 0
            do while ( i.lt. num_data )
               i = i + step_data
c
c....          Put the data in main memory arrays
*              See if only seclected points are to be connected
*                                ! Default to include
               incpt  = .true.
*              if line_type(2) set, then see if point matches
               if( line_type(2).ne.0 .and.
     .             pt_data(1,i).ne.k ) incpt = .false.

*              See if we should add
               if( incpt ) then
                  j = j+1
                  if( .not.map_mode ) then
                     xp(j) = sign_x*x_data(1,i)
                     yp(j) = sign_y*y_data(1,i)
                  else
                     slong = x_data(1,i) + ref_valx
                     slat  = y_data(1,i) + ref_valy
                     call maptrn(slat, slong, xp(j), yp(j))
                  end if
               end if

c****          See if break in data
               line_break = .false.
               if( j.gt.1 ) then
*                                               ! Break in line
                   if( xp(j).lt.xp(j-1) .and.
     .                 line_type(1).lt.0  ) then
                       line_break = .true.
                       j = j - 1
                   end if
               end if

*                                                     ! plot this batch
               if( (j.eq.max_plot .or. i.eq.num_data)
     .             .or. line_break )                then
c
                  call j2ply(j,xp,yp)
c
*                        ! reset to zero
                  j = 0
*                        ! go back one point so that line will be continous
                  i = i - step_data
c
c....             Check to see if we have reached end of data
                  if ( i.eq.num_data-1) i = num_data
c
               end if
c
c....          Check for break
C              if( ifbrk() ) then
C                 call jmcur
C                 call report(' User break: while drawing lines')
*                 Turn off clipping before we return
C                 call gsclip(0)
C                 return
C              end if
            end do
*                 ! Looping over all types of line
        end do
*                 ! line to be drawn
      end if
 
*     Turn off clipping before we return
      call jmcur
      call gsclip(0)
c
c.... Thats all
      return
      end
 
