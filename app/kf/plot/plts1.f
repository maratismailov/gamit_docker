      subroutine plts1

      implicit none 
c
c     This program will set up or turn off the plotting
c     package and/or read a data file depending of the values
c     stored in PEL in the labelled PLOT_CONTROL
c
c     Control of the program though PEL
c     ---------------------------------
c     PEL  Action
c     -1   Just set up GRAPHICS/1000
c     -2   Terminate plotting with GRAPHICS/1000
c     -3   Set up and read a new file (if name has been given)
c     -4   just read the new file
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
*                          ! the ema declaration for pltsl
      include 'plot_ema.h'
 
c.... All we do is check to see if we call setup of fininish up
      if( pel.eq.-1 ) then
         call setup         
      end if
c
c.... See if finish up
      if( pel.eq.-2 ) then
         call finish_up
      end if
c
c.... See if read file
*                                            ! OK read a file (see if setup
      if( pel.eq.-3 .or. pel.eq.-4  ) then
c                                            also needed
         if( pel.eq.-3 ) then
            call setup
         end if
c
         call read_file(ema_data,.false.)
c
      end if
c
c.... Thats all for this segment         
      return
c
      end
 
