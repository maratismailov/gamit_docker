CTITLE TIMR
 
 
      real*4 function timr( mode )

      implicit none
 
*     Function to emulate the LIBHS routine timr.  The mode is 0 to set
*     the clock, and 1 to get a read out of the time.
*     HERE we use we use the Sun dtime routine.
*     Resolution is only good to nearest second.
 
*   mode        - Mode of call, 0 to reset, 1 to get time since
*               - last reset.
 
      integer*4 mode
 
* Local variables
 
*   last_time   - seconds since start of day when mode 0 last used 
*   date(5)     - Current date with hours and minutes
*   sectag      - Seconds tag
 
      integer*4 date(5)
      real*4 last_time
      real*8 sectag
 
      save last_time
 
*                             ! reset timer, amnnd start clock
      if( mode.eq.0 ) then
          call systime( date, sectag)
          last_time = date(4)*3600.d0 + date(5)*60.d0 + sectag
          timr = 0
      else
          call systime( date, sectag )
          timr = date(4)*3600.d0 + date(5)*60.d0 + sectag - last_time
          if( timr.lt.0.0 ) then
              timr = timr + 86400.d0
          end if
      end if
 
****  Thats all
      return
      end
 
