CTITLE plts9
 
      subroutine plts9
c
c     This segment does some of the miscellaneous file outputting
c     operations.  At moment this is mainly the HELP command.
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
*                          ! the main ema declaration
      include 'plot_ema.h'
c
c.... Go to the routine specified by pel
      if( pel.eq.5 ) then
          call read_data( ema_data(ibak_array), ema_data(ix_array),
     .                    ema_data(iy_array), ema_data(ipt_array)  )
      end if
 
      if( pel.eq.30 .or. pel.eq.29 ) then
         call help
      end if
 
      if( pel.eq.32 ) then
          call status
      end if
c
c.... Thats all
*                   ! return to main segment
      return
c
      end
 
