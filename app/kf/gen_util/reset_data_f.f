CTITLE RESET_DATA_FLAG
 
      subroutine reset_data_flag ( data_flag )
 
      implicit none 

 
*     Routine to reset the user parts of the data flag.  (These are marked
*     with '*' in OBS_DATA.FTNI description).  As the contributions and
*     calibrations are applied the bits will be set in need be
*                                  3:08 PM  THU., 16  APR., 1987
*
 
 
      integer*4 data_flag, cand
 
*     Use SBIT to reset
 
*     call sbit(data_flag, 2, 0)  ! Ion bit
*     call sbit(data_flag, 3, 0)  ! WVR antenna temperatures
*     call sbit(data_flag, 4, 0)  ! WVR data base delay values
*     call sbit(data_flag, 6, 0)  ! Cable calibration bad
*     call sbit(data_flag,12, 0)  ! User source delete
*     call sbit(data_flag,13, 0)  ! User time/site dounw weight
*     call sbit(data_flag,15, 0)  ! Data do not close
*     call sbit(data_flag,16, 0)  ! Set to data available
 
*     The above code is achieved using
 
c mod simon 2/22/2002 (layhey compiler compat problem) convert oct 3701 to integer 32768.
c      data_flag = cand( data_flag,  o'3701')
      data_flag = cand( data_flag,  1985)
 
      return
      end
 
