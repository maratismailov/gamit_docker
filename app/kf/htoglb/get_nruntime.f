CTITLE GET_NRUNTIME
 
      subroutine get_nruntime( line, indx, qrun_time)
 
      implicit none

*     This routine will read the run date from the first line of the
*     NGS format and convert it to the runtime used by GLOBK
 
* PASSED VARIABLES
 
*   indx        - Pointer to start of day field
*   qrun_time(7)    - Runstime for the solution
 
      integer*4 indx, qrun_time(7)
 
*   line        - Line read from the input (this is the full
*               - line)
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   yr, mon, day, hr, min   - Date and time of runtime
*   jndx        - Position in string
*   ierr        - IIOSTAT error reading line
*   trimlen     - Length of string
 
      integer*4 yr, mon, day, hr, min, jndx, ierr, trimlen
 
*   sec         - Seconds part of value
 
      real*4 sec
 
*   month_names(12) - Names of the months
*   cmon            - Month as character string
 
      character*3 month_names(12), cmon
 
      data month_names / 'JAN','FEB','MAR','APR','MAY','JUN',
     .                'JUL','AUG','SEP','OCT','NOV','DEC' /
 
****  Read the values from the line
 
      read(line(indx:),100,iostat=ierr) day, cmon, yr, hr, min,
     .                                 sec
 100  format(i2,1x,a3,4(1x,i2))
 
      call report_error('IOSTAT',ierr,'read',line,0,'get_nruntime')
      if( ierr.eq.0 ) then
          jndx = 1
          call get_cmd(cmon, month_names, 12, mon, jndx)
          if( mon.le.0 ) then
              write(*,120) cmon, line(1:trimlen(line))
 120          format(' Month string ',a3,' not found.  Runtime from',
     .                /,a)
              mon = 0
          end if
 
*****     Save the runtime
          qrun_time(1) = yr + 1900
          qrun_time(2) = mon
          qrun_time(3) = day
          qrun_time(4) = hr
          qrun_time(5) = min
          qrun_time(6) = nint(sec)
          qrun_time(7) = (sec-nint(sec))*100
      else
 
*         Save some dummy values
          qrun_time(1) = 1900
          qrun_time(2) = 0
          qrun_time(3) = 0
          qrun_time(4) = 0
          qrun_time(5) = 0
          qrun_time(6) = 0
          qrun_time(7) = 0
      end if
 
***** Thats all
      return
      end
 
 
 
 
 
 
