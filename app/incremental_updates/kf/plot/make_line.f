CTITLE MAKE_LINE
 
      subroutine make_line ( value, field, outline)
 
      implicit none
 
*     Routine to add current values to output line
 
*   date(6)     - Calender date
*   field(3)    - Type, start and errorbar column for field
*   i,j         - Loop counters
*   ierr        - IOSTAT error
*   lenline     - Length of output line
*   trimlen     - HP function for length of string
 
      integer*4 date(6), field(3), i, j, ierr, lenline, trimlen
      integer*4 absf3
 
*   sec_tag     - Seconds tage of epoch
*   value(2)    - Two values to be output
 
      real*8 sectag, value(2)
 
*   outline     - output line to be added to
 
      character*(*) outline
 
*     Get current position in line
      lenline = trimlen(outline)
 
*     See if time field
*                                 ! Time field
* MOD TAH 200331: Added field 3 date type field.
      if( field(1).eq.0 .or. field(1).eq.3 ) then
          call mjd_to_ymdhms( value, date, sectag)
          if( field(3).gt.0 ) then 
              absf3 = field(3)
          else
              absf3 = abs(field(3))+1
          endif
 
          if( absf3.lt.6 ) then
              write(outline(lenline+2:),100)
     .            (date(j),j=1,absf3)
          else
              write(outline(lenline+2:),100)
     .            (date(j),j=1, 5), sectag
 100          format(i4," ",i2," ",i2,1x,i2," ",i2,1x,F6.2)
          endif
*                                 ! Output normal values
      else
          if( abs(value(1)).lt.1.d6 .and. abs(value(1)).gt.1.d-3 ) then
              write(outline(lenline+2:),120, iostat=ierr)
     .             value
 120          format(2(1x,f16.6))
          else 
              write(outline(lenline+2:),140, iostat=ierr)
     .             value
 140          format(2(1x,d16.7))
          end if
      end if
 
****  Thats all
      return
      end
