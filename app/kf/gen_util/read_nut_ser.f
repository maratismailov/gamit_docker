
CTITLE 'read_nut_series'
 
      subroutine read_nut_series( unit, nut_name, plan_name, ierr )
 
 
      implicit none 

*     Routine to read the standard Earth nutation series and the
*     coefficients for the Planetary nutation as well.
 
      include '../includes/nut_eval.h'
 
*         ierr      - IOSTAT error
*         unit      - Unit number to use for reading files.
*         trimlen   - Length of string,
*         i,j       - Loop counters
 
      integer*4 ierr, unit, trimlen, i, j, jerr

      character*(*) nut_name, plan_name

*        period     - Period of the nutation in days. NOT USED but in file.

      real*8 period
 
*        line  - line read from file (non-blank first
*              - character is comment)
 
      character*132 line

****  Read the standard ZMOA-1990 nutation series  
*     Open file
      if( ichar(nut_name(1:1)).eq.0 ) then
          nut_name = ' '
      end if
      if( trimlen(nut_name).eq.0 ) then
          ierr = -1
      else
          open(unit, file=nut_name , iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open', nut_name, 1,
     .                      'read_nut_series')
      end if

      num_series = 0
      j = 0

      line = ' '
 
      do while ( ierr.eq.0 )

          read(unit,'(a)', iostat=ierr, err=100, end=100 ) line
 100      continue
 
*         Decode if no error and not a comment
*                                          ! Process
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 )  then

              j = j + 1

*             Get the fundamental arguments and arguments from the line
              read(line,*) (arg_series(i,j),i=1,5), period, 
     .                     val_nut_coeff(1,j), val_nut_rate(1,j),
     .                     val_nut_coeff(3,j), val_nut_rate(3,j),
     .                     val_nut_coeff(2,j), val_nut_rate(2,j), 
     .                     val_nut_coeff(4,j), val_nut_rate(4,j) 
          end if
*                 ! Looping over file
      end do

*     Save the number of terms
      num_series = j
      series_read = .true.

      close(unit=unit, iostat=ierr )
 
****  Get to see if optional Planetary nutations are to be read
*     Open up the planetary nutation series.
      if( ichar(plan_name(1:1)).eq.0 ) then
          plan_name = ' '
      end if
      if( trimlen(plan_name).eq.0 ) then
          jerr = -1
      else 
          open(unit, file=plan_name, iostat=jerr, status='old')
          call report_error('IOSTAT',jerr,'open', plan_name , 0,
     .                      'read_nut_series')
      end if
 
      num_plan = 0
      j = 0
 
      do while ( jerr.eq.0 )
          read(unit,'(a)', iostat=jerr ) line

*         Decode if no error and not a comment
*                                          ! Process
          if( jerr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 )  then
              j = j + 1
              read(line,*) (arg_plan(i,j),i=1,5), period,
     .              val_plan(1,j), val_plan(3,j),
     .              val_plan(2,j), val_plan(4,j)
          end if
*                     ! Looping over file
      end do
      plan_read = .true.
      num_plan = j
 
***** Thats all
      return
      end
 
