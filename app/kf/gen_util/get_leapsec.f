CTITLE GET_LEAPSEC

      subroutine get_leapsec( jd, ut1_tai, tai_utc )

      implicit none 

*     Routine to get the number of leap seconds at particular
*     Julian date.  The expected value of ut1_tai is passed to
*     routine and if this is less than 15000 ms and is is assumed
*     that the leap seconds have already been applied.  

* PARAMETERS

* max_leap -- Maximum number of leap seconds savable
      integer*4 max_leap
      parameter ( max_leap = 50 )

* PASSED variables

* jd       -- Julian date
* ut1_tai  -- Value of ut1-at  (mas)

      real*8 jd, ut1_tai

* tai_utc -- Number of leap seconds (seconds)

      integer*4 tai_utc

* LOCAL VARIABLES

* i  -- Loop counter
* ierr -- IOSTAT error
* ijd  -- Integer JD (need 0.5 added)
* trimlen -- Length of line

      integer*4 i, ierr, ijd, trimlen 

* leapfile  -- Name for leap second file
* line      -- Line read from file
* home_dir  -- Users home directory (used to get ~/gs/tables/leap.sec)
* MOD TAH 000901: Added gg/tables as check also.

      character*128 leapfile, line, home_dir

* SAVED VARIABLES
* num_leap  -- Number of leap second values in table
* leap_jd(max_leap) -- JD of leap second values
* leap_val(max_leap) -- Values for the leap seconds
* leap_read -- Set true once leap second file has been tried

      integer*4 num_leap, leap_val(max_leap)
      real*8 leap_jd(max_leap)
      logical leap_read

      save num_leap, leap_val, leap_jd, leap_read 

      data leap_read / .false. /

****  Check the magnitude of the ut1_tai
      if( abs(ut1_tai).lt. 15000d0 ) then
          tai_utc = 0
          RETURN
      endif

****  See if we have already read the leap second file
      if( .not.leap_read ) then

****      Try to open the leap second file: Try leap.sec $HOME/gs/tables/leap.sec
      
          leapfile = 'leap.sec'
          open(104,file=leapfile, status='old', iostat=ierr)
          if( ierr.ne.0 ) then
              call getenv('HOME',home_dir)
              leapfile = home_dir(1:max(1,trimlen(home_dir))) // 
     .                   '/gs/tables/leap.sec'
              open(104,file=leapfile, status='old', iostat=ierr)

* MOD TAH 000901: If not found in gs/tables, try gg/tables.
              if( ierr.ne.0 ) then
                 leapfile = home_dir(1:max(1,trimlen(home_dir))) // 
     .                   '/gg/tables/leap.sec'
                 open(104,file=leapfile, status='old', iostat=ierr)
              end if
          end if

          if( ierr.eq.0 ) then    ! Read the leap second file
              num_leap = 1
              leap_jd(1) = 2444786.50d0
              leap_val(1) = -20
              do while ( ierr.eq.0 .and. num_leap.le.max_leap )
                 read(104,'(a)',iostat=ierr) line
                 if( ierr.eq.0 .and. line(1:2).eq.' 2' ) then
                     num_leap = num_leap + 1
                     read(line,'(1x,i7)') ijd
                     leap_jd(num_leap) = ijd + 0.5d0
                     leap_val(num_leap) = leap_val(num_leap-1)-1
                 end if
              end do
          else
              num_leap = 0
              call report_error('IOSTAT',ierr,'open',leapfile,0,
     .                          'get_leapsec')
          endif
          leap_read = .true.
      endif 

****  See if we have anything
      if( num_leap.eq.0 ) then              

*         No Leap-second values, then 'wing-it' from the values
*         we have
          if( ut1_tai.ne.1.d5 ) then
             tai_utc = nint(ut1_tai/15000.d0)
          else
             tai_utc = 0.d0
          end if
      else

*         See where we are in the tables
          tai_utc = 0
          do i = 1, num_leap
C            if( jd.ge.leap_jd(i) .and. jd.lt.leap_jd(i+1) ) 
C    .                 tai_utc = leap_val(i)
* MOD TAH 090101: "Fixed" small rounding error problem when jd
*            is exactly on leap second boundary.  Here we allow
*            ~ 15 second slack because of gamit midpoints.
             if( jd.ge.leap_jd(i)-2.d-4 .and. jd.lt.leap_jd(i+1) ) 
     .                 tai_utc = leap_val(i)
          enddo 
      end if

****  Thats all
      return
      end





