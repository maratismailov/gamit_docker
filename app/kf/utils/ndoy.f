      program ndoy

      implicit none 
 
*     This is generic program which will read a varierty of
*     date types from the runstring and return the other date
*     types. 

* NDOY: is a keyboard input, shortended output version of doy. ZAV 000630.
 
*   date(5)     - Calender date
*   day_of_year - Day of year
*   rcpar       - Get runstring parameters
*   len_run     - Length of string entry
*   ierr        - IOSTAT error flag
*   num_run     - Number of arguments in the runstring
*   i           - Loop counter
*   gps_week    - GPS week number starts with week zero on 1980/1/5 Midnight
*   gps_dow     - GPS doy of week
*   gps_sow     - GPS seconds into week
*   gps_start_mjd - MJD of the start of the GPS week count (80/1/5)
*   indx        - Postition of character in string

c   Add Decimal Year in output, K. Feigl May 97
c   MOD TAH 970520: Changed gps_dow to run from 0-6 (instead of 1-7).
c   MOD SCM 971029: Changed some of the math from REAL*8 to INTEGER*4 for LINUX.
c   MOD SCM 971115: Fixed more problems with mixing REAL and INTEGER arithmatic.
 
      integer*4 date(5), day_of_year,ierr, num_run, i,
     .          gps_week, gps_dow, gps_sow, gps_start_mjd , indx, iyr

      integer*4 irunstr, iswitch, iend, istart
 
*   jd          - Julian date
*   mjd         - Modified Julian date.
*   sectag      - Seconds tag
 
      real*8 jd, mjd, sectag,dyear
 
*   runstring(5)    - Five elments of the string.
 
 
      character*40 runstring(5)
      character inp_str*40

      data gps_start_mjd / 44243 /
 
  5   read (*,10,end=200) inp_str
 10   format (A)
  
      irunstr = 0
      iswitch = 0
      i = 1
      istart = 1

      do  while ( i.le.40 )
         if ( inp_str ( i:i ) .eq. ' ' ) then
            if ( iswitch.ne.0 ) then
               if ( irunstr.ne.5 ) then
                  iend = i - 1
                  iswitch = 0
                  irunstr = irunstr + 1
                  runstring ( irunstr ) = inp_str ( istart:iend )
               end if
            end if
         else
            if ( iswitch.eq.0 ) then
               istart  = i
               iswitch = 1
            end if
         end if
         i = i + 1
      end do

      num_run = irunstr

****  Start looping over the runstring to see how many arguments
*      i = 0
*      len_run = 1
* 
*      do while ( len_run.gt.0 .and. i.lt.5 )
*          len_run = rcpar(i+1,runstring(i+1))
*          if( len_run.gt.0 ) i = i + 1
*      end do

*      num_run = i
 
****  See if we got any
      if( num_run.eq.0 ) then 
          call proper_runstring('doy.hlp','doy',0)

c         time is not GPS (see below)
          indx = 0
       endif
 
****  Based on number, decoade the results
 
*     Set to 0 hr, 0 min
      date(4) = 0
      date(5) = 0
      sectag  = 0

*     Check to see if the first argument has a W in it.  If it
*     does assume that it is GPS week.
      call casefold( runstring(1) )
      indx = index(runstring(1),'W')
      if( indx.gt.0 ) then

*         GPS week number passed.  Replace the W with a blank
*         and decode
          call sub_char(runstring(1),'W',' ')
          read(runstring(1),*,iostat=ierr) gps_week

*         Now see if gps_dow or gps_sow passed as second runstring
          if( num_run.eq.2 ) then
              read(runstring(2),*,iostat=ierr) gps_dow
*             check the size
              if( gps_dow.gt.6 .or. gps_dow.eq.0 ) then
                  gps_sow = gps_dow
                  gps_dow = gps_sow/86400 
              else
                  gps_sow = (gps_dow)*86400
              end if
          else
              gps_dow = 0
              gps_sow = 0
          end if

*         Now compute the other quantities (add 1 because day of week
*         runs from 1 to 7.
          mjd = dble(gps_start_mjd + gps_week*7 + gps_sow/86400 + 1)
          jd  = mjd + 2400000.5d0
          call jd_to_ymdhms(jd, date, sectag)
          call ymd_to_doy(date, day_of_year)

      else
*         Original conversions.

          if( num_run.eq.1 ) then
 
*             Take to Julian date
              read(runstring(1),*,iostat=ierr) jd
              if( jd.lt.100000 ) jd = jd + 2400000.5d0
              call jd_to_ymdhms(jd, date, sectag)
              call ymd_to_doy(date, day_of_year)
 
          else if( num_run.eq.2 ) then
              read(runstring(1),*,iostat=ierr) date(1)
*                         ! January
              date(2) = 1
              read(runstring(2),*,iostat=ierr) date(3)
              call ymdhms_to_jd(date, sectag, jd)
              day_of_year = date(3)
              call jd_to_ymdhms(jd, date, sectag)
          else if( num_run.ge.3 ) then
              read(runstring(1),*,iostat=ierr) date(1)
              read(runstring(2),*,iostat=ierr) date(2)
              read(runstring(3),*,iostat=ierr) date(3)
              if( num_run.ge.4 ) 
     .        read(runstring(4),*,iostat=ierr) date(4)
              if( num_run.ge.5 ) 
     .        read(runstring(5),*,iostat=ierr) date(5)
              
              call ymdhms_to_jd(date, sectag, jd)
              call ymd_to_doy(date, day_of_year)
          else
               write (*,'(a)') '***TODAY*** IS: '
              call systime (date,sectag)
              call ymdhms_to_jd(date, sectag, jd)
              call ymd_to_doy(date, day_of_year)
          end if

*         Now compute the gps date quanities
          mjd = jd - 2400000.5d0
          gps_week = (idint(mjd) - gps_start_mjd - 1)/7
*         This test of if date is before start of gps time.
*         (Usual problem with fortan set int(-0.99) to 0)
          if( mjd-gps_start_mjd-1 .lt.0 ) gps_week = gps_week - 1
          gps_sow  = (idint(mjd) - (gps_start_mjd+gps_week*7+1))*86400
*         This adjustment is also for negative int values.
          if( mjd-gps_start_mjd-1 .lt.0 .and.
     .        gps_sow.ge. 604800 ) gps_sow = gps_sow - 604800
          gps_dow = gps_sow/86400 
      end if

c     calculate decimal year, dealing with leap years
      iyr = date(1)
      if (mod(iyr,4) .eq. 0) then
         dyear = dble(iyr) + dble(day_of_year-1)/366.d0
      else
         dyear = dble(iyr) + dble(day_of_year-1)/365.d0
      endif
           
****  Now write results
      write(*,100) (date(i),i=1,3), day_of_year, jd, mjd, gps_week,
     .             gps_dow, gps_sow
 100  format(i4,' ',i2.2,' ',i2.2,
     .       ' ',i3,' ',F11.2,' ',F9.2,
     .       ' ',i5,' ',i2,' ',i6)
 
      goto 5

 200  continue
      end
 
 

