C********************************************************************
c
c     TIME CONVERSION ROUTINES ARE DOWN HERE
C
C     I tested them for every second of GPS weeks 180 to 365 -kf 871118
c
c     Think before you change these routines.
c
C********************************************************************
      subroutine timcon(itflag,
     .  iwk0,sow0,
     .  iyr0,idn0,ihr0,imn0,sec0,
     .  utcoff)

c     Convert gps time expressed in gps weeks and
c     seconds of week to UTC expressed as year, day_ofyear,
c     hour, minute and second  and vice versa

c     calls:  SECSUM,DAYJUL,DS2HMS
c
c     written by peter morgan for the Apollo, january 1987
c     modified by Kurt Feigl to work june 22, 1987
c     cleaned up and tested by Kurt Feigl November 8, 1987
c     finished beta-test: kurt Feigl: November 16, 1987
c     changed to integer flag: Kurt Feigl April 30, 1989
c
c
c     Tested for every second gps weeks 180 to 365.
c     Change this routine lightly and I'll kill you.

c     Bug found and fixed in igpsdow calculation for the case itflag=4
c     (no conversion).  Not changed lightly.  I've removed much of the
c     comment on the original tests since I found it not instructive.
c     Bob King, 24 Jan 92 .

c
c                      PARAMETER DESCRIPTION
c                      *********************
c
c     itflag:      input:  output:
c     +/-1        GPST    UTC
c     +/-2        UTC     GPST
c     +/-3        UTC     UTC
c     +/-4        GPST    GPST
c
c     if itflag is positive:
c     input:
c              iwk0         GPS week number
c              sow0         second of week
c     output:
c              iyr0         year
c              idn0         day number (day of year)
c              ihr0         hour
c              imn0         minutes
c              sec0         seconds
c
c     if itflag is negative:
c     input
c              iyr0         year
c              idn0         day number (day of year)
c              ihr0         hour
c              imn0         minutes
c              sec0         seconds
c     output:
c              iwk0         GPS week number
c              sow0         second of week
c
c
c     older calling sequence:
c     direction = .true.  ==> itflag = +1
c     direction = .false. ==> itflag = -2
c
c    TEST_TIME.DOC
c    Test the time conversion routine for all gps time.
c
c    The test consists of converting a GPS time in weeks and seconds of week
c    into UTC, and back again.  If the answers are different, print them here.
c    The only errors seem to occur at leap seconds.
c
c    I am going to install this new routine in MAKEX on
c    Wednesday, November 18, 1987   10:25:06 am (EST)
c
      implicit none

      integer*4 iyr0,idn0,ihr0,imn0,
     .          iwk0,igpsday,igpsdow,
     .          jd,jd0,jd1,julday,iyr,idn,ihr,imn,iwkn,
     .          itflag,len,rcpar


      real*8 taiutc,sow0,utcoff,sec0,utcoff2,sod,
     .       sow,sc

      character*80  prog_name
      character*256 message
          

      IF (itflag .gt. 0) THEN
         iwkn = iwk0
         sow  = sow0

c        convert from iwk0, sow0
c                to   iyr0, idn0, ihr0, imn0, sec0
c
c        Day zero of week zero of GPS time is 5 January, 1980
c        which was a Saturday.  Thus, Sunday, 6 January, 1980
c        is day 1 of week 0.

         jd0 = julday(1,5,1980)

c        calculate day of GPS week
         igpsdow = int(sow/86400.d0) + 1

         jd  = jd0 + 7*iwkn + igpsdow

         if (itflag .eq. 1) then          
c           Go from GPST to UTC
c           get the UTC offset according to JD in GPS time
            utcoff = taiutc(jd) - 19.d0
c           subtract the UTC offset
            call secsum(iwkn,sow,-utcoff,iwkn,sow)

c           make sure we have not gone over a week!
c           calculate day of GPS week
            igpsdow = int(sow/86400.d0) + 1

            jd  = jd0 + 7*iwkn + igpsdow

c           now get the UTC offset according to the JD in UTC
            utcoff2 = taiutc(jd) - 19.0D0
            if (utcoff2 .ne. utcoff) then
                call secsum(iwkn,sow,-(utcoff2 - utcoff),
     .                       iwkn,sow)
                igpsdow = int(sow/86400.d0) + 1
                jd  = jd0 + 7*iwkn + igpsdow
                utcoff = utcoff2
c
c               deal with the super - weird case: GPS time maps onto
c               a UTC leap second.  At 00:00:00 UTC a leap second offset
c               is added.  In other words, 00:00:00 becomes 00:00:01, and
c               there is an undefined "hole" in the time scale, because the
c               first minute really contains 61 seconds.
                utcoff2 = taiutc(jd) - 19.0D0
                if (utcoff2 .ne. utcoff) then
c                 get calling program name and m-file name for report_stat
                  len = rcpar(0,prog_name)
                  call report_stat('WARNING',prog_name,'lib/timcon',' '
     .      , 'Time conversion mapped onto a leap second at midnight',0)
                endif
            endif
         else if (itflag .eq. 2) then
c           Go from UTC to GPST               
c           get the GPST offset according to JD in UTC time
            utcoff = taiutc(jd) - 19.d0
c           add the GPST offset
            call secsum(iwkn,sow,utcoff,iwkn,sow)
c           make sure we have not gone over a week!
c           calculate day of GPS week
            igpsdow = int(sow/86400.d0) + 1
            jd  = jd0 + 7*iwkn + igpsdow
c           now get the UTC offset according to the JD in UTC
            utcoff2 = taiutc(jd) - 19.0D0
            if (utcoff2 .ne. utcoff) then
                call secsum(iwkn,sow,(utcoff2 - utcoff),
     .                       iwkn,sow)
                igpsdow = int(sow/86400.d0) + 1
                jd  = jd0 + 7*iwkn + igpsdow
                utcoff = utcoff2
c               deal with the super - weird case: GPS time maps onto
c               a UTC leap second.  At 00:00:00 UTC a leap second offset
c               is added.  In other words, 00:00:00 becomes 00:00:01, and
c               there is an undefined "hole" in the time scale, because the
c               first minute really contains 61 seconds.
                utcoff2 = taiutc(jd) - 19.0D0
                if (utcoff2 .ne. utcoff) then
c                 get calling program name and m-file name for report_stat
                  len = rcpar(0,prog_name)
                  call report_stat('WARNING',prog_name,'lib/timcon',' '
     .      , 'Time conversion mapped onto a leap second at midnight',0)
                endif
            endif
         else if (itflag .eq. 3 .or. itflag .eq. 4) then
c           Do not add leap second
            utcoff = 0.0d0
         else
           len = rcpar(0,prog_name)
           write(message,'(a,i6)') 'Undefined time conversion: ',itflag
           call report_stat('FATAL',prog_name,'lib/timcon',' '
     .                       ,message,0)
         endif

c        now get day of year
         call dayjul(jd,iyr0,idn0)   

c        seconds of day
         sod = sow - 86400.d0*dble(igpsdow-1)

c        get HMS format
         call ds2hms(iyr0,idn0,sod,ihr0,imn0,sec0)

c        do some error checking here
* MOD TAH 210101: Updated 2020 to 2100.
         if ( iyr0  .lt. 1980   .or. iyr0  .gt. 2100 .or.
     .        idn0  .lt. 1    .or. idn0  .gt. 366 .or.
     .        ihr0  .lt. 0    .or. ihr0  .gt. 23  .or.
     .        imn0  .lt. 0    .or. imn0  .gt. 59  .or.
     .        sec0  .lt. 0.d0 .or. sec0  .ge. 60.d0 ) then
           len = rcpar(0,prog_name)
           write(message,'(a,i3,i4,f10.2,2i8,i5,3i3,f6.2)')
     .          'Time conversion error: '
     .          ,itflag,iwkn,sow,jd,jd0,iyr0,idn0,ihr0,imn0,sec0
           call report_stat('FATAL',prog_name,'lib/timcon',' '
     .                       ,message,0)
         endif

c        print some debugging tidbits

c        print *, 'itflag    ',itflag
c        print *, 'iwkn      ',iwkn
c        print *, 'sow       ',sow
c        print *, 'igpsdow   ',igpsdow
c        print *, 'jd        ',jd
c        print *, 'jd0       ',jd0
c        print *, 'iyr0      ',iyr0
c        print *, 'idn0      ',idn0
c        print *, 'ihr0      ',ihr0
c        print *, 'imn0      ',imn0
c        print *, 'sec0      ',sec0
c        print *, 'utcoff    ',utcoff
c        print *, 'sod       ',sod


c--------------------------------------------------------------------

      ELSE if (itflag .lt. 0) then

c        convert  from: iyr0,idn0,imn0,sec0
c                 to    iwk0,sow0

         iyr = iyr0
         idn = idn0
         ihr = ihr0
         imn = imn0
         sc  = sec0

c        day zero of week zero of GPS time is 5 January, 1980
c        which was a Saturday

         jd0 = julday(1,5,1980)

c        julian day of Jan 1st of this year
         jd1 = julday(1,1,iyr)

c        julian date
         jd = idn + jd1 - 1
                             
c        get the UTC offset
         if (itflag .eq. -1) then
            utcoff = -1.0d0 * (taiutc(jd) - 19.0)
         else if (itflag .eq. -2) then
            utcoff = taiutc(jd) - 19.0         
         else if (itflag .eq. -3 .or. itflag .eq. -4) then
            utcoff = 0.0d0
         else
           len = rcpar(0,prog_name)
           write(message,'(a,i6)') 'Undefined time conversion: ',itflag
           call report_stat('FATAL',prog_name,'lib/timcon',' '
     .                      ,message,0)
         endif

c**      This call to (now obsolete) sumday replaced by rwk  95/2/10
c**      add the UTC offset
c**      ignore the case of going over years
c**      call sumday(idn,sc,utcoff,idn,sc)

c        julian date
         jd = idn + jd1 - 1

c**      new rwk 95/2/10
         call timinc(jd,sc,utcoff)

c        days in gps time
         igpsday  = jd - jd0

c        gps week number
         iwk0 = igpsday/7

c        day of gps week
         igpsdow = igpsday - 7*iwk0

         sow0 = dble(igpsdow-1) * 86400.d0 +
     .            dble(ihr)       *  3600.d0 +
     .            dble(imn)       *    60.d0 +
     .            sc

c        make sure we have not gone over a week
         call secsum(iwk0,sow0,0.d0,iwk0,sow0)

c        do some error checking here
         if ( iwk0 .lt. 0    .or. iwk0 .gt. 5000 .or.
     .        sow0  .lt. 0.d0 .or. sow0  .ge. 7.d0*86400.d0) then
           len = rcpar(0,prog_name)
           write(message,'(a,i3,i6,3i4,f6.2,i2,i4,f10.2)')
     .          'Time conversion error: '
     .          ,itflag,iyr,idn,ihr,imn,sc,igpsdow,iwk0,sow0
           call report_stat('FATAL',prog_name,'lib/timcon',' '
     .                      ,message,0)
         endif

c        print some debugging tidbits

c        print *, 'itflag    ',itflag
c        print *, 'jd0       ',jd0
c        print *, 'jd1       ',jd1
c        print *, 'jd        ',jd
c        print *, 'igpsday   ',igpsday
c        print *, 'iwk0      ',iwk0
c        print *, 'sow0      ',sow0
c        print *, 'igpsdow   ',igpsdow
c        print *, 'utcoff    ',utcoff
c        print *, 'iyr       ',iyr
c        print *, 'idn       ',idn

      ENDIF

      return
      end
