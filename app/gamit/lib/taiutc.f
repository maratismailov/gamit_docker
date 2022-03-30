Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      REAL*8 FUNCTION TAIUTC(PJD)
C
C     Compute International Atomic Time (TAI) minus UTC (seconds)
c     note: the leap second is applied at 00:00:00 in UTC time.

c     The input argument is 'PEP' Julian day; i.e., the Julian day number of
c     the true Julian date beginning at noon on the UTC day:
c       PJD = JD  + 0.5  = Modified JD + 2440000 + 1

      implicit none

      logical ifirst

      integer*4 pjd,pjd0,pjdn,pjdnt,pjd_last,ios,len,rcpar

      real*8 taiutc_last

      character*8    leapf,lowerc   
      character*15   stpdt
      character*80   prog_name
      character*256  message

      data pjd0/2444971/ pjd_last/0/,taiutc_last/0/,ifirst/.true./

      save pjd_last,taiutc_last,ifirst,ios

c     simple logic to save time:  open the file only once
c                                 reread the table whenever the day changes

c       print *,'pjd,ifirst,pjd_last',pjd,ifirst,pjd_last

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

      taiutc = 20.0d0

      if ( ifirst ) then
c        use a default name (link) for the table
         leapf = 'leap.sec'
         leapf = lowerc(leapf)
         open (unit=99,file=leapf,status='old',iostat=ios,err=5)
    5    if( ios.ne.0 ) then
            call report_stat('FATAL',prog_name,'lib/taiutc',' '
     .   ,'Cannot open file leap.sec--is it linked in this directory?'
     .       ,0)
         endif
         ifirst = .false.
         endif

      if( pjd .eq. pjd_last ) then
        taiutc = taiutc_last

      else

      pjd_last = pjd
      rewind 99

c     first line is comment
      read (unit=99,fmt=20,iostat=ios,err=100)

c     second line contains last valid date
      read (unit=99,fmt=20,iostat=ios,end=100) pjdnt,stpdt
   20 format(29x,I7,3x,a15)

   10 continue
      read (unit=99,fmt=30,iostat=ios,err=100,end=100) pjdn
   30 format(1x,I7)

  100 continue

      if (ios .ne. int(0)) then
         call report_stat('FATAL',prog_name,'lib/taiutc',leapf
     .                   ,'Error in leap second file:',ios)
      else if (pjd.gt.pjdnt) then
         write(message,998) pjd,stpdt
  998    format ('Date for TAI-UTC (',i8
     .       ,') after stop date in leap.sec: ',a11)
         call report_stat('FATAL',prog_name,'lib/taiutc',' '
     .                    ,message,0)
      else if (pjd .lt. PJD0) then
        write(message,996) pjd
  996    format('Cannot calculate TAI-UTC for date ('
     .          ,i15,')  before Jan 1982')
         call report_stat('FATAL',prog_name,'lib/taiutc',' '
     .                    ,message,0)
      else if (pjd .gt. pjdn) then
         taiutc = taiutc + 1.0d0
         goto 10
      else
         continue
      endif

      endif


      taiutc_last = taiutc

c        print *,'ifirst,pjd,pjd_last,taiutc,taiutc_last'
c     .          ,ifirst,pjd,pjd_last,taiutc,taiutc_last

      return

      end

