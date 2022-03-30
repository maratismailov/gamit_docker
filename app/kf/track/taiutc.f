Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      REAL*8 FUNCTION TAIUTC(PJD)
C
C     Compute International Atomic Time (TAI) minus UTC (seconds)
c     note: the leap second is applied at 00:00:00 in UTC time.

c     The input argument is 'PEP' Julian day; i.e., the Julian day number of
c     the true Julian date beginning at noon on the UTC day:
c       PJD = JD  + 0.5  = Modified JD + 2440000 + 1

      implicit none
* MOD AZ 190305: include Header for shared variable "leapsec_path"
*             All related "open" actions are adjusted from "leapf"
*             to "leapsec_path"
* MOD TAH 200531: Cleaned up AZ code.  Removed unnecessry fotran, 
*     removed GOTO statment, fixed indenting; improved error reporting
*     so that when there is a failure user knows what happened; 
*     Fixed unit number conflict.
* MOD MAF 210624: Changed unit number for leap.sec to 96 from 98 to avoid
*     clash with summary file defined in init_track.f.

      include 'track_com.h'

      logical ifirst
      logical found   ! Set true if JD greater than current value

      integer*4 pjd,pjd0,pjdn,pjdnt,pjd_last,ios

      real*8 taiutc_last

      character*8    leapf,lowerc   
      character*15   stpdt
      character*128  line   ! Line read from leap.sec file.

      data pjd0/2444971/ pjd_last/0/,taiutc_last/0/,ifirst/.true./

      save pjd_last,taiutc_last,ifirst,ios

c     simple logic to save time:  open the file only once
c                                 reread the table whenever the day
c                                 exceeds entry used.

* MOD TAH 200531: See if data too early.
      if (pjd .lt. PJD0) then
          write(*,50) pjd, pjd0
 50       format('Data PepJD ',i7,' before start of leap.sec ',i7)
          stop 'TRACK: Before start of leap.sec'
      endif

*     Set leap seconds at start of file.
      taiutc = 20.0d0

      if ( ifirst ) then
         open (unit=96,file=leapsec_path,status='old',
     .         iostat=ios)
         if( ios.ne.0 ) then

            write(*,153) ios, trim(leapsec_path)
153         format('ERROR: ',I6,' open/reading leap.sec ',a)
            stop 'TRACK: leap.sec out of order'
         endif
         ifirst = .false.
      endif

* MOD TAH 200531: Changed  pjd_last to be last table value found
*     so that we don't keep re-reading file.
      if( pjd .lt. pjd_last ) then
         taiutc = taiutc_last
         RETURN
      else

* MOD TAH 200531: This logic is bad (no need to re-wind file if
*        structured correctly but this version works so keep.
         rewind 96
c        first line is comment
         read (96,'(a)',iostat=ios) line

c        second line contains last valid date
         read (96,'(a)',iostat=ios) line
         read (line,20,iostat=ios) pjdnt,stpdt
   20    format(29x,I7,3x,a15)
* MOD TAH 200531: See if our date is past end of file
         if( ios.ne.0 .or. pjd.gt.pjdnt ) then
            write(*,120) ios, pjd, pjdnt
 120        format('IOS Error ',i5,' or data PepJD ',i7,
     .             ' after end of leap.sec ',i7)
            stop 'TRACK: Past end of leap.sec or read error'
         endif

* MOD TAH:200531: Removed use of GOTO statement.
         found = .false.
         do while ( ios.eq.0 .and. .not. found)          
            read (96,'(a)',iostat=ios) line
            read (line,30,iostat=ios) pjdn
   30       format(1x,I7)

            if (ios .ne. 0 ) then
               write(*,153) ios, line
               stop 'TRACK: Error reading leap.sec'
            elseif (pjd .gt. pjdn) then
               taiutc = taiutc + 1.0d0
            else
* MOD TAH 200531: Save table JD so we can test next call
               pjd_last = pjdn 
               found = .true.
            endif
         enddo  

      endif

*     Save value so can be used next call
      taiutc_last = taiutc

      return

      end

