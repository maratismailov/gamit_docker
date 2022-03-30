      subroutine tagchk( timeset,weekset,ibadtag,epoch
     .                 , iwkntag,sowtag,igpswk,gpssec
     .                 , iwknstart )
C
c     Subroutine to check the validity of time tags and set the
c     appropriate flags.  - R. King February 1990 from MAKEX code

      implicit none
c
      include '../includes/makex.h'

      logical      timeset,weekset

      integer*4    ibadtag,epoch,iwkntag,igpswk,iwknstart

      real*8       sowtag,gpssec,secdif

      character*256 message

c     timeset = true     if the time tag is any good
c                         This is not set to true until we have
c                         found a good epoch.

c     weekset = true     if the week number has not changed or been lost

c     iwkntag, sowtag  =  the current
c     igpswk  , gpssec   =  the most recently read values
c     iwknstart         =  the week number from the file name (use if necessary)

c     ibadtag            =  counter for bad time tags

c     If timeset is false, set the clock from the most recently read value
      if (.not. timeset) then
         iwkntag = igpswk
         sowtag  = gpssec

c     If timeset is true, check to see if the most recently read
c     value is within the requested span; if so set it.
      else if (
     .      secdif(igpswk,gpssec,iwkntag,sowtag).gt.0.d0
     .      .and.
     .      secdif(igpswk,gpssec,iwkntag,sowtag).lt.8.64d4)
     .     then
           iwkntag = igpswk
           sowtag  = gpssec

c     If timeset is true, but the most recently read values are
c     screwed up, set the flag to false
c     and try to set the week number correctly

      else
         timeset = .false.
         ibadtag = ibadtag+1

c        Don't assume the week number has rolled over unless
c        the time tag is nonzero and adding one week still gives
c        a reasonable time tag, less than a day later.
c        If this case actually occurs, we will have to wait
c        for an ephemeris block to straighten us out.

         if ( dabs(gpssec) .gt. 1.0d-20 .and.
     .      secdif(igpswk+1,gpssec,iwkntag,sowtag).lt.8.64d4 ) then
            weekset = .false.
         else if ( dabs( secdif(igpswk,gpssec,iwkntag,sowtag) )
     .            .gt.(7.0d0*8.64d4) ) then
c           if we are off by more than a week, panic
           weekset = .false.
         endif

         if (mod(ibadtag,100) .eq. 0 .or. ibadtag .le. 5) then
           write(uinfor,564)ibadtag,epoch,igpswk,gpssec,iwkntag,sowtag
  564      format (1x,'Bad time tag number ',i5,' at epoch',1x,i4,1x,
     .     'wkn=',1x,i8,1x,'sec=',1x,g22.14,/,
     .     '        Last good one',5x,'         ',1x,4x,1x,
     .     'wkn=',1x,i8,1x,'sec=',1x,g22.14)
           write(message,565)ibadtag,epoch,igpswk,gpssec,iwkntag,sowtag
  565      format ('Bad time tag number ',i5,' at epoch',1x,i4,1x,
     .     'wkn=',1x,i8,1x,'sec=',1x,g22.14,' Last good one',1x,4x,1x,
     .     'wkn=',1x,i8,1x,'sec=',1x,g22.14)
           call report_stat('WARNING','MAKEX','tagchk',' ',message,0)
         endif

c        If the week number has been lost, reset it from the file name
         if (.not. weekset) then
            write (message,566) igpswk,iwknstart
            write (uinfor,566) igpswk,iwknstart
            iwkntag = iwknstart
 566        format (1x,'Week number wrong.  Resetting from ',i5,
     .      ' to ',i5,' but waiting for ephemeris block before',
     .     ' accepting more data.')
            call report_stat('WARNING','MAKEX','tagchk',' ',message,0)
         endif
      endif

      return
      end
