      logical function isgood( ndat,dattyp,phase1,issi1,phase2,issi2
     .                       , pr1,pr2,ibad )

c     Check for minimum signal amplitude and reasonableness
c     of phase and pseudorange observations.  Returns false if
c     either test fails.  
c
c       Written by Kurt Feigl; last modified by R. King 9 July 1990

      include '../includes/dimpar.h'
      include '../includes/makex.h'

c       data types present
      integer*4 ndat,dattyp(maxdat) 

c       phase and pseudorange
      real*8 phase1,phase2,pr1,pr2

c       RINEX signal strength indicators
      integer*4 issi1,issi2

      logical chkdat,ibad,chkrnx,chkpr

      data chkrnx/.false./,chkpr/.false./

      isgood = .true.
      ibad   = .false.

c     We want ISSI = 1 to be marginal, not missing. kurt 901020
c       check for RINEX signal strength greater than worst possible    
c     RWK 070907: Not sure why this was commented out, but make
c     live, but unused code to avoid compiler warning on unused variable
      if( chkrnx) then
        if( chkdat(1,ndat,dattyp) .and. issi1.eq.1 .or.
     .    chkdat(2,ndat,dattyp) .and. issi2.eq.1 ) then
           isgood = .false.
           ibad   = .false.
           return
         endif 
      endif 
cd      print *,'ISGOOD 1 ',isgood



c       check for reasonable phase and pseudorange (non-zero)

      if((chkdat(1,ndat,dattyp) .and. phase1.eq.0.d0) .or.
c yb 12/11/90
     .   (chkdat(2,ndat,dattyp) .and. phase2.eq.0.d0) .or.
     .   (chkdat(3,ndat,dattyp) .and. pr1.eq.0.d0) .or.
c yb 12/29/92 (Do not discard data if P2 is bad)
c     .   (chkdat(4,ndat,dattyp) .and. pr2.eq.0.d0) .or.
     .   (chkdat(5,ndat,dattyp) .and. pr1.eq.0.d0))
     .  then
          isgood = .false.
          ibad   = .true.
      endif                     
cd      print *,'ISGOOD 2 ',isgood
c** rwk/tah 070710: Restore the check of reasonableness, but at a crude
c    level (5 min = 9e10 km) to trap problem with PRs so far off that we
c    go beyond the t-file interval (15 minutes)

      if(                                  
     . ( chkdat(3,ndat,dattyp) .and. dabs(pr1).gt.1.d10 ) .or.
     . ( chkdat(4,ndat,dattyp) .and. dabs(pr2).gt.1.d10 ) .or.
     . ( chkdat(5,ndat,dattyp) .and. dabs(pr1).gt.1.d10 ) ) then
          isgood = .false.
          ibad   = .true.
      endif                     
cd      print *,'ISGOOD 3 ',isgood


c       check for identical L1 and L2 pseudorange, found occasionally
c       with TurboRogues.  This indicates that L2 phase is also identical
c       to L1 except for a bias.

      if (
     .   (chkdat(3,ndat,dattyp) .and. chkdat(4,ndat,dattyp)) .and.
     .    dabs(pr1-pr2).lt.5.d-3 ) then
          isgood = .false.
          ibad   = .true.
      endif                     
cd      print *,'ISGOOD 4 ',isgood

      return
      end

      logical function chkdat ( itype,ndat,dattyp )
c     check to see if data type itype is present in the dattyp array

      include '../includes/dimpar.h'
      integer*4 i,itype,ndat,dattyp(maxdat)

      chkdat = .false.
      do i=1,ndat
      if( dattyp(i).eq.itype ) chkdat=.true.
      enddo
      return
      end
