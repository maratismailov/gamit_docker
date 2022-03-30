      subroutine saveobsa(nbytes, endepoch, nobs, pass_debug, obsbytes)

      implicit none

      include 'trackRTObs.h'
      include 'trackRT.h'

*     F90 routine that interfaces to the c++ routine trackRTComm.cpp
*     to save one data record from the real-time stream.  

* MOD TAH 121231: Updated to allow for variable formats with BNC 2.6 and greater.

      integer*4 nbytes    ! Number of bytes in this records (should be
                          ! multiple of 363 (362+\nl at end of line) or
                          ! 138 with BNC 2.6.
      integer*4 endepoch  ! Set 0 if this is still at same time;
                          ! Set 1 when epoch seems to change.
      integer*4 nobs      ! Number of obs already seen at this time;
                          ! When non-zero value passed, we check times.

      integer*4 pass_debug(10)  !  Debug array passed from trackRT or test_trackRT

      integer*4 recl      ! Record length (363 for BNC 2.5, 138 BNC 2.6)

      integer*4 ver       ! BNC apparent version * 100

      character*(maxbytes) obsbytes

!     Local variables
      character*10 stat_ida
      character*1  sata      ! Satellite system (G,R, S: Geostationary, E Galileo)
      character*1  nl, null        ! String wiht \n  and null in it
      character*4  sysprn    ! Combined system and PRN number (G02 e.g.)

      integer*4 gpswa, 
     .          prna,
     .          sl1c, sl1w, sl2p, sl5x

      real*8  gpssa,               ! GPS seconds of week
     .        c1c, l1c, d1c, s1c,  ! Range, Phase, Doppler, SNR C1C 
     .        c1w, l1w, d1w, s1w,  ! Range, Phase, Doppler, SNR P1 
     .        c2p, l2p, d2p, s2p,  ! Range, Phase, Doppler, SNR P2
     .        c5x, l5x, d5x, s5x   ! Range, Phase, Doppler, SNR C5

*     BNC 2.6 version variables 
      integer*4 mdatype ! Maximumum if FT (1/2/5 + C/P/X) types allowed
      parameter ( mdatype = 10 ) 
      integer*4 slot    ! Glonass slot number 
      integer*4 b(mdatype), bf(2)    ! Band for data 1 L1 2 L2
      integer*4 l(mdatype), lf(2)    ! Slip count
      character*1 t(mdatype), tf(2)  ! Trackning mode (C/P) for GPS
      real*8 r(mdatype), p(mdatype), d(mdatype), s(mdatype) ! Data Range, Phase, Doppler, SNR 
      real*8 rf(2), pf(2), df(2), sf(2) ! Data Range, Phase, Doppler, SNR 
      

      integer*4 ierr, jerr, kerr(6)  ! IOSTAT Error
      integer*4 no     ! Counter 
      integer*4 sti, eni   ! Position in string for decoding current
                       ! Data record.
      integer*4 i,j, k, n
      integer*4 lenw6   ! Length of Word(6) in decoding records.
      integer*4 ind, jf  ! Freqency number 1/2 (or 3 for L5 in future)

* MOD TAH 121231: Added more flexaible decoding based on BNC 2.6 and 
*     greater
      integer*4 iel, indx_sp, indx_nl, obs_recl, indx_end 
* MOD TAH 130716: Add another index to work through start of line where
*     seconds part of time tag may be short
      integer*4 indx_sh    ! Start of GPSW and GPS Seconds of line
 
      character*364  buf_record    ! First line of buffer record with
                      ! site name removed whcih can be different lengths
      character*32 word(6)   ! Word pulled from big_record

      logical dataOK  ! Set true for GPS data.  BNC has bug sometimes
                      ! where PRN is slit so that G23 become 203.
                      ! MOD TAH 110515:
      logical eol     ! End of line (buf_record) indicator
      logical mdaterr ! Set true if we exceeded number of data types



      nl = char(10) 
      null = char(0)

      endepoch = 0
      ver = 0
      mdaterr = .false. 

* MOD TAH 121231: Pull site name from first record to see how long it 
*     is.  Find the \n character
      indx_nl = index(obsbytes,nl) 
      indx_sp = index(obsbytes,' ' )   ! First space to end of name
*     Check for BNC >2.9 format
      if( indx_sp.eq.2 .and. obsbytes(1:1).eq.'>' ) then
*        Time tag from >2.9 version
         ver = 2900
      endif 

* MOD TAH 170329: See if version specified.
      if( BNC_Vers.gt.0 ) then
          ver = BNC_Vers
      endif

      if ( pass_debug(10).gt.0 ) then
          write(*,105) nbytes, indx_nl, indx_sp, ver
 105      format('DEBUG: NBYTES ',i5,' \n Index ',I6,' sp index ',i6,
     .           ' Ver ',i4)
      endif 


      if( indx_sp.gt.2 .and. ver.eq.0 ) then     ! Pull off name
          stat_ida = obsbytes(1:indx_sp-1)
*         Compute record length
          obs_recl = indx_nl-indx_sp
*         For 2.5 version od trackRT we expect this to 352-348 depending
*         on whether 1-2 digits in the slip counts.
          if( obs_recl.ge.348 .and. obs_recl.le.352 ) then
              ver = 2500
* MOD Set conditions for 2.6 version 
          elseif( obs_recl.ge.138 .and. obs_recl.le.142 ) then  
              ver = 2600
*         else     ! Set default to be 2900
              ver = 2900
          end if
      elseif( ver.eq.0 ) then    ! Short record; report and return
          write(*,110) nbytes, indx_nl, indx_sp, obsbytes(1:80) 
 110      format('UNKNOWN BNC record: Length ',i6,' New Line at ',i6,
     .           ' Space at ',i2,'. First 80-characters are:',/,a)
          nbytes = 0
          return
      endif

      if( pass_debug(10).gt.0 )
     .write(*,120) nbytes, ver/1000., obs_recl
 120  format('NBYTES ',I6,' BNC Ver ',f5.2,' OBS_RECL ',I6)

!     Loop over this ascii block, extracting the data by reading the
!     lines.  Due to errors in transferrs not all lines are complete.
!     Work through the lines
!     MOD TAH 120710: Try a version where we assume lines are complete
!     This seems to be the case recently (Error will now be printed 
!     because version will not compute correctly.)
      sti = 1
      eni = indx_nl

* MOD TAH 170329: See if this is a time tag. 
* MOD TAH 170402: Record can be shorter at start of the week
      if( eni.le.22 .and. ver.eq.2900 ) then
          read(obsbytes(2:22),*,iostat=jerr) gpswa, gpssa
          if( jerr.ne.0 ) then
             write(*,160) ierr, obsbytes(1:22)
 160         format('IOSTAT Error ',i6,' reading GSPW GPSS from ',a)
          elseif ( pass_debug(10).gt.0 ) then
             write(*,165) gpswa, gpssa
 165         format('FOUND TimeTag GSPW ',I5,' GPSS ',F14.7)
          endif
          sti = indx_nl+1
          eni = index(obsbytes(sti:),nl)+sti-2
          if ( pass_debug(10).gt.0 ) then 
             write(*,170) sti,indx_nl, obsbytes(sti:eni)
 170         format('Update Start INDX ',2I4,' Rec: |',a)
          endif  

      endif

      no  = nobs
      endepoch = 0
!     print *,'ENI etc ',indx_sp, indx_nl, eni, nbytes, nobs
      indx_end = min( index(obsbytes,null), len_trim(obsbytes))
!      print *,'obsbytes String ', indx_end
!      print *,obsbytes(1:indx_end)

!     See if we see a \n\n at the end of obs_recl.  If present then
!     this is end of epoch
      if( obsbytes(nbytes-1:nbytes).eq. nl // nl ) then
          endepoch = 1
          if( pass_debug(10).gt.0 )
     .    write(*,180) nbytes, nobs 
 180      format('END OF RECORD found Length ',i6,' #Obs so far ',i3)
      end if

      do while ( eni.le.nbytes )

!         Decode differently based on version
!         Pull off next records line
          if ( ver.lt.2900 ) then
             buf_record = obsbytes(sti:eni)
          else
             buf_record = obsbytes(sti:eni+1)
          endif 
!         print *,'buf_record ',sti,eni, trim(buf_record)
          indx_sp = index(obsbytes,' ' )   ! First space to end of name
!         Get site name and strip out next data record
          if( indx_sp.gt.2 .and. ver.lt.2900 ) then
              stat_ida = buf_record(1:indx_sp-1) 
          elseif ( ver.lt.2900 ) then 
              write(*,210) ver/1000.0, nbytes, trim(buf_record)
 210          format('NO STAT_ID in BNC Ver ',F5.2,' Record with ',
     .                I6,' bytes. Record: ',/,a)
              exit    ! Exit out of loop 
          endif    
         

          if( ver .eq. 2500 ) then 
!            Read the line between start in stop
             read(buf_record(indx_sp:),*,iostat=ierr) 
     .       	 gpswa, gpssa, sysprn, 
     .       	 c1c, l1c, sl1c, d1c, s1c,
     .       	 c1w, l1w, sl1w, d1w, s1w,
     .       	 c2p, l2p, sl2p, d2p, s2p,
     .       	 c5x, l5x, sl5x, d5x, s5x
  220         format(A10,1x,I4,1x,F14.7,1x,a1,I2.2,
     .           4(F14.3,F16.3,I3,F15.3,F16.3,2x))
             sata = sysprn(1:1)
             read(sysprn(2:),*,iostat=ierr) prna  

             if( pass_debug(10).gt.0 )
     .       write(*,220) stat_ida(1:4), gpswa, gpssa, sata, prna, 
     .       	 c1c, l1c, sl1c, d1c, s1c,
     .       	 c1w, l1w, sl1w, d1w, s1w,
     .       	 c2p, l2p, sl2p, d2p, s2p,
     .       	 c5x, l5x, sl5x, d5x, s5x

!            Only process the record if there was no error reading the line
!            and the data type is GPS
!            Test to see if sata really is letter
             call check_num(sata, jerr)
             if( jerr.eq.0 ) then    ! SATA contain number not letter
                 read(sata,*,iostat=jerr) k
                 if( pass_debug(10).gt.0 )
     .           write(*,280) sata, k, prna, k*10+prna, jerr
 280             format('BAD SAT/PRN: Sat ',a,1x,i3,1x,' Low PRN ',i2,
     .                  ' Final PRN ',i2,' Err ',I5)
                 if( k.le.3 .and. jerr.eq.0 ) then
                    sata = 'G'
                    prna = k*10+prna
                 end if
             end if


             if( ierr.eq.0 .and. sata(1:1) .eq. 'G' ) then

!               Increment count of number of observations
                no = no + 1

! MOD TAH 130716: 

!               See if times have changed
                if( no.gt.1 ) then   ! Check times
!                   if( gpssa.gt. RT_GPSweeks(no-1) ) then
! MOD TAH 111229: Only change epoch if time diffferent by > 0.1 seconds
                    if( abs(gpssa-RT_GPSweeks(no-1)).gt.0.1d0 ) then
!                       Epoch has changed; Set endepoch and return.  Block will 
!                       be processed later
                        endepoch = 1
                        no = no - 1 
                        if( pass_debug(10).gt.0 )
     .                  write(*,310) no, RT_GPSweeks(no-1), gpssa
 310                    format('NEW Epoch # ',i5,' Times ',2F16.7)
!                       See if this was not first observation in this block in
!                       which case, we will copy the data down for use next
!                       time
                        if( sti.gt.1 ) then
                            if( pass_debug(10).gt.0 )
     .                      write(*,320) nbytes-sti+1, sti,nbytes
 320                        format('MID-Block epoch change: copy  1:',
     .                             I6.6,' from ',I6.6,':',I6.6) 
                            obsbytes(1:nbytes-sti+1) = 
     .                               obsbytes(sti:nbytes)
                            obsbytes(nbytes-sti+2:nbytes) = ' '   ! Clear rest of line
                            nbytes = nbytes-sti+1      ! Save number of bytes now
                        endif
                        RETURN
                    endif
                end if

                nobs = no 
                RT_flags(nobs) =        0        
                RT_StatID(nobs) =       stat_ida    ! Station ID
                RT_satSys(nobs) =       sata        ! Satellite System ('G' or 'R')
                RT_satNum(nobs) =       prna        ! Satellite Number (PRN for GPS NAVSTAR)
                RT_slot(nobs) =         0           ! Slot Number (for Glonass)
                RT_GPSWeek(nobs) =      gpswa       ! Week of GPS-Time
                RT_GPSweeks(nobs) =     gpssa       ! Second of Week (GPS-Time)
                RT_C1(nobs) =           c1c         ! CA-code pseudorange (meters)
                RT_C2(nobs) =           0           ! CA-code pseudorange (meters)
                RT_P1(nobs) =           c1w         ! P1-code pseudorange (meters)
                RT_P2(nobs) =           c2p         ! P2-code pseudorange (meters)
                RT_L1(nobs) =           l1c         ! L1 carrier phase (cycles)
                RT_L2 (nobs) =          l2p         ! L2 carrier phase (cycles)
                RT_slip_cnt_L1(nobs) =  sl1c        ! L1 cumulative loss of continuity indicator (negative value = undefined)
                RT_slip_cnt_L2(nobs) =  sl2p        ! L2 cumulative loss of continuity indicator (negative value = undefined)
                RT_lock_timei_L1(nobs) = 0          ! L1 last lock time indicator                (negative value = undefined)
                RT_lock_timei_L2(nobs) = 0          ! L2 last lock time indicator                (negative value = undefined)
                RT_S1(nobs) =           s1c         ! L1 signal-to noise ratio
                RT_S2(nobs) =           s2p         ! L2 signal-to noise ratio
                RT_SNR1(nobs) =         nint(s1c)   ! L1 signal-to noise ratio (mapped to integer)
                RT_SNR2(nobs) =         nint(s2p)   ! L2 signal-to noise ratio (mapped to integer)

!               Compute the MLD of this obs
                RT_MJD_obs(nobs) =  44244 + RT_GPSWeek(nobs)*7 + 
     .                              RT_GPSweeks(nobs)/86400 
!               Set the initial error flags (Maybe RT_flag could be used too?)
                RT_errflag(nobs) = 0 

*               Now make sure we have 'P1' and 'P2' ranges
                if( RT_P1(nobs).eq.0 ) then
                    RT_P1(nobs) = RT_C1(nobs)
                    call sbit(RT_flags(nobs),16,1)
                endif
                if( RT_P2(nobs).eq.0 ) then
                    RT_P2(nobs) = RT_C2(nobs)
                    call sbit(RT_flags(nobs),17,1)
                endif

     
                num_rtobs = no

             endif   ! Data read OK

!            Now goto next records
             sti = eni+1
             if( sti+1.ge.nbytes ) then
                eni = sti
             else
                eni = sti+index(obsbytes(sti:),nl)-1
                indx_end = min(index(obsbytes,null), len_trim(obsbytes))
!               print *,'ENI+STI ',eni, sti, indx_end, nbytes
!     .,                 (obsbytes(sti:indx_end))
             end if
             if ( eni.eq.sti-1 .or. eni.eq.sti ) eni = nbytes+1
             if( pass_debug(10).gt.0 )
     .       write(*,420) no, ierr, sti, eni, nbytes-eni, endepoch
 420         format('Obs ',i4,' IOSTAT ',i6,' STI, ENI ',2I8,
     .              ' Remaining ',I8,' EndEpoch ',i2)

***************************************************************
          ELSEIF ( ver.eq.2600 ) then     ! Version 2.6
!            Version 2.6 Format:
!            StationID | GPSWeek | GPSWeekSec | PRN, G=GPS, R=GLO | 
!            SlotNumber (if GLO) | Band/Frequency & trackingMode | 
!            Code | Phase | Doppler | SNR | SlipCount | ....
!P498_RTCM3 1695 392234.0000000 G32     1C   21778302.736  114445787.184 0.0   47.000 -1  2P   21778293.076   89178543.508 0.0   34.750 -1
!P497_RTCM3 1748 0.0000000 G17     1C   22990270.368  120814756.710 0.0   45.750 -1  1P 0.0 0.0 0.0 0.0 -1  2P   22990260.548   94141385.600 0.0   34.250 -1  5C 0.0 0.0 0.0 0.0 -1 
C     integer*4 slot    ! Glonass slot number 
C     integer*4 b(2)    ! Band for data 1 L1 2 L2
C     integer*4 l(2)    ! Slip count
C     character*1 t(2)  ! Trackning mode (C/P) for GPS
C     real*8 r(2), p(2), d(2), s(2) ! Data Range, Phase, Doppler, SNR 

!            With version 2.6 and greater we need to decode by blocks.
!            First get the time and PRN values
!            read(buf_record(indx_sp+1:),620,iostat=ierr) 
!    .            gpswa, gpssa, sata, prna, slot
!620         format(I4,1x,F14.7,1x,a1,I2.2,1x,I2) 
!
! MOD TAH 130716: Problem is this record is not fixed format so can't be read
!     as above.  Pull the line to peices element by element
             indx_sh = indx_sp+1
             call GetWord(buf_record,word(1), indx_sh)   ! GPS Week
             read(word(1),'(I4)',iostat=ierr) gpswa
             if( ierr.ne.0 ) then 
                 write(*,622) 'GPS Week', indx_sh,indx_sp,
     .                 trim(word(1)), trim(buf_record(indx_sp+1:))
 622             format('ERROR Decoding ',a,' Indices ',2i5,
     .                  ' Word ',a,' Buffer ',a)
             endif

             call GetWord(buf_record,word(2), indx_sh)   ! GPS seconds
             read(word(2),'(F14.7)',iostat=ierr) gpssa
             if( ierr.ne.0 ) then 
                 write(*,622) 'GPS second', indx_sh,indx_sp, 
     .                 trim(word(2)), trim(buf_record(indx_sp+1:))
             endif
!            Rest of line should be readable under format
             read(buf_record(indx_sh+1:),625,iostat=ierr) 
     .             sata, prna, slot
 625         format(a1,I2.2,1x,I2) 

****         Now loop over the reset of the record to see what we find.
             eol = .false.
             iel = indx_sh+ 6    ! Mod'd to account for gps second variable format
             n = 0
             do while ( .not.eol )
!                 Pull of next 6 values and see what they are
                  n = n + 1
                  if( n.gt.mdatype ) mdaterr = .true.
                  if( n.ge.mdatype) then
                       eol = .true.
                       write(*,630) n, trim(buf_record)
 630                   format('Too Many Data types: Found max ',i3,
     .                        ' types. Buffer: ',/,a)
                       write(*,635) trim(word(5)),
     .                         word(6)(lenw6-1:lenw6+1),
     .                         (ichar(word(6)(lenw6+k:lenw6+k)),
     .                            k = -1,1)
 635                   format('Words W5 ',a,'| W6 |',a,' ICHAR ',3i6)
                  end if
                  do i = 1, 6
                      call GetWord(buf_record, word(i), iel)
                  end do
!                 Read word(1) -- Freq + letter type
                  read(word(1),'(i1,a1)',iostat=kerr(1)) b(n),t(n)
                  read(word(2),*,iostat=kerr(2)) r(n)
                  read(word(3),*,iostat=kerr(3)) p(n)
                  read(word(4),*,iostat=kerr(4)) d(n)
                  read(word(5),*,iostat=kerr(5)) s(n)
                  read(word(6),*,iostat=kerr(6)) l(n)
!                 See if next character is nl.  If so this is end of line
                  lenw6 = len_trim(word(6))
                  if( word(6)(lenw6:lenw6).eq.nl ) then
!                     end of this buffer, so save and continue
                      eol = .true.
                  end if
             end do

             if ( n.lt.mdatype  ) mdaterr = .false.

             if( pass_debug(10).gt.0 .or. mdaterr )
     .       write(*,640) stat_ida, gpswa, gpssa, sata, prna, slot,
     .           (b(i), t(i), r(i), p(i), d(i), s(i), l(i),i=1,n)
 
 640         format('Record ',A10,1x,I4,1x,F14.7,1x,a1,I2.2,1x,I2,
     .              20(2x,i1,a1,F15.3,F15.3,F4.1,F9.3,1x,i2))

!            Save values to L1 and L2 frequencies (skip L5 for the moment).
             ierr = 0
             do i = 1, n
                jf = b(i)
                if( jf.eq.1 .or. jf.eq.2 ) then
                   ierr = ierr + abs(kerr(i))
                   if( r(i).ne.0.0 ) then
                      rf(jf) = r(i)
                      tf(jf) = t(i)
                   endif
                   if( p(i).ne.0.0 ) pf(jf) = p(i)
                   if( s(i).ne.0.0 ) sf(jf) = s(i)
                   if( l(i).ne.0.0 ) lf(jf) = s(i)
                endif
             end do

             if( ierr.eq.0 .and. sata .eq. 'G' ) then

!               Increment count of number of observations
                no = no + 1

!               See if times have changed
                if( no.gt.1 ) then   ! Check times
!                   if( gpssa.gt. RT_GPSweeks(no-1) ) then
! MOD TAH 111229: Only change epoch if time diffferent by > 0.1 seconds
                    if( abs(gpssa-RT_GPSweeks(no-1)).gt.0.1d0 ) then
!                       Epoch has changed; Set endepoch and return.  Block will 
!                       be processed later
                        endepoch = 1
                        no = no - 1 
                        if( pass_debug(10).gt.0 )
     .                  write(*,650) no, RT_GPSweeks(no), gpssa
 650                    format('NEW Epoch # ',i5,' in block GPS Sec ',
     .                         2F16.7)
!                       See if this was not first observation in this block in
!                       which case, we will copy the data down for use next
!                       time
                        if( sti.gt.1 ) then
                            if( pass_debug(10).gt.0 )
     .                      write(*,660) nbytes-sti+1, sti,nbytes
 660                        format('MID-Block epoch change: copy  1:',
     .                             I6.6,' from ',I6.6,':',I6.6) 
                            obsbytes(1:nbytes-sti+1) = 
     .                               obsbytes(sti:nbytes)
                            obsbytes(nbytes-sti+2:nbytes) = ' '   ! Clear rest of line
                            nbytes = nbytes-sti+1      ! Save number of bytes now
                        endif
                        RETURN
                    endif
                end if

****            Get index for L1 and L2
                nobs = no 
                RT_flags(nobs) =        0        
                RT_StatID(nobs) =       stat_ida    ! Station ID
                RT_satSys(nobs) =       sata        ! Satellite System ('G' or 'R')
                RT_satNum(nobs) =       prna        ! Satellite Number (PRN for GPS NAVSTAR)
                RT_slot(nobs) =         0           ! Slot Number (for Glonass)
                RT_GPSWeek(nobs) =      gpswa       ! Week of GPS-Time
                RT_GPSweeks(nobs) =     gpssa       ! Second of Week (GPS-Time)
                RT_C1(nobs) =           rf(1)       ! CA-code pseudorange (meters)
                RT_C2(nobs) =           rf(2)       ! CA-code pseudorange (meters)
                RT_P1(nobs) =           rf(1)       ! P1-code pseudorange (meters)
                RT_P2(nobs) =           rf(2)       ! P2-code pseudorange (meters)
                RT_L1(nobs) =           pf(1)       ! L1 carrier phase (cycles)
                RT_L2 (nobs) =          pf(2)       ! L2 carrier phase (cycles)
                RT_slip_cnt_L1(nobs) =  lf(1)       ! L1 cumulative loss of continuity indicator (negative value = undefined)
                RT_slip_cnt_L2(nobs) =  lf(2)       ! L2 cumulative loss of continuity indicator (negative value = undefined)
                RT_lock_timei_L1(nobs) = 0          ! L1 last lock time indicator                (negative value = undefined)
                RT_lock_timei_L2(nobs) = 0          ! L2 last lock time indicator                (negative value = undefined)
                RT_S1(nobs) =           sf(1)       ! L1 signal-to noise ratio
                RT_S2(nobs) =           sf(2)       ! L2 signal-to noise ratio
                RT_SNR1(nobs) =         nint(sf(1)) ! L1 signal-to noise ratio (mapped to integer)
                RT_SNR2(nobs) =         nint(sf(2)) ! L2 signal-to noise ratio (mapped to integer)

!               Compute the MLD of this obs
                RT_MJD_obs(nobs) =  44244 + RT_GPSWeek(nobs)*7 + 
     .                              RT_GPSweeks(nobs)/86400 
!               Set the initial error flags (Maybe RT_flag could be used too?)
                RT_errflag(nobs) = 0 

*               Now make sure we have 'P1' and 'P2' ranges
                if( tf(1).eq.'C' ) then
                    RT_P1(nobs) = RT_C1(nobs)
                    call sbit(RT_flags(nobs),16,1)
                endif
                if( tf(2).eq.'C' ) then
                    RT_P2(nobs) = RT_C2(nobs)
                    call sbit(RT_flags(nobs),17,1)
                endif
     
                num_rtobs = no
             endif

***************************************************************
          ELSE      ! Version 2.9
!            Version 2.9 Format:
!            StationID | PRN, G=GPS, R=GLO | <Obs> <Value> | <Obs> <Value> | ...
!            
!P513_RTCM3 G01 C1C   21226688.820 L1C  111546993.239  127 S1C   50.000 C2W   21226690.320 L2W   86919734.635  127 S2W   48.000
C     integer*4 slot    ! Glonass slot number 
C     integer*4 b(2)    ! Band for data 1 L1 2 L2
C     integer*4 l(2)    ! Slip count
C     character*1 t(2)  ! Trackning mode (C/P) for GPS
C     real*8 r(2), p(2), d(2), s(2) ! Data Range, Phase, Doppler, SNR 
!
! MOD TAH 130716: Problem is this record is not fixed format so can't be read
!     as above.  Pull the line to peices element by element
!            Get the station ID
             indx_sp = index(buf_record,' ');
             stat_ida = buf_record(1:indx_sp-1)

             indx_sh = indx_sp+1
 
!            Get the Type and PRN 
             read(buf_record(indx_sh:),725,iostat=ierr) 
     .             sata, prna
 725         format(a1,I2.2) 
!            write(*,726) gpswa, gpssa, sata, prna, indx_sh+4
!726         format('DATA ',I4,1x,F14.7,1x,a,1x,i2,1x,' Index ',i4)


****         Now loop over the reset of the record to see what we find.
             eol = .false.
             iel = indx_sh+4    ! Add count for Gnn (PRN number)
             n = 0
             do while ( .not.eol )
                  n = n + 1
!                 Pull off Data type
                  call GetWord(buf_record,word(1),iel)
!                 Read the value
                  call GetWord(buf_record,word(2),iel) 
              
                  if( n.gt. mdatype ) mdaterr = .true.
                  if( n.ge. mdatype) then
                       eol = .true.
                       write(*,730) n, trim(buf_record)
 730                   format('Too Many Data types: Found max ',i3,
     .                        ' types. Buffer: ',/,a)
                  end if

!                 Read word(1) -- Freq + letter type
                  read(word(1),'(a1,I1)',iostat=kerr(1)) t(n),b(n)
                  read(word(2),*,iostat=kerr(2)) r(n)
!                 If t(n) is 'L' then phase measurement which has a flag
                  if( t(n)(1:1).eq.'L' ) then
                     call GetWord(buf_record,word(3),iel) 
                     read(word(3),*,iostat=kerr(3)) l(n)
                  endif

*                 Start assigning values to be saved in the data arrays
                  if( t(n)(1:1).eq.'L' ) then
                      if( b(n).eq.1 ) l1c = r(n)
                      if( b(n).eq.2 ) l2p = r(n)
                      if( b(n).eq.1 ) sl1c = l(n)
                      if( b(n).eq.2 ) sl2p = l(n)
                  elseif ( t(n)(1:1).eq.'C' ) then
                      if( b(n).eq.1 ) c1c = r(n)
                      if( b(n).eq.2 ) c2p = r(n)
                  elseif ( t(n)(1:1).eq.'P' ) then
                      if( b(n).eq.1 ) c1w = r(n)
                      if( b(n).eq.2 ) c2p = r(n)
                  elseif ( t(n)(1:1).eq.'S' ) then
                      if( b(n).eq.1 ) s1c = r(n)
                      if( b(n).eq.2 ) s2p = r(n)
                  endif 
 
!                 See if next character is nl.  If so this is end of line
                  if( buf_record(iel-1:iel-1).eq.nl ) then
!                     end of this buffer, so save and continue
                      eol = .true.
                  end if
             end do

             if ( n.lt.mdatype ) mdaterr = .false.

             if( ierr.eq.0 .and. sata .eq. 'G' ) then

!               Increment count of number of observations
                no = no + 1

!               See if times have changed
                if( no.gt.1 ) then   ! Check times
!                   if( gpssa.gt. RT_GPSweeks(no-1) ) then
! MOD TAH 111229: Only change epoch if time diffferent by > 0.1 seconds
                    if( abs(gpssa-RT_GPSweeks(no-1)).gt.0.1d0 ) then
!                       Epoch has changed; Set endepoch and return.  Block will 
!                       be processed later
                        endepoch = 1
                        no = no - 1 
                        if( pass_debug(10).gt.0 )
     .                  write(*,750) no, RT_GPSweeks(no), gpssa
 750                    format('NEW Epoch # ',i5,' in block GPS Sec ',
     .                         2F16.7)
!                       See if this was not first observation in this block in
!                       which case, we will copy the data down for use next
!                       time
                        if( sti.gt.1 ) then
                            if( pass_debug(10).gt.0 )
     .                      write(*,760) nbytes-sti+1, sti,nbytes
 760                        format('MID-Block epoch change: copy  1:',
     .                             I6.6,' from ',I6.6,':',I6.6) 
                            obsbytes(1:nbytes-sti+1) = 
     .                               obsbytes(sti:nbytes)
                            obsbytes(nbytes-sti+2:nbytes) = ' '   ! Clear rest of line
                            nbytes = nbytes-sti+1      ! Save number of bytes now
                        endif
                        RETURN
                    endif
                end if

****            Get index for L1 and L2
                nobs = no 
                RT_flags(nobs) =        0        
                RT_StatID(nobs) =       stat_ida    ! Station ID
                RT_satSys(nobs) =       sata        ! Satellite System ('G' or 'R')
                RT_satNum(nobs) =       prna        ! Satellite Number (PRN for GPS NAVSTAR)
                RT_slot(nobs) =         0           ! Slot Number (for Glonass)
                RT_GPSWeek(nobs) =      gpswa       ! Week of GPS-Time
                RT_GPSweeks(nobs) =     gpssa       ! Second of Week (GPS-Time)
                RT_C1(nobs) =           c1c         ! CA-code pseudorange (meters)
                RT_C2(nobs) =           0           ! CA-code pseudorange (meters)
                RT_P1(nobs) =           c1w         ! P1-code pseudorange (meters)
                RT_P2(nobs) =           c2p         ! P2-code pseudorange (meters)
                RT_L1(nobs) =           l1c         ! L1 carrier phase (cycles)
                RT_L2 (nobs) =          l2p         ! L2 carrier phase (cycles)
                RT_slip_cnt_L1(nobs) =  sl1c        ! L1 cumulative loss of continuity indicator (negative value = undefined)
                RT_slip_cnt_L2(nobs) =  sl2p        ! L2 cumulative loss of continuity indicator (negative value = undefined)
                RT_lock_timei_L1(nobs) = 0          ! L1 last lock time indicator                (negative value = undefined)
                RT_lock_timei_L2(nobs) = 0          ! L2 last lock time indicator                (negative value = undefined)
                RT_S1(nobs) =           s1c         ! L1 signal-to noise ratio
                RT_S2(nobs) =           s2p         ! L2 signal-to noise ratio
                RT_SNR1(nobs) =         nint(s1c)   ! L1 signal-to noise ratio (mapped to integer)
                RT_SNR2(nobs) =         nint(s2p)   ! L2 signal-to noise ratio (mapped to integer)

!               Compute the MLD of this obs
                RT_MJD_obs(nobs) =  44244 + RT_GPSWeek(nobs)*7 + 
     .                              RT_GPSweeks(nobs)/86400 
!               Set the initial error flags (Maybe RT_flag could be used too?)
                RT_errflag(nobs) = 0 

*               Now make sure we have 'P1' and 'P2' ranges
                if( RT_P1(nobs).eq.0  ) then
                    RT_P1(nobs) = RT_C1(nobs)
                    call sbit(RT_flags(nobs),16,1)
                endif
                if( RT_P2(nobs).eq.0  ) then
                    RT_P2(nobs) = RT_C2(nobs)
                    call sbit(RT_flags(nobs),17,1)
                endif

     
                num_rtobs = no

             endif   ! Data read OK

!            Now goto next records
             sti = eni+1
             if( sti+1.ge.nbytes ) then
                eni = sti
             else
                eni = sti+index(obsbytes(sti:),nl)-1
                indx_end = min(index(obsbytes,null), len_trim(obsbytes))
!               print *,'ENI+STI ',eni, sti, indx_end, nbytes
!     .                 (obsbytes(sti:indx_end))
             end if
             if ( eni.eq.sti-1 .or. eni.eq.sti ) eni = nbytes+1
!            if( pass_debug(10).gt.0 )
!    .       write(*,780) no, ierr, sti, eni, nbytes-eni, endepoch
!780         format('Obs ',i4,' IOSTAT ',i6,' STI, ENI ',2I8,
!    .              ' Remaining ',I8,' EndEpoch ',i2)
         endif

      end do

      nobs = no
      if( pass_debug(10).gt.0 ) call flush

      return

      end

      subroutine reportobs( nobs )

      implicit none

      include 'trackRTObs.h'
      include 'trackRT.h'

      integer*4 nobs      ! Number of obs already seen at this time;
                          ! When non-zero value passed, we check times.
      integer*4 i

****  Routine to print data that will be used in trackRT
      write(*,110) nobs, num_rtobs
 110  format('SAVED DATA: Number of obs passed ',i4,' Saved ',I4)
      do i = 1,nobs
           write(*,140) i, RT_StatID(i), RT_GPSWeek(i), 
     .                     RT_GPSweeks(i),RT_satSys(i), 
     .                     RT_satNum(i), RT_slot(i),
     .                'L1',RT_P1(i), RT_L1(i), RT_S1(i), 
     .                'L2',RT_P2(i), RT_L2(i), RT_S2(i)
 140         format('RT_OBS ',i5,1x,A10,1x,I4,1x,F14.7,1x,a1,I2.2,1x,I2,
     .            2x, 20(a,1x,F15.3,F15.3,1x,F5.2,1x))
      end do
      write(*,160) nobs, num_rtobs
 160  format('AT END DATA: Number of obs passed ',i4,' Saved ',I4)

      return
      end 
      


