      subroutine saveobsa( nbytes, endepoch, nobs, obsbytes )

      implicit none

      include 'trackRTObs.h'
      include 'trackRT.h'

*     F90 routine that interfaces to the c++ routine trackRTComm.cpp
*     to save one data record from the real-time stream.  

      integer*4 nbytes    ! Number of bytes in this records (should be
                          ! multiple of 363 (362+\nl at end of line).
      integer*4 endepoch  ! Set 0 if this is still at same time;
                          ! Set 1 when epoch seems to change.
      integer*4 nobs      ! Number of obs already seen at this time;
                          ! When non-zero value passed, we check times.

      character*(maxbytes) obsbytes

!     Local variables
      character*10 stat_ida
      character*1  sata      ! Satellite system (G,R, S: Geostationary, E Galileo)
      character*1  nl        ! String wiht \n in it

      integer*4 gpswa, 
     .          prna,
     .          sl1c, sl1w, sl2p, sl5x

      real*8  gpssa,               ! GPS seconds of week
     .        c1c, l1c, d1c, s1c,  ! Range, Phase, Doppler, SNR C1C 
     .        c1w, l1w, d1w, s1w,  ! Range, Phase, Doppler, SNR P1 
     .        c2p, l2p, d2p, s2p,  ! Range, Phase, Doppler, SNR P2
     .        c5x, l5x, d5x, s5x   ! Range, Phase, Doppler, SNR C5

      integer*4 ierr, jerr   ! IOSTAT Error
      integer*4 no     ! Counter 
      integer*4 sti, eni   ! Position in string for decoding current
                       ! Data record.
      integer*4 k
      integer*4 ind

      logical dataOK  ! Set true for GPS data.  BNC has bug sometimes
                      ! where PRN is slit so that G23 become 203.
                      ! MOD TAH 110515:

      nl = char(10) 
      endepoch = 0

      if( debug(10).gt.0 )
     .write(*,120) nbytes, int(nbytes/363), nbytes-int(nbytes/363)*363
 120  format('NBYTES ',I6,' Blocks',I5,' Remainder ',i5)

!     Loop over this ascii block, extracting the data by reading the
!     lines.  Due to errors in transferrs not all lines are complete.
!     Work through the lines
      sti = 1
      eni = index(obsbytes,nl)

      no  = nobs
      endepoch = 0
      do while ( eni.le.nbytes )
!         Read the line between start in stop
          read(obsbytes(sti:eni),220,iostat=ierr) 
     .        stat_ida, gpswa, gpssa, sata, prna, 
     .        c1c, l1c, sl1c, d1c, s1c,
     .        c1w, l1w, sl1w, d1w, s1w,
     .        c2p, l2p, sl2p, d2p, s2p,
     .        c5x, l5x, sl5x, d5x, s5x
 220      format(A10,1x,I4,1x,F14.7,1x,a1,I2.2,
     .           4(F14.3,F16.3,I3,F15.3,F16.3,2x))

          if( debug(10).gt.0 )
     .    write(*,220) stat_ida(1:4), gpswa, gpssa, sata, prna, 
     .        c1c, l1c, sl1c, d1c, s1c,
     .        c1w, l1w, sl1w, d1w, s1w,
     .        c2p, l2p, sl2p, d2p, s2p,
     .        c5x, l5x, sl5x, d5x, s5x

!        Only process the record if there was no error reading the line
!        and the data type is GPS
!        Test to see if sata really is letter
         call check_num(sata, jerr)
         if( jerr.eq.0 ) then    ! SATA contain number not letter
             read(sata,*,iostat=jerr) k
             if( debug(10).gt.0 )
     .       write(*,280) sata, k, prna, k*10+prna, jerr
 280         format('BAD SAT/PRN: Sat ',a,1x,i3,1x,' Low PRN ',i2,
     .              ' Final PRN ',i2,' Err ',I5)
             if( k.le.3 .and. jerr.eq.0 ) then
                sata = 'G'
                prna = k*10+prna
             end if
         end if


         if( ierr.eq.0 .and. sata .eq. 'G' ) then

!           Increment count of number of observations
            no = no + 1

!           See if times have changed
            if( no.gt.1 ) then   ! Check times
!               if( gpssa.gt. RT_GPSweeks(no-1) ) then
! MOD TAH 111229: Only change epoch if time diffferent by > 0.1 seconds
                if( abs(gpssa-RT_GPSweeks(no-1)).gt.0.1d0 ) then
!                   Epoch has changed; Set endepoch and return.  Block will 
!                   be processed later
                    endepoch = 1
                    no = no - 1 
                    if( debug(10).gt.0 )
     .              write(*,310) no, RT_GPSweeks(no-1), gpssa
 310                format('NEW Epoch # ',i5,' Times ',2F16.7)
!                   See if this was not first observation in this block in
!                   which case, we will copy the data down for use next
!                   time
                    if( sti.gt.1 ) then
                        if( debug(10).gt.0 )
     .                  write(*,320) nbytes-sti+1, sti,nbytes
 320                    format('MID-Block epoch change: copy  1:',I6.6,
     .                         ' from ',I6.6,':',I6.6) 
                        obsbytes(1:nbytes-sti+1) = obsbytes(sti:nbytes)
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

!           Compute the MLD of this obs
            RT_MJD_obs(nobs) =  44244 + RT_GPSWeek(nobs)*7 + 
     .                          RT_GPSweeks(nobs)/86400 

!           Set the initial error flags (Maybe RT_flag could be used too?)
            RT_errflag(nobs) = 0 

*           Now make sure we have 'P1' and 'P2' ranges
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

!        Now goto next records
         sti = eni+1
         eni = sti+index(obsbytes(sti:),nl)-1
         if ( eni.eq.sti-1 .or. eni.eq.sti ) eni = nbytes+1
         if( debug(10).gt.0 )
     .   write(*,420) no, ierr, sti, eni, nbytes-eni
 420     format('Obs ',i4,' IOSTAT ',i6,' STI, ENI ',2I6,
     .          ' Remaining ',I5)

      end do

      return

      end

       

      


