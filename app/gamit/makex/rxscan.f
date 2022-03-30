      program RXSCAN

c     scan a RINEX file for useful information

c     calling arguments: RINEX-file-name  gnss 

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

c     passed values

c     passed values
C     gps sec of week for epoch
      real*8      sow
C     L1, L2 doppler phase in cycles
      real*8      dofl1(maxchn), dofl2(maxchn)
C     L1, L2  pseudorange
      real*8      prgl1(maxchn),prgl2(maxchn)
C     snr and loss-of-lock for phase and pseudorange observables
      integer*4    issi(maxchn,4),illi(maxchn,4)


C     .true. at end of file
      logical     fend
C     .true. on error
      logical     ferr
C     0 for time scan only; 1 for selecting observations
      integer*4   iflag
C     for broadcast ephem
      integer*4   nprn
C     sat id num PRN 
      character*1 asvid(maxchn)
      integer*4   isvid(maxchn)
C     gps week number for epoch
      integer*4  iwkn
c     flag to debug
      logical debug/.false./

c     I/O Error code
      integer    ioerr

c     strings for buffering data
      character*10 buff10
      character*80 aline

c     RINEX defined items
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      integer irxcom
      character*60 rxcom(maxlin)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiever serial number, type and SW version
      character*20 rcvnum
      character*20 rctype
      character*20 rcvers
c     antenna serial number and type
      character*20 antnum
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
c     observation types
      integer nobtyp
      character*3 rxobtyp(maxobt)
c     data interval in seconds
      real*8 rxint
c     time type (always 'GPS' for mixed files)
      character*3 rxtime
c      data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1
c     seconds of day
      real*8 sod
  
c     indices in rxobtyp for observables to be selected
      integer*4 iobtypx(6)

c     array of ALL PRN numbers, not just at one epoch
      integer kprn(maxsat)

c     first and last time tags:
      integer*4 iwkn0,iwkn1
      real*8    sow0,sow1,secdif

c     for gap analysis (L1, L2 and Combined)
      integer k1old(maxsat),k2old(maxsat),kcold(maxsat)
      integer k1new(maxsat),k2new(maxsat),kcnew(maxsat)
      integer ngood1(maxsat),ngood2(maxsat),ngoodc(maxsat)
      logical l1gap(maxsat),l2gap(maxsat),lcgap(maxsat)
      integer n1gap(maxsat),n2gap(maxsat),ncgap(maxsat)
      logical lfound

c     pointer in PRN (sat id) array
      integer iprn
         
c     GNSS id
      character*1 gnss

c     miscellaneous variables
      integer i,j

c     number of epochs
      integer*4 nepoch,iepoch

c     total number of sats
      integer*4 nsat

c     bargraph count
      character*96  bargph(maxsat)
      integer       kount(maxsat)
      character*1   letter
      integer       jj(24),kbin,nbin,nwide

c     function to return day of year
      integer*4 idoy

c     function to return number of non-blank characters in a string
      integer*4 nblen

c     report_stat message variable
      character*256 message

      call report_stat('STATUS','RXSCAN','rxscan',' ',
     .                 'Started RXSCAN ',0)

c      use one bin/15 minutes for 24 hours = 96 bins
       nbin = 96

      do i = 1,maxsat
         write (bargph(i),13)
  13     format (96('.'))
      enddo


      iprn = 0
      nsat = 0
      nwide = 0
      do 5 i = 1,maxsat
         kprn(i)  = 0
         k1old(i) = 0
         k2old(i) = 0
         kcold(i) = 0
         k1new(i) = 0
         k2new(i) = 0
         kcnew(i) = 0
         l1gap(i) = .false.
         l2gap(i) = .false.
         lcgap(i) = .false.
         ngood1(i) = 0
         ngood2(i) = 0
         ngoodc(i) = 0
         kount(i) = 0
 5    continue

      uscren = 6                               
      gnss = ' '                    
     
      write (uscren,*) 'RINEX file name?'
      read  (*,'(2a)',iostat=ioerr) frinex,gnss
      if( gnss.eq.' ' ) gnss = 'G'

      urinex = 10
      open (unit = urinex,
     .      file = frinex,
     .      iostat = ioerr,
     .      status = 'old')

      if (ioerr .eq. 0) then
        call report_stat('STATUS','RXSCAN','rxscan',frinex,
     .  'Opened file: ',0)
      else
        call report_stat('FATAL','RXSCAN','rxscan',frinex,
     .  'Error opening file: ',0)
      endif
                    
      call rrxhed ( debug,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1)
                                  
c     select the phase and pseudorange observables to be checked
      call sel_obtyp(gnss,nobtyp,rxobtyp,iobtypx)  

c     check file name
      j = index(frinex,'.')
      write (buff10(1:2),'(i2.2)') mod(irxyr0,100)
      if (buff10(1:2) .ne. frinex(j+1:j+2)) then
         write (message,*) frinex,'Start year (',buff10(1:2),
     .   ') does not match file name.'
         call report_stat('WARNING','RXSCAN','rxscan',frinex,message,0)
      endif
      write (buff10(1:3),'(i3.3)') idoy(irxyr0,irxmo0,irxdy0)
      if (buff10(1:3) .ne. frinex(j-4:j-2)) then
         write (message,*) frinex,'Start day of year (',
     .   buff10(1:3),') does not match file name.'
         call report_stat('WARNING','RXSCAN','rxscan',frinex,message,0)
      endif                                     

c     read until end of file
      nepoch = 0
      debug = .false.                                 
c     iobtypx not defined for rxscan, set to zero
      do i=1,6
       iobtypx(i) = 0
      enddo
      do while(.not.fend)                              
        call rrinex( debug,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .             , nprn,isvid,rxtime,iwkn,sow,nepoch
     .             , dofl1,dofl2,prgl1,prgl2,illi,issi
     .             , anth,ante,antn,fend,ferr )
        if( .not.fend ) then 
          if (nepoch .eq. 1) then
            sow0 = sow
            iwkn0 = iwkn
            sow1 = sow
            iwkn1 = iwkn
            iepoch = 1
            nwide = 100
          else
            sow1 = sow
            iwkn1 = iwkn
            if (nepoch .eq. 2 .and. rxint .eq. 0.d0 ) then
              rxint = secdif (iwkn1,sow1,iwkn0,sow0)
            endif
          endif

c         figure out the width of the bins
          if (nepoch .eq. 2) then
            if (rxint .gt. 0.d0 ) then
               nwide = int(900./rxint)
            else
              nwide = 60
             endif
          endif   

c         seconds since start of GPS day
          sod = dmod(sow,86400.0d0)
c         aproximate epoch number
          iepoch = int(sod/rxint) + 1   

c         figure out which sat we have
          do j =1,nprn
            lfound = .false. 
            do  i=1,nsat
              if (isvid(j) .eq. kprn(i)) then
                iprn = i
                lfound = .true.
              endif
            enddo 

c           This is the first time we have seen this satelite.
            if (.not. lfound) then
               nsat = nsat + 1
               iprn = nsat
               kprn(iprn) = isvid(j)
            endif                 
         
C           Do gap analysis

c           1 if L1 is good
            if (issi(j,iobtypx(1)) .eq. 0 .or. 
     .         issi(j,iobtypx(1)) .ge. 2) then
              k1new(iprn) = 1
              ngood1(iprn) = ngood1(iprn) + 1
            else
              k1new(iprn) = 0
            endif                                

c           1 if L2 is good
            if (issi(j,iobtypx(2)) .eq. 0 .or. 
     .          issi(j,iobtypx(2)) .ge. 2) then
              k2new(iprn) = 1
              ngood2(iprn) = ngood2(iprn) + 1
            else
               k2new(iprn) = 0
            endif        

c           1 if LC (both L1 and L2) are good
            if (k1new(iprn) .eq. 1 .and. k2new(iprn) .eq. 1) then
              kcnew(iprn) = 1
              ngoodc(iprn) = ngoodc(iprn) + 1
              kount(iprn) = kount(iprn) + 1
            else
              kcnew(iprn) = 0
            endif   

            if (k1old(iprn).eq.1 .and. k1new(iprn) .eq. 0) then
               l1gap(iprn) = .true.
            endif
            if (k2old(iprn).eq.1 .and. k2new(iprn) .eq. 0) then
               l2gap(iprn) = .true.
            endif
            if (kcold(iprn).eq.1 .and. kcnew(iprn) .eq. 0) then
               lcgap(iprn) = .true.
            endif

            if (l1gap(iprn)) then
              if (k1old(iprn).eq.0 .and. k1new(iprn).eq.1) then
                  n1gap(iprn) = n1gap(iprn) + 1
                  l1gap(iprn) = .false.
              endif
            endif

            if (l2gap(iprn)) then
              if (k2old(iprn).eq.0 .and. k2new(iprn).eq.1) then
                 n2gap(iprn) = n2gap(iprn) + 1
                 l2gap(iprn) = .false.
              endif
            endif

             if (lcgap(iprn)) then
               if (kcold(iprn).eq.0 .and. kcnew(iprn).eq.1) then
                  ncgap(iprn) = ncgap(iprn) + 1
                  lcgap(iprn) = .false.
               endif
            endif

c           update the values for previous epoch
            k1old(iprn) = k1new(iprn)
            k2old(iprn) = k2new(iprn)
            kcold(iprn) = kcnew(iprn)

c           make bargraph
            if (mod(iepoch,nwide) .eq. 0) then
c              convert count in bin to tenths.
               kount(iprn) = 10*kount(iprn)/nwide
               if (kount(iprn) .eq. 0) then
                  write (letter,'(a)') '.'
               else if (kount(iprn) .le. 9) then
                  write (letter,'(i1)') kount(iprn)
               else
                 write (letter,'(a)')
     .            char(kount(iprn)-10+ichar('A'))
               endif
               kbin = iepoch/nwide
               bargph(iprn)(kbin:kbin) = letter
               kount(iprn) = 0
            endif
          enddo
        endif 
      enddo

      write(uscren,'(1x,2a,i6,a,f10.3,a)') frinex(1:nblen(frinex))
     .     ,' Data span found   : ',nepoch,' epochs  '
     .     ,secdif(iwkn1,sow1,iwkn0,sow0)/3600.,' hours'
      write(uscren,'(1x,2a,i6,a,f10.3,a)') frinex(1:nblen(frinex))
     .     ,' Data span expected: '
     .     ,nint(secdif(iwkn1,sow1,iwkn0,sow0)/rxint)
     . ,' epochs  ',secdif(iwkn1,sow1,iwkn0,sow0)/3600.,' hours'
      write (aline,'(1x,a,1x,a)') frinex(1:nblen(frinex)), 'Start'
      call wtime (uscren,iwkn0,sow0,'GPS',aline)
      write (aline,'(1x,a,1x,a)') frinex(1:nblen(frinex)), 'Stop '
      call wtime (uscren,iwkn1,sow1,'GPS',aline)
      write (uscren,*) ' ',frinex(1:nblen(frinex)),
     .'    PRN     Good observations     Number of Gaps  '
      write (uscren,*) ' ',frinex(1:nblen(frinex)),
     .'            L1     L2     LC     L1     L2     LC'
      do 800 i = 1,nsat
         write (uscren,'(1x,a,1x,7(i6,1x))') frinex(1:nblen(frinex)),
     .   kprn(i),ngood1(i),ngood2(i),ngoodc(i),
     .   n1gap(i),n2gap(i),ncgap(i)
 800  continue

c     time line
      do i =1,24
         jj(i) = i-1
      enddo
      write (uscren,220) frinex(1:nblen(frinex)),(jj(i),i=1,24)
  220 format (1x,'HHHH',1x,2x,1x,a,1x,24('-',i2.2,'-'))

c     bar graph
      do i = 1,nsat
         write (uscren,225) kprn(i),frinex(1:nblen(frinex)),bargph(i)
  225    format (1x,'PRN:',1x,i2,1x,a,1x,a96)
      enddo

      call report_stat('STATUS','RXSCAN','rxscan',' ',
     .'Normal end to RXSCAN',0)

      stop
      end











