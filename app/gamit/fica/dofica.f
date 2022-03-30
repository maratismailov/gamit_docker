      subroutine dofica (
     .                debug,iflag,fend,ferr,nprn,svid,tmode,igpswk
     .              , gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .              , bcephem,bcclock,subfr1,iwknfile,weekset
     .              , iwknstart,sowstart)

c
c     Read a FICA format data file open on logical unit lu
c       Written originally by K.Feigl
c       Broken into subroutines by R. King  3 Octobe 88
c
c         iflag = 1   Record read contained phase & pseudorange data
c               = 2   Record read contained ephemeris data
c
c     Note that not all the returned values will be valid, depending
c     on the value of iflag. For example, bcephem and bcclock
c     will be bogus if iflag = 1, and the phase data will be
c     meaningless if iflag = 2. Iflag is returned to the main program
c     according to what kind of record was read in.
c
c     Block ID numbers currently allowed:
c
c     For TI 4100 : 6, 55, 401, 9, 101, 70, 80
c     For MACROMETER II: 670, 1080
c     For Minimac : 1101, 1180
c     For Trmible : 1201, 1280
c     For Rogue   : 1301, 1380

c     Note: block 70 was defined at MIT, but not approved by Clynch
c     at ARL.  The ARL approved block is numbered 80,  not 70.
c     This code will continue to read the MIT flavor of block 70,
c     but provide a warning.

c     For TI 4100 all times in GPS time

      implicit none

      include '../includes/makex.h'

c     passed values
C     .true. at end of file
      logical     fend,
C     .true. to print things
     &            debug,
C     .true. on error
     &            ferr,
C     .true. once we have a week num
     &            weekset

C     1 for data, 2 for ephemeris
      integer*4   iflag,
C     for broadcast ephem
     &            nprn,
C     sat id num PRN
     &            svid(maxchn),
C     tracking mode
     &            tmode(maxchn),
C     gps week number for epoch
     &            igpswk(maxchn),
C     week number of sceno start
     &            iwknstart,
C     quality vector L1,L2
     &            iqvec(maxchn,2)
C     gps sec of week for epoch
      real*8      gpssec(maxchn),
C     L1 doppler phase (DATA!)
     &            dofl1(maxchn),
C     L2 doppler phase (DATA!)
     &            dofl2(maxchn),
C     L1 pseudorange
     &            prgl1(maxchn),
C     L2 pseudorange
     &            prgl2(maxchn),
C     signal to noise ratio L1,L2
     &            denrat(maxchn,2),
C     second of week for sceno start
     &            sowstart,
C     broadcast ephemeris params
     &            bcephem(16),
C     broadcast clock params
     &            bcclock(6),
C     additional SV quantites from sub-frame 1
     .            subfr1(8)

c     other values
c
c     values from FICA header files
      character*16 auser
     .,            asite
C     reciever serial number
     .,            arcvr
C     operator input antenna ht.
     .,            antht
C     atmospheric pressure
      real*8       apress
C     atmospheric temperature
     .,            atemp
C     atmospheric pressure
     .,            ahumid
c     collection interval
     .,            asampl
C     software version
     .,            swveru
C     L1 phase center wrt monument
     .,            offsl1(3)
C     L2 phase center wrt monument
     .,            offsl2(3)


C     error count
      integer*4    errct
      integer*4    ioerr


c     things for FICA
C     number of elements in array
      integer nf,ni,nc
C     number of the block
      integer iblkid
c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     week number, from block 9
      integer*4 iwknfile
c     error message for report_stat routine
      character*256 message

      integer i

C     count the number of block 70s & 670s
      integer i70,i670

      save i70,i670
      data i70/0/,i670/0/


      iflag = -1
      errct = 0           
c     --st to true on error
      ferr = .false. 
c     --set to true on EOF
      fend = .false.       

c     initialize SV prn array
      do i=1,maxchn
       svid(i) = 0
      enddo

 100  call rfica (uficaf,uinfor,iblkid,ff,fi,fc,nf,ni,nc,ioerr)
c      print *,'DOFICA RFICA iblkid ioerr iflag weekset iwknfile '
c     .       ,iblkid,ioerr,iflag,weekset,iwknfile
      

C     end of file
      if (ioerr .eq. -1) then
         fend = .true.
         return
      else if (ioerr. gt. 0) then
         weekset = .false.
         write(message,'(a)') 'Error reading FICA file'
         call report_stat('WARNING','MAKEX','dofica',fficaf
     .                   ,'Error reading FICA file ',ioerr)
         write (uinfor,'(a)') message
         ferr = .true.
         iflag = -1
         return
      endif


c     TI 4100 Blocks

      if ( iblkid .eq.6 .or. iblkid.eq.55 ) then
        call blk6( debug,iflag,svid,tmode,igpswk,gpssec
     .          , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .          , iwknfile,weekset,iwknstart,sowstart
     .          , ff,fi,fc,nf,ni,nc )  

      else if ( iblkid .eq.70 ) then
c       issue a warning
        i70 = i70+1
         write (message,'(2a)') 'Found block 70--'
     .     ,'if the FICA files are not MIT in-house, something wrong'
         call report_stat('WARNING','MAKEX','dofica',' ',message,0)
         write(uinfor,'(a)') message
         call blk70( debug,iflag,svid,tmode,igpswk,gpssec
     .             , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .             , iwknfile,weekset,iwknstart,sowstart
     .             , ff,fi,fc,nf,ni,nc)

      else if ( iblkid .eq.80 ) then
         call blk80 ( debug,iflag,svid,tmode,igpswk,gpssec
     .              , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .              , iwknfile,weekset,iwknstart,sowstart
     .              , ff,fi,fc,nf,ni,nc) 

      else if ( iblkid .eq. 101) then
         call blk101( debug,iflag,apress,atemp,ahumid,asampl,swveru
     .              , auser,asite,arcvr,antht
     .              , ff,fi,fc,nf,ni,nc )

      else if ( iblkid .eq.9  ) then
         call blk9 ( debug,iflag,nprn,bcephem,bcclock,subfr1
     .             , iwknfile,weekset
     .             , ff,fi,fc,nf,ni,nc )
          igpswk(1)= iwknfile

      else if ( iblkid .eq.401 ) then
         call blk401 ( debug,iflag,svid,tmode,igpswk,gpssec
     .               , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .               , iwknfile,weekset,ff,fi,fc,nf,ni,nc ) 


c   MACROMETER II Blocks

      else if ( iblkid. eq.1001 ) then
         call blk1001
     .     (debug,iflag,offsl1,offsl2,arcvr,asite,ff,fi,fc,nf,ni,nc)

      else if ( iblkid .eq.670 ) then
c        issue a warning
         i670 = i670+1
         if (i670 .eq. 1) then
           write (message,'(2a)') 'Found block 670--'
     .       ,'if the FICA files are not MIT in-house, something wrong'
           call report_stat('WARNING','MAKEX','dofica',' ',message,0)
           write(uinfor,'(a)') message
         endif
         call blk670 ( debug,iflag,svid,tmode,igpswk,gpssec
     .               , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .               , iwknfile,weekset,ff,fi,fc,nf,ni,nc )  

      else if ( iblkid .eq.1080 ) then
         call blk1080 ( debug,iflag,svid,tmode,igpswk,gpssec
     .                , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .                , iwknfile,weekset,ff,fi,fc,nf,ni,nc )  

c    Provisional Minimac blocks

      else if ( iblkid .eq.1180 ) then
         call blk1180 ( debug,iflag,svid,tmode,igpswk,gpssec
     .             , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .             , iwknfile,weekset,ff,fi,fc,nf,ni,nc )

c    Provisional Trimble CIGNET blocks

      else if ( iblkid .eq.1280 ) then
         call blk1280 ( debug,iflag,svid,tmode,igpswk,gpssec
     .                , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .                , iwknfile,weekset,ff,fi,fc,nf,ni,nc )

c    Provisional Rogue CIGNET blocks

      else if ( iblkid .eq.1380 ) then
         call blk1380 ( debug,iflag,svid,tmode,igpswk,gpssec
     .                , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .                , iwknfile,weekset,ff,fi,fc,nf,ni,nc )   

c     If unwanted block id, go read another record

      else
         goto 100
      endif

      if (debug) then
         call wtime (uscren,iwknfile,gpssec(1),'gps','DOFICA:')
         write (uscren,*) 'DOFICA: PRNs ',(svid(i),i=1,maxchn)
      endif
           
      return
      end
