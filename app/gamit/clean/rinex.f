      subroutine rinex (lrinex,fname,isite,nobs,nsat)

C     read a RINEX file into the CVIEW arrays
c
c     input:
c            lrinex logical file unit
c            fname RINEX file name
c            isite  Index number of site
c     output
c            nepoch total number of epochs (including empties at start)
c            nobs   number of observations
c            nsat   total number of satellites

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../../libraries/includes/freq_def.h'

      logical debug/.false./

C     VALUES PASSED FROM RRINEX
c     File unit for RINEX file
      character*(*) fname
      integer     lrinex
C     Value set to .true. at end of file
      logical fend
C     Value set to .true. on error
      logical ferr
c     number and types of observables for selected GNSS
      integer*4 nobtyp                                   
      character*3 rxobtyp(maxobt)
      integer*4 iobtypx(6)
c     type of satellites
      character*1 asvid(maxchn)
c     number of satellites
      integer*4 nprn
C     sat id numbers (PRN)
      integer*4 isvid(maxchn)
c     GPS week number for this epoch
      integer*4 iwkn
C     second of GPS week for this epoch
      real*8  sow
C     L1, L2 doppler phase in cycles
      real*8  dofl1(maxchn), dofl2(maxchn)
C     L1, L2  pseudorange in meters
      real*8  prgl1(maxchn), prgl2(maxchn)
C     signal-to-noise ratio for L1, L2 phase
      integer*4 isnr(maxchn,4)
c     loss-of-lock indicator
      integer*4 illi(maxchn,4)
c     number of this epoch
      integer     nepoch
c     epoch counters
      integer*4 kepoch,kepoch0
      integer*2 kepoch2

      integer*4 ihr,imn,l
      real*8 sec

c     logical units for each file
      integer*4  lscren

c     I/O Error code
      integer    ioerr

c     read a RINEX header on logical unit lrinex, assumed open

c     ITEMS IN RINEX HEADER PASSED FROM RRXHED
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
c     receiver unit number
c** old integer definintion integer ircvr
      character*40 rcvnum
c     receiever type and SW version
      character*20 rctype
      character*20 rcvers
c     antenna serial number and type
c** old integer definintion integer iant  
      character*20 antnum
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
      real*8 rxint
c     data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1
c     time from origin time
      real*8 dsec
c     large number
      integer noffepc
      save noffepc                   
c     time type (always 'GPS' for mixed RINEX files)
      character*3 rxtime

c     first and last time tags:
      integer*4 iwkn0
      real*8    sow0,secdif

      logical lfound

c     pointer in PRN (sat id) array
      integer isat

c     miscellaneous variables
      integer*4 i,j,k,isite,iflag,itflag,iday,ierfl
      real*8 utcoff

c     total number of sats
      integer*4 nsat

c     number of observations
      integer*4 nobs

c     variables needed to read svnav.dat for Glonass frequency
      integer*4 frqchn(maxsat),isvn,svnstart(5),svnstop(5)
      real*8 dumr8
      character*20 dumc20
                               
c     function to return day of year
      integer*4 idoy

c     function to figure out error flags
      logical lmarg,lgood

c     report_stat message variable
      character*256 message

c     variables needed to get prog_name
      character*80 prog_name
      integer*4 len,rcpar

c     functions to handle I*2 min and max 
      integer*2 min02,max02

c     use for initial assignment of frequencies
      logical first/.true./

      data lscren/6/

c     get the main program name
      len = rcpar(0,prog_name)

c     initialize
      nepoch = 0 
      isat = 0
      do 20 i = 1,maxsat
         kk0(i,isite) = maxepc
         kk1(i,isite) = 1
 20   continue

c     open the file
      open (unit = lrinex,
     .      file = fname,
     .      iostat = ioerr,
     .      status = 'old')

      if (ioerr .eq. 0) then
         write (lscren,*) 'Opened ',fname
      else
         call report_stat('FATAL',prog_name,'rinex',fname,
     .   'Error opening rinex file: ',ioerr)
      endif

      urinex = lrinex
      call rrxhed ( debug,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1)
           
c     Put the times into the common block,
c     but only if they are not already there!
c     This will be the origin time.
      if (isite .eq. 1) then
         iit0(1) = irxmo0
         iit0(2) = irxdy0
         iit0(3) = irxyr0
c        while we are at it, initialize sat array.
         nsat = 0
         do i = 1,maxsat
            isprn(i) = -1
         enddo
      else if (iit0(1) .ne. irxmo0 .or.
     .         iit0(2) .ne. irxdy0 .or.
     .         iit0(3) .ne. irxyr0) then
         write(message,25)irxyr0,irxmo0,irxdy0
25       format('Wrong date in rinex file: ',i3,1x,i3,1x,i4,
     .   ' Skipping file: ')
         call report_stat('WARNING',prog_name,'rinex',fname,message,0)
         return
      endif

c     Put the interval into the common block, but only if not
c     already there.
      if (inter .eq. 0 .and. rxint .gt. 0.d0 ) then
         inter = int(rxint)
      else if (inter .eq. 0 .and. rxint .lt. 0.001d0 ) then
c        RINEX file claims interval of zero.
c        This occurs when RINEX translator selects all data.
c        Rather than guess, we will ask the operator.
         print *,'RINEX: The file claims 0 sec between epochs.'
         print *,'RINEX: Please enter your value.'
         read (5,*) inter
      else if ( int(rxint) .ne. inter) then
         write(message,'(a,i3,a,f7.2,a)') 
     .      'Data sampling mismatch; inter='
     .     ,inter,'  rxint=',rxint
         call report_stat('WARNING',prog_name,'rinex',' ',message,0)
      endif

c     Set origin time (t0)
c     week number = iwkn0
c     second of week = sow0
c     Choose this to be at midnight, ie, hr,min,sec=0
      iday = idoy(iit0(3),iit0(1),iit0(2))
      itflag= -4
      ihr = 0
      imn = 0
      sec = 0.0d0   
      call timcon(itflag,iwkn0,sow0,iit0(3),iday,ihr,imn,sec,utcoff)

      kepoch  = 0
      kepoch0 = 0

c     call wtime (6,iwkn0,sow0,'gps','Origin time')

c     determines the number of epochs between good data points
      inext = 1
      jnext(isite) = inext              
   
c     select the observables two phase and two pseudorange observables to be used 
 
      iflag = 1 
      call sel_obtyp(gnss,nobtyp,rxobtyp,iobtypx)


c     read until end of file
 10   continue 
         call rrinex( .false.,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .              , nprn,isvid,rxtime,iwkn,sow,nepoch
     .              , dofl1,dofl2,prgl1,prgl2,illi,isnr
     .              , anth,ante,antn,fend,ferr )  
c        Set a flag if an antenna event record enc    
c        do this for phase data only
         if (iflag .eq. 1) then

c           calculate seconds after origin time
            dsec = secdif(iwkn,sow,iwkn0,sow0)
            i = nint(dsec/dble(inter))

            if (isite .eq. 1 .and. kepoch .eq. 0) then
c              Note that this number is SAVED for posterity!
               noffepc = i
            endif

            kepoch = i - noffepc + 1

c           check to avoid disaster
c           check that time is not flowing backwards!
            if (kepoch .lt. 1) then
               call wtime (6,iwkn,sow,'gps','Time < t0.  Skip: ')
               kepoch = 1
            else if (kepoch .lt. kepoch0) then
                call wtime (6,iwkn,sow,'gps',
     .          'Time < t_last  Skip: ')
               kepoch = kepoch0
            else if (kepoch .gt. maxepc) then
               call wtime (6,iwkn,sow,'gps','Time > t_end.  Skip: ')
               kepoch = maxepc
            endif

c           write (buff20,*) 'epoch ',kepoch
c           call wtime (6,iwkn,sow,'gps',buff20)
c           print *, 'epoch ',kepoch,dsec,kepoch0

c           set intervening epochs to zero
c           and assign a time tag to avoid confusion.
            do i = 1,maxsat
               do l=kepoch0+1,kepoch
                  ierr(l,i,isite) = ignone
                  tag(l,isite) = (l + noffepc)  * inter
               enddo
            enddo

            tag(kepoch,isite) = dsec

c           figure out which sat we have
            do j =1,nprn
               lfound = .false.
               do i=1,nsat
                  if (isvid(j) .eq. isprn(i)) then
                     isat = i
                     lfound = .true.
                  endif
               enddo

c              This is the first time we have seen this satelite.
               if (.not. lfound) then
                  nsat = nsat + 1
                  isat = nsat
                  isprn(isat) = isvid(j)
               endif

               if (lambds(isite,isat,1) .eq. 0) then
c                 figure out the observable type
                  do k=1,nobtyp
                     if( rxobtyp(k)(1:2).eq.'L1') then
                        if( nwave1.eq.1 ) lambds(isite,isat,1) = -1
                        if( nwave1.eq.2 ) lambds(isite,isat,1) = -2
                     endif
                     if( rxobtyp(k)(1:2).eq.'L2') then
                        if( nwave2.eq.1 ) lambds(isite,isat,2) = -1
                        if( nwave2.eq.2 ) lambds(isite,isat,2) = -2
                     endif
                     if( rxobtyp(k)(1:2).eq.'P1') then
                        if( nwave1.eq.1 ) lambds(isite,isat,3) = +1
                        if( nwave1.eq.2 ) lambds(isite,isat,3) = +2
                     endif
                     if( rxobtyp(k)(1:2).eq.'C1') then
                        if( nwave1.eq.1 ) lambds(isite,isat,3) = +1
                        if( nwave1.eq.2 ) lambds(isite,isat,3) = +2
                     endif
                     if( rxobtyp(k)(1:2).eq.'P2') then
                        if( nwave2.eq.1 ) lambds(isite,isat,4) = +1
                        if( nwave2.eq.2 ) lambds(isite,isat,4) = +2
                     endif
c                 end loop on SVs
                  enddo
               endif

c              Cannot know receiver clock offset
c              until we modify lib/rinex to read RINEX 2 format.
               clk(kepoch,isat,isite) = 0.


c              Assign the frequencies based on GNSS 
               if( first ) then
                 if( gnss.eq.'G' )then
                   fL1(isat) = gps_f1
	                fL2(isat) = gps_f2
                 elseif( gnss.eq.'R') then 
* MOD TAH 190702: Added place holder for antpwr to snav_read call
                   call svnav_read(-1,iit0(3),iday,ihr,imn,gnss
     .                            , isprn(isat),isvn,frqchn(i)
     .                            , dumc20,dumr8,dumr8,dumr8, dumr8
     .                            , svnstart,svnstop )
                   fL1(isat) = glonass_f1 + frqchn(isat)*glonass_df1
                   fL2(isat) = glonass_f2 + frqchn(isat)*glonass_df2
                 elseif( gnss.eq.'C') then
                   fL1(isat) = beidou_f2
                   fL2(isat) = beidou_f7     
                 elseif( gnss.eq.'E' ) then
                   fL1(isat) = galileo_f1
                   fL2(isat) = galileo_f1
                 elseif( gnss.eq.'J' ) then
                   call report_stat('FATAL','CLEAN','rinex',' '
     .                 ,'QZSS not yet supported',0 )
                 elseif( gnss.eq.'I' ) then 
                   fL1(isat) = irnss_f9
                   fL2(isat) = irnss_f5 
                 else
                   call report_stat('FATAL','MODEL','model',' '
     .                  ,'GNSS not recognized',0)
                 endif  
                 first = .false.
               endif 
                
c              Phase observables in cycles
c              RRINEX returns DOPPLER convention,
c              CVIEW plots PSEUDORANGE convention
               yl1(kepoch,isat,isite) = -dofl1(j)
               yl2(kepoch,isat,isite) = -dofl2(j)
    
c              Pseudoranges in cycles
               pr1(kepoch,isat,isite) = prgl1(j)*fL1(isat)/clight
               pr2(kepoch,isat,isite) = prgl2(j)*fL2(isat)/clight
c              error code (ad hoc.)
c              Keep in mind that isnr = 0 is of uknown quality.
               if (illi(j,iobtypx(1)).eq.1 .or. illi(j,2).eq.1) then
                  ierfl = igbias
               else if (isnr(j,iobtypx(1)).gt.1 .and. 
     .                  isnr(j,iobtypx(2)).gt.1) then
                  ierfl = iggood
               else if (isnr(j,iobtypx(1)).eq.1 .or. 
     .                  isnr(j,iobtypx(2)).eq.1) then
                  ierfl = iglamp
               else if (isnr(j,iobtypx(1)).eq.0 .or.  
     .                  isnr(j,iobtypx(2)).eq.0) then
                  ierfl = iggood
               else
                  ierfl = ignone
               endif

               ierr(kepoch,isat,isite) = ierfl

c              update start and end pointers
               if (lmarg(ierfl).or.lgood(ierfl)) then  
                  kepoch2 = kepoch
                  kk0(isat,isite) = min02(kk0(isat,isite),kepoch2)
                  kk1(isat,isite) = max02(kk1(isat,isite),kepoch2)
               endif

c              hold last good epoch number
               kepoch0 = kepoch
            enddo
         else
c           if no data, then flag this epoch as such
            do  isat =1,maxsat
               ierr(kepoch+1,isat,isite) = ignone
            enddo
         endif

      if (.not. fend .and. kepoch .lt. maxepc) then
         goto 10
      else if (fend .and. kepoch .lt. maxepc) then
c         fill out empty epochs.
          do i = 1,maxsat
             do l = kepoch+1,maxepc
                ierr(l,i,isite) = ignone
                tag(l,isite) = (l + noffepc)  * inter
             enddo
          enddo
      else
         close (lrinex)
         if (kepoch .ge. maxepc) then
           write (message,'(a,i5,a)') 'Too many epochs, Read first :'
     .         ,maxepc,' epochs.'
           call report_stat('WARNING',prog_name,'rinex',' ',message,0)
         endif
      endif    

      write (6,'(//,5(1x,a7))') 'Channel','PRN','First','Last','Ndata'
      do i = 1,nsat
        write (6,'(5(1x,i7))')  i,isprn(i),kk0(i,isite),kk1(i,isite),
     .                          max(kk1(i,isite)-kk0(i,isite)+1,0)
        nobs = max(int(kk1(i,isite)),nobs)
      enddo
           
      return
      end




