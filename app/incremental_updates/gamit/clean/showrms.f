      program showrms
c
c purpose:      present rms in a table format in terms of percent
c
c input:      rms file name (default=rms.srt)
c            rms file (output of scandd or scanrms)
c
c output:      screen display of the rms distribution
c
c author:      Peng Fang, pfang@ucsd.edu
c
c date:            10/23/91
c update:      08/28/95
c
      include '../includes/dimpar.h'
c     allow an extra row and column
      integer maxstn,maxchn
      parameter (maxstn=maxsit+1, maxchn=maxsat+1)
      integer iflg,nchn,luin,nsit,iop
      real    rms(maxchn,maxchn,maxstn,maxstn),tbl(maxchn,maxstn),grand
      character*8 filnam

      common /tables/ rms,tbl

      data luin/10/,nchn/0/,nsit/0/,iflg/2/
        
      do while (iop.lt.1.or.iop.gt.3)

      write(*,'("Enter 1 for rms.qui 2 for rms.ful 3 for rms.tot")')
      read(*,*) iop
      filnam='rms.qui'
      if (iop.eq.2) filnam='rms.ful'
      if (iop.eq.3) filnam='rms.tot'
      enddo
c
      open(unit=luin,file=filnam,status='old')
c
      call rdrms(luin,nchn,nsit,iop)
      do while (iflg.ge.1.and.iflg.le.2)
      call getsum(nchn,nsit,grand,iflg)
      call showtb(nchn,nsit,grand)
      call rdrms(luin,nchn,nsit,iop)
      call remove(luin,nchn,nsit,grand,iflg)
      enddo
      end

c ---------------------------------
      subroutine remove(luin,nchn,nsit,grand,iflg)

      include '../includes/dimpar.h'
c     allow an extra row and column
      integer maxstn,maxchn
      parameter (maxstn=maxsit+1, maxchn=maxsat+1)
      integer luin,nchn_rm,nsit_rm,iflg,nchn,nsit
     .      , i,j,ii,jj,i1,i2,j1,j2
      integer chn_rm(maxchn), sit_rm(maxstn)
      real    rms(maxchn,maxchn,maxstn,maxstn),tbl(maxchn,maxstn),grand
      character*120 tmptxt

      common /tables/ rms,tbl

      iflg=0
      do i=1,nchn
      chn_rm(i)=0
      enddo
      do i=1,nsit
      sit_rm(i)=0
      enddo
      do while (iflg.ne.1.and.iflg.ne.2)
      write(*,'("Enter 1 fix, 2 new, 0 quit, otherwise help")')
      read *,iflg
      if (iflg.eq.0) stop
      if (iflg.ne.1.and.iflg.ne.2) call help
      enddo
      write(*,'("Enter channel #s to be removed (- for all sites)")')
      read(*,'(a)') tmptxt
      read(tmptxt,*,end=60) chn_rm
60      write(*,'("Enter site #s to be removed (- for all channels)")')
      read(*,'(a)') tmptxt
      read(tmptxt,*,end=80) sit_rm
c
80      i=1
      do while (chn_rm(i).ne.0)
            i=i+1
      enddo
      nchn_rm=i-1
       i=1
      do while (sit_rm(i).ne.0)
            i=i+1
      enddo
      nsit_rm=i-1
c
      do j1=1,nsit_rm
c remove channel related on all sites
      if (sit_rm(j1).lt.0) then
            do i=1,nchn
            do ii=1,nchn
            do j=1,nsit
            rms(i,ii,j,iabs(sit_rm(j1)))=0
            rms(i,ii,iabs(sit_rm(j1)),j)=0
            enddo
            enddo
            enddo
      endif
      do i1=1,nchn_rm
c remove channel and site pairs
      if (sit_rm(j1).gt.0.and.chn_rm(i1).gt.0) then
            do i=1,nchn
            do j=1,nsit
            rms(i,chn_rm(i1),j,sit_rm(j1))=0
            rms(chn_rm(i1),i,j,sit_rm(j1))=0
            rms(i,chn_rm(i1),sit_rm(j1),j)=0
            rms(chn_rm(i1),i,sit_rm(j1),j)=0
            enddo
            enddo
      endif
      do j2=j+1,nsit_rm
      do i2=i+1,nchn_rm
c remove channel_1 channel_2 site_1 site_2 combination
      if (chn_rm(i1).gt.0.and.chn_rm(i2).gt.0.and.
     *            sit_rm(j1).gt.0.and.sit_rm(j2).gt.0) then
            rms(chn_rm(i1),chn_rm(i2),sit_rm(j1),sit_rm(j2))=0
      endif
      enddo
      enddo
      enddo
      enddo
c
      do i1=1,nchn_rm
c remove site related on all channels
      if (chn_rm(i1).lt.0) then
            do i=1,nchn
            do j=1,nsit
            do jj=1,nsit
            rms(iabs(chn_rm(i1)),i,j,jj)=0
            rms(i,iabs(chn_rm(i1)),j,jj)=0
            enddo
            enddo
            enddo
      endif
      enddo
      return
      end
c ---------------------------------
      subroutine showtb(nchn,nsit,grand)
c
      include '../includes/dimpar.h'
c     allow an extra row and column
      integer maxstn,maxchn
      parameter (maxstn=maxsit+1, maxchn=maxsat+1)
      integer nchn,nsit,i,j,k
      real    rms(maxchn,maxchn,maxstn,maxstn),tbl(maxchn,maxstn),grand

      common /tables/ rms,tbl

      write(*,
     .'("RMS DISTRIBUTION  (total rms = 1000)  example: 234 = 23.4%")')
      write(*,'("stn|sum|",62(1h-))')
c fill table
      do j=1,nsit
      do i=1,nchn
c get % rms
      tbl(i,j)=tbl(i,j)*1000/grand
c sum up all sites for one channel
      tbl(i,nsit+1)=tbl(i,nsit+1)+tbl(i,j)
c sum up all channels for one site
      tbl(nchn+1,j)=tbl(nchn+1,j)+tbl(i,j)
      enddo
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
      write(*,900) j,int(tbl(nchn+1,j)),(int(tbl(k,j)),k=1,nchn)
900      format(i2,1h:,i4,1h:,50i4)
      enddo
      write(*,'(70(1h-))')
      write(*,'("sum --> ",50i4)') (int(tbl(k,nsit+1)),k=1,nchn)
      write(*,'("chn --> ",50i4)') (k,k=1,nchn)
      return
      end
c ---------------------------------
      subroutine getsum(nchn,nsit,grand,iflg)
c
      include '../includes/dimpar.h'
c     allow an extra row and column
      integer maxstn,maxchn
      parameter (maxstn=maxsit+1, maxchn=maxsat+1)
      integer nchn,nsit,iflg,i1,i2,j1,j2,i,j
      real    rms(maxchn,maxchn,maxstn,maxstn),tbl(maxchn,maxstn),grand

      common /tables/ rms,tbl

c clean table
      do i=1,nchn+1
      do j=1,nsit+1
      tbl(i,j)=0
      enddo
      enddo
      if (iflg.eq.2) grand=0
c
      do i1=1,nchn
      do i2=1,nchn
      do j1=1,nsit
      do j2=1,nsit
      tbl(i1,j1)=tbl(i1,j1)+rms(i1,i2,j1,j2)
      tbl(i2,j1)=tbl(i2,j1)+rms(i1,i2,j1,j2)
      tbl(i1,j2)=tbl(i1,j2)+rms(i1,i2,j1,j2)
      tbl(i2,j2)=tbl(i2,j2)+rms(i1,i2,j1,j2)
      if (iflg.eq.2) grand=grand+rms(i1,i2,j1,j2)
      enddo
      enddo
      enddo
      enddo
      if (iflg.eq.2) grand=grand*4
      if (grand.eq.0) then
            write(*,'("All RMSs = 0  Please check input rms file")')
            stop
      endif
      return
      end
c ---------------------------------
      subroutine rdrms(luin,nchn,nsit,iop)
c
      include '../includes/dimpar.h'
c     allow an extra row and column
      integer maxstn,maxchn,iop,ioerr
      parameter (maxstn=maxsit+1, maxchn=maxsat+1)
      integer nchn,ist1,ist2,luin,nsit,ntot,i1,i2,j1,j2,ich1,ich2
      real rms(maxchn,maxchn,maxstn,maxstn),tbl(maxchn,maxstn),
     *     rmstmp(3)

      common /tables/ rms,tbl

c clean rms
      do i1=1,maxchn
      do i2=1,maxchn
      do j1=1,maxstn
      do j2=1,maxstn
      rms(i1,i2,j1,j2)=0
      enddo
      enddo
      enddo
      enddo
c
      ntot=0
      do while (ntot.ge.0)
        read(luin,900,iostat=ioerr) rmstmp,ich1,ich2,ist1,ist2
        if( ioerr.ne.0) exit
900      format(t7,f9.2,5x,f9.2,4x,f9.2,t52,4i8)
c store the rms
      rms(ich1,ich2,ist1,ist2)=rmstmp(iop)
c get highest numbers of channel and site
      nchn=max0(ich1,ich2,nchn)
      nsit=max0(ist1,ist2,nsit)
      ntot=ntot+1
      enddo
      rewind (luin)
      if (ntot.gt.2000) print 920
920      format(' WARNING: the number of rms is greater then 2000',/
     *    ,'It might be better if a higher threshold is used in sorter')
      return
      end
c ---------------------------------

      subroutine help

c **  character*80 tmptxt
c **  for now, embed the help text in the code
c **  open(unit=64,file='showrms.hlp',status='old')
c **  i=1
c **  do while (i.ge.1)
c **  read(64,'(a)',end=90) tmptxt
c **  print *, tmptxt
c      if (mod(i,24).eq.0) then
c            print *,('any key to continue')
c            read *,tmptxt
c      endif
c **  i=i+1
c **  enddo
      character*60 help_txt(18)

      data help_txt/
     .  '                                                            ',
     .  'fix   means the total rms remains unchanged even after      ',
     .  '      a simulated removal of channel/site(s) is performed.  ',
     .  '      This mode is good for comparing effects ofdifferent   ',
     .  '       removal trials.                                      ',
     .  '                                                            ',
     .  'new   means a new total rms will be computed after a        ',
     .  '      simulated removal of channel/site(s) is performed.    ',
     .  '      This mode is good for inspecting rms distribution     ',
     .  '      after removing known bad channel/site(s).             ',
     .  '                                                            ',
     .  '  To remove one channel at one site, use positive entries.  ',
     .  '  To remove a channel (site) at all sites (channel), use    ',
     .  '  negative entries.  Multiple entries are allowed: e.g.     ',
     .  '     sit:  4                                                ',
     .  '     chn:  3 -6                                             ',
     .  '  to exclude channel 3 at site 4 and channel 6 at all sites.',
     .  '                                                            '/

      write(*,'(a60)') help_txt

90    return
      end

