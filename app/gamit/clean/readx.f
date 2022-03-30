      subroutine readx
     .  (luin,filnam,imode,isite,nepoch,nsat,nband,nwave)
c
c     read an X-file called filnam on logical unit luin
c
c     input:
c            luin   logical file unit
c            filnam file name
c            imode  = 0 pre-fit, no partials
c     output
c            nepoch total number of epochs   
c            nsat   total number of sattelites
c            isite  the index number of this site
c            nband  =  1 for monofrequency
c                      2 for bifrequency
c            nwave  =  1 for full-cycle data
c                   =  2 for half-cycle data
c            isprn  = sat PRN numbers
c
c     cleaned up by K. Feigl May 89
c     roots:
c        obsred - R.W. King
c        xhdred - R.W. King

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../../libraries/includes/freq_def.h'


      INTEGER DATTYP(MAXDAT),LAMBDA(MAXSAT,MAXDAT),NTEXT
     .      , ihr0,min0,latd,latm,lond,lonm,ndat,isat,nsat,id0
     .      , ierx,iy0,im0,ischan,iepoch,iamp1,iamp2,icount,jdoy0
     .      , mtime,isessn,i,j,k,m,len,rcpar

      REAL*8  sec0,offarp,height,seclat,seclon
     .      , data1,data2,data3,data4

      CHARACTER*80 TEXT(MAXTXT),buff80,prog_name
      character*20 rctype,rcvnum,anttyp,antnum
      CHARACTER*16 SITNAM
      DIMENSION offarp(3),ISCHAN(MAXSAT)
      character*1 latflag,lonflag
      character*4 rcvrsw
      character*6 antcod             
      character*16 satnam(maxsat)
      character*3 rxobtyp(maxdat)
      real*4 swver

c     variables needed to read svnav.dat for Glonass frequency
      integer*4 frqchn(maxsat),isvn,svnstart(5),svnstop(5)
      real*8 dumr8       
      character*20 dumc20
                           
      character*16 filnam
      integer*4 ioerr
      integer*4 idoy,jdoy,jm,jy,jh
      real*8 sec,clkdif

      integer luin, isite, imode, nepoch, nband, nwave
      logical lmarg,lgood,lbias

CD     print *,'Wavelengths ',clight/frq1,clight/frq2

c     get program name calling this routine
      len = rcpar(0,prog_name)

c     x-file is by definition pre-fit
      imode = 0

c     initialize start and end indexes
      do 3 i = 1,maxsat
         kk0(i,isite) = maxepc
         kk1(i,isite) = 1
  3   continue

      write (6,5) filnam
  5   format (1x,'Opening: ',a16)

c     open current X-file
      open (unit   = luin,
     .      file   = filnam,
     .      form   = 'formatted',
     .      status = 'old',
     .      iostat = ioerr,
     .      err    = 1000)

c     Read and display the X-file header.
      call XHDRED ( luin,6,6
     1,                   NEPOCH,INTER,ircint,isessn
     2,                   MTIME,IY0,IM0,ID0,IHR0,MIN0,SEC0
     3,                   nsat,ISCHAN,satnam,NDAT,DATTYP,rxobtyp,LAMBDA
     4,                   offarp,SITNAM,rcvrsw,swver,antcod
     5,                   rctype,rcvnum,anttyp,antnum
     6,                   LATFLAG,LATD,LATM,SECLAT
     7,                   LONFLAG,LOND,LONM,SECLON,HEIGHT
     8,                   NTEXT,TEXT,gnss)
    
c     get the GNSS 
      gnss = satnam(1)(1:1)
      if( gnss.eq.'N' ) gnss = 'G' 
c     put the times and sampling epoch into the common block
      iit0(1) = im0
      iit0(2) = id0
      iit0(3) = iy0
      tt00(1) = ihr0
      tt00(2) = min0
      tt00(3) = sec0
      jdoy0 = idoy(iy0,im0,id0)
c     determines the number of epoches between good data points
      if((inter .ge. ircint) .or. (inter .le. 0)) then
           inext = 1
      else
           inext = ircint/inter
      endif
      jnext(isite) = inext

c     copy lambda array and PRNs
      do 135 i=1,nsat
         isprn(i) = ischan(i)
         do 133 j=1,ndat
            lambds(isite,i,j) = lambda(i,j)
 133     continue
 135  continue

      do 200 j=1,nepoch
         do 205 k=1,nsat
            ierr(j,k,isite) = 1
            yl1 (j,k,isite) = 0.0d0
            yl2 (j,k,isite) = 0.0d0
            pr1 (j,k,isite) = 0.0d0
            pr2 (j,k,isite) = 0.0d0
  205    continue

         read (luin,'(/,a34)') buff80
         if (buff80(8:8) .eq. '0') then
            read (buff80,150,iostat=ioerr) iepoch,icount
  150       format (2i4)
            jdoy = 0
            jy = 0
            jh = 0
            jm = 0
            sec = 0.0d0
         else
            read(buff80,151,iostat=ioerr)
     .      iepoch,icount,jy,jdoy,jh,jm,sec
  151       format (2I4,I5,I4,2I3,F11.7)
         endif

         if (ioerr .ne. 0) then
            call report_stat('FATAL',prog_name,'readx',filnam,
     .     'Error reading X-file: ',ioerr)
         endif

  155    format (1x,i4,1x,i4,1x,f2.0)
CD        print 155,iepoch,icount

         if (iepoch .ne. j) then
            call report_stat('WARNING',prog_name,'readx',filnam,
     .     'Epoch mismatch reading X-file: ',0)
         endif

c        Get seconds after first epoch.
c        Note!  This is NOT the same as seconds of day.
         if (jdoy .gt. 0) then
           tag(j,isite) = clkdif(jdoy,jh,jm,sec,jdoy0,ihr0,min0,sec0)
         else
           tag(j,isite) = dble(iepoch - 1) * dble(inter)
         endif

         do 170 m=1,icount
  160       format (10x,2i2,2(1x,d22.15,1x,i3),2(2x,d22.15))
            if (ndat .eq. 4) then
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2,data3,data4
            else if (ndat .eq. 3) then
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2,data3
            else if (ndat .eq. 2) then
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1,data2,iamp2
            else
               read (luin,160,err=1000,end=1020,iostat=ioerr)
     .         ierx,isat,data1,iamp1
            endif

c           update start and end pointers
            if (lgood(ierx) .or. lmarg(ierx) .or. lbias(ierx)) then
               if (j.lt. kk0(isat,isite)) kk0(isat,isite) = j
               if (j.gt. kk1(isat,isite)) kk1(isat,isite) = j
            endif

c           error flag
            ierr(j,isat,isite)=ierx

c           factor for o's - when mixing receiver types
            if (ndat .ge. 1) yl1(j,isat,isite) = -data1
            if (ndat .ge. 2) yl2(j,isat,isite) = -data2
c           pseudoranges in cycles
            if (ndat .ge. 3) pr1(j,isat,isite) = data3*fL1(isat)/clight
            if (ndat .ge. 4) pr2(j,isat,isite) = data4*fL2(isat)/clight

CD           print *,'Phases',j,isat,isite,
CD    .      yl1(j,isat,isite),yl2(j,isat,isite)
CD           print *,'Ranges',j,isat,isite,
CD    .      pr1(j,isat,isite),pr2(j,isat,isite)
  170    continue
  200 continue   

c     assign the frequencies based on GNSS 
      do i=1,nsat
        if( gnss.eq.'G' )then
          fL1(i) = gps_f1
	       fL2(i) = gps_f2
        elseif( gnss.eq.'R') then 
          gnss = satnam(i)(1:1)
* MOD TAH 190702: Added place holder for antpwr to snav_read call
          call svnav_read(-1,iy0,jdoy0,ihr0,min0,gnss,isprn(i),isvn
     .                   , frqchn(i),dumc20,dumr8,dumr8,dumr8, dumr8
     .                   , svnstart,svnstop )
          fL1(i) = glonass_f1 + frqchn(i)*glonass_df1
          fL2(i) = glonass_f2 + frqchn(i)*glonass_df2
        elseif( gnss.eq.'C') then
          fL1(i) = beidou_f2
          fL2(i) = beidou_f7     
       elseif( gnss.eq.'E' ) then
          fL1(i) = galileo_f1
          fL2(i) = galileo_f5
       elseif( gnss.eq.'J' ) then
          call report_stat('FATAL','CLEAN','readx',' '
     .       ,'QZSS not yet supported',0 )
       elseif( gnss.eq.'I' ) then               
          fL1(i) = irnss_f9
          fL2(i) = irnss_f5
      else
         call report_stat('FATAL','MODEL','model',' '
     .                  ,'GNSS not recognized',0)
        endif   
      enddo

      write (6,'(//,5(1x,a7))') 'Channel','PRN','First','Last','Ndata'
      do 305 i = 1,nsat
        write (6,'(5(1x,i7))')  i,isprn(i),kk0(i,isite),kk1(i,isite),
     .                          max(int(kk1(i,isite)-kk0(i,isite)),0)
 305  continue

      close (luin)
 1000 continue
      if (ioerr .ne. 0) then
         call report_stat('FATAL',prog_name,'readx',filnam,
     .   'Error reading X-file: ',ioerr)
      endif
 1020 continue
      return
      end


