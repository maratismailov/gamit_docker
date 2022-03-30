      PROGRAM TTOICS
C
C     Written by Yehuda Bock
c     Mods for new GSATEL rwk 12/28/91
C
C     Dump a T-file to an ASCII file with full state vector elements
C      (positions and velocities)
C
      implicit none

      include '../includes/dimpar.h'

      CHARACTER*1 COMT(80,3)
      CHARACTER*2 BUF1, BUF2
      CHARACTER*16 SATNAM(maxsat)
      character*38 text
      character*80 head,tmpnam
      character*120 version
      character*256 message

      integer*4    iscrn,iterm,iprnt,iutin,itstp,ioerr,iclarg
     .        ,    iepoch,nepchs,jlast,jil,ji0,jnow,iy1,iy2,iendf
     .        ,    isat,nsat,nics,ite,itb
     .        ,    nintrs,nepcht
     .        ,    jd,iday,imin,imon,idoy,jds,ihr,iyr
     .        ,    nstrtt,ii,ksat,julday,i,j

      real*8       xpole,ypole,rvec,trun,ut1,ytrp,yy,satcor,tstp,tee
     .        ,    sdelt,satics,tbb,ts,sec,sod
C
      DIMENSION RVEC(MAXYT2)
      DIMENSION YTRP(5,2,MAXYT2,MAXSAT),YY(10,MAXYTP,MAXSAT)
      dimension ITE(3),TEE(3),ITB(3),TBB(3),ITSTP(3),TSTP(3)
      DIMENSION SATCOR(6,maxsat),SATICS(maxorb,maxsat)
C
      iterm=5
      iscrn=6
      iprnt=7
      iutin=8

c     Print the version and machine

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5)version
    5 format('Started TTOICS ',a120)
      call report_stat('STATUS','TTOICS','orbits/ttoics',' ',message,0)
c
C     Open the T-file and output ASCII file

c     get input file name from command line
c     ask if it is not there.
      ii = iclarg(1,tmpnam)
      if (ii .le. 0) then
         WRITE(ISCRN,11)
11       FORMAT(/,' Enter tabular-ephemeris (T-) file name : ')
         READ(ITERM,'(A)') TMPNAM
      endif

c     open the file
      call topens(tmpnam,'old',iutin,ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','TTOICS','orbits/ttoics',tmpnam,
     .  'Error opening T-file: ',ioerr)
      endif

c     output file name
      tmpnam(1:1)='a'

      OPEN (UNIT=IPRNT,FILE=TMPNAM,STATUS='UNKNOWN',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','TTOICS','orbits/ttoics',tmpnam,
     .  'Error opening ASCII ICs file: ',ioerr)
      endif
c
      write(iscrn,'(a)')' '
      WRITE(IPRNT,21) 't'//TMPNAM(2:80)
      WRITE(message,21) 't'//TMPNAM(2:80)
21    FORMAT('Dump by TTOICS of Tabular Ephemeris T-file: ',A16)
      call report_stat('STATUS','TTOICS','orbits/ttoics',' ',message,0)

C Read and Write the first header record of the T-file

      READ(IUTIN)     HEAD,ITE,TEE,ITB,TBB,ITSTP,TSTP,SDELT,NEPCHT,
     1            UT1,XPOLE,YPOLE
      WRITE(iprnt,70) HEAD,ITE,TEE,ITB,TBB,ITSTP,TSTP,SDELT,NEPCHT,
     1            UT1,XPOLE,YPOLE

70    FORMAT(//,' T-File Header Information:',//,1X,A80,/,
     $ /,3I5,2F4.0,F7.3,2X,'Epoch of initial conditions',
     $ /,3I5,2F4.0,F7.3,2X,'Ephemeris start',
     $ /,3I5,2F4.0,F7.3,2X,'Ephemeris end',
     $ /,F10.3,'  Tabular interval (sec) ;',3X,I5,' T-file epochs'
     $ ,/,1X,'UT1,XPOLE,YPOLE : ',3F8.3)

c Read the second header record of the T-file

      WRITE(BUF2,'(A2)') HEAD(79:80)
      READ (BUF2,'(I2)') NICS
      IF (NICS.le.0) then 
         write(message,'(a,i4)') 
     .       'Old T-file; number of ICs unknown: ',nics
         call report_stat('WARNING','TTOICS','orbits/ttoics',' '
     .                   ,message,0)
          NICS=maxorb
      endif

c**   Hardwire for NOAA comparison
      if( nics.eq.8 ) nics = 9

      READ(IUTIN) COMT,nsat,NINTRS,
     .   (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)
      write(IPRNT,75) ((COMT(i,j),i=1,80),j=1,3),NSAT,NINTRS
   75 format(' T-file comment :',/,3(80a1,/),
     1       ' nsat,nintrs :',2i5,/,
     2       ' Satellite name and ics :')
c     text='(__(/,1X,a6,11x,3d12.5/,18x,_d12.5))'
      text='(__(/,1X,a6,11x,3f14.5/,18x,__f14.9))'
      WRITE (BUF2,'(I2)') NSAT
      READ  (BUF2,'(A2)') TEXT(2:3)
      WRITE (BUF1,'(I2)') NICS-3
      READ  (BUF1,'(A2)') TEXT(29:30)
      write(iprnt,fmt=text)
     $ (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)

c Initialize the interpolation pointers

      JI0= 0
      JIL= 0
      IY1= 0
      IY2= 0
      JLAST= 0
      IENDF= 0
      trun=0.d0

c Convert the start T-file time to sec-of-day
           
      ts = tbb(1)*3600.d0 + tbb(2)*60.d0 + tbb(3)
      jds = julday(itb(1),itb(2),itb(3))

c For correct interpolation of velocities, must skip the first
c and last 5 records of the T-file

       nstrtt = 6
       nepchs = nepcht - 10
       trun = (nstrtt-2) * sdelt

c Sequence through the T-file, interpolating position and velocity

      write(iscrn,80) nepchs
   80 format(/,1x,'Evaluating position and velocity for ',i4
     .        ,' epochs')
      write(iprnt,81)
   81 format(//,1x,'Position and velocity beginning at epoch 6'
     .      ,//,1x,'Times are UTC',/)

      do 100 iepoch=1,nepchs

      trun = trun + sdelt

c        Compute and write the time of the epoch
         jd= jds
         sod= ts
         call timinc( jd,sod,trun )
         call dayjul( jd,iyr,idoy )
         call monday( idoy,imon,iday,iyr )
c        avoid roundoff by assuming that the epochs are even seconds
         sod = anint(sod)
         call ds2hms( iyr,idoy,sod,ihr,imin,sec )
         write(iprnt,60) iyr,imon,iday,ihr,imin,sec
   60    format(//,1x,i4,4(1x,i2),1x,f10.7)

      DO 55 KSAT=1,NSAT
          call gsatel( 2,trun,iutin,ksat,rvec,ytrp,yy,nsat,sdelt
     .               , nintrs,ji0,jil,iy1,iy2,jlast,jnow,iendf,nepcht )
          DO 55 I=1,6
   55     SATCOR(I,KSAT)= RVEC(I)

c        Write the coordinates for each satellite
         WRITE(iprnt,110) ((SATCOR(I,ISAT),I=1,6),ISAT=1,NSAT)
110      FORMAT(/(6F20.9))

100   CONTINUE

c Close the files and display the output file name

      close(unit=iutin)
      close(unit=iprnt)

      call report_stat('STATUS','TTOICS','orbits/ttoics',tmpnam,
     .'ASCII ICs dump on file : ',0)
      call report_stat('STATUS','TTOICS','orbits/ttoics',' ',
     .'Normal end in TTOICS ',0)

      stop
      end

