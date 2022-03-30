      PROGRAM TTOASC
C
C     Written by Yehuda Bock and Robert King
C
C     Dump a T-file to an ASCII file
C
      implicit none

      include '../includes/dimpar.h'

      CHARACTER*1 COMT(80,3)
      CHARACTER*2 BUF2
      CHARACTER*16 SATNAM(maxsat)
      character*49 text
      character*80 HEAD,tmpnam
      character*120 version
      character*256 message

      integer*4    iscrn,iterm,iprnt,iutin,ioerr,itb,ite,itstp
     .           , nics,isat,nsat,nepchs,nintrs,iys,idns,ihs,iwkn
     .           , iclarg,ii,i,j,iepoch,idoy,mins

      real*8       xpole,ypole,ut1,satics,tstp,tbb,tee,x,sdelt,secs
     .           , sow,utcoff

c     number of terms per sat:
c     X,Y,Z:                     3
c     6 elements                 6
c     3 radiation pressures      9
c     total:                    18

      dimension ITE(3),TEE(3),ITB(3),TBB(3),ITSTP(3),TSTP(3)
      DIMENSION X(MAXYT2,MAXSAT),SATICS(maxorb,MAXSAT)

      iterm=5
      iscrn=6
      iprnt=7
      iutin=8

c     Print the version and machine

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5) version
    5 format('Started TTOASC ',a120)
      call report_stat('STATUS','TTOASC','orbits/ttoasc',' ',message,0)
C
C     Version 1.00     Sent to NGS - 12/04/1987 YB
c     version 2.1      Installed on Apollo SCCS kurt 880119
c     Version 6.2      Compiled INTEGER*4 by rwk     890822
c     version 7.2      Library dimpar.fti, maxtrm included in main MHM 890825
c     version 8.1      Longer file name and command line argument.  Kurt 910703
C
C     Open the T-file and output ASCII file
C
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
        call report_stat('FATAL','TTOASC','orbits/ttoasc',tmpnam,
     .  'Error opening T-file: ',ioerr)
      endif

c     output file name
      tmpnam(1:1)='A'

      OPEN (UNIT=IPRNT,FILE=TMPNAM,STATUS='UNKNOWN',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','TTOASC','orbits/ttoasc',tmpnam,
     .  'Error opening ASCII T-file dump: ',ioerr)
      endif
c
      write(iscrn,'(a)')' '
      WRITE(IPRNT,21) 't'//TMPNAM(2:80)
      WRITE(message,21) 't'//TMPNAM(2:80)
21    FORMAT('Dump by TTOASC of Tabular Ephemeris T-file: ',A16)
      call report_stat('STATUS','TTOASC','orbits/ttoasc',' ',message,0)

C     Read and Write the first header record of the T-file
      READ(IUTIN)     HEAD,ITE,TEE,ITB,TBB,ITSTP,TSTP,SDELT,NEPCHS,
     1            UT1,XPOLE,YPOLE
      WRITE(iprnt,70) HEAD,ITE,TEE,ITB,TBB,ITSTP,TSTP,SDELT,NEPCHS,
     1            UT1,XPOLE,YPOLE

70    FORMAT(//,' T-File Header Information:',//,1X,A80,/,
     $ /,3I5,2F4.0,F7.3,2X,'Epoch of initial conditions',
     $ /,3I5,2F4.0,F7.3,2X,'Ephemeris start',
     $ /,3I5,2F4.0,F7.3,2X,'Ephemeris end',
     $ /,F10.3,'  Tabular interval (sec) ;',3X,I5,' T-file epochs'
     $ ,/,1X,'UT1,XPOLE,YPOLE : ',3F8.3)

c     Read the second header record of the T-file
      WRITE(BUF2,'(A2)') HEAD(79:80)
      READ (BUF2,'(I2)') NICS
      IF (NICS.le.0) then
      write(message,'(a,i3)') 'Old T-file, number of ICs unknown: ',nics  
      call report_stat('WARNING','TTOASC','orbits/ttoasc',' ',message,0)
      NICS=maxorb
      endif

c**   Hardwire for NOAA comparison
      if( nics.eq.8 ) nics = 9
c
      READ(IUTIN) COMT,nsat,NINTRS,
     .   (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)
crwk  trap nsat=0 which will create a synatx error in the formts
      if( nsat.le.0 )   call report_stat('FATAL','TTOASC'
     .   ,'orbits/ttoasc',' ',' nsat=0 on t-file header ',0)
      write(IPRNT,30) ((COMT(i,j),i=1,80),j=1,3),NSAT,NINTRS
   30 format(' T-file comment :',/,3(80a1,/),
     1       ' nsat,nintrs :',2i5,/,
     2       ' Satellite name and ics :')
      if (nics .gt. 9) then
c       text='(__(/,1X,a6,11x,3d12.5/,18x,_d12.5))'
c       text='(__(/,1X,a6,11x,3f14.5/,18x,__f14.9))'
c       text='(__(/,1X,a6,11x,3f14.5/,18x,6f14.9/,18x,_f14.9))'
c       text='(__(/,1X,a16,1x,3f14.5/,18x,6f14.9/,18x,_f14.9))'
        text='(__(/,1X,a16,1x,3f14.5/,18x,6f14.9/,18x,__f14.9))'
        WRITE (BUF2,'(I2)') NSAT
        READ  (BUF2,'(A2)') TEXT(2:3)
        WRITE (BUF2,'(I2)') NICS-9
        READ  (BUF2,'(A2)') TEXT(41:42)
      else
        text='(__(/,1X,a16,1x,3f14.5/,18x,_f14.9/))'
        WRITE (BUF2,'(I2)') NSAT
        READ  (BUF2,'(A2)') TEXT(2:3)
        WRITE (BUF2,'(I2)') NICS-3
        READ  (BUF2,'(A2)') TEXT(29:29)
      endif
      write(iprnt,fmt=text)
     $ (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)

C     loop over epochs

      do iepoch = 1,nepchs

         if( iepoch.eq.1 ) then
c           convert start (GPST) epoch to GPST week number and seconds
            iys = itb(3)
            idns = idoy( iys,itb(1),itb(2) )
            ihs = tbb(1)
            mins = tbb(2)
            secs = tbb(3)
            call timcon( -4,iwkn,sow,iys,idns,ihs,mins,secs,utcoff )
         else
c           add delta seconds
            call secsum(iwkn,sow,sdelt,iwkn,sow)
         end if
c        write the epoch number and time to the print file
         call timcon(4,iwkn,sow,iys,idns,ihs,mins,secs,utcoff )
         READ(IUTIN,END=200) ((X(I,ISAT),I=1,NINTRS),ISAT=1,NSAT)
         WRITE(iprnt,110) iepoch,iys,idns,ihs,mins,secs,iwkn,sow
     .                  , ((X(I,ISAT),I=1,NINTRS),ISAT=1,NSAT)
110      FORMAT( /,1x,'Epoch >: ',i3,3x,i4,i4,2i3,f4.0,i6,f10.2
     .         , /,(3F20.9))
      enddo

200   continue
      close(unit=iutin)
      close(unit=iprnt)

      call report_stat('STATUS','TTOASC','orbits/ttoasc',tmpnam,
     .'ASCII dump on file : ',0)
      call report_stat('STATUS','TTOASC','orbits/ttoasc',' ',
     .'Normal end in TTOASC ',0)

      stop
      end
