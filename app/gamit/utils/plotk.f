      PROGRAM PLOTK
C
C     K. Feigl

C     Plot a K-File

      implicit none

      include '../includes/dimpar.h'

      CHARACTER*1 latflag,lonflag,nojumpfit,jumpfit,clkft,gnss
      CHARACTER*3 rcvrsw,rxobtyp(maxdat)
      CHARACTER*4 SNAME
      character*6 antcod
      CHARACTER*16 sitnam,satnam(maxsat)
      CHARACTER*16 XFILE,KFILE
      character*20  rctype,rcvnum,anttyp,antnum
      character*80 TEXT(MAXTXT),XFILES(10),list

      logical fcheck

      INTEGER DATTYP(MAXDAT),LAMBDA(MAXSAT,MAXDAT),NTEXT
     .      , iobs,ircint,isessn,iy,im,id,ihr,min,iclk
     .      , latd,latm,lond,lonm,iscrn,ioerr,inter,idoy,idyoyr
     .      , isprn,nsat,it0,mtime,ndat,nepoch,iprnt,n

      real*4 swver

      REAL*8  t00,offarp,height,sec,seclat,seclon
     .     ,  epoch,rate,accel,cubic

      DIMENSION offarp(3)
      DIMENSION ISPRN(maxsat)
      DIMENSION IT0(3),T00(3)

      iscrn = 6

      write (iscrn,*) ' Enter name of K-file '
      read (*,'(a)') kfile

      xfile = kfile
      write(xfile,'(a16)') kfile
      xfile(1:1) = 'x'
      call lowers(xfile)

c     cannot find X-file to match the K-file name.
c     Go look for one in the directory
      if (.not. fcheck(xfile)) then
         list = xfile
         list(6:6) = '?'
         call getdir (list,10,xfiles,n)
         if (n .eq. 0) then
           write(iscrn,'(a,a16,a)') '**Cannot find X-file ('
     .        ,xfile,') to match K-file.  Stop in PLOTK'
           stop
         endif
         xfile = xfiles(1)(1:16)
      endif


      write (iscrn,*) ' Using X-file: ',xfile

c     read an X-file to get start and stop time
      iobs = 14
      OPEN (UNIT=iobs,FILE=XFILE,FORM='FORMATTED',
     .iostat=ioerr,ERR=70,STATUS='OLD')
      iprnt = 8
      OPEN (UNIT= iprnt,FORM='FORMATTED',iostat=ioerr,
     .ERR=70,STATUS='SCRATCH')

      CALL XHDRED ( iobs,iprnt,iscrn
     .            , NEPOCH,INTER,ircint,isessn
     .            , MTIME,IY,IM,ID,IHR,MIN,SEC
     .            , NSAT,ISPRN,satnam
     .            , NDAT,DATTYP,rxobtyp,LAMBDA
     .            , offarp,SITNAM,rcvrsw,swver,antcod
     .            , rctype,rcvnum,anttyp,antnum
     .            , LATFLAG,LATD,LATM,SECLAT
     .            , LONFLAG,LOND,LONM,SECLON,HEIGHT
     .            , NTEXT,TEXT,gnss )

c     if X times are UTC, convert to GPST (expected by CLKERA)
      if( mtime.eq.1 ) then
        idyoyr = idoy( iy,im,id )
        call utc2gps( idyoyr,iy,ihr,min,sec )
        call monday( idyoyr,im,id )
      elseif( mtime.ne.2 ) then
        call suicid('IMAKEF: mtime neither 1 (UTC) nor 2 (GPST)')
      endif
c     put the times into the arrays
      it0(1) = im
      it0(2) = id
      it0(3) = iy
      t00(1) = ihr
      t00(2) = min
      t00(3) = sec
      sname = xfile(2:5)
      call uppers(sname)
      write(*,'(a)') 'Enter polynomial order for no-jump fit (L, Q, C)'
      read(*,'(a)') nojumpfit
      write(*,'(a)') 'Enter polynomial order for with jump fit'
      read(*,'(a)') jumpfit
      clkft = 'I'
c     no print file from PLOTK (6th arg = 0)
      CALL CLKERA( KFILE,SNAME,nojumpfit,jumpfit,clkft,IT0,T00,0
     .           , ICLK,EPOCH,RATE,ACCEL,CUBIC)

      stop

  70  continue
      print *,'Open error on ',xfile
      if (ioerr .ne. 0) then
         call ferror (ioerr,6)
      endif

      STOP 1
      END
