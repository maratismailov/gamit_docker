      Program TTOT
C
C Modification of TTOASC to convert integer*2 T-files to integer*4
C   N.B.  MUST be compiled as -i*2
C
C  Modified by PT 950620 to allow option to transform from inertial
c  to earth fixed and vice versa

      implicit none

      include '../includes/dimpar.h'

      character*16 Tin,tout
      character*1 upperc,ans,gnss
      CHARACTER*2 BUF2
      character*4 icsnam(maxorb)
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      CHARACTER*16 SATNAM(maxsat)
      CHARACTER*80 HEAD,comt(3)
      character*120 version
      CHARACTER*256 message

      integer*2 ite(3),itb(3),itstp(3),nepchs,nsat,nintrs

      INTEGER*4 ITE4(3),ITB4(3),ITSTP4(3),nepch4,nsat4,nintr4
     .        , maxtrm,iterm,iscrn,iutin,iutout,isat,i,sel
     .        , iprnt,itsat(maxsat),jdb,jds,nepcht,jde
     .        , nics,idir,ioerr

      parameter(maxtrm=(maxorb+3)*maxsat)

      real*8 tee(3),tbb(3),tstp(3),satics(maxorb,maxsat)
     .      , x(maxtrm,maxsat),tb,ts,sdelt,te,ut1,xpole,ypole
C
C Open the Input T-file and output T-file
C
      iterm=5
      iscrn=6
      iutin=8
      iutout=9
      iprnt = 0
      ioerr = 0
c
c          Print the version and machine

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5) version
    5 format('Started TTOT ',a120)
      call report_stat('STATUS','TTOT','orbits/ttot',' ',message,0)
c
c Begin
      write(iscrn,1000)
1000  format(//,10x,' TTOT - program to convert tfiles ... ',//
     .       ,5x,'1.   From I*2 to  I*4  ',/
     .       ,5x,'2.   Efixed   -> Inertial transformation',/
     .       ,5x,'3.   Inertial ->   Efixed transformation'/)

      write(*,'(a,$)')'  Enter selection : '
      read(5,*) sel

3     WRITE(6,'(a,$)') ' Input T-File Name :> '
      READ(5,101) tin
  101 FORMAT(A16)
      OPEN(UNIT=iutin,FILE=Tin,STATUS='OLD',FORM='UNFORMATTED'
     .     ,iostat= ioerr)
      if(ioerr .ne. 0) then
        call report_stat('WARNING','TTOT','orbits/ttot',tin,
     .'Error opening input T-file, try again:',ioerr)
        ioerr = 0
        goto 3
      endif

10    CONTINUE
      WRITE(6,'(a,$)') ' Output T-File Name :> '
      READ(5,101) tout
      OPEN(UNIT=IUTOUT,FILE=tout,STATUS='NEW',FORM='UNFORMATTED'
     .    ,iostat=ioerr)
      if(ioerr.ne.0)then
        write(iscrn,'(a,a,$)')tout,' exists. Overwrite ? (y): '
        read(iterm,'(a1)')ans
        if(ans.eq.'y'.or.ans.eq.'Y'.or.ans.eq.' ')then
          OPEN(UNIT=IUTOUT,FILE=tout,STATUS='unknown',FORM='UNFORMATTED'
     .       ,iostat=ioerr)
          if(ioerr .ne. 0) then
            call report_stat('FATAL','TTOT','orbits/ttot',tout,
     .      'Error opening output T-file: ',ioerr)
          endif
        else
          goto 10
        endif
      endif

c PT 950620: now branch program according to which option is required (ie sel = 1,2,3)

      if(sel.eq.1)then
C Read and Write the first header record of the T-file
c
        READ(IUTIN) HEAD,ITE,TEE,ITB,TBB,ITSTP,TSTP,SDELT,NEPCHS,
     1            UT1,XPOLE,YPOLE
        DO I = 1, 3
           ITE4(I) = ITE(I)
           ITB4(I) = ITB(I)
           ITSTP4(I) = ITSTP(I)
        ENDDO
        nepch4=nepchs
        WRITE(IUTOUT) HEAD,ITE4,TEE,ITB4,TBB,ITSTP4,TSTP,SDELT,
     1                nepch4,UT1,XPOLE,YPOLE
c
c Read the second header record of the T-file
c
        WRITE(BUF2,'(A2)') HEAD(79:80)
        READ (BUF2,'(I2)') NICS
        IF (NICS.le.0) then  
          write(message,'(a,i3)') 'Old T-file; number of ICs unknown: '
     .        ,nics
          call report_stat('FATAL','TTOT','orbits/ttot',' ',message,0)
        endif

        READ(IUTIN) COMT,NSAT,NINTRS,
     .     (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)
        WRITE(IUTOUT) COMT,nsat4,nintr4,
     .     (SATNAM(ISAT),(SATICS(I,ISAT),I=1,NICS),ISAT=1,NSAT)

C
C loop over epochs
C
100     READ(IUTIN,END=200) ((X(I,ISAT),I=1,NINTRS),ISAT=1,NSAT)
        WRITE(IUTOUT) ((X(I,ISAT),I=1,NINTRS),ISAT=1,NSAT)
        GOTO 100
200     continue
        close(unit=iutin)
        close(unit=IUTOUT)
        call report_stat('STATUS','TTOT','orbits/ttot',tout,
     .  'Succesfully converted T-file from I*2 to I*4. Output T-file: '
     .  ,0)
C
      elseif(sel.eq.2.or.sel.eq.3)then

c  transform efixed <-> inertial
c  determine direction to transform
        if(sel.eq.2)then
          idir = 1
71        write(iscrn,'(a,$)')' B1950 or J2000 inertial frame? : '
          read(*,'(a5)')frame
          frame(1:1) = upperc(frame(1:1))

          if(frame.ne.'B1950'.and.frame.ne.'J2000')then
            write(*,'(a,1x,a5)')' Invalid inertial frame: ',frame
            goto 71
          endif
        endif
        if(sel.eq.3)then
          idir = -1
          frame = 'INERT'
c  read tfile header to determine ref frame of input tfile
          call thdred( iutin,iscrn,iprnt,nsat4,gnss,itsat,satnam
     .               , jdb,tb,jds,ts,sdelt,nepcht,jde,te
     .               , nics,satics,nintr4,icsnam
     .               , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .               , antradmod )
        endif


c  close the tfiles
        close(iutin)
        close(unit=IUTOUT)

c transform input tfile into output tfile
        call trot ( tin,tout,idir,frame )
      else    
        write(message,'(a,i3)') 'Option not supported, sel = ',sel
        call report_stat('FATAL','TTOT','orbits/ttot',' ', message,0)
      endif

      call report_stat('STATUS','TTOT','orbits/ttot',' ',
     .'Normal end in TTOT. ',0)
      stop

      end
