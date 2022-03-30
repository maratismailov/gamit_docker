      Program TMERGE
c
c     program to take two identical length tfiles and merge the postions
c     from the first with the partials of the second, the header of the
c     second file is kept.
c
      implicit none

      include '../includes/dimpar.h'

      character*7 status                         
      character*1 gnss
      character*4 icsnam(maxorb) 
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*16 tname(2),pname,satnam1(maxsat),satnam2(maxsat)
        character*120 version
      character*256 message

      integer*4 iutin1,iutin2,iscrn,iprnt,nsat(2),itsat1(maxsat),
     .          itsat2(maxsat),jdb(2),jds(2),nepchs(2),jde(2),nics(2),
     .          nintrs(2),iutout,isat,i,maxtrm,iterm,ioerr

      parameter ( maxtrm = (maxorb+3)*3 )

      real*8 tb(2),ts(2),sdelt(2),te(2)
     .      , tsatic1(maxorb,maxsat),tsatic2(maxorb,maxsat)
     .      , x(2,maxtrm,maxsat),y(maxtrm,maxsat)

      iterm=5
      iscrn=6
      iprnt=6
      iutin1=8
      iutin2=9
      iutout=10   
      frame = 'UNKWN'

c     Print the version and machine

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5) version
    5 format('Started TMERGE ',a120)
      call report_stat('STATUS','TMERGE','orbits/tmerge',' ',message,0)
c
c Open the Input T-files and output T-file
c
      write(iscrn,50)
50    format('Program to take two identical length T-Files and merge',/,
     .       'the postions from the first with the partials of the',/,
     .       'second. The header of the second T-file is used.',/)

      WRITE(iscrn,'(A)') 'Name of input T-File without partials (A) :> '
      write(iscrn,*) 'X,Y,Z from this T-File copied to Partials T-File '
      READ(iterm,101) TNAME(1)
  101 FORMAT(A16)
      status = "old    "
      call topens(tname(1),status,iutin1,ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('STATUS','TMERGE','orbits/tmerge',tname(1),
     .  'Error opening input T-file: ',ioerr)
      endif
c
      WRITE(iscrn,'(A)') 'Name of input T-File with partials :> '
      write(iscrn,*) 'X,Y,Z in this T-File replaced by X,Y,Z from T-File
     . (A)'
      READ(iterm,102) TNAME(2)
  102 FORMAT(A16)
      status = "old    "
      call topens(tname(2),status,iutin2,ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('STATUS','TMERGE','orbits/tmerge',tname(2),
     .  'Error opening input T-file: ',ioerr)
      endif
c
      WRITE(iscrn,'(A)') ' Output T-File Name :> '
      READ(iterm,101) PNAME
      status = "unknown"
      call topens(pname,status,iutout,ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('STATUS','TMERGE','orbits/tmerge',pname,
     .  'Error opening output T-file: ',ioerr)
      endif
c
c Read the header of the first T-File
c
      call thdred(iutin1,iscrn,iprnt,nsat(1),gnss,itsat1,satnam1
     .      , jdb(1),tb(1),jds(1),ts(1),sdelt(1),nepchs(1),jde(1),te(1)
     .      , nics(1),tsatic1,nintrs(1),icsnam
     .      , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .      , antradmod )
c
c Read the header of the second T-File
c
      call thdred(iutin2,iscrn,iprnt,nsat(2),gnss,itsat2,satnam2
     .      , jdb(2),tb(2),jds(2),ts(2),sdelt(2),nepchs(2),jde(2),te(2)
     .      , nics(2),tsatic2,nintrs(2),icsnam
     .      , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .      , antradmod )
c
c Write the header of the output T-File
c
      call thdrit( iutout,jde(2),te(2),jdb(2),tb(2),jds(2),ts(2)
     .           , sdelt(2),nepchs(2),nintrs(2)
     .    , nsat(2),gnss,itsat2,satnam2,tsatic2,nics(2),tname(2),icsnam
     .    , precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod )
c
c loop over epochs and write out T-File data records.
c
100   READ(IUTIN1,END=200) ((X(1,I,ISAT),I=1,NINTRS(1)),ISAT=1,NSAT(1))
      READ(IUTIN2,END=200) ((X(2,I,ISAT),I=1,NINTRS(2)),ISAT=1,NSAT(2))
      do ISAT=1,NSAT(2)
        do I=1,NINTRS(2)
          if(i.le.3) then
c       print*,' coordinate t1 = ',x(1,i,isat)
c       print*,'            t2 = ',x(2,i,isat)
            Y(I,ISAT) = X(1,I,ISAT)
          else
            Y(I,ISAT) = X(2,I,ISAT)
          endif
        enddo
      enddo
      WRITE(IUTOUT) ((Y(I,ISAT),I=1,NINTRS(2)),ISAT=1,NSAT(2))
cd    WRITE(iscrn,*)  ((Y(I,ISAT),I=1,NINTRS(2)),ISAT=1,NSAT(2))
      GOTO 100
200   continue
c
      close(unit=iutin1)
      close(unit=iutin2)
      close(unit=IUTOUT)
c
      call report_stat('STATUS','TMERGE','orbits/tmerge',' ',
     .'Normal end in TMERGE: ',0)
c
      stop
      end
