      Subroutine xcsats( xfile, xorc, lxfil, nchan, ischan
     .                 , itb, tbb, itstp, tstp )


c       Read the satellites and start, stop times from the X- or C-file

c       R. King  920911 from S. Shimada SETMRD (removed)

C     arguments
C        xfile  : X- or C-file name                              input
C        xorc   : First character of input X-,C-file             input
C        lxfil  : unit number of X- or C-file                    input
C        lcfil  : unit number of C-file                          input
c        nchan  : number of satellites                           output
c        ischan : array of PRN numbers                           output
C        itb    : begin time ( month,day,year)                   output
C        tbb    : begin time (hour,minute,second)                output
C        itstp  : end time   (month,day,year)                    output
C        tstp   : end time   (hour,minute,second)                output

      implicit  none

      include '../includes/dimpar.h'

      character* 1  xorc,latflag,lonflag,gnss
      character* 3  rcvrsw,rxobtyp(maxdat)
      character*6   antcod
      character*16  xfile,sitnam16,satnam(maxsat)
      character*20   rctype,rcvnum,anttyp,antnum
      character*32  sitnam
      character*80  text(maxtxt)
C
      integer*4     iscrn,iprnt,ntext,ndat,lambda(maxsat,maxdat)
     .            , itb(3),itstp(3),ischan(maxsat)
     .            , dattyp(maxdat),lxfil,nslip,inter,ircint,isessn
     .            , latd,lond,latm,lonm,idoy,jdoy
     .            , iyr,jyr,imon,jmon,iday,jday,ihr,jhr,imin,jmin
     .            , isec,jsec,nchan,nepoch,mtime,ioerr

      integer*2     islip(maxcsb),islpst(maxcsb)

      real*4      swver

      real*8        tbb(3),tstp(3),offarp(3),offsl1(3),offsl2(3)
     .            , seclat,seclon,height,sec

      parameter (iscrn=6,iprnt=8)

c     Open a scratch file for the header print information

         open( iprnt, status='scratch' )


c     Read the header information from the X- or C-file

      if( xorc .eq. 'x' )  then
C
         open ( lxfil, file=xfile, status='old', err=90 )
         call xhdred ( lxfil,iprnt,iscrn,nepoch,inter,ircint,isessn
     .               , mtime,iyr,imon,iday,ihr,imin,sec
     .               , nchan,ischan,satnam
     .               , ndat,dattyp,rxobtyp,lambda
     .               , offarp,sitnam16,rcvrsw,swver,antcod
     .               , rctype,rcvnum,anttyp,antnum
     .               , latflag,latd,latm,seclat
     .               , lonflag,lond,lonm,seclon,height
     .               , ntext,text,gnss )
         close ( lxfil )
c        print *,'nepoch,inter,ircint,isessn',nepoch,inter,ircint,isessn
c        print *,'mtime,iyr,imon,iday,ihr,imin,sec'
c    .          , mtime,iyr,imon,iday,ihr,imin,sec
c        print *,'nchan,ischan,ndat,dattyp,lambda'
c    .          , nchan,ischan,ndat,dattyp,lambda
c        print *,'offsl1,offsl2,sitnam16,rcvrsw,swver'
c    .          , offsl1,offsl2,sitnam16,rcvrsw,swver
c        print *,'rctype,rcvnum,anttyp,antnum '
c    .            rctype,rcvnum,anttyp,antnum
c        print *,'latflag,latd,latm,seclat',latflag,latd,latm,seclat
c        print *,'lonflag,lond,lonm,seclon',lonflag,lond,lonm,seclon
c        print *,'ntext,text,gnss'
c    .          , ntext,text,,gnss

      elseif ( xorc.eq.'c' ) then

         call copens ( xfile, 'old', lxfil, ioerr )
         call chdred ( lxfil, iprnt, nepoch, inter, mtime
     .               , iyr, imon, iday, ihr, imin, sec
     .               , nchan, ischan, ndat, dattyp, lambda
     .               , offarp,offsl1, offsl2, sitnam, rcvrsw, swver
     .               , rctype,rcvnum,anttyp,antnum
     .               , ircint,isessn,ntext, text
     .               , nslip, islip, islpst )
c        print *,'nepoch,inter,mtime',nepoch,inter,mtime
c        print *,'iyr,imon,iday,ihr,imin,tbb',iyr,imon,iday,ihr,imin,tbb
c        print *,'nchan,ischan,ndat,dattyp,lambda'
c    .          , nchan,ischan,ndat,dattyp,lambda
c        print *,'offarp,offsl1,offsl2,sitnam,rcvrsw,swver'
c    .          , offarp,offsl1,offsl2,sitnam,rcvrsw,swver
c        print *,'rctype,rcvnum,anttyp,antnum '
c    .            rctype,rcvnum,anttyp,antnum
c        print *,'ircint,isessn,ntext,text'
c    .          , ircint,isessn,ntext,text
c        print *,'nslip,islip,islpst',nslip,islip,islpst
         close ( lxfil )

      else
         call report_stat('FATAL','FIXDRV','xcsats',' '
     .          ,'No input X- or C-file',0)

      endif

      close (iprnt)


C     Calculate the stop time from the start time, interval, and number of epochs

c     FIXDRV, MODEL, and ARC now use GPST - may need to convert X- or C-file values
      if( mtime.eq.1 ) then
         jdoy = idoy(iyr,imon,iday)
         call utc2gps(jdoy,iyr,ihr,imin,sec )
         call monday(jdoy,imon,iday,iyr)
      endif
      isec = sec + 0.01D0
      call nextt1( iyr, imon, iday, ihr, imin, isec, inter*(nepoch-1)
     .           , jyr, jmon, jday, jhr, jmin, jsec )


c     Fill the start and stop arrays for BMAKE

      itb(3)   = iyr
      itb(1)   = imon
      itb(2)   = iday
      tbb(1)   = ihr
      tbb(2)   = imin
      tbb(3)   = sec
      itstp(3) = jyr
      itstp(1) = jmon
      itstp(2) = jday
      tstp(1)  = jhr
      tstp(2)  = jmin
      tstp(3)  = jsec + mod(tbb(3),1.0D0)

      return

  90  call report_stat('FATAL','FIXDRV','xcsats',xfile
     .  ,              'Error opening X-file: ',0)

      end



