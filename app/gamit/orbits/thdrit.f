Copyright (c) Massachusetts Institute of Technology and the Unviersity of
California at San Diego, 1994.  All rights reserved.

      Subroutine thdrit( iut,jde,te,jds,ts,jdf,tf,delt,nepoch,nintrs
     .         , nsat,gnss,itsat,satnam,satics,nics,source_file,icsnam
     .         , precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod )
C
c     R.W. King   March 1988   from R.Abbot WRTHED (Nov 84) as modified
C                              by Y. Bock for BCTOT (Jan 88)
c   mods:
c   Mar 1990 by pch - Apollo doesn't like long sequential unformatted writes
c                     so I broke the tfile (unit 16) write into 2 pieces
c   July 1990 by ss - delete FIRST from argument of THDRIT and NSNPRN
c   July 1994 by rwk - modify for joint use by NGSTOT and TTONGS; generalize
c                      source file (g-, t-, or NGS), add efixed to argument list
c
c   Apr 1995 by PT - also write out the inertial frame precession,nutation and gravity models
c                    used in the tfile (precmod,nutmod,gravmod)
c   Mar 1997 by PT  - remove efixed from argument list and replace with use of "frame"
c   Nov 2015 by rwk - Mods for GNSS

      implicit none

      include '../includes/dimpar.h'

      character*1 gnss 
      CHARACTER*2 BUF2
      CHARACTER*4 BUF4,BLANK,icsnam(maxorb)
      CHARACTER*16 SNAME,SATNAM(MAXSAT),source_file
      CHARACTER*80 HEAD,comt(3)
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod

      INTEGER*4 iut,JDE,JDS,JDF,idoy
     .        , itsat,nics,isat,nsat,nintrs,nepoch
     .        , iyear,imonth,iday,ihr,imn
     .        , iy0,im0,id0,ih0,imin0,ihnsec,isec
     .        , iye,iyf,ime,imf,ide,idf,ihe,ihf,imine,iminf
     .        , i,k


      real*8 xpole,ypole,ut1s,delt,satics
     .     , h0,he,hf,te,tf,ts,sec0,sece,secf,xmin0,xmine,xminf

c     for debug
c      logical dump
c      parameter( dump=.false. )

      DIMENSION SATICS(MAXORB,MAXSAT),ITSAT(MAXSAT)
 
cd      print *,'THDRIT nsat nics ',nsat,nics
cd      do i=1,nsat
cd        print *,itsat(i),satnam(i),(satics(k,i),k=1,nics) 
cd      enddo

      BLANK='    '
      HEAD(1:40)= '                                        '
      HEAD(41:80)='                                        '
      do i=1,3
       comt(i)=' '
      enddo

c        Set time and reference frame - time now always GPST internally

      comt(1)(11:14) = 'GPST'
      if( frame.eq.'EFIXD' ) then
        comt(1)(21:40) = 'EARTH-FIXED     '
      else
        comt(1)(21:29) = 'INERTIAL'
        comt(1)(36:40) = frame
      endif

c  set the precession,nutation and gravity models
      comt(1)(42:46) = precmod
      comt(1)(48:52) = nutmod
      comt(1)(54:58) = gravmod
      comt(1)(60:64) = srpmod
      comt(1)(66:70) = eradmod
      comt(1)(72:76) = antradmod

c        Set the IC parameters in comment line 2
      k = 1
      do i = 1,nics
        comt(2)(k:k+3) = icsnam(i)
        k = k+5
      enddo

C        Set up T-File Header

      HEAD(1:39)=' Tabular ephemeris file generated from '
      HEAD(40:55)=source_file(:)
      CALL GETDAT(IYEAR,IMONTH,IDAY)
      CALL GETTIM(IHR,IMN,ISEC,IHNSEC)
      WRITE (BUF2,'(I2)') IDAY
      READ  (BUF2,'(A2)') HEAD(57:58)
      HEAD(59:59)='-'
      WRITE (BUF2,'(I2)') IMONTH
      READ  (BUF2,'(A2)') HEAD(60:61)
      HEAD(62:62)='-'
      WRITE (BUF4,'(I4)') IYEAR
      READ  (BUF4,'(A4)') HEAD(63:66)
      WRITE (BUF2,'(I2)') IHR
      READ  (BUF2,'(A2)') HEAD(68:69)
      HEAD(70:70)=':'
      WRITE (BUF2,'(I2)') IMN
      READ  (BUF2,'(A2)') HEAD(71:72)
      HEAD(73:73)=':'
      WRITE (BUF2,'(I2)') ISEC
      READ  (BUF2,'(A2)') HEAD(74:75)
      WRITE (BUF2,'(I2)') NICS
      READ  (BUF2,'(A2)') HEAD(79:80)
C
C
C        Convert times of IC epoch, start and stop times of
C        the integration to calender day, hr, min, sec
C
c      print *,'THDRIT jde jds jdf ',jde,jds,jdf
c*      CALL MDYJUL( IME,IDE,IYE,JDE )
      call dayjul( jde,iye,idoy )
      call monday( idoy,ime,ide,iye )
      CALL ds2hms( iye,ide,te,IHE,IMINE,SECE )
c      print *,'iye idoy ime ide ',iye,idoy,ime,ide
      HE= IHE
      XMINE= IMINE
c*      CALL MDYJUL( IM0,ID0,IY0,JDS )
      call dayjul( jds,iy0,idoy)
      call monday( idoy,im0,id0,iy0 )
      CALL ds2hms( iy0,id0,ts,IH0,IMIN0,SEC0 )
c      print *,'iy0 idoy im0 id0 ',iy0,idoy,im0,id0
      H0= IH0
      XMIN0= IMIN0
c*      CALL MDYJUL( IMF,IDF,IYF,JDF ) 
      call dayjul( jdf,iyf,idoy )
      call monday( idoy,imf,idf,iyf )
      CALL ds2hms( iyf,idf,tf,IHF,IMINF,SECF ) 
c      print *,'iyf idoy imf idf ',iyf,idoy,imf,idf
      HF= IHF
      XMINF= IMINF

C        Compute the number of epochs for the T-file header

      NEPOCH= ( (JDF-JDS)*86400.D0 + (TF-TS) ) / DELT + .000001D0 + 1


C       Write the T-file headers

      UT1S=0.D0
      XPOLE=0.D0
      YPOLE=0.D0
      WRITE(iut) HEAD,IME,IDE,IYE,HE,XMINE,SECE
     1              ,IM0,ID0,IY0,H0,XMIN0,SEC0
     2              ,IMF,IDF,IYF,HF,XMINF,SECF
     3              ,DELT,NEPOCH,UT1S,XPOLE,YPOLE

      write( iut) comt,nsat,nintrs
     .          , (satnam(isat),(satics(i,isat),i=1,nics),isat=1,nsat)


c      IF( DUMP ) THEN
c        WRITE( 6,50) head,iye,ime,ide,ihe,imine,sece
c     .              , iy0+1900,im0,id0,ih0,imin0,sec0
c     .              , iyf+1900,imf,idf,ihf,iminf,secf
c     .              , delt,nepoch,comt,nsat,nintrs,comt(2)
c   50   format(//,
c     .        '-------------------------------------------------------',
c     .        /,
c     .        ' THDRIT: T-File Header Information  :',1x,a80,/,
c     .        ' Epoch of initial conditions:',I4,4I3,F10.6,/,
c     .        ' Ephemeris start            :',I4,4I3,F10.6,/,
c     .        ' Ephemeris end              :',I4,4I3,F10.6,/,
c     .        ' Tabular interval (sec)     :',f10.3,/,
c     .        ' T-file epochs              :',i5,/,
c     .        ' Comments                   :',/,3(/,1x,a80),/,
c     .        ' Number of satellites       :',1x,i2,/,
c     .        ' NINTRS                     :',i3,/,
c     .        ' Orbital parameter types    :',1x,a80,//,
c     .        ' Satellite name and ics     :')
c        do isat = 1,nsat
c          write(6,60) satnam(isat),(satics(i,isat),i=1,nics)
c        enddo
c   60   format(1x,a16,1x,3f14.6,3f8.4,9f9.4)
c      endif
C

      RETURN
      END
