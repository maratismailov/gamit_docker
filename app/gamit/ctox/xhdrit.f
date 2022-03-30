Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.
C
      Subroutine XHDRIT ( iux,cfname,shift )

C     Write the X-File header records

      implicit none

      include '../includes/dimpar.h'
      include '../includes/global.h'
      include '../includes/model.h'

      CHARACTER*1 LATFLAG,LONFLAG,UPPERC,yawbias
      CHARACTER*2 BUF2
      CHARACTER*4 BUF4,x_time_flag
      character*5 radome
      character*9 exptyp
      CHARACTER*16 CFNAME,satnam
      character*20 antbody
      character*22 exfmt
      CHARACTER*80 HEAD
C
      INTEGER INSN(MAXSAT),iyr,idoy1,ihour,min
     2      , IHNSEC,IYEAR,IMONTH,IDAY,IHR,IMN,ISEC
     3      , LATD,LATM,LOND,LONM,IUX,iy2
     4      , IWKN,itflag,nsn,frqchn,svnstart(5),svnstop(5),I,J

      REAL*8 OFFSL1(3),OFFSL2(3),DOUT(3),TMPLAT,TMPLON,TMPRAD
     1     , SEC,SECLAT,SECLON,UTCOFF,sbmass,yawrate
     2     , sod,shift,sow
      real*8 antpwr   ! Tranmission power (W)

C
C        Get the run Time and Date and Write the X-File Information Records
C
      head=' '
      HEAD(1:40)='CTOX    :  X-File written from C-File : '
      HEAD(40:55)=CFNAME(:)
      CALL GETDAT(IYEAR,IMONTH,IDAY)
      CALL GETTIM(IHR,IMN,ISEC,IHNSEC)
      WRITE (BUF2,'(I2)') IDAY
      READ  (BUF2,'(A2)') HEAD(56:57)
      HEAD(58:58)='-'
      WRITE (BUF2,'(I2)') IMONTH
      READ  (BUF2,'(A2)') HEAD(59:60)
      HEAD(61:61)='-'
      WRITE (BUF4,'(I4)') IYEAR
      READ  (BUF4,'(A4)') HEAD(62:65)
      WRITE (BUF2,'(I2)') IHR
      READ  (BUF2,'(A2)') HEAD(68:69)
      HEAD(70:70)=':'
      WRITE (BUF2,'(I2)') IMN
      READ  (BUF2,'(A2)') HEAD(71:72)
      HEAD(73:73)=':'
      WRITE (BUF2,'(I2)') ISEC
      READ  (BUF2,'(A2)') HEAD(74:75)
      HEAD(76:80)='   '
      DO 10 I=1,NTEXT
c     replace any old-code nulls with blanks in C-file header
      call sub_char(text(i),'\0',' ')
      WRITE(IUX,15) TEXT(I)
10    CONTINUE
      WRITE(IUX,15) HEAD
15    FORMAT(A80)
      WRITE(IUX,20)
20    FORMAT(/,'END')
C
C
C        Write the Station Description information
C            
      CALL RADDEG( latr_sph,DOUT )
      LATD= DOUT(1)
      LATM= DOUT(2)
      SECLAT= DOUT(3)
      CALL RADDEG( lonr,DOUT )
      LOND= DOUT(1)
      LONM= DOUT(2)
      SECLON= DOUT(3)


C     Use 'N'/'S' and 'E'/'W' instead of '+'/'-'
      if(latd.lt.0.or.latm.lt.0.or.seclat.lt.0.d0) then
         latflag=upperc('S')
         latd=iabs(latd)
         latm=iabs(latm)
         seclat=dabs(seclat)
      else
         latflag=upperc('N')
      end if
      if(lond.lt.0.or.lonm.lt.0.or.seclon.lt.0.d0) then
         lonflag=upperc('W')
         lond=iabs(lond)
         lonm=iabs(lonm)
         seclon=dabs(seclon)
      else
         lonflag=upperc('E')
      endif
      write(iux,'(20x,a)') 'SPHERICAL COORDINATES' 
      WRITE(IUX,31)
31    FORMAT(1X,'STATION_NAME___      LATITUDE         LONGITUDE      RA
     .DIUS     RVCR SWVER RCVCOD'
     .,/,18X,'_DG MN SS.SSSSS __DG MN SS.SSSSS ____(M).___')
      call read_rcvant(2,2,antcod,anttyp,radome,rcvcod,rctype,pcncod)
c     c-file does not carry the 3-character MAKEX code
      rcvrswx = ' '
cd      print *,'DEBUG XHDRIT rcvtype rcvcod rcvrswx,swver ircint '
cd     .     ,rctype,rcvcod,rcvrswx,swver,ircint
      WRITE(IUX,32)
     . SITNAM,LATFLAG,LATD,LATM,SECLAT,LONFLAG,LOND,LONM,SECLON,radius
     . ,rcvrswx,swver,rcvcod
32    FORMAT
     . (1X,A12,5X,A1,I2,1X,I2,1X,F8.5,1X,A1,I3,1X,I2,1X,F8.5,F13.4
     .  ,2x,a3,1x,f5.2,1x,a6) 
      write(iux,33)
33    format(' ANT ARP OFFSETS (M)   UP     NORTH   EAST ')   
cd      print *,'XHDRIT anttyp ',anttyp
      call read_rcvant(2,1,antcod,anttyp,radome,rcvcod,rctype,pcncod)           
      write(iux,34) antcod,offarp
34    format(1x,a6,12x,3f8.4,/) 

C     Write the Satellite Information Lines
      call dayjul( jd0,iyr,idoy1 )
      DO I=1,nchan
* MOD TAH 190702: Added place holder for antpwr to snav_read call
         call svnav_read( -1,iyr,idoy1,ihour,min,gnss,ischan(i),
     .          insn(i),frqchn,antbody,sbmass,yawbias,yawrate, antpwr,
     .          svnstart,svnstop )
      enddo
      WRITE(IUX,41) nchan,NDAT,(DATTYP(J),J=1,NDAT)
41    FORMAT(/,1X,I2,' SATELLITES',12X,I2,' DATA TYPES:',7I3)
      do i=1,nchan  
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read(-1,iyr,idoy1,ihour,min,gnss,ischan(i),nsn,
     .           frqchn,antbody,sbmass,yawbias,yawrate, antpwr,
     .           svnstart,svnstop ) 
        write(satnam,'(a1,i2,1x,i3,1x,a8)') 
     .           gnss,ischan(i),nsn,antbody(7:14)
         write(iux,'(a,i2,a,a16,a,7i3)')  ' CHANNEL ',i,'  SV  '
     .         ,satnam,'LAMBDA ',(lambda(i,j),j=1,ndat)
c old         WRITE(IUX,46) I,ISCHAN(I),gnss_name,INSN(I)
c     .     ,(LAMBDA(I,J),J=1,NDAT)
c 46    FORMAT(' CHANNEL ',I2,'  PRN# ',I2,2X,a8,I2,2X,'LAMBDA'
c     1      ,  7I3 )
      enddo

C     Write start time and data interval 

      t0 = t0 + shift 
      call timinc(idoy1,t0,0.0d0)
      call ds2hms(iyr,idoy1,t0,ihour,min,sec)
c     get GPST-UTC for the X-file header (output iwk, sow not used)
      itflag = -2
      call timcon(itflag,iwkn,sow,iyr,idoy1,ihour,min,sec,utcoff)
c     internal and X-file time now GPST
      x_time_flag = 'GPST'
      iy2 = mod(iyr,100)
      write(iux,51) x_time_flag,utcoff,IY2,IDOY1,IHOUR,MIN,SEC,INTER
     .            , ircint,isessn
   51    format(/,' FIRST EPOCH (',a4,')    GPST-UTC= ',f4.1,/
     . ,' YR DAY HR MN SECS   INTERVAL(SECS)   DATA INTERVAL  SESSION'
     . ,/,1X,I2,1X,I3,1X,I2,1X,I2,1X,F6.3,1X,I6,11x,i6,11x,i2,/)


C     Write Number of Epochs
      WRITE(IUX,60) NEPOCH
60    FORMAT(I4,' EPOCHS')


C     Old Band and Frequency Descriptors are now blank
      WRITE(IUX,70)
70    FORMAT(21X)

c     Write the Experiment Descriptor  
c    RWK 160420: This no used since kinematic not supported, but could be a 
c                format descriptor if needed
c      if( upperc(skd).eq.' ' ) then
      exfmt = '                      '
      exptyp= '         '
c      elseif( upperc(skd).eq.'S' ) then
c        exfmt = ' EXPANDED DATA FORMAT '
c        exptyp = 'STATIC   '
c      elseif( upperc(skd).ne.'S' ) then
c        exfmt = ' EXPANDED DATA FORMAT '
c        if( upperc(skd).eq.'K' ) exptyp='KINEMATIC'
c        if( upperc(skd).eq.'D' ) exptyp='DYNAMIC  '
c      endif
C     Write the Data Header Line
      WRITE(IUX,80) exfmt,exptyp
80    FORMAT(/,a22,a9,/,' EPOCH # FLG CH     L1 PHASE',11X,'AMP'
c     jfg end change
     1      ,7X,'L2 PHASE',11X,'AMP'
     2      ,4X,'L1 PSEUDORANGE',11X,'L2 PSEUDORANGE')
C
      RETURN
      END
