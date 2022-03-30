Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.
C
      Subroutine XHDRED ( IOBS,IPRNT,ISCRN
     1,                   NEPOCH,INTER,ircint,isessn
     2,                   MTIME,IY,IM,ID,IHR,MIN,SEC
     3,                   NCHAN,ISCHAN,satnam
     4,                   NDAT,DATTYP,rxobtyp,LAMBDA
     5,                   offarp,sitnam16,rcvrsw,swver,antcod
     6,                   rctype,rcvnum,anttyp,antnum
     7,                   LATFLAG,LATD,LATM,SECLAT,
     8                    LONFLAG,LOND,LONM,SECLON,HEIGHT
     9,                   NTEXT,TEXT,gnss )    


C     Read an X-File Header
C     R. King  April 1987

      implicit none

      include '../includes/dimpar.h'

      INTEGER*4 DATTYP(MAXDAT),LAMBDA(MAXSAT,MAXDAT),NTEXT
     .        , iscrn,inter,npart,iprnt,id,im,ihr,min,ifac,iy,latd,ndat
     .        , iobs,lond,latm,lonm,ischan,idyoyr,nchan,mtime,nepoch
     .        , ircint,isessn,nblen,i,j

      REAL*8  offarp,height,seclat,seclon,sec

      CHARACTER*1 LATFLAG,LONFLAG,UPPERC,skd,gnss 
      CHARACTER*3 rcvrsw,TXTEND,BUF3,x_time_flag,rxobtyp(maxdat)
      character*6 antcod
      CHARACTER KINHED*4
      character*8 anthead  
      CHARACTER*12 IDUAL,FRQFAC
      character*16 sitnam16,satnam(maxsat)
      character*20 rctype,rcvnum,anttyp,antnum,temp20,temp20x
      CHARACTER*80 SCARR,TEXT(MAXTXT)
      
      real*4 swver

c     variables for status reporting
      integer*4 len,rcpar,ioerr
      character*6 source_type
      character*16 xfname
      character*80 prog_name
      character*256 message


      DIMENSION offarp(3),ISCHAN(MAXSAT)

      DATA TXTEND/'END'/
    


c     get calling program name and X-file name for report_stat
      len = rcpar(0,prog_name)
      inquire( unit=iobs, name=xfname, iostat=ioerr )
      if( ioerr.ne.0 )
     .  call report_stat('FATAL',prog_name,'lib/xhdred',xfname
     .           ,'Cannot find X-file name for error reporting',ioerr)


C     Case Insensitive Variables
      call lowers(txtend)
                
c     Dummy to avoid compiler warning for unused 'iscrn'
      if( iscrn.ne.0 ) then
        continue
      endif


C     Read and Print the Header, which is terminated by END
      if( iprnt.gt.0 ) then
        WRITE(IPRNT,20)
   20   FORMAT('----------------------------------------------------',/,
     1       1x,'** X-FILE HEADER INFORMATION **',/)
      endif
      NTEXT = 0     
      rctype = ' '
      rcvnum = ' '            
      anttyp = ' '
      antnum = ' '
      do i = 1,maxsat
        ischan(i)=0
      enddo
      source_type = 'unknwn'
   50 READ(IOBS,'(A80)',iostat=ioerr) SCARR   
      if( ioerr.eq.-1) call report_stat('WARNING',prog_name,'lib/xhdred'
     .  ,xfname,'Unexpected end reading X-file header',ioerr)
c     extract the RINEX receiver and antenna information from the header comments  
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading header comments',ioerr)
      if( scarr(14:24).eq.'SOURCE FICA' ) source_type = 'fica  '
      if( scarr(14:25).eq.'SOURCE RINEX' ) source_type = 'rinex ' 
      if( source_type.eq.'fica  ') then
        if( scarr(17:31).eq.'receiver serial') then  
          temp20 = ' '
          read(scarr,'(37x,a20)',iostat=ioerr) temp20 
          rcvnum = temp20(1:nblen(temp20)) 
        endif
      elseif( source_type.eq.'rinex ') then
        if( scarr(2:11).eq.'Receiver  ') then
          temp20 = ' '
          read(scarr,'(20x,a20)',iostat=ioerr) temp20
          rctype = temp20(1:nblen(temp20)) 
        endif
        if( scarr(2:11).eq.'Receiver s') then
          temp20 = ' '
          read(scarr,'(20x,a20)',iostat=ioerr) temp20 
          rcvnum = temp20(1:nblen(temp20)) 
        endif
        if( scarr(2:11).eq.'Antenna   ') then      
          temp20 = ' ' 
          temp20x = ' '
          read(scarr,'(20x,a20,11x,a20)',iostat=ioerr) temp20,temp20x  
          anttyp = temp20(1:nblen(temp20)) 
          antnum = temp20x(1:nblen(temp20x)) 
        endif
        if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading rcvr/ant info',ioerr)
      endif
      if( iprnt.gt.0 ) WRITE(IPRNT,'(A80)') SCARR
      BUF3 =SCARR(1:3)
      call lowers(buf3)
      IF (BUF3.EQ.TXTEND.and.nblen(scarr).eq.3) GO TO 65
      IF( NTEXT+1.LE.MAXTXT) THEN
         NTEXT = NTEXT + 1
         TEXT(NTEXT) = SCARR
      ELSE  
         write(message,'(a,i4,a,i4,a)') 'ntext (',ntext,') > maxtxt ('
     .        ,maxtxt,')--ignorning remaining text'
         call report_stat('WARNING',prog_name,'lib/xhdred',xfname
     .        ,message,0)  
         goto 65
      ENDIF
      GO TO 50


C     READ SITE NAME, COORDINATES, and Receiver, software id
C     New coordinate convention - Yehuda Bock 4/25/90
   65 rcvrsw='   '
      swver= 0.0
      READ(IOBS,70,iostat=ioerr) sitnam16
     .     ,LATFLAG,LATD,LATM,SECLAT,LONFLAG,LOND,LONM,SECLON,HEIGHT
     . , rcvrsw,swver
   70 FORMAT
     .  (///,1X,A16,1X,A1,I2,1X,I2,1X,F8.5,1X,A1,I3,1X,I2,1X,F8.5,F13.4
     .      ,2x,a3,1x,f5.2) 
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading coordinates',ioerr)

C     old conventions
      IF(LATFLAG.EQ.' ') LATFLAG='N'
      IF(LONFLAG.EQ.' ') LONFLAG='W'
C     If longitude is negative assume that
C     it refers to a left-handed coordinate system (earlier GAMIT asumption)
      if(lonflag.eq.'-' .or. lond.lt.0) lonflag='W'
      if(latflag.eq.'-' .or. latd.lt.0) latflag='S'

      LATFLAG=UPPERC(LATFLAG)
      LONFLAG=UPPERC(LONFLAG)

C     READ ANTENNA OFFSETS: Changed rwk 050928
c       old-style: Up, North, and East , L1 AND L2  
c       new_sytle: Up, North, and East, ARP  
c     skip blank line
      read(iobs,'(1x)') 
c     read the antenna header
      read(iobs,'(a8)',iostat=ioerr) anthead
      if( anthead.eq.' ANT ARP' ) then
c       new-style, read 2 values                    
        read(iobs,'(1x,a6,12x,3f8.4)',iostat=ioerr) antcod,offarp
c*        print *,'DEBUG XHDRED new style'
      else    
c       old-style: cannot reliably convert to ARP offsets since don't know PCV model   
        call report_stat( 'WARNING',prog_name,'lib/xhdred',xfname
     .   , 'Cannot convert x-file L1 L2 ant offsets to ARP, set = 0.',0)
        read(iobs,'(1x,a6)',iostat=ioerr) antcod
        do i=1,3
          offarp(i)= 0.d0
        enddo    
c*        READ(IOBS,71,iostat=ioerr) antcod,OFFSL1(1),OFFSL1(2),OFFSL1(3),
c*     $            OFFSL2(1),OFFSL2(2),OFFSL2(3)
c*   71 FORMAT (1x,a6,12x,3F8.4,3X,3F8.4)   
      endif

c     Read the satellite information

      read(iobs,'(/,a80)',iostat=ioerr) scarr
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading SV line',ioerr)                        
c     for now hard-wire the observable types to 4 since that's all GAMIT supports
      read( scarr,'(1x,i2,23x,i2)',iostat=ioerr) nchan,ndat
      if( nchan.gt.maxsat) then
         write(message,'(2a,i2)') 'Number of satellites on X-file '
     .   ,' greater than MAXSAT = ',maxsat
         call report_stat('FATAL',prog_name,'lib/xhdred',xfname,message
     .                    ,0)  
      endif
      if( ndat.gt.4) call report_stat('FATAL',prog_name,'lib/xhdred'
     .             , xfname,'Number of observables on X-file > 4 ',0)
      read( scarr,'(1x,i2,23x,i2,12x,4i3,2x,4(1x,a3))',iostat=ioerr) 
     .            nchan,ndat,(dattyp(j),j=1,4),(rxobtyp(j),j=1,4)
      if(ioerr.ne.0) call report_stat('FATAL',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading nchan ndat dattyp',ioerr)
      do i=1,nchan   
        read(iobs,'(a80)',iostat=ioerr) scarr
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .     ,'lib/xhdred',xfname,'Error reading SV line',ioerr) 
        if( ndat.eq.0 )  then           
c         this read for very old formats which appear no longer to exist  
c         read(scarr,'(18x,i2)')) ischan(i)
          call report_stat('FATAL',prog_name,'lib/xhdred', ' ' 
     .                    ,'Obsolete x-file format',0)                
        elseif( scarr(23:25).eq.'NAV' ) then 
          read(scarr,'(18x,i2,20x,7i3)',iostat=ioerr) 
     .          ischan(i),(lambda(i,j),j=1,ndat) 
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .              ,'lib/xhdred',xfname
     .              ,'Error reading SV line for pre-GNSS format',ioerr)
          satnam(i) = scarr(23:32)//'      '
        else
          read(scarr,'(17x,a16,7x,7i3)',iostat=ioerr) 
     .               satnam(i),(lambda(i,j),j=1,ndat)
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .               ,'lib/xhdred',xfname
     .               ,'Error reading SV line for GNSS format',ioerr)
          read(satnam(i)(1:1),'(a1)') gnss 
          read(satnam(i)(2:3),'(i2)') ischan(i)    
        endif
c       trap zero L1 and L2 wavelength factors (lib/rrxhed bug with RINEX 2)
        if( lambda(i,1).eq.0.and.lambda(i,2).eq.0 ) 
     .      call report_stat('FATAL',prog_name,'lib/xhdred',' '
     .                       ,'Both phase LAMBDAs zero: bad X-file',0)
      enddo        

C     Skip 1 line
      READ(IOBS,'(A)') SCARR

c     Read time type (UTC or GPST)
      read(iobs,'(14x,a3)',iostat=ioerr) x_time_flag   
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading time flag',ioerr)
      if( x_time_flag.eq.'   '.or.x_time_flag.eq.'UTC') then
        mtime = 1
      elseif (x_time_flag.eq.'GPS') then
        mtime = 2
      else
        call report_stat('FATAL',prog_name,'lib/xhdred',' '
     .                  ,'Unknown time flag',0)
      endif

c     Skip 1 line
      read(iobs,'(a)') scarr

C     Read the start time, interval, and session number
      READ(IOBS,120,iostat=ioerr) 
     .        IY,IDYOYR,IHR,MIN,SEC,INTER,ircint,isessn
  120 FORMAT (1X,I2,1X, I3,1X, I2,1X,I2,1X,F6.3,1X,I6,11x,i6,11x,i2)
cd      print *,'XHDRED start iy idyoyr ',iy,idyoyr
      call fix_y2k(iy)
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading start time',ioerr)
C     Convert day of year to month and day of month
      CALL MONDAY(IDYOYR,IM,ID,IY)
c     If the session number is missing, set = 1
      if( isessn.eq.0 ) isessn=1


C     Read the number of epochs
      READ (IOBS,130,iostat=ioerr) NEPOCH
  130 FORMAT (/,I4)                        
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading number epochs',ioerr)
c      write(*,*) 'nepoch',nepoch

C      For old X-files, set up lambda from data and frequency factors
c      (line is blank for new X-files)

C     See if have single band or dual band data and read frequency factor
c        set up lambda array containing wavelength factors
c           lambda = 1 for unambiguous, undoubled values
c           lambda =-1 for   ambiguous, undoubled values
c           lambda = 2 for unambiguous,  doubled values
c           lambda =-2 for   ambiguous,  doubled values

      READ(IOBS,140,iostat=ioerr) IDUAL,FRQFAC
  140 FORMAT(1X,2A12)      
      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
     .    ,xfname,'Error reading header frequency factor',ioerr)
c      write(*,*) 'idual,frqfac',idual,frqfac
      IF (NDAT.EQ.0) THEN
C        Backwards compatibility for old X-files
         IF (FRQFAC.EQ.'HALF-CYCLES ') THEN
            IFAC = 2
         ELSE
            IFAC = 1
         ENDIF
         IF (IDUAL.EQ.'DUAL-BAND   ') THEN
            NDAT = 4
            DATTYP(1) = 1
            DATTYP(2) = 2
            DATTYP(3) = 3
            DATTYP(4) = 4
            DO 142 I = 1, NCHAN
               LAMBDA(I,1) = -IFAC
               LAMBDA(I,2) = -IFAC
               LAMBDA(I,3) =  IFAC
               LAMBDA(I,4) =  IFAC
 142        CONTINUE
         ELSE
            NDAT = 2
            DATTYP(1) = 1
            DATTYP(2) = 3
            DO 143 I = 1, NCHAN
               LAMBDA(I,1) = -IFAC
               LAMBDA(I,2) =  IFAC
 143        CONTINUE
         ENDIF
      ENDIF

C     Check for expanded data format
C     SKD = 'S'TATIC, 'K'INEMATIC, or 'D'YNAMIC
      READ(IOBS,160,iostat=ioerr) KINHED,SKD
  160 FORMAT (/,A4,18X,A1)                    
c   RWK 160420: Kinematic/dynamic no longer supported, so line may be blank (ignore)
c      if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/xhdred'
c     .    ,xfname,'Error reading S/K/D ',ioerr)

C     Uninitialized variable: NPART (determined from Tfile in SETUP)
      NPART = 0

C     Print the X-File header information
      if( iprnt.gt.0 ) then
      WRITE(IPRNT,150) sitnam16,LATD,LATM,SECLAT,LOND,LONM,SECLON,HEIGHT
     A              , rcvrsw,swver,antcod,offarp
     1              , IY,IM,ID,IHR,MIN,SEC,x_time_flag
     .              , NCHAN,NEPOCH,INTER,ircint,isessn
     3              , NPART,IDUAL,FRQFAC
  150 FORMAT(/,1X,'Site: ',A16,/
     2        ,1X,'Coordinates (lat,lon,rad):',2(2X,I3,1X,I2,1X,F8.5)
     3        ,      2x,F14.6,/
     A        ,1X,'Receiver software:',2x,a3,2x,f5.2,'  Antenna: ',a6,/
     4        ,1X,'Antenna offsets (Up North East) ',3F7.4,'  meters',/
     6        ,1X,'Start time:',I5,2I3,2X,2I3,F7.3,3x,a,//
     7        ,1X,'Channels=',I2,3X,'Epochs= ',I4,3X,'Interval=',I4,//
     8        ,1x,'Original data interval=',i4,3x,'Session=',i2,//
     8        ,1X,'Partials=',I2,3X,'IDUAL=',A12,3X,'FRQFAC=',A12)
      write(iprnt,'(/,1x,a,6i3)') ' Channel  PRN   Data type: '
     .                   ,(dattyp(i),i=1,ndat)
      write(iprnt,'(2a)') '                             LAMBDA'
     .           ,'        RXOBTYP '
      do i=1,nchan
        write(iprnt,'(4x,i2,6x,i2,14x,4i3,2x,4(1x,a3))')
     .                    i,ischan(i),(lambda(i,j),j=1,ndat)
     .             ,      (rxobtyp(j),j=1,ndat) 
      enddo
      endif

C     Skip 1 line and be ready to read data
      READ(IOBS,170) scarr
  170 FORMAT (A80)
c      write(*,*) 'scarr',scarr      
      return

  900 write(message,'(2a,i2)') 'Number of satellites on X-file '
     .   ,' greater than MAXSAT = ',maxsat
      call report_stat('FATAL',prog_name,'lib/xhdred',' ',message,0)

      END

