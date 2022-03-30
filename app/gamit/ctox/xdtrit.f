Copyright (c) Massachusetts Institute of Technology,1987. All rights reserved.

      Subroutine XDTRIT ( iux,iep,bias )
c** rwk 160420 , latflag,dlat,mlat,slat,lonflag,dlon,mlon,slon,bias )
C
C     Write a data record from the X-File - R. King  13 April 1987
C     Mod:  write data with all flags except 1 -- MHM 880331
C     include '../includes/dimpar.h'date more sats and pseudo-ranges
c     Mod:  add logic to avoid extraneous data with error flags=2
c           880726 mhm
c     Mod:  No longer write the kinematic information on the data records
C
      implicit none
C
      include '../includes/dimpar.h'
      include '../includes/model.h'

      LOGICAL BIAS(MAXSAT)    
  
      character*1 latflag,lonflag

      INTEGER IUX,iep,IAMP1,IAMP2
     2      , ICOUNT
     3      , IYR,iidoy,IHR,MIN,I
c     4      , dlat,mlat,dlon,mlon
c    the kinematic flag now unused
     5      , kflag      

      REAL*8 DATA1,DATA2,DATA3,DATA4,SOD,SEC
c    .     , slat,slon


c        Convert JD, seconds-of-day from model.h to yr, doy, hr, min, sec

      call dayjul( jdobs,iyr,iidoy ) 
      call ds2hms( iyr,iidoy,tobs,ihr,min,sec )

C        Determine the number of non-empty channels

      ICOUNT= 0
      DO I=1,NCHAN
        IF( IER(I).NE.1 ) ICOUNT = ICOUNT + 1
      enddo
c     RWK 160420: No longer write the kinematic information
c      write(iux,592) iepoch,icount,iyr,iidoy,ihr,min,sec
c     .     , kflag,sitecd,latflag,dlat,mlat,slat
c     .     , lonflag,dlon,mlon,slon,radius,offarp
c  592    format(/,2I4,I5,I4,2I3,F11.7,1x,i2,1x,a4,1x,
c     .          a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4,
c     .          1x,3F8.4)
      write(iux,'(/,2i4,i5,i4,2i3,f11.7)') 
     .     iep,icount,iyr,iidoy,ihr,min,sec


      DO 40 I=1,NCHAN

      IF( IER(I).NE.1 ) THEN

C**     Note sign reversal in phase from C-File
        DATA1= -OBS(1,I)
        DATA2= -OBS(2,I)
        DATA3=  OBS(3,I)
        DATA4=  OBS(4,I)
        IAMP1= ISNR(1,I)
        IAMP2= ISNR(2,I)
        IF( IER(I).LT.-1 .OR. IER(I).GT. 98) THEN
          call report_stat('WARNING','CTOX','xdtrit',ier(i),
     .    'Problem with flags: should be between -1 and 98 ',0)
        ENDIF
        IF ( BIAS(I) ) THEN
          IER(I)=10
          BIAS(I)=.FALSE.
          ENDIF
        IF( LAMBDA(I,1).EQ.0 )
     1       WRITE(IUX,35) IER(I),I,DATA1,IAMP1
        IF( LAMBDA(I,1).NE.0 )
     1    WRITE(IUX,35) IER(I),I,DATA1,IAMP1,DATA2,IAMP2,DATA3,DATA4
35      FORMAT (10X,2I2,2(1X,D22.15,1X,I3),2(2X,D22.15))

      ELSE
        CONTINUE
      ENDIF

40    CONTINUE

C
      RETURN
      END
