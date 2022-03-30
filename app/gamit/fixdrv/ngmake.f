      Subroutine  NGMAKE( bfil2, tfile, yr, doy, lsess, ierr_sestbl )

C     generate a batch file to run TTONGS
C     (convert t-file to NGS SP1 or SP3 precise ephemeris format)
C
C     arguments
C        BFIL2  : Batch file name
C        TFILE  : T-file name
C        YR     : year(2 character)
C        DOY    : day of year
C        LSESS  : Logical unit number of SESTBL. file
C        ILL    : return code
C
      implicit none

      include '../includes/dimpar.h'
C
      CHARACTER*2   YR
      CHARACTER*3   ORG,DOY,ORFM,upperc
      CHARACTER*4   ORID
      CHARACTER*5   ORBRF
      CHARACTER*16  TFILE,BFIL2
      CHARACTER*12  NFILE

      logical reqd

      integer*4 lsess, ill, ierr_sestbl
      integer*4 gpsweek,gpsday,idoy,iyr

      NFILE = '            '

      WRITE( 17, '(A,A16)' )  'ttongs < ', BFIL2
C
      OPEN( 21, FILE=BFIL2, FORM='FORMATTED',status='unknown' )
C
C     write t-file name
C
         WRITE( 21, '(A16)' )  TFILE
C
C     NGS file name - old garner/pgga format
      reqd = .false.
      ORID = '    '
      CALL  RDSEST( 8, 'Orbit id', 4, ORID, LSESS
     .            , reqd, ILL )
      if( ill.ne.0 ) ierr_sestbl = ill
      IF( ORID .EQ. '   ')  THEN
c       File name and organization required for export orbits
        reqd = .true.
C       NGS file name - new IGS format
        ORFM = '   '
        CALL  RDSEST( 12, 'Orbit Format', 3, ORFM, LSESS
     .            , reqd, ILL )
        if( ill.ne.0 ) ierr_sestbl = ill
        IF(ORFM.NE.'SP1' .AND. ORFM.NE.'SP3') then
           call report_stat('WARNING','FIXDRV','ngmake',' '
     .          , 'Unrecognized Orbit Format in sestbl.',0)
           ierr_sestbl = -1
         endif
      ENDIF

C     read organization from session table
C
      CALL  RDSEST( 18, 'Orbit organization', 3, ORG, LSESS,
     .              reqd, ILL )
      if( ill.ne.0 ) ierr_sestbl = ill
      IF( ORG.EQ.'   ' ) THEN
          call report_stat('WARNING','FIXDRV','ngmake',' '
     .          , 'No Orbit Organization in sestbl.',0)
           ierr_sestbl = -1
      endif

C Old format
      if(ORID.NE.'    ') then
      NFILE(1:1)='n' 
      call lowers(orid)
      NFILE(2:5)= orid
      NFILE(6:7)=yr
      NFILE(8:8)='.'
      NFILE(9:11)=doy

      WRITE( 21, '(A12)' )  NFILE
c F  (Type of orbit for header; F=fitted, E=extrapolated, B=broadcast)
c 90 (Coordinate system; 84=WGS84, 90=ITRF90)
      WRITE( 21, '(A3,/,A1,/,A2)' )  ORG,'F','90'

      elseif(ORFM.NE.'   ') then

C     Reference system for orbit
      CALL  RDSEST( 26, 'Reference System for Orbit', 5, ORBRF, LSESS
     .            , reqd, ILL )
      if( ill.ne. 0 ) then
         call report_stat('WARNING','FIXDRV','ngmake',' '
     .          , 'No Reference System for Orbit in sestbl.',0)
         ierr_sestbl = ill
       endif

C Output orbit format
      WRITE( 21, '(A3,/,A3)' ) upperc(ORFM),upperc(ORG)

C Convert yr and doy to gpsweek and gpsday
      read(yr,'(i2)') iyr
      call fix_y2k(iyr)
      read(doy,'(i3)') idoy
      call doygwk(idoy,iyr,gpsweek,gpsday)
c      write(*,*) gpsweek,gpsday
      NFILE(1:3)=ORG
      write(NFILE(4:7),'(i4)') gpsweek
	if (nfile(4:4).eq.' ') nfile(4:4)='0'
      write(NFILE(8:8),'(i1)') gpsday
      NFILE(9:9)='.'
      NFILE(10:12)=ORFM
      call lowers(nfile)
      WRITE( 21, '(A12)' )  nfile

c F  (Type of orbit for header; F=fitted, E=extrapolated, B=broadcast)
c  & Reference Coordinate system for orbit
      WRITE( 21, '(A3,/,A5)' )  'FIT',ORBRF

      ENDIF
C
C
      CLOSE( 21 )
      RETURN
      END



