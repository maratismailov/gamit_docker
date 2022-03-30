Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE OPENB
     1   ( ITERM,ISCRN,IPRNT,IUTIN,IUTOUT,IUBC,IUX,INUT,IUT1,IPOLE
     2   , tfin,tfout,bcfile,xfile,printfile )
C
C      Open all the input and output files for subroutine TROT
C        and also its calling programs NGSTOT, TTONGS, BCTOT
C      R. King   17 November 1987
C      Modified for BCTOT - YB 1/9/88
c      mod to remove the calls to tables - pch mar 1990
C
      implicit none

      CHARACTER*16 tfin,tfout,bcfile,xfile,printfile
      CHARACTER*80 UT1F,POLEF,NUTF,prog_name

      integer*4 iscrn,iterm,iutin,iprnt,iut1,inut,iutout,ipole
     .        , iubc,iux,ioerr,len,rcpar
      logical fcheck 
c
      len = rcpar(0,prog_name)
c
      UT1F='ut1.'
c      CALL TABLES(UT1F)
      POLEF='pole.'
c      CALL TABLES(POLEF)
      NUTF='nutabl.'
c      CALL TABLES(NUTF)
c
      ITERM =   5
      ISCRN =   6
      IUX   =  14
      IUBC  =  15
      IUTIN =  16
      IUTOUT=  17
      INUT  =  20
      IUT1  =  21
      IPOLE =  22
      IPRNT =  30
C
C        Open the print files
C
      OPEN( UNIT=IPRNT,FILE=printfile,FORM='FORMATTED'
     .      ,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',printfile,
     .    'Error opening print file: ',ioerr)
      else
        call report_stat('STATUS',prog_name,'orbits/openb',printfile,
     .    'Opened print file: ',ioerr)
      endif
C
C        Open the output T-file
C
      OPEN( UNIT=IUTOUT,FILE=tfout,FORM='UNFORMATTED',STATUS='UNKNOWN'
     .    , iostat=ioerr )
      if (ioerr .ne. 0 ) then
         call report_stat('FATAL',prog_name,'orbits/openb',tfout,
     .    'Error opening output T-file: ',ioerr)
      else
        call report_stat('STATUS',prog_name,'orbits/openb',tfout,
     .    'Opened output T-file: ',ioerr)
      endif
      WRITE(IPRNT,22) tfout
22    FORMAT(/,' Output Ephemeris (T-) File :  ',A16)


c        Open the input T-file if it is required (TROT only)

      if( tfin.ne.'               ') then
      OPEN ( UNIT=IUTIN,FILE=tfin,FORM='UNFORMATTED',STATUS='OLD'
     .     ,iostat=ioerr )
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',tfin,
     .    'Error opening input T-file: ',ioerr)
      else
        call report_stat('STATUS',prog_name,'orbits/openb',tfin,
     .    'Opened input T-file: ',ioerr)
      endif

      WRITE(IPRNT,21) tfin
21    FORMAT(/,' Input Ephemeris (T-) File :   ',A16)
      endif

c         Open the E-file if it is required (BCTOT only)

      if( bcfile.ne.'                ') then
        WRITE(IPRNT,31) bcfile
31      FORMAT(/,' Input Navigation File : ',A16)
        open(unit=iubc,file=bcfile,status='old',form='formatted'
     1      ,iostat=ioerr)
        if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',bcfile,
     .    'Error opening input E-file: ',ioerr)
        else
          call report_stat('STATUS',prog_name,'orbits/openb',bcfile,
     .    'Opened input navigation file: ',ioerr)
        endif
      endif

c         Open the X-file if it is required

      if( xfile(1:1).ne.' ') then
      WRITE(IPRNT,41) xfile
41    FORMAT(/,' Input X-file :   ',A16)
      open(unit=iux,file=xfile,status='old',form='formatted',
     1     iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',xfile,
     .    'Error opening X-file: ',ioerr)
      endif
      endif

C
C        Open TAI-UT1 Table File
C
      OPEN (UNIT=IUT1,FILE=UT1F,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',ut1f,
     .    'Error opening UT1 table: ',ioerr)
      endif
C
C        Open Nutation Table File
C             
      if( .not.fcheck('nbody') ) then 
        OPEN (UNIT=INUT,FILE=NUTF,STATUS='OLD',iostat=ioerr)
        if (ioerr .ne. 0 ) then
            call report_stat('FATAL',prog_name,'orbits/openb',nutf,
     .      'Error opening nutation table: ',ioerr)
        endif
      endif 

C        Open Pole Table File
C
      OPEN (UNIT=IPOLE,FILE=POLEF,STATUS='OLD',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/openb',polef,
     .    'Error opening pole table: ',ioerr)
      endif
C
      RETURN
      END
