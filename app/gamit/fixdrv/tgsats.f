      Subroutine  tgsats ( tfile, gfile, ltfil, lgfil, iscrn
     .                   ,  ntsat, itsat, satnam,  nintrs, nics
     .                   , itb, tbb, itstp, tstp )

c      R. King  920911

c      Read the T-(G-)file to get satellites and start, stop times
c      Check the consistency of T-file partials and experiment type

C     input
c        tfile  : T-file name
c        gfile  : G-file name
C        itfil  : unit number for T-file
c        lgfil  : unit number for G-file
C        iscrn  : unit number for screen
c        efixed : .true. if expecting earth-fixed orbits

C     output
C        ntsat  : number of satellites on T-file
c        itsat  : array for PRN numbers 
c        satnam : array for SV names 
C        nics   : number of satellite parameters (9 if 3 rad parms; 15 if 9 rad parms)
C        itb    : begin time
C        tbb    : begin time
C        itstp  : end time
C        tstp   : end time

C     local
C        jdb    : start (PEP) Julian Day of ephemeris
C        tb     : start seconds of day of ephemeris
C        jds    : stop JD
C        ts     : stop seconds of day
C        sdelt  : tabular interval of ephemeris in seconds


      implicit none

      include '../includes/dimpar.h'

      character*1  gnss(maxsat)
      character*4  icsnam(maxorb),time_type
      character*5  precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*16 tfile, gfile, satnam(maxsat), gsatnam(maxsat)

      integer * 4  ltfil,lgfil,iscrn,iprnt,jdb,jds,jde,nintrs
     .          ,  nepcht,nics,ntsat,itsat(maxsat),ngsat,igsat(maxsat)
     .          ,  iyr,idoy,ihr,min
     .          ,  ite(3),itb(3),itstp(3),ioerr,i

      real * 8     te,tb,ts,tee(3),tbb(3),tstp(3),sec
     .          ,  satics(maxorb,maxsat),sdelt 

      logical fcheck

      parameter   ( iprnt=8 )

c       Set default number of orbital parameters to 9 (3 rad parameters)
c       overriden by value on T-file or, later in BMAKE, by sestbl. input

      nics = 9


c       Determine the SVs from the G-file

      if( fcheck(gfile) ) then  
         call read_gfile( gfile,lgfil,jde,te,frame,precmod,nutmod
     .                  , gravmod,srpmod,eradmod,antradmod,time_type
     .                  , ngsat,nics,icsnam,gsatnam,satics ) 
cd         print *,'TGSATS ngsat gsatnam(1)',ngsat,gsatnam(1)
         do i=1,ngsat
           read(gsatnam(i)(2:3),'(i2)') igsat(i)
         enddo
       endif

c       Deterime the SVs and start/stop times from the T-file

      if( fcheck(tfile) ) then

        call lowers( tfile )    
        call topens( tfile, 'old', ltfil, ioerr )
        if( ioerr.ne.0 )
     .      call report_stat('FATAL','FIXDRV','tgsats',tfile
     .                      , 'Error opening T-file: ',0)
        open( iprnt, status='scratch' )    
        frame = 'UNKWN'  
        call thdred( ltfil, iscrn, iprnt, ntsat, gnss, itsat, satnam
     .             , jdb, tb, jds, ts, sdelt, nepcht, jde, te
     .             , nics, satics, nintrs,  icsnam
     .             , precmod, nutmod, gravmod, frame, srpmod
     .             , eradmod, antradmod )
        close ( ltfil )
        close ( iprnt )
        call dayjul( jdb,iyr,idoy )
        call ds2hms( iyr, idoy, tb, ihr, min, sec )
        tbb(1) = ihr
        tbb(2) = min
        tbb(3) = sec
        itb(3) = iyr
        call monday( idoy, itb(1), itb(2), itb(3) )
        call dayjul( jds,iyr,idoy )
        call ds2hms( iyr, idoy, ts, ihr, min, sec )
        tstp(1) = ihr
        tstp(2) = min
        tstp(3) = sec
        itstp(3) = iyr
        call monday( idoy, itstp(1), itstp(2), itstp(3) )
        call dayjul( jde,iyr,idoy )
        call ds2hms( iyr, idoy, te, ihr, min, sec )
        tee(1) = ihr
        tee(2) = min
        tee(3) = sec
        ite(3) = iyr
        call monday( idoy, ite(1), ite(2), ite(3) )


C        Debug: Print the header information
c            write( iscrn, '(//,A,//,1X,/,3(/,3I5,2F4.0,F7.3,A),
c     .                   /,F10.3,A,I8,A,/,i3,a,/)' )
c     .        ' T-File Header Information:',
c     .        ITE,   TEE,  '  Epoch of initial conditions (GPST)',
c     .        ITB,   TBB,  '  Ephemeris start',
c     .        ITSTP, TSTP, '  Ephemeris end',
c     .        SDELT,       '  Tabular interval (sec) ;',
c     .        nepcht,      '  T-file epochs',
c     .        nics,        '  Number of orbit parameters'
C
c            WRITE( ISCRN, '(A,6X,3(A4,9X),A,/,26X,3(A4,9X),A)' )
c     .        '   PRN #   (km) / (km/s) / (dimensionless)'
c            do  j= 1, ntsat
c               write(iscrn, '(1X,I3,1X,3F13.5,/,5x,3F13.9,
c     .                      /5x,9F8.5)' )
c     .                itsat(j), (satics(i,j),i=1,nics)
c             enddo

      endif   

c     compare the SV lists from the g and t files  

      if( fcheck(gfile) ) then
        if( fcheck(tfile) ) then
           if( ngsat.lt.ntsat ) 
     .         call report_stat('WARNING','FIXDRV','tgsats',' '
     .    ,'More SVs on T-file than G-file; G-file list used for ARC',0)
        endif  
        ntsat = ngsat
        do i=1,ngsat
          itsat(i) = igsat(i)                    
          satnam(i) = gsatnam(i)
        enddo
      endif                                                              
cd      print *,'end of TGSATS ngsat ',ntsat
cd      do i=1,ntsat
cd       print *,itsat(i),gsatnam(i)
cd      enddo

c     warn if neither t nor g file available

      if ( .not.fcheck(gfile) .and. .not.fcheck(tfile) )
     .   call report_stat('WARNING','FIXDRV','tgsats',' '
     .          , 'No T- or G-file available',0)

      return
      end
