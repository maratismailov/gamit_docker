      SUBROUTINE  ARMAKE( BFIL2, TFILE, GFILE, ISAT, itsat, satnam
     .                  , ITB, TBB, ITSTP, TSTP, LSESS, ierr_sestbl
     .                  , lexp, EXPO, norb )

C     Generate a batch file to run ARC

C     Input:
C        TFILE  : T-file name
C        GFILE  : gfile name from header line
C        ISAT   : number of satellites on T-file
C        itsat  : prn numbers of satellites on T-file
c        satnam : satellite name (shortened form of svnav.dat entries0
C        ITB    : begin time (m/d/y)
C        TBB    : begin time (h/m/s)
C        ITSTP  : end time   (m/d/y)
C        TSTP   : end time   (h/m/s)
C        EFIXED : .true. if we are expecting earth-fixed
C        LSESS  : unit number of session file
C        IERR_SESTBL : return code for error reading sestbl.
c        LEXP   : type of experiment (1=BASELINE, 2=RELAX, 3=ORBIT, 4=KINEN)
C        EXPO   : Export Orbit code (Y/N)

c     Output:
c        norb   : Number of orbit parameters (needed for SOLVE), determined from
c                   radiation-pressure model
C
      implicit none
C
      include '../includes/dimpar.h'

      integer*4 nblen,idoyb,lsess,itstp(3),itb(3),ierr_sestbl,ill,i
     .      , idoyst,lexp,isat,idoy,i1,ihup,minup,itsat(maxsat),norb
     .      , idbprt,iug

      real*8 tbb(3),tstp(3),secup,tabint,stepsize

      CHARACTER* 1 EXPO,buf1
      character* 2 buf2,gravdeg,etidedeg,otidedeg
* MOD TAH Added string for lbody (0 or 1 currently, but may increase latter
*     if more than just Jupiter and Venus (JV) included.
      character*(2) lbody_str   ! String to contain 0 or 1 for JV perturbations
      character* 4 time_type
      CHARACTER* 5 buf5,radmod,REFSYS,frame5,precession
     .            ,eradmod,antradmod
      character* 6 buf6
      character* 8 buf8
      CHARACTER*16  TFILE, GFILE, outfile, SATNAM(MAXSAT), BFIL2
      character*19 frame
      character*256 message

      logical reqd
C
      WRITE( 17, '(A,A16)' )  'arc    < ', BFIL2
C
      OPEN( 21, FILE=BFIL2, FORM='FORMATTED',status='unknown' )
C

c     determine the satellite names from the PRN numbers

      write(buf2,'(i2)')

c Write the satellite names
             
      do i=1,isat
         write( 21, '(a16)' ) satnam(i)
      enddo

* MOD TAH 200304: Set defaults before we start.  This way if entries are
*     not in sestbl. We will have the defaults.  (This is for new
*     Inertial Reference System).
      precession = 'IAU0A'
      
C Read the reference system from the sestbl

      reqd = .false.
* MOD TAH 200304: Changed default to EGR08. Relatisitic verison.
*     refsys = 'EGM08'
      refsys = 'EGR08'
      CALL  RDSEST( 24, 'Reference System for ARC', 5, buf5, LSESS,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.ne.'     ') refsys = buf5
      if( refsys.NE.'WGS84' .and. refsys.ne.'WGS72' .and.
     .    refsys.ne.'MERIT' .and. refsys.ne.'IGS92' .and.
     .    refsys.ne.'EGM96' .and. refsys.ne.'EGM08' .and.
     .    refsys.ne.'EGR08' )  then
          write(message,'(a,a5)')
     .       'Improper Reference System for ARC: ',refsys
          call report_stat('WARNING','FIXDRV','armake',' ',message,0)
          ierr_sestbl = -1
      endif
                                             
c Read the degree of gravity field, solid-Earth tides, and ocean-tides to be used

      reqd = .false.
      gravdeg = '12'
      call rdsest(11,  'ARC gravdeg', 2, buf2, lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf2.ne.'  ' ) gravdeg = buf2
      etidedeg = ' 4'
      call rdsest(12,  'ARC etidedeg', 2, buf2, lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf2.ne.'  ' ) etidedeg = buf2
      otidedeg = '12'
      call rdsest(12,  'ARC otidedeg', 2, buf2, lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf2.ne.'  ' ) otidedeg = buf2

* MOD TAH 190622: Added reading planetary terms (JV) to be included.
      lbody_str = '0'
      call rdsest(11,  'ARC planets', 2, buf2, lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf2.ne.'  ' ) lbody_str = buf2

c Read the radiation-pressure models from the sestbl

      reqd = .false.
      radmod = 'BERNE'
      call  rdsest( 23, 'Radiation Model for ARC', 5, buf5, LSESS,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.eq.'     ') then
         radmod = 'BERNE'
      else
         radmod = buf5
      endif
      if( radmod.eq.'BERNE'.or.radmod.eq.'ECOM1'.or.
     .    radmod.eq.'ECOM2'.or.
     .    radmod.eq.'UCLR1'.or.radmod.eq.'UCLR2' ) then
         norb = 15  
       elseif( radmod.eq.'ECOMC' ) then 
         norb = 19 
       else
         norb = 9
       endif
      reqd = .false.
      eradmod = 'NONE '
      call  rdsest( 21, 'Earth radiation model', 5, buf5, LSESS,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.eq.'     ') then
         eradmod = 'NONE '
      else
         eradmod = buf5
      endif
      reqd = .false.
      antradmod = 'NONE '
      call  rdsest( 20, 'Antenna thrust model', 5, buf5, LSESS,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.eq.'     ') then
         antradmod = 'NONE '
      else
         antradmod = buf5
      endif


c Read the integration controls for ARC

c     get 'Tabular interval for ARC'
      call  rdsest( 7, 'Tabular',6,buf6,lsess, reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf6.eq.'      ' ) then
          tabint = 900.
      else
          read(buf6,'(f6.0)') tabint
      endif
c     get 'Stepsize for ARC'
      call  rdsest( 8, 'Stepsize',6,buf8,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf8.eq.'        ' ) then
          stepsize = 75.
      else
          read(buf8,'(f6.0)') stepsize
      endif


c  Hardwire GPST but input inertial reference frame

      reqd = .false.
      frame5 = 'J2000'
      call  rdsest( 14, 'Inertial frame', 5, buf5, lsess,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.ne.'     ' )  frame5 = buf5
      if( frame5.ne.'B1950' .and. frame5.ne.'J2000' ) then
          write(message,'(a,a5,a)') 'Inertial frame (=',frame5
     .     ,') must be B1950 or J2000'
          call report_stat('WARNING','FIXDRV','armake',' '
     .          ,message,0)
          ierr_sestbl = -1
      endif
      if( frame5.eq.'B1950') then
        frame = 'INERTIAL      B1950'

      elseif  (frame5.eq.'J2000') then
        frame = 'INERTIAL      J2000'
      endif
      
      reqd = .false.
      call  rdsest( 14, 'Inertial Reference System', 5, buf5, lsess,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5.ne.'     ' ) then
        precession = buf5
      else
        write(message,'(a)')  'Inertial Reference System not defined &
     .  & sestbl. set to old defaults: IAU68 or IAU76'
        call report_stat('WARNING','FIXDRV','armake',' '
     .          ,message,0)
        if(frame5.eq.'B1950') precession = 'IAU68'
* MOD TAH 190918: Changed default to be IAU0A to be consistent with sp3tot.f
C       if(frame5.eq.'J2000') precession = 'IAU76'
        if(frame5.eq.'J2000') precession = 'IAU0A'

      endif     

      time_type = 'GPST'   


c  Set flags to turn off models (not normally used) to turn off solid-Earth tides (no sestbl. entry)
c  and to turn off lunar eclipses
c  RWK 120523: These options no longer written explicitly with the arcout file name, but rather
c              controlled by the binary-coded idbprt flag

      call  rdsest( 14,'Arc debug flag',6,buf8,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf8.eq.'        ' ) then
          idbprt = 0 
      else
          read(buf8,*) idbprt
      endif


c  Write the batch-file entries

      write(21,'(a3)') 'END'
* MOD TAH 190622: Added lbody to output line.
      write(21,'(a5,1x,a5,1x,f6.1,2x,f6.2,3x,a4,2x,a19,1x,a5,1x
     .     ,a5,1x,a5,1x,a2,1x,a2,1x,a2,1x,a)') 
     .          refsys,radmod,tabint,stepsize,time_type,frame,precession
     .        , eradmod,antradmod
     .        , gravdeg,etidedeg,otidedeg, lbody_str
c        list file and G-file names
      i1 = nblen( gfile )
      outfile = 'arcout.'//gfile(i1-2:i1)    
c     if( the debug flag is not zero, write it into the end of the archive file name
      if( idbprt.ne.0 ) write(outfile(13:16),'(i4)') idbprt
      write( 21, '(a16,/,a16,/)' )  outfile, gfile


c Check for validity of integration interval (may be undefined if no
c X-, C-, or T-file available in directory--this shouldn't normally happen)
       if( itb(3).lt.80 .or. itb(3).gt.2100 .or. itstp(3).lt.80  .or.
     .     itstp(3).lt.80 .or. itstp(3).gt.2100 ) then
         write(message,'(a,6i5,a)') 'Invalid T-file interval = '
     .     , itb,itstp
     .     ,' --need either T-file or one X-file to run FIXDRV'
         call report_stat('FATAL','FIXDRV','armake',' ',message,0)
       endif

C Expand integration interval if Export Orbits
      IF(EXPO.EQ.'Y') THEN
c      start time
       ihup=-12
       minup=0
       secup=0.d0
       call timeup(itb,tbb,ihup,minup,secup)
c      end time
       ihup=12
       minup=0
       secup=0.d0
       call timeup(itstp,tstp,ihup,minup,secup)
      ENDIF

C     write start, stop time in yr,doy,hr,mn,sec format
C
      idoyb = idoy( itb(3),itb(1),itb(2) )
      write( 21, '(i4,1x,i3,1x,2(i2,1x),f8.5)' )
     .      itb(3), idoyb,
     .      INT(TBB(1)), INT(TBB(2)), TBB(3)
c      WRITE( 21, '(5(I2,1X),F8.5)' )
c     .     MOD(ITB(3),100), ITB(1), ITB(2),

      idoyst = idoy( itstp(3),itstp(1),itstp(2) )
      write( 21, '(i4,1x,i3,1x,2(i2,1x),f8.5)' )
     .      itstp(3), idoyst,
     .      INT(TSTP(1)), INT(TSTP(2)), TSTP(3)
c      WRITE( 21, '(5(I2,1X),F8.5)' )
c     .     MOD(ITSTP(3),100), ITSTP(1), ITSTP(2),

C
      if ( lexp.eq.1 .or. lexp.eq.4 ) then
         WRITE( 21, '(A1,/,A16)' )  'N', TFILE
      else
         WRITE( 21, '(A1,/,A16)' )  'Y', TFILE
      endif
C
      CLOSE( 21 )
      RETURN
      END



