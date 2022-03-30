      Subroutine TROT ( tfin, tfout, idir, frame  )
c
C
C Convert earth-fixed ephemeris to 1950.0 or J2000 inertial ephemeris
C  or vice versa
C
C Written by Yehuda Bock and Robert King
c mods:
c    Mar 1990 by pch - Apollo couldn't handle the large unformatted reads,
c                      so I broke it into several reads - also see subroutine thdrit
C
c    Nov 1991 by rwk - change call arguments to be more explicit
c
c    Jul 1994 by rwk - use library version of THDRED which converts time
c                      to GPST
c    May 1995 by pt  - pass in the inertial frame required. This needs to be written
c                      on the tfile header when going from efixed to inertial     
c    Mar 1997 by pt  -  replace efixed variable with better use of frame variable
c    Oct 2015 by lei - read iut1pol from sestbl.
c    May 2020 by tah - Updated the read iut1pol from sestbl. to use standard
c                       subroutine get_iut1pol.f (based on Lei's addition)

c
c     tfin :  input T-file name  (assigned unit = 16 by OPENB)
c     tfout:  output T-file name (assigned unit = 17 by OPENB)
c     idir :  direction of rotation  (checked against file header)
c             =  1  earth-fixed to inertial
c             = -1  inertial to earth-fixed   

c     Note:  OPENB also assigns units 5,6,14,15,20,21,22, and 30.

      implicit none

      include '../includes/dimpar.h'

      integer*4 maxtrm
      parameter(maxtrm=(3+maxorb)*maxsat)


      character*1  gnss
      character*2  aerp
      character*4  icsnam(maxorb)
      character*5  precmod,nutmod,gravmod,frame,tin_frame,srpmod
     .          ,  eradmod,antradmod
      character*16 tfin,tfout,bcfile,xfilef,printfile,satnam(maxsat)
      character*80 prog_name
      character*256 message
c
      real*8 tb,te,tf,tdtgpst
     .     , stveci,stveco,t,fract,x,satics
     .     , rot,rotdot,sdelt,sidtm,xpole,ypole
c
      integer*4 jd,jde,jdb,jdf
     .        , iscrn,iterm,iutin,iutout,iprnt,iut1,inut,ipole,iubc,iux
     .        , ifrstn,ifrstp,ifrstu,nintrs,itsat(maxsat)
     .        , ip,idir,ipar,isat,nvel,nsat,nics
     .        , nepchs,iut1pol,i,j,len,rcpar
C
C
      dimension stveci(maxtrm),stveco(maxtrm),x(maxtrm,maxsat)
     .        , satics(maxorb,maxsat),rot(3,3),rotdot(3,3)
c   
      logical reqd, kbit, debug/.false./, first/.true./
      integer ill, ioerr
C
c
C Flags for Tables
      ifrstn=0
      ifrstp=0
      ifrstu=0

c
c     Get the program name calling trot
      len = rcpar(0,prog_name)
C
C Open all the files

      bcfile='                '
      xfilef='                '
      printfile = 'trot.out'
      call openb( iterm,iscrn,iprnt,iutin,iutout,iubc,iux
     .          , inut,iut1,ipole,tfin,tfout,bcfile,xfilef,printfile )


c Set the frame indicator to the expected value from the direction flag
c  --thdred will verify that the input T-file is correctd. 

      if( idir.gt.0 ) then
         tin_frame = 'EFIXD'              
         call report_stat('STATUS',prog_name,'orbits/trot',' '
     .     ,'Input T-file is Earth-fixed, converting to inertial',0)
      elseif( idir.lt.0 ) then
         tin_frame = 'INERT'   
         call report_stat('STATUS',prog_name,'orbits/trot',' '
     .     ,'Input T-file is inertial, converting to Earth-fixed',0)
      else                                              
         write(message,'(a,i5)')  'idir is undefined. idir=',idir
         call report_stat('FATAL',prog_name,'orbits/trot',' ',message,0)
      endif

c Read the input T-file header

c       iprnt = 6
       call thdred( iutin,iscrn,iprnt,nsat,gnss,itsat,satnam
     .            , jdb,tb,jdf,tf,sdelt,nepchs,jde,te
     .             ,nics,satics,nintrs,icsnam
     .            , precmod,nutmod,gravmod,tin_frame,srpmod
     .            , eradmod,antradmod )
              
c Set the "frame" variable for the correct inertial frame - either from the
c tfile header (inert-efixed) or from the passed in value (efixed-inert)

c if input tfile is earth-fixed, set it according to the frame required (as passed
c in to this subroutine
      if(tin_frame.eq.'EFIXD')then
        if(frame.eq.'B1950')precmod = 'IAU68'
c        if(frame.eq.'J2000')precmod = 'IAU76' 
        if(frame.eq.'J2000') then
c Should put some logic in here to get the precession/nutation 
c model from the sestbl. rather than have it hardwired based on frame?
c Current possibilities are
c         precmod = 'IAU76'  ! Old GAMIT IAU76/IAU00 model (GAST)
c         precmod = 'IAU0A'  ! SOFA IAU00A/IAU00A equinox based (GAST)
c         precmod = 'IAU0C'  ! SOFA IAU00A/IAU00A CIO based (ERA)
c         precmod = 'IAU06'  ! SOFA IAU06/IAU00A  CIO based (ERA)
c         precmod = 'IAU6A'  ! SOFA IAU06  CIO based x/y series (ERA)
          tin_frame = frame
        endif
      endif
c Convert the initial conditions

      tdtgpst = 32.184d0 + 19.d0
c      **  Hardwire no diurnal pole/ut1 terms for now
c      iut1pol = 0
c      ** change to use rwk/tah 990310
c      iut1pol = 3  
c      ** change to use new (Ray) model  tah/rwk 990319
c      iut1pol = 7   
c      ** change to use IERS2010 model tah/rwk 150325 
c MOD TAH 200504: Add the libration terms as default as well. 
       iut1pol = 11 + 16  ! Bits 1,2,4 and 5.  Ideally sestbl. will be used.
c      ** change to read iut1pol from sestbl. lei 151002
      open( unit=31, file='sestbl.',status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
          call report_stat('WARNING',prog_name,'orbits/trot','sestbl.',
     .    'Errors in opening session table,assuming iut1pol=11',ioerr)
      else
* NID TAH 200505: Use the standard get_iut1pol subroutine.  Also implements
*     the Desai and Subious model.
        call get_iut1pol( iut1pol ) 
      endif
         
c------         
      call rotsnp( idir,jde,te,tdtgpst,iut1pol
     .            , tin_frame,precmod,iut1,ipole,inut
     .            , rot,rotdot,sidtm,xpole,ypole )
        
      do j=1,nsat
        do i=1,6
          stveci(i) = satics(i,j)
        enddo                         
        call rotcrd( rot,rotdot,stveci,stveco,1,1 )
        do i=1,6
           satics(i,j) = stveco(i)
        enddo
      enddo

c Write the output T-file header

      if( idir.lt.0 ) then
         frame = 'EFIXD'
      endif
          
      call thdrit( iutout,jde,te,jdb,tb,jdf,tf,sdelt,nepchs,nintrs
     .            , nsat,gnss,itsat,satnam,satics,nics,tfin,icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .            , antradmod )


c Set the velocity and partials counters

      nvel = 0
C                      NVEL= 0   No velocities on T-file
C                          = 1   Velocity coordinates on T-file
C
      ipar = nintrs/3
C                      IPAR= 1   Coordinates only on T-file
C                          > 1   Coordinates + (IPAR-1) partials


C        Begin loop over epoch records of T-file
C
      jd= jdb
      t=  tb

C Re-initialize flags for tables
      ifrstn=0
      ifrstp=0
      ifrstu=0

100   read(iutin,end=200) ((x(i,isat),i=1,nintrs),isat=1,nsat)

      fract= t/86400.D0
c     write(*,'(1x,i8,f15.6,f15.12)') jd,t,fract
c      print *,'jd t tdtgpst ',jd,t,tdtgpst
      call rotsnp( idir,jd,t,tdtgpst,iut1pol
     .           , tin_frame,precmod,iut1,ipole,inut
     .           , rot,rotdot,sidtm,xpole,ypole )
      if( debug.and.first ) then
         print *,'idir,jd,t,tdtgpst,iut1pol '
     .      ,idir,jd,t,tdtgpst,iut1pol
         print *,'rot ',rot
         print *,'rotdot ',rotdot
         print *,'sidtm xpole,',sidtm,xpole,ypole
         write(*,'(/,3f20.9)') ((x(i,isat),i=1,nintrs),isat=1,nsat)
         first = .false.
      endif
      do 130 isat=1,nsat
        do ip=1,nintrs
           stveci(ip)= x(ip,isat)
        enddo 
        call rotcrd( rot,rotdot,stveci,stveco,nvel,ipar )
c       if( t.gt.1300.d0 .and. t.lt.1600.d0 .and. isat.eq.2) then
c         if( isat.eq.2 .and. t.lt.70000.d0 ) then
c         print *,'jd,fract,sat 2 stveci,stveco: '
c     .        , jd,fract,(stveci(i),i=1,3),(stveco(i),i=1,3)
c       stop
c       endif
        do ip=1,nintrs
          x(ip,isat)= stveco(ip)
        enddo
130   continue
c     write(*,'(/,3f20.9)') ((x(i,isat),i=1,nintrs),isat=1,nsat)
      write(iutout) ((x(i,isat),i=1,nintrs),isat=1,nsat)

c     increment the time counter
      call timinc (jd,t,sdelt )
      goto 100

c Close the files so that other programs may run afterward

200   call closeb( iutin,iutout,inut,iut1,ipole )

c Parting words
      if( idir.gt.0) then
        call report_stat('STATUS',prog_name,'orbits/trot',tfout
     .  ,'Successfully wrote Inertial T-file: ',0)
      elseif( idir.lt.0) then
        call report_stat('STATUS',prog_name,'orbits/trot',tfout
     .  ,'Successfully wrote earth-fixed T-file: ',0)
      else                          
        write(message,'(a,i3)') 
     .    'Error, unknown direction of rotation idir=',idir
        call report_stat('FATAL',prog_name,'orbits/trot',' ',message,0)
      endif

      return
      end
