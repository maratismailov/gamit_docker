      Subroutine  SETUP 

      implicit none

      include '../includes/dimpar.h'
      include '../includes/units.h'  
      include '../includes/global.h'
      include '../includes/model.h'
c      include '../../libraries/includes/const_param.h' 
      include '../../libraries/includes/freq_def.h'

      character*1 upper1,latflag,lonflag
     .          , lower1,yawbia,yawbias,gnssx,gnsst
      character*3 rxobtyp(maxdat)
      character*4 none,site1
     .          , prject,orbit
     .          , cvalue,metsrc4(5),gptnam
      character*5 htcod,datum,nutmod1
      character*6 etid6,antcodx,antcods
      character*8 asite2,metsrcmod
      character*9 slprt
      character*16 sitnam16,amodel,ytfile,gpt_filnam,satnam(maxsat)
      character*20 rctypex,anttypx,rcvrsnx,antsnx
      character*80 nut_header
      character*256 message,line
      logical stnfo_exist,first,found,kbit,found_svant(maxsat)
     .        , old_stinf,fcheck,warnings,eof

      INTEGER*4 julday,iyr,nindim,iarray,iend,istart
     .      , latd,latm,lond,lonm
     .      , nintrp,icall
     .      , isat,iday,imo,ihr,imin,idoy,icheck
     .      , IDUM(MAXSAT)
     .      , isvn(maxsat),nyr,ndoy,nsod
     .      , ioerr,ksod_starts,ksod_stops
     .      , nyversn,nyinter,nyepoch,nysat,ysat(maxsat),ysvn(maxsat)
     .      , yblk(maxsat),ispan,dlat,dlon,mlat,mlon,frqchn(maxsat)
     .      , svnstart(5),svnstop(5),minelv,jdstop,indx,ivalue,ifreq
     .      , istr,mchkey,i,j,jj
     .      , intdeg 
      

      real*8  sec,radiusx,seclat,seclon
     .     ,  testt,tb,ts,fract,radian,taiutc,tdtgpst
     .     ,  span,xjd
     .     ,  humid0,e0,fjd,oblq,dum
     .     ,  xl,ascm,d,f
     .     ,  TSATIC(MAXORB,MAXSAT)
     .     ,  shftdot(3),rot(3),rotdot(3),scale,scaledot,gepoch
     .     ,  anth,antn,ante,rnut(3,3)
     .     ,  utcoff
     .     ,  twopi,atpart,sbmass,el,az
     .     ,  svoffl1(3),svoffl2(3),svcorrl1,svcorrl2,sign(2) 
     .     ,  yaw_s,yaw_e,elvcutdummy
     .     ,  t0_yrs,pos(3),site_epoch,decyrs
     .     ,  slat,slon
     .     ,  rvalue,undu,dmjd,docmc(3)
     .     ,  dt,dxp,dyp                 
     .     ,  tstop,sec_to_rename,kpos2(3),kvel2(3)
     .     ,  kepoch2,yawrate,nadangd,yatt 

* RWK MOD 201216: Added from GPT3 call but not yet used
      real*8  Tm,la,Gn_h,Ge_h,Gn_w,Ge_w
    
      real*8 antpwr  ! Antenna Transmit power (W)

c     Function
       integer*4 nblen 
       real*8 timdif

      DATA TWOPI/6.283185307179586D0/

      REAL*4 swverx      
                         
      logical debug/.false./
      data none /'none'/


      RADIAN=DATAN(1.D0)/45.D0
                         
      ioerr = 0
      simulation = .false.                    

c     print warnings for missing antenna info
      warnings = .true.

c---------------------------------------------------------------------------
              
c  Read and Print the Observation Header Information
                     
      IF(upper1(obfiln(1:1)).eq.upper1('C') ) THEN

         call chdred( iyr,imo,iday,ihr,imin,sec)
         call check_y2k(iyr)         

c
c     need to modify xhdred to read the full 16 char station name from the xfile. McClusky 951111
      ELSEIF(upper1( obfiln(1:1)).EQ.upper1('X') ) THEN
         CALL XHDRED(  IOBS,IPRNT,ISCRN
     1,                NEPOCH,INTER,ircint,isessn
     2,                MTIME,iyr,imo,iday,ihr,imin,sec
     3,                NCHAN,ISCHAN,satnam
     4,                NDAT,DATTYP,rxobtyp,LAMBDA
     4,                offarp,sitnam16,rcvrswx,swverx,antcodx
     5,                rctypex,rcvrsnx,anttypx,antsnx
     6,                LATFLAG,LATD,LATM,SECLAT
     7,                LONFLAG,LOND,LONM,SECLON
     8,                radiusx,NTEXT,TEXT,gnssx )    
  
cd         print *,'SETUP xhdred iyr imo iday ',iyr,imo,iday
         call check_y2k(iyr) 

c        Check GNSS code from X-file against MODEL batch input
         if( lower1(gnssx).ne.lower1(gnss) ) then
            write(message,'(a,a1,a,a1,a)') 'GNSS code on X-file (',gnssx
     .           ,') does not match MODEL batch file (',gnss,')'
           call report_stat('WARNING','MODEL','setup',' ',message,0)
         endif
c
         nslip=1
         islip(1) = 0
         islpst(1)= 0
c        extra c-file variables not currently used so need to be initialized for output
         niextra = 1
         iextra(1) = 0
         nrextra = 1
         rextra(1) = 0.d0
         ncextra = 1
         cextra(1) = ' '

      ELSEIF( upper1(obfiln(1:1)).eq.upper1('S') ) then   
         call simred( iyr,imo,iday,ihr,imin,sec )                          
         simulation = .true.
         rcvrswx = ' '
         rcvrsn = ' '
         antsn = ' ' 
         niextra = 1
         iextra(1) = 0
         nrextra = 1
         rextra(1) = 0.d0
         ncextra = 1
         cextra(1) = ' '
         nslip=1
         islip(1) = 0
         islpst(1)= 0         

      ELSE
         call report_stat('FATAL','MODEL','setup',obfiln(1:1),
     1   'Observation file type unknown',0)
      ENDIF

c
c**   Temporary for backwards compatibility
      if( ircint.eq.0 ) ircint=inter

c     Convert the start time to Julian Day and seconds of day for processing
c     and to decimal year for reading an apr file with velocities
c     Compute the stop time
      
      jd0 = julday(imo,iday,iyr)
      t0  = ihr*3600.d0 + imin*60.d0 + sec
cd      print *,'SETUP jd0 t0 ',jd0,t0 
      if( mtime.eq.1 ) then
c          X- or C-file is UTC: convert to GPST
           utcoff = taiutc(jd0) - 19.d0
           call timinc(jd0,t0,utcoff)
      elseif (mtime.ne.2 ) then  
        write(message,'(a,i3)') 'mtime neither 1 (UTC) nor 2 (GPST)'
     .                        , mtime
        call report_stat('FATAL','MODEL','setup',' ',message,0)
      endif   
      call jd_to_decyrs( dfloat(jd0) + t0/86400.d0, t0_yrs )
      call dayjul( jd0,iyr,idoy )
      call monday( idoy,imo,iday,iyr)
      jdend = jd0
      tend  = t0                   
      ispan = inter*(nepoch-1) 
      span = dble(ispan)
      call timinc(jdend,tend,span)   

c     Check the X- or C-file header to see if the RINEX file was run
c     through FIXASH to alter its time tags
              
      fixash = .false.
      do i=1,ntext    
c        test each line of the comments for the FIXASH id
         istr = mchkey(text(i),'offset fixed',80,12)    
         if( istr.gt.0 ) then
            fixash = .true.
         endif
      enddo    
      if( fixash ) then
         call report_stat('WARNING','MODEL','setup',' '
     .    ,'ASHL12 RINEX files has clock corrections already applied',0)
         call report_stat('WARNING','MODEL','setup',' '
     .     ,' --do not reapply in MODEL',0) 
      endif


c----------------------------------------------------------------------

c  Read the Header Information from the Ephemeris (T-) File
               
      write(iprnt,'(/)')                   
c     set the input frame to 'INERT' so that THDRED will assign it
      frame = 'INERT' 
      call thdred ( iut,iscrn,iprnt,ntsat,gnsst,itsat,satnam
     .            , jdbt,tbt,jdft,tft,sdelt,nepcht,jdet,tet
     .            , norbprm,tsatic,nintrs,icsnamt
     .            , precmod,nutmod,gravmod,frame,srpmod
     .            , eradmod,antradmod )    
cd       print *,'SETUP aft THDRED frame ',frame  
      if( lower1(gnsst).ne.lower1(gnss) ) then
         write(message,'(a,a1,a,a1,a)') 'GNSS code on t-file (',gnsst
     .      ,'differs from MODEL input (',gnss,')' 
          call report_stat('WARNING','MODEL','setup',' ',message,0)
      endif      
      do i=1,nchan               
         sv_missing(i) = .false.
         ICHECK = IARRAY(ISCHAN(I),ITSAT,NTSAT)
         IF( ICHECK.EQ.0 ) THEN
c           WRITE(ISCRN,30) I,ISCHAN(I),(ITSAT(J),J=1,NTSAT)
           WRITE(IPRNT,30) I,ISCHAN(I),(ITSAT(J),J=1,NTSAT)
           WRITE(message,31) I,ISCHAN(I),(ITSAT(J),J=1,NTSAT)
30         FORMAT(/,1X,'Satellite in channel ',I2,' (PRN=',I2,') not on'
     1        ,' T-File',/,1X,'T-File satellites are PRN numbers ',45I3)
31         FORMAT('Satellite in channel ',I2,' (PRN=',I2,') not on'
     1        ,' T-File; T-File satellites are PRN numbers ',45I3)
           call report_stat('WARNING','MODEL','setup',' ',message,0)
           sv_missing(i) = .true.
         ENDIF
      enddo          

      do j=1,nchan
         ISAT= IARRAY(ISCHAN(J),ITSAT,NTSAT)
         do i=1,norbprm
            saticst(i,j) = tsatic(i,isat)
         enddo
      enddo                                            
c     if using, check the nutation used for MODEL versus the t-file
      if( .not.fcheck('nbody') ) then 
        read(inut,'(a)',iostat=ioerr) nut_header
        if( ioerr.ne.0) call report_stat('FATAL','MODEL','setup',' '
     .   ,'Error reading first line of nutation file to get model',0)
        rewind(inut)
* MOD TAH 200303: Changed default IAU0A which is fixdrv default. (Override 
*        sestbl. if needed
* MOD TAH 200304: If nutabl. then IAU00
        if(mchkey(nut_header,'IAU20',80,5).gt.0 ) then
          nutmod1='IAU00'
        else
          nutmod1='IAU80'
        endif  
      else
*         MOD TAH/SCM: option with IAUIAU0A, IAU0C and IAU06 precmod
          if ( precmod .eq. 'IAU76' ) then 
* MOD TAH 200303: Changed default IAU0A which is fixdrv default.
* MOD TAH 200304: Change back to IAU00 for IAU76 precession model.
           nutmod1 = 'IAU00'
          else 
            nutmod1 = precmod
          endif
      endif     
      if( nutmod1.ne.nutmod ) then   
        write(message,'(a,a5,a,a5,a)') 'Available nutation model ('
     .     ,nutmod1,') differs from t-file (',nutmod,'('  
        call report_stat('FATAL','MODEL','setup',' ',message,0)
        nutmod = nutmod1 
      endif
                     

  
c----------------------------------------------------------------------------

c  Print the header for the models and controls for this run

          write(iprnt,33)
33    format(///,'---------------------------------------------------',/
     1      ,1x,'** INPUT MODELING INFORMATION **',/) 

      if( simulation ) then
         write(iprnt,'(/,a,4f5.2)') 
     .    'Observaton noise level for simulation (L1 L2 P1 P2 in cyc): '
     .    ,(noise(i),i=1,ndat)
         write(iprnt,'(a,4f12.4)')
     .    'Station displacement for simulation (N E U in m): ',dispneu
      endif
                            

c---------------------------------------------------------------------------------

c  Read the receiver-clock (I-)to get the clock polynomials

      if( iui.gt.0 ) then 
cd        print *,'calling READI iyr idoy isessn jd0 t0 sitecd '
cd     .         ,               iyr,idoy,isessn,jd0,t0,sitecd
        call readi( iui,iscrn,iprnt,iyr,idoy,isessn,sitecd,jd0,t0
     .            ,  clkepc,clkrat,clkacc,clkcub )
      else
        clkepc = 0.d0
        clkrat = 0.d0
        clkacc = 0.d0 
        clkcub = 0.d0
      endif
c**   temporary until C-file format updated  (so also extra(2) below) 
c***  this removed rwk 050129: not needed since 'clock' added, 1999?
c**     store cubic part of clock term in extra array for autcln
c**      extra(1) = clkcub        
                                 


c---------------------------------------------------------------------------------

c   Read the station.info and antmod.dat files to get the antenna offsets

      inquire( file='station.info', exist=stnfo_exist )
      if( stnfo_exist ) then     
        call even_minute(iyr,idoy,t0,nyr,ndoy,nsod)  
* MOD TAH 200203: Added AntDAZ to list of values from station.info
        call rstnfo( istnfo, sitecd, nyr, ndoy,nsod,ispan
     .             , sitnam16, anth, antn, ante, antdaz, rcvcod, antcods     
     .             , htcod, radome_in, swver, rcvers, rcvrsn, antsn
     .             , kstarts, kstops )
cd        print *,'SETUP aft rtnfo  kstarts kstops ',kstarts,kstops
c       Abort if missing receiver or antenna (warnings in other programs)
        message = ' ' 
        if( rcvcod.eq.'      ')  write(message,'(a,a4,i5,i4)')
     .      'Station.info missing receiver for ',sitecd,nyr,ndoy
        if( antcods.eq.'      ')  write(message,'(a,a4,i5,i4)') 
     .      'Station.info missing antenna for ',sitecd,nyr,ndoy
        if(message(1:2).ne.'  ' )  then   
           write(iprnt,'(/,a)') message
           call report_stat('FATAL','MODEL','setup',' ',message,0)
           write(iprnt,'(/,a)') message
        endif   
c       rwk 100909: this dangerous for getting PCVs if there is a change mid-session, so leave blank
c        if( antsn.eq.'                    ') antsn = antsnx
c        if( rcvrsn.eq.'                    ') rcvrsn = rcvrsnx

c       Overwrite xfile or cfile full station name with name from station.info file...
        sitnam(1:16) = sitnam16
        sitnam(17:32) = ' ' 
                              

c       Check compatibility of receiver software versions
        if( .not.simulation ) then
          if( abs(swver-swverx).gt.0.0001 ) then
            write(message,'(a,f5.2,a,f5.2,2a)') 
     .           'Rcvr firmware from X or C-file (',swverx
     .          ,') does not match station.info (',swver,')---'
     .          ,'station.info value written on C- and H-file'
            call report_stat('WARNING','MODEL','setup',' ',message,0)
            write(iprnt,'(/,a)') message
          endif 
        endif
c       If blank, set the receiver software 3-character code from station.info
        if( rcvrswx.eq.'   ' ) then
            rcvrswx = rcvcod(1:3)
          if( rcvrswx.eq.'TI4' ) then
            if( abs(swver-1.11).le..05 ) then
                   rcvrswx = 'ROM'
            elseif( swver.ge.-1. .and. swver.le.4.0 ) then
                   rcvrswx = 'GES'
            elseif( swver.gt.4.0 .and. swver.le.6.0 ) then
                   rcvrswx = 'COR'
            endif
          endif
        endif
                                                      
c     Translate the station.info antenna code to a unique GAMIT code 
c       (done also in hisub and read_rcvant)
        call ant_alias(antcods,antcod)

c       Convert the raw antenna offset position to antenna reference point (ARP)
         call hisub( ihi,anth,antn,ante,antcod
     .               , htcod,sitecd,iyr,idoy,isessn,offarp,warnings) 
      else
        write(iprnt,'(//,1x,2a)')'MODEL warning: No station.info file;'
     .                ,' antenna offsets taken from the X-file'
        call report_stat('WARNING','MODEL','setup',' ',
     .          'No station.info file;  antenna offsets from X-file',0)
        
      endif


c     Get the official receiver and antenna names from rcvant.dat translation table.... 

      call read_rcvant(1,1,antcod,anttyp,radome_in,rcvcod,rctype,pcncod)
      call read_rcvant(1,2,antcod,anttyp,radome_in,rcvcod,rctype,pcncod)
      if( debug ) 
     .  print *,'after READ_RCVANT radome_in anttyp ',radome_in,anttyp  

      
c-----------------------------------------------------------------------

c  Write the geodetic datum, tide, EOP, and loading models to the p-file

c     set datum to be geocentric, but get WGS84 parameters for atmospheric correction  
      datum = '     '
      call read_gdatum( iud,datum,semi,finv,shft
     .                , shftdot,rot,rotdot,scale,scaledot,gepoch ) 
      write(iprnt,'(/,a,f12.3,f14.9)') 
     .    'Ellipsoid parameters tropo corrections (a 1/f): ',semi,finv  
      atidemod = '        '
* MDO TAH 201019: Added 64 to ietide for IERS2020 secular (mean) pole.
      write(iprnt,125)ietide,isptide
  125 format(/,'Tidal corrections applied are:',i4,/,
     a'The value is decoded as:'/,
     1'                      1 = solid earth tides',/,
     2'                      2 = frequency dependant K1 model'/,
     3'                      4 = Pole tide applied to zero mean pole',/,
     4'                      8 = ocean tides',/,      
     4'                     16 = Pole tide applied to IERS2010 ',
     .                           'mean pole',/,
     4'                     32 = atmospheric tides',/,
     4'                     64 = Pole tide applied to IERS2020 ',
     .                           'mean pole',/,
     5'The corrections are applied as per the binary option'/,
     6'That is a value of 31 implies all corrections',//,
     7'Short period earth orientation corrections applied are: ',i4,/,
     8'The value is decoded as:',/,
     9'                      0 = no short period corrections',/,
     b'                      1 = short period pole corrections',/,
     c'                      2 = short period ut1 corrections',/,
     d'                      4 = Ray model (IERS96)',/,
     e'                      8 = IERS 2010 model (IERS10)',/,
     .'                     16 = UT1 Libration term',/,
     .'                     32 = Gipson VLBI model',/,
     .'                     64 = Desai and Sibois Repro3 model',/,
     f'The corrections are applied as per the binary option') 
* MOD TAH 200505: Updated the model designation for new models
      if( kbit(isptide,7) ) then 
        speopmod = 'D&S2016'
      elseif( kbit(isptide,6) ) then 
        speopmod = 'GIPSON17'
      elseif( kbit(isptide,4) ) then 
        speopmod = 'IERS10'
      else
        speopmod = 'IERS96'  
      endif

      if(etidemod.eq.'        ') then
         write(iprnt,'(/,a)') 
     .     'Earth tide model not specified, using old IERS 1992 default'   
            call report_stat('WARNING','MODEL','setup',' ',
     .'Earth tide model not in batch file, using old IERS92 default',0)
         etidemod = 'IERS92  '
      else
         write(iprnt,'(/,a,a6)') 'Solid Earth tide model is '
     .            ,etidemod
      endif    
* MOD TAH 202018: Replaced with code to use actual values for mean pole
*     Work backwards though bits so that latest selected model is used.                   
C     if( kbit(ietide,3) ) then
C       if( kbit(ietide,5) ) then   
C         dt = (dfloat(jd0)-2451544.5d0)/365.25d0   ! time in years
c         fraction of day ignored since change is 1e-5 arcsec/day
c         Linear trend from IERS Conventions 2000.
C         dxp = 0.054d0+0.00083d0*dt
C         dyp = 0.357d0+0.00395d0*dt
C         write(iprnt,'(/a,f7.4,a,f7.4,a)') 
C    .     'Mean pole tide removed is Xp ',dxp,'  Yp ',dyp,' arcsec' 
C       else
C         write(iprnt,'(/a)') 'Mean pole tide not removed'
C       endif
C     else
C       write(iprnt,'(/a)') 'Pole tide not applied'
C     endif  
*     Go from latest to earlier models: Added 0,5 to PEP JD to make
*     real JD and convert integer to real*8.
      dxp = 0
      dyp = 0
      if( kbit(ietide,7) ) then  ! IERS 2020 model on
         call mean_pole(jd0+0.5d0,'IERS20',dxp, dyp)
         write(iprnt,210) 'IERS20',dxp, dyp
 210     format('Pole tide applied to ',a,' MXP/MYP ',2F8.3,' mas')
      elseif( kbit(ietide,5) ) then  ! IERS 2020 model on
         call mean_pole(jd0+0.5d0,'IERS10',dxp, dyp)
         write(iprnt,210) 'IERS10',dxp, dyp
      elseif ( kbit(ietide,3) ) then  ! Zero pole
         write(iprnt,210) 'Zero',dxp, dyp
      else
        write(iprnt,'(/a)') 'Pole tide not applied'
      endif 
*     Now convert back to arc seconds that GANIT wants
      dxp = dxp / 1000.d0
      dyp = dyp / 1000.d0

c--------------------------------------------------------------------------------------
 
c Read the u-file to get empirical models for ocean tidal loading, atmospheric loading,
c atmospheric tidal loading, met values at the site, and mapping function coefficients
             
      if( iuu.gt.0 ) call readu   
cd          print *,'SETUP lotl latml latl lmet lmap '
cd     .          ,lotl,latml,latl,lmet,lmap
cd           print *,' ot bit atmlflg at bit '
cd     .         ,kbit(ietide,4),atmlflg,kbit(ietide,6)  
      if( kbit(ietide,4) ) then
        if( lotl ) then       
          write(iprnt,'(/,a,a8,1x,54(1x,a3))') 
     .        'Ocean tidal loading model is '
     .        ,otidemod,(otlwaves(i),i=1,notl)
          if( otidemod(8:8).eq.'E') then
            write(iprnt,'(a)') '  CM - CE correction to be applied'  
c           read in the coefficients (saved in otlcmc; docmc not calculated)
            call otlcmc( jd0,t0,otidemod,notl,1,docmc )  
          endif
        else
         call report_stat('FATAL','MODEL','setup'
     .      ,' ', 'Ocean tides requested but not on  u-file',0) 
        endif
      else
c       if tides on u-file but not applied, blank out the mode name for p- and h-file
        if( otidemod.ne.'        ') then
           otidemod = ' '   
           write(iprnt,'(/,a)') 
     .         'Ocean tidal loading model read but not used'
           call report_stat('WARNING','MODEL','setup'
     .      ,' ', 'Ocean tidal loading model read but not used',0)
         endif
      endif
      if( atmlflg.eq.'Y' ) then
        if( latml ) then
          write(iprnt,'(/,a,a8)') 'Atmospheric loading model is '
     .            ,atmlmod       
        else
          call report_stat('FATAL','MODEL','setup'
     .      ,' ', 'Atmospheric loading requested but not on  u-file',0)
         endif 
      else
        atmlmod = '        '
      endif
      if( kbit(ietide,6) ) then
         if( latl ) then
           write(iprnt,'(/,a,a8)') 'Atmospheric tidal model is '
     .            ,atidemod       
         else
           call report_stat('FATAL','MODEL','setup'
     .     ,' ', 'Atmospheric tides requested but not on u-file',0)   
         endif
      else
        atidemod = '       '
      endif
      hydrolmod = '        '
      if(hydrlflg.eq.'Y') then
          call report_stat('WARNING','MODEL','setup',' '
     .     , 'Hydrological loading requested but not yet coded',0) 
          write(iprnt,'(/,a)') 
     .     'Hydrological loading requested but not yet coded '
      endif

c  Temporary: zero out the average atm loading values (N E U)
      do i=1,3
        atmlavg(i) = 0.d0
      enddo

          
c  Temporary: zero out the average hydrological loading values (N E U)
      do i=1,3
        hydrolavg(i) = 0.d0
      enddo

c
c-----------------------------------------------------------------------------

c  Get the SV antenna/body-types from svnav.dat 
      do i=1,nchan 
* MOD TAH 190702: Added antpwr to snav_read call
        call svnav_read( -1,iyr,idoy,ihr,imin,gnss,ischan(i),isvn(i), 
     .       frqchn(i),svantbody(i),sbmass,yawbias,yawrate, antpwr,
     .       svnstart,svnstop )

* MOD TAH 180310: Added warning if frequency not found.
        if( frqchn(i).eq.99 ) then
           write(*,320) i, ischan(i),isvn(i), svantbody(i)
 320       format('**WARNING** No GLONASS Frequency for SV ',I2,
     .            ' PRN ',I2,' SVN ',I3,' Body ',a)
        end if
      enddo
 
c-------------------------------------------------------------------------------
c MOD TAH 180309: Moved block from above block above becuase GLONASS frequencies
c     not read yet.
          
c  Set the (two, for now) satellite frequencies (same for all channels except GLONASS)
c    --check to make sure they match the observables on the x- or c-file
                             
cd      print *,'DEBUG gnss rxobtyp ',gnss,rxobtyp 
* MOD TAH 200512: Added selections here based on -lfreq selection of lower GNSS
*     frequency 
      do i=1,nchan
        if( gnss.eq.'G' )then
          fL1(i) = gps_f1
          atxfrq(1) = 'G01'
* MOD TAH 200512: Check rxobtyp(2) to see what F2 to use
          if( rxobtyp(2)(1:2).eq.'L2'.or.rxobtyp(2)(1:2).eq.'  ') then
             fL2(i) = gps_f2
             atxfrq(2) = 'G02'
             atxfrq(3) = 'G02'
          else if( rxobtyp(2)(1:2).eq.'L5' ) then
             fL2(i) = gps_f5
             atxfrq(2) = 'G05'
             atxfrq(3) = 'G02'
          else
             call report_stat('WARNING','MODEL','setup',' '
     .          ,'GPS frequencies on x-/c-file file not L1,L2/L5',0 )
          endif
        elseif( gnss.eq.'R') then 
          fL1(i) = glonass_f1 + frqchn(i)*glonass_df1
          fL2(i) = glonass_f2 + frqchn(i)*glonass_df2
          if( rxobtyp(1)(1:2).ne.'L1'.or.rxobtyp(2)(1:2).ne.'L2' ) 
     .      call report_stat('WARNING','MODEL','setup',' '
     .       ,'Glonass frequencies on x-/c-file file not L1 and L2',0 )  
          atxfrq(1) = 'R01'
          atxfrq(2) = 'R02'
        elseif( gnss.eq.'C') then
* MOD TAH/RWK: 200223: Modified frequency assignments with preference to L6/C6
*         Allow different frequencies depending on available observables
*         (Re-ordered logic as well).
*         Change L1 to L2 for old RINEX 2 files that used L1 for L2.
*         (sel_obtyp will take L2 over L1 for BDS 3 satellites transmitting
*         true L1.
          if( rxobtyp(1)(1:2).eq.'L1') rxobtyp(1)(1:2) = 'L2'
*         Now see what we have: L2/C2 should always be present
          fL1(i) = beidou_f2
          atxfrq(1) = 'C02'
*         Assign second frequnecy
          if( rxobtyp(2)(1:2).eq.'L6' ) then
             fL2(i) = beidou_f6      
             atxfrq(2) = 'C06'
             atxfrq(3) = 'C07'   ! Most common (after C01?)
* MOD MAF 210729: Added L5 (B2a; B2b = L7) lower frequency choice
          elseif( rxobtyp(2)(1:2).eq.'L5' )    then
             fL2(i) = beidou_f5      
             atxfrq(2) = 'C05'
             atxfrq(3) = 'C07'
          elseif( rxobtyp(2)(1:2).eq.'L7' )    then
             fL2(i) = beidou_f7      
             atxfrq(2) = 'C07'
             atxfrq(3) = 'C06'
          else   ! No match; warn user
* MOD TAH 200620: Updated message to be more accurate
             call report_stat('WARNING','MODEL','setup',' '
     .       ,'Beidou missing lower frequency (try -lfreq option)',0 )

          endif
       elseif( gnss.eq.'E' ) then
          fL1(i) = galileo_f1
          atxfrq(1) = 'E01'
* MOD TAH 200512: Now check frequency choices.
          if( rxobtyp(2)(1:2).eq.'L5' ) then
             fL2(i) = galileo_f5
             atxfrq(2) = 'E05'
             atxfrq(3) = 'E07'
          elseif( rxobtyp(2)(1:2).eq.'L6' ) then
             fL2(i) = galileo_f6
             atxfrq(2) = 'E06'
             atxfrq(3) = 'E05'
          elseif( rxobtyp(2)(1:2).eq.'L7' ) then
             fL2(i) = galileo_f7
             atxfrq(2) = 'E07'
             atxfrq(3) = 'E05'
          elseif( rxobtyp(2)(1:2).eq.'L8' ) then
             fL2(i) = galileo_f8
             atxfrq(2) = 'E08'
             atxfrq(3) = 'E07'
          else
            call report_stat('WARNING','MODEL','setup',' ',
     .      'Galileo missing lower frequency (try -lfreq option)',0)
          endif

       elseif( gnss.eq.'J' ) then
          call report_stat('FATAL','MODEL','setup',' '
     .       ,'QZSS not yet supported',0 )
       elseif( gnss.eq.'I' ) then  
          fL1(i) = irnss_f9
          fL2(i) = irnss_f5
          if( rxobtyp(1)(1:2).ne.'L9'.or.rxobtyp(2)(1:2).ne.'L5' ) 
     .      call report_stat('WARNING','MODEL','setup',' '
     .       ,'IRNSS frequencies on x-/c-file file not I9 and E5',0 )
          atxfrq(1) = 'I09'
          atxfrq(2) = 'I05'
       else
         call report_stat('FATAL','MODEL','setup',' '
     .                  ,'GNSS not recognized',0)
        endif   
      enddo                
      write(iprnt,'(/,a,2(2x,a3,4pd16.7))') 'Observation frequencies: '
     .            ,atxfrq(1),fL1(1),atxfrq(2),fL2(1)

c----------------------------------------------------------------------------             

c Get the the phase-center offset  and variations from the ANTEX file
               
      if(debug) print *,'SETUP calling get_svantpcv jd0 yatt',jd0,yatt
      do i=1,nchan                      
        first = .true.
        found_svant(i) = .false.
        j=0    
        do jj = 0,14,1
          nadangd = float(jj)
          j=j+1              
          call get_svantpcv( jd0,i,nadangd,yatt,first
     .                     , found_svant(i),svoffl1,svoffl2
     .                     , svcorrl1,svcorrl2 )
        enddo     
c       convert the offsets from mm to m
        do j=1,3
          svantdx(j,1,i) = svoffl1(j)/1.d3
          svantdx(j,2,i) = svoffl2(j)/1.d3
        enddo 
      enddo                                      
c     convert the offsets from mm to m
     
                                       

c------------------------------------------------------------------------------
        
c  Print the SV antenna models and offsets
    
      write(iprnt,'(/,a)') 
     .   'SV antenna models '
      write(iprnt,'(2a)') 
     .  ' PRN  SVN  AntBody                Requested  Used '
     .  ,' Offsets L1  DX DY DZ  (m)  Offsets L2             Source'  
      do i=1,nchan   
        if( found_svant(i) ) then 
cd           print *,'svantmod(1) ',svantmod(1)
           write(iprnt,'(2i4,3x,a20,5x,a4,5x,a4,2x,6f8.4,2x,a)') 
     .         ischan(i),isvn(i),svantbody(i),svantmod_in,svantmod(i)
     .       , (svantdx(j,1,i),j=1,3),(svantdx(j,2,i),j=1,3)
     .       , svantmod_snx(i)
        else  
          write(message,'(a,a1,i3.3,a)') 'SV antenna offsets for SVN '
     .             ,gnss,isvn(i),' not found in antmod.dat '
          write(iprnt,'(/,a)') message
          call report_stat('FATAL','MODEL','setup',' ',message,0)
        endif
c       change units to km for processing : rwk 050208: no, keep in meters
c       for the c-file, make change to km only in svant.f.
cd       print *,'prn svantbody svantdx '
cc     .        ,ischan(i),svantbody(i),(svantdx(j,i),j=1,3) 
      enddo 
                        
 
c-----------------------------------------------------------------------------------

c  Get the receiving antenna phase center model

      write(iprnt,'(/,a,a6,/,a,a6,a,a15,a,a20,a,a6)')
     .    'Antenna type from station.info: ',antcods
     .   ,'Standard code: ',antcod,'  Full name: ',anttyp(1:15)
     .   ,'  Serial No: ',antsn
     .   ,'  Radome: ',radome_in    
      write(iprnt,'(a,a5,7f8.4)')
     .'Offset from monument (code, U,N,E, DHARP): '
     .    ,htcod,anth,antn,ante,offarp(1)   
      if( radome_in(1:4).eq.'UNKN') then    
         call report_stat('WARNING','MODEL','setup',' '
     .    ,'Radome status unknown, use antenna alone',0)
         radome_in = 'NONE '
      endif            
      call get_antinfo(debug) 
                  
c--------------------------------------------------------------------------------

c  Read the L-file and compute the Earth-fixed coordinates
           
c     set defaults (in model.h)  on valid times for site from eq/rename file 
      kstartr(1) = 1900
      kstopr(1) = 2100   
      kstopr(2) = 1
      do i=3,5
        kstartr(i) = 0
        kstopr(i) = 0
      enddo      
c     get the midpt times for the site-coordinate epoch
      jdmidpt = jd0          
      tmidpt = t0
      call timinc(jdmidpt,tmidpt,span/2.d0)

c     call lread at the beginning of the session and check for renames or EQs within the session
      xjd = dfloat(jd0) - 0.5d0 + t0/86400.d0      
      call jd_to_decyrs( xjd,site_epoch)   
cd      print *,'SETUP calling lread for start ',xjd,site_epoch
      site1 = sitecd
      call lread( site1,site_epoch ) 
cd      print *,'SETUP lread returned start ',asite,kpos
      nsod = kstopr(3)*3600 + kstopr(4)*60 + kstopr(5)
      tstop = int(nsod)
      call yds_to_jd( kstopr(1),kstopr(2),nsod,xjd )
      jdstop = int(xjd+0.5d0)
      sec_to_rename = timdif(jdstop,tstop,jd0,t0)   
cd      print *,'kstopr jdstop tstop sec_to_rename '
cd     .       , kstopr,jdstop,tstop,sec_to_rename       
c     call lread at the midpoint of the session (will always be right for the longer half)
      xjd = dfloat(jdmidpt) - 0.5d0 + tmidpt/86400.d0
      call jd_to_decyrs( xjd,site_epoch )   
cd      print *,'SETUP calling lread at midpt ',xjd,site_epoch
      call lread ( site1,site_epoch ) 
cd      print *,'SETUP lread returned midpt asite,kstopr,kpos '
cd     .   ,asite,kstopr,kpos
      do i=1,3
        pos(i) = kpos(i) + kvel(i)*(site_epoch-kepoch0)
      enddo  
c     see which data to delete
      if( sec_to_rename.gt.(span-60.d0) ) then
c       l-file stop epoch is beyond end of session: use all the data
        nepoch_use(1) = 1
        nepoch_use(2) = nepoch  
cd        print *,'DEBUG no break site_epoch asite pos,nepoch_use '
cd     .     ,site_epoch,asite,pos,nepoch_use
      elseif( sec_to_rename.le.(span/2.d0) ) then
c       l-file stop epoch is in the 1st half of the session: delete before the break
        nepoch_use(1) = sec_to_rename/dfloat(inter) + 1
        nepoch_use(2) = nepoch
cd        print *,'DEBUG break before midpt site_epoch asite pos ',
cd     .     ' nepoch_use ',site_epoch,asite,pos,nepoch_use
      elseif( sec_to_rename.gt.(span/2.d0) ) then
c       l-file stop epoch is in the 2nd half of the session: delete after the break    
        nepoch_use(1) = 1
        nepoch_use(2) = sec_to_rename/dfloat(inter)  
cd        print *,'DEBUG break after midpt site_epoch asite pos '
cd     .   ,' nepoch_use ',site_epoch,asite,pos,nepoch_use
      endif
      call xyz2sph(pos,latr_sph,lonr,radius)     
cd      print *,'SETUP kpos kvel pos ',kpos,kvel,pos
cd      print *,'  latr_sph lonr radius ',latr_sph,lonr,radius
      
c     L-file values stored in model.h  (see documentation in update_coords.f) 
      write(iprnt,'(/,a)') 'Initial coordinates'
      if( kfflg.eq.0 ) then
c       L-file was old-style, spherical  
        call raddms( latr_sph,latflag,dlat,mlat,slat )
        if( latflag.eq.'-') then
          latflag = 'S'
        else
          latflag = 'N'
        endif
        call raddms( lonr,lonflag,dlon,mlon,slon )  

        if( lonflag.eq.'-' ) then
          lonflag = 'W'
        else
          lonflag = 'E'
        endif      
cd        print *,'From XHDRED latd,latm,seclat,lond,lonm,seclon ' 
cd     .         ,             latd,latm,seclat,lond,lonm,seclon
cd        print *,'  dlat,mlat,slat,dlon,mlon,slon '
cd     .         ,   dlat,mlat,slat,dlon,mlon,slon
        write(iprnt,'(a,a1,i4,i3,f9.5,4x,a1,i4,i3,f9.5,f15.4)')
     .      ' L-file (lat lon rad) : '
     .      ,latflag,dlat,mlat,slat,lonflag,dlon,mlon,slon,radius
      else
c       L-file was new-style, Cartesian
        write(iprnt,'(a,a8,3f14.4,2x,3f9.4,f11.4)')
     .        ' L-file (m  m/yr)     : '
     .        ,asite,kpos,kvel,kepoch0  
        if( nepoch_use(1).ne.1.or.nepoch_use(2).ne.nepoch ) then
          write(iprnt,'(2a,2i6)') 
     .        '** WARNING: earthquake or rename during session,'
     .        ,' using only epochs ',nepoch_use(1),nepoch_use(2)
          call report_stat('WARNING','MODEL','setup',' '
     .        ,'L-file coordinate change during session: see p-file',0)
        endif
      endif
        write(iprnt,'(a,2f14.10,f14.4,2x,3f14.4)') 
     .      ' Converted (radians,m): '
     .      ,latr_sph,lonr,radius,(pos(i),i=1,3)   
      call cortyp( pos,offl1,offl2,semi,finv,shft
     .           , evec0,latr,height,sitepart0 )  
cd       print *,'SETUP aft CORTYP '
cd     .    , ' pos l1l2 semi shift evec0 latr height sitepart0 ',    
cd     . pos,offl1,offl2,semi,shft,evec0,latr,height,sitepart0
      write(iprnt,'(a,/,a,3f14.6,/,a,3f14.6)') 
     .    'Initial Earth-fixed antenna phase-center coordinates (km): '
     .   ,'  L1 ',(evec0(i,1),i=1,3)
     .   ,'  L2 ',(evec0(i,2),i=1,3)
cd      write(*,'(a,3f20.6)') 'DEBUG L1 ',(evec0(i,1)*1.d3,i=1,3)
cd      write(*,'(a,3f20.6)') 'DEBUG L2 ',(evec0(i,2)*1.d3,i=1,3)
            

c     compute cartesian coordinates (SIMVEC0) for the simulated displacement
      if( simulation ) call corsim( pos,offl1,offl2,semi,finv,shft
     .                            , dispneu,simvec0 )    
      
           
c     set a flag to increment positions with velocity during the session if the 
c     velocity exceeds 1 mm/day
      kvflg = 0
      if ( dsqrt(kvel(1)**2+kvel(2)**2+kvel(3)**2)/365.d0 .gt. .001d00 )
     .    kvflg =1 
      if( kvflg.le.0 ) then
         write(iprnt,'(a)') 
     .       'Velocity < 1 mm/d, no updates applied during sessions'
      else
         write(iprnt,'(a)') 
     .       'Velocity > 1 mm/d, coordinates updated during session'
      endif
           
                      

c--------------------------------------------------------------------------

c  Get the DCBs and write the receiver and firmware to the P-file

      write(iprnt,'(/,a,a6,a,a20,a,a20)') 
     .    'Receiver type from station.info: '
     .     ,rcvcod,'  Full name: ',rctype,'  Serial No: ',rcvrsn
      write(iprnt,'(a,f5.2)')  'Firmware: ',swver  
      write(iprnt,'(a,a3,1x,f5.2,a)')  '( MAKEX ',rcvrswx,swverx,' )'
      if( rcvcod.eq.'ASHL12' ) then                              
        if( fixash ) then
          write(iprnt,'(a)') 
     .  ' Codeless Ashtech, clock corrected in RINEX by 45 microseconds'   
        else
          write(iprnt,'(a)') 
     .  ' Codeless Ashtech, clock corrected in MODEL by 45 microseconds' 
        endif
      endif
      if( pcncod.eq.'P' ) then
        write(iprnt,'(2a)') 'Differential code biases (dcb) to be '
     .                    , 'applied to C1 and P2 '
      elseif( pcncod.eq.'C' ) then
        write(iprnt,'(2a)') 'Differential code biases (dcb) to be '
     .                    , 'applied to C1 ' 
      elseif (pcncod.eq.'N' ) then
        write(iprnt,'(2a)') 'No differential code biases (dcb) to be '
     .                    ,   'applied  ' 
      endif
      if( pcncod.eq.'P'.or.pcncod.eq.'C' ) then 
c       read either an old-style or new-sytle (Version 2.0) dcb.dat file
        read( iudcb,'(a)',iostat=ioerr ) line 
        if( line(11:21).eq.'Version 2.0' ) then 
c         call the routine one SV at a time             
          do i=1,nchan
            isat = ischan(i) 
c           the prndcb array will not needed after version 1 removed since
c           we'll use the same list as for the observations
            prndcb(i) = isat
            rewind(iudcb)
            call get_dcb2( iudcb,jd0,gnss,isvn(i),dcb(i) )
          enddo
          numdcb = nchan       
          do i=1,nchan
            write(iprnt,'(2i4,f9.3)') prndcb(i),isvn(i),dcb(i)
          enddo
        else
          call get_dcb( iudcb,jd0,numdcb,prndcb,dcb )  
          write(iprnt,'(a)') ' PRN  DCB (ns)'
          do i=1,numdcb
            write(iprnt,'(i4,f9.3)') prndcb(i),dcb(i)
          enddo
        endif
      endif


c-----------------------------------------------------------------------

c  How to handle clocks and yaw?
      
      if(klock.lt.1 .or. klock.gt.3) then 
        write(message,'(a,i3,a)') 
     .   'Klock parameter (=',klock,') must be between 1 and 3'
        call report_stat('FATAL','MODEL','setup',' ',message,0)
      endif
C       Check to see if the ephemeris span is long enough
c         1) for no yaw modelling case - need 5 epochs for interpolation
c         2) for new yaw modelling case - need an additional hour
      if(yfiln(1:5).eq.'     ')then
        testt= t0-tbt+(jd0-jdbt)*86400.D0-5.D0*sdelt
      else
        testt= t0-tbt+(jd0-jdbt)*86400.D0-5.d0*sdelt-3600.d0
      endif
      if (testt.LE.0.0D0) then
        call report_stat('FATAL','MODEL','setup',' ',
     .  'T-file tabular ephemeris starts too late.',0)
      endif
      if(yfiln(1:5).eq.'     ')then
        testt= tft-tend+(jdft-jdend)*86400.D0-5.D0*sdelt
      else
        testt= tft-tend+(jdft-jdend)*86400.D0-5.d0*sdelt-3600.d0
      endif
      if (testt.LE.0.0D0) then
        call report_stat('FATAL','MODEL','setup',' ',
     .  'T-file tabular ephemeris ends too early.',0)
      endif
        
c  Read the header information for the Y (yaw) file
        
c     Determine whether a yaw filename has been supplied to model. If so, read the
c     file. If not, set the unit number to zero to indicate no yaw modelling
      if (yfiln.eq.'          ') then
        iuy = 0 
        write(message,'(a)')  
     .  'No y-file in batch file: yaw modelling will not be implemented'
        call report_stat('STATUS','MODEL','setup',' ',message,0)
        write(iprnt,'(a)') message
      else
c       open yaw attitude file 
        open(iuy,file=yfiln,form='unformatted'
     .     ,access='sequential',status='old',iostat=ioerr)
        if(ioerr.ne.0)then 
          write(message,'(a,a)')
     .    'Error opening yaw table: ',yfiln
          call report_stat('FATAL','MODEL','setup',' ',message,ioerr)
        endif  
        write(message,'(a)')'Yaw modelling is implemented'
        call report_stat('STATUS','MODEL','setup',' ',message,0)
c       Need to read the header record gingerly to determine if this is an
c       First check on the naming convention
        if(yfiln(6:6).ne.'t') call report_stat('WARNING','MODEL','setup'
     .     ,yfiln,'Y-file 6th character not t, may be old-style file',0)
c       In a new-sytle file, the first variable will be an integer version number
        read(iuy,iostat=ioerr) nyversn
        if( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','setup',yfiln
     .      ,'Error reading version from y-file header--old-style file?'
     .      ,ioerr)  
        else
          if(nyversn.lt.1051.or.nyversn.gt.1061  ) then
            write(message,'(a,i4,a)') 'Incompatible y-file version ('
     .      ,nyversn,')  File: '
            call report_stat('FATAL','MODEL','setup',yfiln,message,0)     
          endif
          rewind(iuy)
          if( nyversn.eq.1051 ) then
            read(iuy) nyversn,ytfile,yaw_s,yaw_e,nyinter,nyepoch
     .          , nysat,(ysat(isat),ysvn(isat),yblk(isat),isat=1,nysat)       
          elseif( nyversn.eq.1061 ) then
            read(iuy) nyversn,ytfile,yaw_s,yaw_e,nyinter,nyepoch,nysat
     .           , (ysat(isat),ysvn(isat),svantbody(isat),isat=1,nysat)
          endif
                        
c  PT990505: check that the yfile interval is the same as the x-file interval.
c         If not, for now just fatally stop MODEL.   
          if(nyinter.ne.inter)then
            write(message,'(a,i4,a,a,i4,a)')
     .        'Interval on yfile (',nyinter,')'
     .        ,' is different to interval on observation file ('
     .        ,inter,').'
            call report_stat('WARNING','MODEL','setup',yfiln,message,0)
              write(message,'(a,a,a,a,i4,a)')
     .        ' Change session.info to match obs file interval and'
     .        ,' rerun FIXDRV or recreate ',yfiln,' with interval '
     .        ,inter,' seconds.'
              call report_stat('FATAL','MODEL','setup',yfiln,message,0)
           endif      
c  RWK140125: check that the y-file satellites are the same as the t-file
          if( nysat.ne.ntsat ) then
            write(message,'(a,i2,a,i2,a)') '# SVs on t-file (',ntsat
     .         ,') differs from y-file (',nysat,')' 
            call report_stat('FATAL','MODEL','setup',yfiln,message,0)
          endif
          do i=1,nysat
            if( ysat(i).ne.itsat(i) ) 
     .         call report_stat('FATAL','MODEL','setup',yfiln
     .                   , 't-file and y-file PRN #s do not match ',0)
          enddo 
c            possibly write the yaw-file header info to the p-file here
        endif
      endif  
                               

       
c   Read the satellite clock coefficients into storage from the J-file
c     secret back door for debugging:
c     Can avoid use of J-file by specifying NONE for name
      if(index(jfiln,none) .gt. 0) then
         write (iscrn,*) 'J-file is ignored'
         write (iprnt,*) 'J-file is ignored'
      else
         ICALL = 0     
         CALL READJ ( IUJ,ISCHAN,NCHAN,IDUM,IDUM,DUM,DUM,ICALL,
     .            IDUM,DUM,DUM,DUM,DUM,DUM)
      endif

c-------------------------------------------------------------------------------

C Determine the number of partials to be written on the C-file 

c     The T-file header has the number of orbit parameters for the force
c     model (norbprm) and the number of equations integrated (nintrs),
c     from which we can determine how many orbital parameter partials
c     were integrated (norbpart).  With position only, there are 3 equations
c     for the motion and 3 for each of the partials of motion wrt a parameter,
c     so the total number of equations is nintrs = 3*(1+norbpart).  If 
c     velocities are included, then the number doubles.   Note that only 
c     norbpart, the # of orbit partials is written onto the C-file since
c     CFMRG and SOLVE do care about the number of parameters if they are
c     not estimated.

c     The C-file is allowed to have either 1) station partials only 
c     (3 coordinates, 1 zenith delay, and 1 clock epoch), or 2) these
c     5 partials + for each satellite 'norbpart' orbit partials and 
c     3 SV antenna offsets + 6 Earth orientation parameters 

c     Case 1: no orbit partials
      if (nintrs.eq.3) then  
        norbpart = 0
        npart=5 
             
c     Case 2: 9 orbit partials (6 ICs + 3 non-gravitational force parameters) 
      elseif (nintrs.eq.30) then  
        if( norbprm.ne.9 ) call report_stat('FATAL','MODEL','setup',' '
     .     ,'# orbit parameters inconsistent with # T-file equations',0)
        norbpart = 9
        npart=23     

c     Case 3: 15 orbit partials (6 ICs + 9 prms of ECOM1 or ECOM2 models)
      elseif (nintrs.eq.48.or.nintrs.eq.96) then  
        if( norbprm.ne.15 ) call report_stat('FATAL','MODEL','setup',' '
     .     ,'# orbit parameters inconsistent with # T-file equations',0)
        norbpart = 15
        npart=29
             
c     Case 4: 19 orbit partials (6 ICS + 13 prms of ECOMC model)
      elseif( nintrs.eq.60.or.nintrs.eq.120 ) then 
        norbpart = 19
        npart = 33 
      else                                       
        write(message,'(a,i2,a)') 'Unknown number of partials ('
     .                             ,nintrs,') on T-file'
        call report_stat('FATAL','MODEL','setup',' ',message,0)
      endif

  
c-------------------------------------------------------------------------------

C        Initialize quantities for interpolation setup

C     NINDIM= 3*(1 + no. partials)*NTSAT
      NINDIM= 3*(1 + MAXORB)*MAXSAT
      NINTRP=NINTRS*NTSAT
      if (NINTRP.GT.NINDIM) then
        call report_stat('FATAL','MODEL','setup','nintrp',
     .  'Number of T-file quantities too large.',0)
      endif
 

c--------------------------------------------------------------------------------

c Select the atmospheric models          
            
c  Priority set by the 'metopts' values from the sestbl and written  on the b-file :
c      RNX : RINEX met file if available
c      UFL : U-file (global circulation model) if availible
c      GPT : Vienna harmonic model (subroutine gpt2 and gpt3 supersede gpt )
c      STP : Standard temperature and pressure (constant)
c      PTH : Pressure, temperature, humidity from sittbl. 
c      Note: The last 'metopts' value may be numerical, for humidity; if
c            negative values from the RINEX file, u-file, or GPT will
c            override it.
c                    
c     Initialize source flags (overall primary and individual defaults)
      metsrc = '   '
      zsrc = '   '
      psrc = '   '
      tsrc = '   '   
      wsrc = '   '
      lsrc = '   '
c     set sensor height (RINEX met) to a large number to indicate not available
      sensor_ht = 1.01d5
      write(iprnt,'(/,a)') 'Atmospheric models'   
      write(iprnt,'(a,f8.3,a)') ' Geodetic height = ',height*1.d3,' (m)'
c     check geodetic height for reasonableness to avoid pressure overflow
      if( height.lt.-0.5d0 .or. height.gt.50.d0 ) 
     .   call report_stat('FATAL','MODEL','setup',' '
     .  , 'Geodetic height unreasonable: check p- and l-files ',0)
     
c     Get initial values from GPT even if it's not requested since these will
c     always be available and are needed for lapse rate, for temperature if 
c     VMF1 grids are read, and for mapping funtions if VMF1 not used

      if( fcheck('gpt.grid') ) then         
        dmjd = (dfloat(jd0)+t0/86400.d0) - 2400001.d0  

* RWK MODS 201216: Open the GPT  grid file here instead of within the subroutine
*       to allow selection of GPT2 or GPT3 and get the grid spacing for GPT3 
*       from the file name  
        open(41,file='gpt.grid',iostat=ioerr)
        if( ioerr.ne.0 )  call report_stat('FATAL','MODEL','setup'
     .    ,'gpt.grid','Error opening GPT grid file: ',ioerr)
        read(41,'(a)') line
        if( line(3:4).ne.'gp' ) then
          call report_stat('FATAL','MODEL','setup','gpt.grid'
     .      ,'File name missing from first line of GPT grid file',0)
        else
          gpt_filnam = line(3:18)
        endif
        gptnam(1:2) = 'GP'
        gptnam(3:3) = gpt_filnam(4:4)
        gptnam(4:4) = gpt_filnam(6:6) 
        read(gptnam(4:4),'(i1)',iostat=ioerr) intdeg
        if(ioerr.ne.0) then
         write(message,'(a,a8)') 'Error reading GPT version from line: '
     .      ,line(1:8)
          call report_stat('FATAL','MODEL','setup'
     .      ,gpt_filnam(1:nblen(gpt_filnam)),message,ioerr) 
        elseif(intdeg.ne.1.and.intdeg.ne.5) then
          write(message,'(3a,i1,a)') 'Grid spacing from '
     .       ,gpt_filnam(1:nblen(gpt_filnam)),' (',intdeg,') ne 1 or 5 '
          call report_stat('FATAL','MODEL','setup'
     .      ,gpt_filnam(1:nblen(gpt_filnam)),message,ioerr) 
        endif                 
        if(debug) print *
     .    ,'SETUP gptnam dmjd latr lonr height pres0 temp0 lapse e0 ' 
     .    ,       gptnam, dmjd,latr,lonr,height,pres0,temp0,lapse,e0 
        if( gptnam(3:3).eq.'2') then 
          call gpt2( dmjd,latr,lonr,height*1.d3,1,0
     .             , pres0,temp0,lapse,e0,ah,aw,undu) 
        elseif( gptnam(3:3).eq.'3') then 
* MOD MAF 20210308: Multiplied height by 1000 to convert to m
c         call gpt3( 41,intdeg,dmjd,latr,lonr,height,0
          call gpt3( 41,intdeg,dmjd,latr,lonr,height*1.d3,0
     .             , pres0,temp0,lapse,tm,e0,ah,aw,la,undu
     .             , Gn_h,Ge_h,Gn_w,Ge_w )
        else
           call report_stat('FATAL','MODEL','setup'
     .      ,gpt_filnam(1:nblen(gpt_filnam))
     .      ,'GPT grid spacing neither 1 nor 5',0)
        endif
        if(debug) 
     .     print *,'GET_MET_SOURCE: '
     .     , 'gptnam,intdeg,dmjd latr lonr height pres0 temp0 lapse e0 '
     .     , gptnam,intdeg,dmjd,latr,lonr,height,pres0,temp0,lapse,e0  
        dryzen = gptnam
        wetzen = gptnam  
        call wpress(-1,wetvar0,e0,temp0) 
        lapse = -lapse
        lsrc = 'GPT'   
        close(unit=41)
      else
        call gpt(dmjd,latr,lonr,height*1.d3,pres0,temp0,undu )  
        lapse = 6.5d0
        lsrc = 'FIX'
      endif
               
c     Echo the requested input heirarchy
      write(iprnt,'(1x,a,a20)') 'Requested priority: ',metopts

c     See what's available   
c     try reading numerical values to see if P T H is given
      read(metopts,*,iostat=ioerr) pres0,temp0,humid0
      if( ioerr.eq.0 ) then
c        numerical values successfully read
         metsrc = 'PTH'
         psrc = 'PTH'
         tsrc = 'PTH'
         wsrc = 'PTH'
         wetvar0 = humid0/100.d0
         write(iprnt,'(a)') ' Using met values from sittbl.'  
         sealoc = 'S'
      else 
c       assume that the line contains an alphameric hierarchy
        read(metopts,'(5a4)',iostat=ioerr) (metsrc4(i),i=1,5)
        if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','setup',' '
     .     ,'Error reading met tokens from b-file',ioerr)   
        call get_met_source( metsrc4 )
cd        print *,'SETUP zsrc psrc tsrc sealoc wsrc '
cd     .         ,       zsrc,psrc,tsrc,sealoc,wsrc 
      endif
      if( metsrc.eq.'STP' ) then
        pres0 = 1013.25d0
        temp0 = 20.d0 
        humid0 = 0.5 
        lapse = 6.5d0
      endif 
                    
c     Write the selected values to the p-file
      write(iprnt,'(1x,a)') 'Source and values selected at start epoch'
c     ZHD or pressure        
      metsrcmod = ' ' 
      if( zsrc.ne.'   '.and.psrc.eq.'   ') then
        if( zsrc.eq.'UFL' ) then
          metsrcmod = mapmod
        elseif( zsrc.eq.'GPT' ) then
          metsrcmod = gptnam//'    '
        endif  
        write(iprnt,'(1x,a,2x,a3,1x,a8,f9.1,a)') 
     .          'ZHD      ',zsrc,metsrcmod,zhd0,' mm'
      elseif( zsrc.eq.'   '.and.psrc.ne.'   ') then 
        if( psrc.eq.'GPT' ) then
          metsrcmod = gptnam//'    '
        endif    
        if( sealoc.eq.'L' ) then
          write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .          'P local   ',psrc,metsrcmod,pres0,' hPa'
        elseif( sealoc.eq.'S' ) then
          write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .          'P sealevel',psrc,metsrcmod,pres0,' hPa'
        else      
          call report_stat('FATAL','MODEL','setup',' '
     .           ,'sealoc not set for pressure',0)  
        endif
      else
        call report_stat('FATAL','MODEL','setup',' '
     .         ,'Neither ZHD nor pressure given for dry delay',0)  
      endif                       
c     Temperature       
      metsrcmod = ' ' 
      if (tsrc.eq.'GPT' ) then
        metsrcmod = gptnam//'    '
      else
c       currently we have no temperature on the VMF1 grid files
        metsrcmod = ' '                                   
      endif
      if( sealoc.eq.'L' ) then   
        write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .     'T local   ',tsrc,metsrcmod,temp0,' C'
      elseif( sealoc.eq.'S' ) then
        write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .     'T sealevel',tsrc,metsrcmod,temp0,' C'
      endif   
c     Water vapor       
      metsrcmod = ' ' 
      if( wsrc.eq.'INP' ) then
        write(iprnt,'(1x,a,1x,a,f9.1,a)') 
     .       'Humidity  ','sestbl input',wetvar0*100.,' %'
      elseif( wsrc.eq.'RNX' ) then
        write(iprnt,'(1x,a,1x,a3,a,f9.1,a)')
     .     'Humidity  ',wsrc,'         ',wetvar0*100.,' %'
      elseif( wsrc.eq.'GPT' ) then
        metsrcmod = gptnam//'    '
        write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .     'WV press  ',wsrc,metsrcmod,e0,' hPa'
       elseif( wsrc.eq.'STP' ) then
         write(iprnt,'(1x,a,1x,a3,a,f9.1,a)')
     .     'Humidity  ',wsrc,'         ',wetvar0*100.,' %'
      endif          
c     Lapse rate           
      metsrcmod = ' '     
      if( lsrc.eq.'GPT' ) metsrcmod = gptnam//'    '
      write(iprnt,'(1x,a,1x,a3,1x,a8,f9.1,a)')
     .    'Lapse rate',lsrc,metsrcmod,lapse,' C/km'      
           
c rwk 160426: Remove the code that sets default PTH values if they are missing
c             for an epoch, which can occur only for a RINEX met file input.  
c             Rather, we now just issue a warning and extrapolate the last values.
                   
c     Set the 4-character zenith delay and mapping function values for the 
c     and computation c-file from the chosen 3-character sources and b-file 
c     input
      if( psrc.eq.'INP' ) then
        dryzen = 'INP '
      elseif( psrc.eq.'RNX' ) then
        dryzen = 'RNX '
      elseif( zsrc.eq.'UFL' ) then
        if( metmod(1:2).ne.'  ' ) then
          dryzen = metmod(1:4)
        elseif( mapmod(1:2).ne.'  ' ) then
          dryzen = mapmod(1:4) 
        endif
      elseif( psrc.eq.'GPT'.or.zsrc.eq.'GPT' ) then
        dryzen = gptnam
      endif
      if( wsrc.eq.'INP' ) then
        wetzen = 'INP '
      elseif( wsrc.eq.'RNX' ) then
        wetzen = 'RNX '
      elseif( wsrc.eq.'UFL' ) then
        if( metmod(1:2).ne.'  ' ) then
          wetzen = metmod(1:4)
        elseif( mapmod(1:2).ne.'  ' ) then
          wetzen = mapmod(1:4) 
        endif
      elseif( wsrc.eq.'GPT' ) then
        wetzen = gptnam
      endif                                        
c     drymap and wetmap will be ok as read from the b-file for VMF1 but
c     should be modified slightly for 'GPT' to reflect the particular model
      if( drymap(1:2).eq.'GP  ') then
        drymap = gptnam 
      endif
      if ( wetmap(1:2) .eq. 'GP') then
        wetmap = gptnam
      endif           
      write(iprnt,'(1x,a,8x,a4)') 
     .    'Dry mapping function',drymap
      write(iprnt,'(1x,a,8x,a4)')
     .    'Wet mapping function',wetmap       

c     Compute the delay at zenith for storing in PREVAL. 
c     Note: this will be the value at the first epoch; actual a priori values
c     at each epoch are now on Record 4 of the C-file (10.60)
      elev = twopi/4.d0                                   
      call get_atmdel( jd0,t0,0,0.d0,elev,zendel0,atpart )  
c     set flag for met-file availability (recode this soon)
      avlmet = 7

c     The elevation cutoff is now used only for simulation mode; otherwise model all obs
      if( .not.simulation ) elvcut = 0.d0
                          
c--------------------------------------------------------------------------------

c     Report the ionospheric model
    
c      Source read from the b-file in open.f
c      Currently the only options are blank (no model) and an f- (IONEX) file

      write(iprnt,'(/,a)') 'Ionospheric models'   
c        'I' means old-style batch file (Inertial frame)   
      if( ionsrc.eq.'NONE') then 
        write(iprnt,'(a)') '---no model applied'
      else
        if( iuf.eq.0 ) then
          call report_stat('FATAL','MODEL','setup',' '
     .     ,'2rd & 3rd order ion requested but no f-file',0)
        endif
        write(iprnt,'(2a)') 
     .     '  2d & 3rd order ionospheric terms computed from'
     .    ,' IONEX map (f-file)'   
        call read_ionex
        ionsrc = 'GMAP'                        
        if( magfield(1:6).eq.'DIPOLE' .or. magfield(1:6).ne.'IGRF10'
     .       .or. magfield(1:6).eq.'IGRF11'
     .       .or. magfield(1:6).eq.'IGRF12' 
     .       .or. magfield(1:6).eq.'IGRF13' ) then
           write(iprnt,'(2a)') 'Magnetic field model: ',magfield
        else
         call report_stat('FATAL','MODEL','setup',magfield
     .      , 'Unrecognized magnetic field model in batch file:',ioerr)
        endif
      endif

c---------------------------------------------------------------------------------

c     Determine the earth rotation values at the span midpoint

      jdmidpt = jd0  
      tmidpt  = t0
      call timinc (jdmidpt,tmidpt,span/2.d0)
      fract= tmidpt/86400.d0
      call ut1red ( iut1,jdmidpt,fract,ut1,ut1dot,iuttyp )
c     Convert TAI-UT1 to UT1-TAI
      ut1 = -ut1
      ut1dot = -ut1dot
C     Add tidal correction if UT1-UTC is regularized (UT1R)
      if( iuttyp.eq.2 ) then
        fjd = jdmidpt + fract - 0.5d0
        call funarg ( fjd,xl,f,d,ascm )
        call ut1tid( ut1,xl,f,d,ascm )
      else if (iuttyp.eq.0 ) then   
        call report_stat('FATAL','MODEL','setup',' '
     .     ,'UT1 type = 0, set = 2 or 4 in UT1. table',0)
      endif
      call polred ( ipole,jdmidpt,fract,xp,yp,xpdot,ypdot )
      tdtgpst = 32.184d0 + 19.d0
      fjd= dble(jdmidpt) + fract + tdtgpst/86400.d0
      call nuttab ( inut,fjd,oblq,psi,eps,rnut )
      psi=  psi / radian * 3600.d0
      eps=  eps / radian * 3600.d0

c     Nutation rates not yet computed

      psidot=.0d0
      epsdot=.0d0       
                                       
c--------------------------------------------------------------------------

c       Write Headers on Output C-File
                                         
      call wrthed     

      write(iprnt,'(a,/,a,/)') 
     . '---------------------------------------------------'
     .,' ** END HEADER INFO -- BEGIN PROCESSING MESSAGES'

cd      if( debug ) stop          
      RETURN
      END


