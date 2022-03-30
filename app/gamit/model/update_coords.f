      Subroutine update_coords

c     Update the Earth-fixed coordinates with new position (eq/rename)
c     velocity,  or antenna offsets (station.info)

c     R. King 8 May 2003 from (now-obsolete) routine update_offs (RWK 14 Feb 97)

c     Input
          
c         iepoch  I*4  : epoch of observation  (model.h)
c         jdobs   I*4  : PEP JD of observation time (model.h)
c         tobs    R*8  : GPST seconds-of-day of observation time (model.h) 
c         semi    R*8  : semi-major axis of geodetic datum (m) (in model.h)
c         finv    R*8  : inverse flattening of datum           (in model.h)
c         shft(3) R*8  : center-of-figure offset of datum (m)  (in model.h)     

c      Output 
c         newant     L*4 : True if new station.info entry found (in commmon /antcom/ in model.h
c         evec0(3,2) R*8 : Updated cartesian coordinates for L1 and L2

c     Antenna information in common /antcom/ in model.h  
c       kstarts(5) kstops(5) : yr doy hr min sec of start of current station.info entry
c       kstartr(5) kstopr(5) : yr doy hr min sec of start/stop of current eq/rename   
c       anttyp               : 20-character antenna name including radome
c       antsn                : 20-character antenna serial number
c       antmod1              : requested antenna model (AZEL ELEV or NONE)
c       antmod               : antenna model available and used (AZEL ELEV or NONE)
c       offarp(3) offapr(3)  : offset of antenna ARP from monument (U N E) (m)
c       offl1(3) offl2(3)    : mean offset of antenna phase center from monument (U N E) (m)
c       antdaz               : Alignment from True N (deg).  TAH 2020203

c     New L-file coordinates through common /kinpar1/ in model.h
c       sitecd  : 4-character site code (used by GAMIT)
c       asite   : 8-character site code (used to document lread return)
c       sname   : 12-character site name 
c       kepoch0 : apr file epoch ( decimal yrs) 
c       kpos    : Cartesian position from the L-file (m)
c       kvel    : Cartesian velocity from the (apr-style) L-file (m/yr) 
c       iul l   : unit number of the L-file     
c       iueqrn  : unit number of the eq-rename file (= 0 if not available)
c       kfflg   : flag for type of L-file (spherical = 0, Cartesian = 1)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/units.h'
      include '../includes/model.h'

      integer*4 iy,idoy,nyr,ndoy,nsod,itime(5),i,j
      integer*8 itimdif
                                                                      
      real*8 anth,antn,ante,dhpab,sec,pos(3),decyrs,tobs_yrs
     .     , pcoff_l1(3),pcoff_l2(3)
                      
      character*1 upper1
      character*5 htcod
      character*6 antcod1 
      character*16 stanam,amodel 
      character*20 rcvtyp 
      character*256 message 

* MOD TAH 200527: Saved values so that we can see if antenna really changed.
      real*8 anth_save, antn_save, ante_save  ! Saved values for testing
      character*5  htcod_save
      character*20 anttyp_save    ! Check size in model.h (also declared in 
                                  ! utils/rinex.h
      
      logical newcrd,warnings,radome_sub,fcheck,debug  
                           
      data debug/.false./

      save anth_save, antn_save, ante_save
      save anttyp_save, htcod_save
             
c  Variables and logic for stations coordinates:
c
c     L-file values read as spherical (old-style) or Cartesian (GLOBK-style) 
c         but now stored only as Cartesian in the commons of modkin.h.
c
c     Earth-fixed values corrected for antenna offsets and velocity (kvel)
c         are evec0(3,2) (km)
c
c     Inertial values are evec(6,2) (km)
c
c     The L-file values are read initially in SETUP and written into the
c     C-file header.  They can be updated during a session (UPDATE_COORDS
c     called by MODEL), however, if there is an earthquake or other rename 
c     (eq_rename entry) or change in antenna type or offset (new station.info 
c     entry), or if velocities are included (GLOBK apr-file) and exceed 1 mm/day.  
c     For kinematic or dynamic measurements, they are updated at each epoch.
c     Whenever the L-file values are updated, the Earth-fixed and inertial
c     values are recomputed.  
c      
c  Keys for updates in modkin.h commons:
c      kstarts, kstops : current stnfo entry     defaults 
c      kstartr, kstopr:  current eq/rename entry   default 1900, 2100
c      kvflg : true if velocities to be used 

  
c  Initialize flags for new initial coordinates, offsets, or position 

      newcrd = .false.

C  Convert JD, sod to yr,doy,hr,min,sec [itime(1-5)] for station.info/eq-rename checks
c  Get time in years for coordinate updates
                              
        call dayjul(jdobs,itime(1),itime(2))   
        if(debug)  print *,'UPDATE_COORDS jdobs yr day '
     .            , jdobs,itime(1),itime(2)
        call ds2hms( itime(1),itime(2),tobs,itime(3),itime(4),sec) 
        itime(5) = int(sec)   
        if(debug) print *,  'itime kstart kstopr ',itime,kstartr,kstopr
        tobs_yrs = decyrs( itime(1),itime(2),tobs )

              
c  If outside the bounds of the current station.info entry, read another one

      if(debug)  print *,'itime kstarts kstops ',itime,kstarts,kstops
      if( itimdif(itime,kstarts).lt.-60 ) then
c       obs time more than 1 min earlier than station.info time, entry is missing
        write(message,'(a,a4,a,i5,i4,3i3,a,i5,i4,3i3,a)') 
     .   'Obs time for ',sitecd,' (',itime
     .   ,') too early for station.info start (',kstarts,')'
        call report_stat('FATAL','MODEL','update_coords',' ',message,0) 
      elseif( itimdif(itime,kstops).gt.5 ) then
        newant = .true.
c       obs time more than 5s later than station.info stop time, 
c       read another entry and update coordinates  
        if(debug) print *
     .    ,'updating station.info entry sitecd itime,kstops '
     .         ,sitecd,itime,kstops
        call dayjul(jdobs,iy,idoy)
        if(debug) print *,'calling even_minute iy idoy tobs '
     .       , iy,idoy,tobs        
c *** rwk 100331  replace this by a rounding to the nearest second
c        call even_minute(iy,idoy,tobs,nyr,ndoy,nsod)  
        nyr = iy
        ndoy = idoy
        nsod = idint(tobs) 
        if(debug)  print *,'called even_minute ',nyr,ndoy,nsod
* MOD TAH 200203: Added AntDAZ to list of values from station.info
        call rstnfo( istnfo,  sitecd, nyr, ndoy, nsod, 0
     .              , stanam, anth, antn, ante, antdaz, rcvcod, antcod1
     .              , htcod, radome_in, swver, rcvers, rcvrsn, antsn
     .              , kstarts, kstops )
        if( debug )   print *,'new anth ',anth  

* MOD TAH 200527: Only try to update antenna information if it is 
*       different from previous values
        if( anttyp.ne.anttyp_save .or. htcod.ne.htcod_save .or.
     .      anth  .ne.anth_save  .or. antn .ne.antn_save  .or.
     .      ante  .ne.ante_save ) then 
           call ant_alias(antcod1,antcod)
c          Convert the raw antenna offset position to phase-center-above-mark
c          The quantities offl1 and offl2 are offsets of L1 and L2 phase centres
c          from the monument, computed by hisub from the offsets measured in the
c          field (anth,antn,ante) and the L1/L2 APR to phase centre offsets
c          read from antmod.dat.  
           warnings = .true.    
           call hisub( ihi,anth,antn,ante,antcod
     .               , htcod,sitecd,iy,idoy,isessn,offarp,warnings )
c          Get the official receiver and antenna names from rcvant.dat translation table.... 
           call read_rcvant(1,1,antcod,anttyp,radome_in,rcvcod,rcvtyp
     .                     ,pcncod)
           call read_rcvant(1,2,antcod,anttyp,radome_in,rcvcod,rcvtyp
     .                     ,pcncod)
           write(iprnt,'(a,i5,2x,i4,1x,i3,1x,3i3)')
     .          'Station.info update at epoch '
     .          ,iepoch,(itime(i),i=1,5)
           write(iprnt,'(a,1x,a6,1x,a20,1x,a5,3f8.4)') 
     .       'New values (rcvcod anttyp htcod U N E): '
     .           ,rcvcod,anttyp,htcod,anth,antn,ante
           write(iprnt,'(a,6f8.4)') 
     .'              (L1/L2 U N E): ',offl1,offl2 
           call get_antinfo(debug)
c          Note in the P-file that we've done an update
           write(message,'(a,i5,i4,2i3,a )') 'Station.info update at '
     .                               ,(itime(i),i=1,4),' --see P-file'
           call report_stat('STATUS','MODEL','update_coords',' ',
     .                       message,0)
           call report_stat('WARNING','MODEL','update_coords',' '
     .                     , message,0)

* MOD TAH 200527: Save values so we can compare later
           anttyp_save = anttyp
           htcod_save = htcod
           anth_save  = anth 
           antn_save  = antn 
           ante_save  = ante 
        else
* MOD TAH 200527: Record in p-file that we are checking
           newant = .false.
           write(iprnt,'(a,i5,2x,i4,1x,i3,1x,3i3)')
     .          'Check station.info at epoch '
     .          ,iepoch,(itime(i),i=1,5)
        endif
      endif
                     

c     Apply antenna offsets and compute the geodetic and cartesian 
c     coordinates and partials 
             
      if( newcrd .or. newant .or.kvflg.gt.0 ) then 
        do i=1,3
          pos(i) = kpos(i) + kvel(i)*(tobs_yrs-kepoch0)
        enddo    
      
       if( debug)  print *,'UPDATE_COORDS bef CORTYP '
     .     , ' pos l1l2 semi shift evec0 latr height sitepart0 ',
     . pos,offl1,offl2,semi,shft,evec0,latr,height,sitepart0
        call cortyp( pos,offl1,offl2,semi,finv,shft
     .             , evec0,latr,height,sitepart0 ) 
       if(debug)    print *,'UPDATE_COORDS aft CORTYP '
     .     , ' pos l1l2 semi shift evec0 latr height sitepart0 ',
     .       pos,offl1,offl2,semi,shft,evec0,latr,height,sitepart0
cd        write(*,'(a,3f20.6)') 'DEBUG L1 ',(evec0(i,1)*1.d3,i=1,3)
cd        write(*,'(a,3f20.6)') 'DEBUG L2 ',(evec0(i,2)*1.d3,i=1,3)    
        write(iprnt,'(a,3f14.4)')  '  New coordinates ',(pos(i),i=1,3)  
        if( newant.or.newcrd ) then
          write(iprnt,'(a,/,a,3f14.6,/,a,3f14.6)')         
     .      '  New Earth-fixed antenna phase-center coordinates (km): '
     .     ,'   L1 ',(evec0(i,1),i=1,3)
     .     ,'   L2 ',(evec0(i,2),i=1,3)  
         endif
      endif 

      return
      end

