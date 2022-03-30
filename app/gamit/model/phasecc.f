Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1994.   All rights reserved.

      Subroutine phasecc( ichan,yatt,delaycor)

c Purpose    :  To calculate delay corection to L1 and L2 carriers
c               for elevation and azmith dependant vertical antenna phase
c               centre errors.   Coded originally S. McClusky for receiver
c               antennas; modified by R. King to include satellite 
c               antennas.    Phase center patterns are read from file antmod.dat
c               which can be either GAMIT format (rcvr antennas only) or 
c               ANTEX format.  
c
c Parameters :            
c         In :  ichan     : channel index (satellite)                           i*4  
c               jdobs     : Julian day of observation  (model.h)                i*4  
c               elev      : elevation angle of sat (radian) (model.h            r*8
c               azim      : azmith of satellite (radian) (model.h               r*8 
c               nadang    : nadir angle of satellite antenna (model.h) (rad)    r*8
c               iprn      : PRN # of SV  (ischan(ichan) model.h)                i*4
c               svantbody : SV ant/body type (BLOCK II-R, BEIDOU etc            c*20
c               iepoch    : session epoch number (1=first call)                 i*4

c                   in common /lunits/ of model.h
c                ipcv     : unit number for antmod.dat file          
c                   in common /lfcom/    of model.h:
c               sitecd    : 4-character site code                              
c                   in common /rcvantcom/ of model.h:     
c               antcod               : 6-character GAMIT antenna code 
c               radome               : 5-character radome code
c               antmod1              : requested antenna model (AZEL ELEV or NONE)
c               antmod               : antenna model available and used (AZEL ELEV or NONE)
c               offarp(3)            : offset of antenna ARP from monument (U N E) (m)
c               offl1(3) offl2(3)    : mean offset of antenna phase center from monument (U N E) (m)   
c               pcvminelev           : minimum elevation in PCV model from ANTEX file (deg)

c                   in  common /svantcom/ of model.h: 
c               svantmod_in          : requested antenna model (AZEL ELEV or NONE)
c               svantmod(maxsat)     : antenna model available (AZEL ELEV or NONE)
c               svantmod_snx(maxsat) : 10-character SINEX code for antenna model   

c                    In commmon /antcom/  of model.h
c               newant    : true if new antenna, read antmod.dat                logical
                     
c                    In common /obscon/ of model.h
c               atxfrq(2)  : 3-character code for ANTEX frequency (e.g. E01)    c*3 

c        Out :  delaycor : correction to l2 and l2 time delays caused by
c               azmith and elevation dependant ground and SV antenna
c               phase centre errors (sec)                                 r*8(2,ichan)
c
c MOD TAH 200205: Changed pcoffl1 and pcoffl2 from get_antpcv to reflect the 
c       effects of rotated antenna.  No change to this routine needed. antdaz
c       in model.h and read from station.info (0 if no AntDAZ column in station.info
c  
c subroutines called : report_stat, get_antpcv, get_svantpcv
c
c Created    :  93/09/17              Last modified :  rwk 2003/5/2
c
      Implicit none

      include '../includes/dimpar.h' 
      include '../../libraries/includes/const_param.h'
      include '../includes/model.h'

      logical      first,found,warnings,debug
      integer*4    iprn,ichan
      real*8       elevd,azimd,nadangd
      real*8       pcoffl1(3),pcoffl2(3),delaycor(2,maxsat)
     .           , svoffl1(3),svoffl2(3),corrl1,corrl2
      real*8       yatt(maxsat)

      data debug/.false./
        
                   
      if( ipcv.eq.0 ) then 
         call report_stat('FATAL','MODEL','phasecc','antmod.dat',
     .   'no antenna table antmod.dat',0) 
      endif

      iprn = ischan(ichan)     
      delaycor(1,ichan) = 0.d0
      delaycor(2,ichan) = 0.d0    
      if ( antmod(1:2).ne."NO".and.elev(ichan).gt.0.d0 ) then
c       convert azimuth, elevation, and nadir angle to degrees.
        elevd = elev(ichan)*180.0d0/pi
        azimd = azim(ichan)*180.0d0/pi    
        nadangd = nadang(ichan)*180.0d0/pi                                         
        if( debug ) then  
          print *,'PHASECC calling GET_ANTPCV iepoch newant antcod '
     .             ,iepoch,newant,antcod
          print *,'PHASECC newant antcod ',newant,antcod    
          print *,'iepoch ichan elevd azimd nadangd '
     .        ,iepoch,elevd,azimd,nadangd
        endif     
cd        if( iepoch.eq.894.or.iepoch.eq.1037.or.iepoch.eq.1080 ) then     
cd         debug = .true.    
        warnings = .true.  
        call get_antpcv( jdobs,elevd,azimd,debug,pcoffl1,pcoffl2
     .                 , corrl1,corrl2 )
        if(debug) then
          print *,' antmod antmod_snx ',antmod,antmod_snx
          print *,' elev corrl1 corrl2 ',elevd,corrl1,corrl2  
          print *,' pcoffl1 pcoffl2 ',pcoffl1,pcoffl2    
        endif
c       convert antenna variation in mm to sec (vlight in km/s)
        delaycor(1,ichan) = delaycor(1,ichan) + corrl1/vel_light_km/1.d6
        delaycor(2,ichan) = delaycor(2,ichan) + corrl2/vel_light_km/1.d6
        if(debug)  print *,'Site PCV sat elev corrl1 corrl2 (mm) '
     .     ,ichan,elevd,corrl1,corrl2 
      endif

      if ( svantmod(ichan)(1:2).ne."NO" ) then
                    
        first= .false.   
        if( debug ) then  
          print *,
     .     'PHASECC calling GET_SVANTPCV svantbody prn jdobs nadangd'
     .             ,svantbody(ichan),iprn,jdobs,nadangd
        endif     
        call get_svantpcv( jdobs,ichan,nadangd,yatt(ichan),first
     .                    ,found,svoffl1,svoffl2,corrl1,corrl2 )
c       convert antenna variation in mm to sec (vlight in km/s)  
        delaycor(1,ichan) = delaycor(1,ichan) + corrl1/vel_light_km/1.d6
        delaycor(2,ichan) = delaycor(2,ichan) + corrl2/vel_light_km/1.d6 
        if(debug) then
          print *,' PHASECC nadangd yatt corrl1 corrl2 '
     .      ,nadangd,yatt(ichan),corrl1,corrl2  
        endif
        if( debug ) then
          print *,'total delay correction L1/L2r '
     .     ,delaycor(1,ichan),delaycor(2,ichan)
        endif
cd      debug = .false.
      endif

      return
      end
