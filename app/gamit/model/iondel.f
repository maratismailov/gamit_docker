      Subroutine IONDEL( ichan,prn,satcrd, pnsmat, rhat
     .                 , grpdel, phsdel )

                            
c       Calculate the higher order ionospheric delay using one of several models.
c         R. King Apr 2007  
C
C        EJP 15 Jan 2008
C        EJP 20 Feb 2008 - rotated maps for earlier periods, different mapping functions for slant_tec.
c        EJP 11 Dec 2008
c        EJP 29 Apr 2009
c        RWK  4 May 2009 - put EJP code from parallel version into a subroutine 'mag_dipole'
C	 EJP 18 Jul 2009 edited mapping function (see also slant_tec) and pass tec as well as STEC to iono_corr
c*** Temporary: There are redundant inputs here to allow some checking of
c               calculations; later use only cartesian (Fritsche et al ) or
c               trig (Kedar et al)  rwk 070515
c
c      Input
c
c        ionsrc :  Source of TEC values (must be IONEX file in current code) (in common in model.h)
c        jdobs  :  PEP Julian Day of observation    
c        tobs   :  Seconds of day observation  
c        lonr   :  Longitude in radians (ejp - of site) (model.h)
c        latr   :  Geodetic latitude in radians (ejp -  of site) (model.h)
c        height :  Site geodetic height (km) (in model.h)
c        prn    :  Satellite PRN number (for warning) 
c        sitpos :  Station coordinates (xyz, km)
c        evec0  :  sitecoords (xyz, km) to check (model.h)
c        rhat   :  Unit vector in direction of satellite    
c        pnsmat :  Rotation matrix from E-fixed to inertial
c        ---input variables available from common/obscon? in ../includes/model.h
c        fl1(ichan),fl2(ichan): Frequencies in Hz of the L1 and L2 signals    
c        ---input variable available from common/obsrecords/ in ./includes/model.h
c        elev(ichan)     :  Elevation angle of satellite
c        azim(ichan)     :  Azimuth of satellite               
c        ---input variables avaiable from  common/ioncom/ in ../includes/model.h    
c        ilon1  :  First longitude value of ion maps (deg)
c        ilat1  :  First latitude value of ion maps (deg)  
c        dilat  :  Latitude interval (deg)
c        dilon  :  Longitude interval (deg)
c        nilon  :  Number of longitude values per map
c        nilat  :  Number of latitude values per map 
c        ntion  :  Number of ion maps (times)
c        ion_time(ntion)  : Times of the maps (day-of-year)
c        ion_val(nilon,nilat,ntion) :  Gridded maps   
c         --input variables available from commons in ../includes/model.h
c        cfiln       : MODEL output filename - avoid doubled output    
c        sitecd      : 4-character code of site (for debugging)
c      Output

c        grpdel(2) : Group delay in meters at L1 and L2
c        phsdel(2) : Phase delay in meters at L1 and L2

c     References;
c       Bassiri, S., and G. A. Hajj, Higher-order ionospheric effects on the globak positioning
c        system observables and means of modeling them, Man. Geod., 18, 28--289, 1993.   
c       Kedar, S., G.A. Hajj, B. D. Wilson, and M. B. Heflin, The effect of the second order
c         GPS ionospheric correction on receiver positions, Geophy. Res. Lett., 30, 1829, doi:
c         10.1029/2003GL017639, 2003. 
c       Fritsche, M, R. Dietrich, C. Knofel, A. Rulke, and S. Vey, Geophys. Res., Lett., 32, 
c        L23322, doi:10.1029/2005/GL024342, 2005.


      implicit none         

      include '../includes/dimpar.h' 
      include '../../libraries/includes/const_param.h'  
      include '../includes/model.h'   

      integer*4 doy,iyr,prn,latindx,lonindx,it,i,j,IRC2,ichan
    
      real*4 lond,latd

      real*8 londp
          
      real*8 sitpos(3),rhat(3) 
     .     , pnsmat(3,3),epstol,htion,re
     .     , nmax,radian,obstime,p,q,dt,tec1,tec2,tec
     .     , theta,thetam,thetamp,phim,sbeta,cbeta
     .     , cdelta,sdelta,stheta,ctheta,sphi,cphi,sphim,cphim
     .     , dcap,freq(2)
     .     , grpdel(2),phsdel(2)
     .     , bcoff 
C     .     ,grpdel3(2),phsdel3(2)

c     for debug only 
      integer*4 debug
  
      character*256 message  

      logical first,found,start_warning,end_warning,tec1null,tec2null

C---------------------------------------------------------------
C     2nd order group delay/m, 3rdorder group delay/m, 2nd order phase delay/m 3rd order phase delay/m
      Real*8 grpdel2(2),grpdel3(2),phsdel2(2),phsdel3(2)

C     Inertial satellite coords + vel?
      Real*8 satcrd(6)
C     Inertial cartesian satellite coords/km
      Real*8 satpos(3)
C     Earthfixed cartesian satellite coords/km
      Real*8 satfix(3)
C     Unit vector satellite to site (cartesian geocentric)
      Real*8 unvec(3)           
C     geocentric lat, lon/degrees of ionospheric pierce point
      Real*8 slat, slon
C     site position efixed Cartesian/km (evec0(i,1), ie L1 phase centre)
c     Real* 8 evecpos(3)
C     Rotation matrix, inertial to earthfixed coordinates
      Real*8 xrotinv(3,3)
C     Epoch in MJD time
      Real*8 MJDt           

C Variable for magnetic field calculation (2nd order ionospheric correction)
C Date of epoch in decimal years AD,leap year flag, year of epoch AD, day of year of epoch
      Real*8 datedec, leap
C Components(along XYZ) of magnetic vector, magnitude of magnetic vector
      Real*8 mag_vec(3), magsize

C Mapping function F(z) = 1/cos(z'), where sin(z') =R(sin(alpha*z)/(R+MLSM_H)
C z is zenith angle, R is earth radius
C Variables for slant TEC calculation
        real*8 zenith, MSLMhght, STEC,alpha,shellht1,shellht2
C Earth radius for MLSM mapping function/km
      Real*8 radiM
C Height of thin layer ionosphere
C 	Height pre 1998 doy 87
      	  Parameter (shellht1 = 400.0d0)
C  	Height 1998 doy 87 and after
      	  Parameter (shellht2 = 450.0d0)
C       MSLM model best approx height (Schaer 1999)
      Parameter (MSLMhght = 506.7d0)
C Values as best fit of Modified Single Layer Model to JPL ESM model. 
C This is what was used to make the CODE TEC maps after doy 251 2001.
C Correction factor
      Parameter (alpha = 0.9782d0)
      Parameter (radiM = 6371.d0)
C---------------------------------------------------------------

                          
c     epstol is a small number tolerance to allow use of the end values
c     currently set to be 90s if times are in days
      data epstol/.001/


c     Initialize 'first' for setup call
      data first /.true./,start_warning/.false./,end_warning/.false./      

      save first, radian, thetam, start_warning, end_warning

C     EJP Debug
C      Print*, 'Beginning of model/iondel'
C      Print*, 'sitpos',sitpos
C      Print*, 'satcrd'
C      Print*, satcrd
C      Print*, 'rhat',rhat
C      Print*, 'evec0'
C      Print*, evec0
C      Print*, 'elev',elev(ichan),'azim',azim(ichan)
C      Print*, 'prn',prn
c      Print*, 'Beginning of model/iondel pnsmat'
c      Print*, pnsmat
c..   If no correction available, return zero

      if( ionsrc.eq.'NONE') then
        do i=1,2
          grpdel(i) = 0.0d0
          phsdel(i) = 0.0d0    
        enddo
        return  
      elseif (ionsrc.eq.'GMAP' ) then
        continue
      else
        call report_stat('FATAL','MODEL','iondel',ionsrc
     .      , 'Unrecognized source for ionospheric corrections:',0)
      endif

C Set up debugging value  - currently so some sites will output some data
C this will print no values if 1, basic values if 2 and nearly all debugging values if 3
C Have checked and iondel is only called when ipass is 2. EJP Dec 2007

C Debug - output data sample eg at GOLD and DARW when year ends in 3.
C      If (cfiln(6:6).eq.'3') then
C       If (sitecd.eq.'GOLD'.or.sitecd.eq.'DARW') then
C              debug = 2
C       Else
C              debug = 1
C       Endif
C      Else
C              debug = 1
C      Endif     
   
      If (cfiln(6:6).eq.'3') then
       If (sitecd.eq.'DRAO'.or.sitecd.eq.'GOLD'.or.sitecd.eq.'MKEA'.or.
     .  sitecd.eq.'THU3'.or.sitecd.eq.'ASC1'.or.sitecd.eq.'TOW2'.or.
     .  sitecd.eq.'UNSA'.or.sitecd.eq.'SANT'.or.sitecd.eq.'OHI2'.or.
     .  sitecd.eq.'AREQ'.or.sitecd.eq.'TIXI'.or.sitecd.eq.'BILI'.or.
     .  sitecd.eq.'POL2'.or.sitecd.eq.'MCM4'.or.sitecd.eq.'PERT'.or.
     .  sitecd.eq.'MALI'.or.
     .  sitecd.eq.'SYOG'.or.sitecd.eq.'CRO1'.or.sitecd.eq.'DARW') then
              debug = 2
       Else
              debug = 1
       Endif
      Else
              debug = 1
      Endif           
C Debug, test output
       debug = 1
C      debug = 2
c      debug = 3
C      debug = 4
C   -------------------------------------------------------------    

      if( first ) then  
* MOD TAH 200125: Commented out print statements. 
        !  Print*, 'MODEL\iondel first call'
        !  Print*,'MODEL\iondel - ',magfield,' magnetic field'
        if(debug.ge.4) then
          print *,'IONDEL DEBUG ionsrc first ',ionsrc,first
          print *,'nilon nilat ntion ',nilon,nilat,ntion
          write(*,'(a,f7.3,5(/,20f5.1))') 
     .       't1 lat1 lon ',ion_time(1),(ion_val(i,1,1),i=1,nilon) 
          write(*,'(/,a,f7.3,5(/,20f5.1))') 
     .      't5 lat5 lon ',ion_time(5),(ion_val(i,5,5),i=1,nilon)
          write(*,'(/,a,f7.3,5(/,20f5.1))') 
     .      'tlast latlast lon ',ion_time(ntion)
     .       ,(ion_val(i,nilat,ntion),i=1,nilon)
         endif       
c       degrees to radians
        radian = datan(1.d0)/45.d0
c Make sure height of ionospheric thin shell is as wished. 
C htion 450km is as for Fritsche et al.(2005) and Hernandes-Pajares et al. (2007).
C Value is passed to call_mag.f
        re = 6378.13655d0      
C        htion = 400.d0
        htion = 450.d0
C        htion = 500.d0
C        htion = 300.d0
C        htion = 350.d0


        bcoff = 2.61d-18*(re/(re+htion))**3              
        first = .false.
      if(debug.ge.2) then
        Print*,'Model/iondel debug pre return'
      endif
C        return
      endif

C----------------------------------------------------------------------
C Debug - makes it easier to read output 
      if(debug.ge.2) then
C      Print*, 'pi', pi, 'radian', radian
      Print*,'------------------------------------------------------'
      Print*, 'MODEL/iondel SITE: ',sitecd, '  T:',tobs,' prn',prn
      Print*,'------------------------------------------------------'
      End if
c.. Initialize return variables to zero
      phsdel(1) = 0.d0
      phsdel(2) = 0.d0
      grpdel(1) = 0.d0
      grpdel(2) = 0.d0
    
c.. Get the observation time

      call dayjul( jdobs,iyr,doy)         
      obstime = dfloat(doy) + tobs/86400.d0 
    
C Deal with end of year problem when have 13 maps in ionex file 
      if (ntion.eq.13 .AND.doy.ge.365 
     .  .AND.ion_time(ntion).eq.1.0000)   then
           ion_time(ntion)=1.0000+ion_time(1)
      endif
      if( debug.ge.4 ) then
        print *,'MODEL/iondel:obstime ion_time1 ion_timelast '
     .     ,obstime,ion_time(1),ion_time(ntion)
      endif    
                             
C debug
C      Print*, 'Model/iondel: doy',doy,'iyr',iyr
Calculate date in decimal years

      if (mod(iyr,400).EQ. 0) then 
          leap=1.d0
        elseif (mod(iyr,100).eq. 0) then 
          leap=0.d0
        elseif (mod(iyr,4).eq. 0) then 
          leap=1.d0
        else 
          leap=0.d0
      endif
      datedec=real(iyr)+((real(doy))/(365.d0+leap))
C Debug
        if(debug.ge.3) then
      Print*,'MODEL\iondel: jdobs',jdobs,'yearad',iyr,'doy',doy
      Print*,'MODEL\iondel: datedec', datedec
        end if
c.. See if the requested time is outside the available array of values
c However, as 'rotating' the maps to follow the sun, can still get a value out.

      if( obstime.lt.(ion_time(1)-epstol) ) then
        if( .not.start_warning ) then
          write(message,'(a,f8.4,a,f8.4)') 'Obs epoch '
     .        ,obstime,' before start of IONEX file ',ion_time(1)  
          call report_stat('WARNING','MODEL','iondel',' ',message,0)
          start_warning = .true.
        endif
C        return
      elseif( obstime.gt.(ion_time(ntion)+epstol) ) then
        if( .not.end_warning ) then 
          write(message,'(a,f8.4,a,f8.4)') 'Obs epoch '
     .     ,obstime,' after end of IONEX file ',ion_time(ntion)  
          call report_stat('WARNING','MODEL','iondel',' ',message,0) 
          end_warning = .true.
        endif
C        return
      endif
           
      if (debug.ge.4 ) then
        print *,'jdobs tobs ',jdobs,tobs   
        print *,'latr lonr height sitpos ',latr,lonr,height,sitpos
        print *,'  rhat elev/deg azim/deg'
     .            ,rhat,elev(ichan)/radian,azim(ichan)/radian
        print *,'  az radians ',azim(ichan)
        print*, 'pnsmat'
        Print*,pnsmat
C        print*, 'prn',prn
      endif

c.. Find the times of IONEX maps bracketing the observation

c     see if the requested time is outside the range of the array   
      it = 0     
      found = .false.
      do while (.not.found.and.it.lt.ntion)       
        it = it + 1
        if( dble(ion_time(it+1)).ge.obstime ) found = .true.     
c        print *,'obstime it ion_time found '
c     .    ,obstime,it,ion_time(it),found
      enddo  
C      if( .not.found ) then
C         write(message,'(a,f8.4,a,2f8.4,a)') 'Obs epoch '
C     .     ,obstime,' not found in IONEX file ('
C     .     ,ion_time(1),ion_time(ntion),')'
C          call report_stat('WARNING','MODEL','iondel',' ',message,0) 
C          return
C      endif
C------------------------------------------------------------        
c.. Get lat and lon in degrees to interpolate from the ion maps: 
C EJP - Calculate rotation matrix for inertial to earthfixed coordinates      
      call transp(pnsmat,xrotinv,3,3)
C      Print*,'MODEL\iondel  xrotinv'
C      Print*, xrotinv
C      Print*,'MODEL\iondel  pnsmat'
C      Print*, pnsmat
C EJP - convert satellite inertial coordinates to earthfixed  
c RWK 141028: remove checks and make assignment such that evec0 can be in common
      Do i=1,3
c RWK    satpos(i)=satcrd(i)
C        evecpos(i)=evec0(i,1)
        satpos(i) = evec0(i,1) 
c end RWK
      Enddo
C      Print*,'MODEL\iondel  satcrd'
C      Print*, satcrd 
C      Print*,'MODEL\iondel  satpos'
C      Print*, satpos
      Call matmpy (xrotinv,satpos,satfix,3,3,1)
C Debug
        if(debug.ge.4) then
      Print*,'MODEL\iondel satfix: earthfixed in km'
      Print*, satfix
C      Print*,'MODEL\iondel  evec0'
C      Print*,evec0 
C      Print*,'MODEL\iondel evecpos (receiver coords in km)'
C      Print*,evecpos
        endif
C EJP -  calculate the piercepoint of GPS signal and layer (i.e. latitude and longitude)
C Satfix is satellite coords, sitpos is receiver coords, latlon is pierce point coords 

        call callppt(satfix,sitpos,unvec,slat,slon,debug,htion)
        if(debug.ge.2) then
      Print*,'MODEL\iondel ppt slat,slon/deg',slat, slon
       Print*, 'MODEL\iondel unit vector : ', unvec
         endif
      
c------------------------------------------------------------     
c     Shift the longitude by the time since the first of the two bracketing
c     maps to compensate for the correlation between the ionsophere and the
c     Sun's position [Schaer, Gurtner, and Feltens, IONEX Version 1, Proceedings
c     of the IGS AC Workshop, Darmstadt, 1998]
         
C debug - test values
C      slat=-20.5
C      slon=-30.5
c      print *,'latr latd lonr lond ',latr,latd,lonr,lond

      londp=0.d0 
      latindx=0
      lonindx=0
       londp = slon + (obstime -ion_time(it))*360.d0
      if( londp.gt.360.d0 ) londp = londp - 360.d0
      if( londp.gt.180.d0 ) londp = londp - 360.d0
      if( londp.lt.-180.d0) londp = londp + 360.d0

c slat is pierce point latitude in -90 to 90 deg. The CODE IONEX files only go from -87.5 to 87.5 deg
C 8 Dec 2008  This is not a proper fix, but only happens at very high latitude sites e.g. MCM4 and THU3, for a few points.
      if (slat.gt.87.5d0 ) slat = 87.499999d0
      if (slat.lt.-87.5d0 ) slat = -87.499999d0


      if( debug.ge.4 ) then       
        print *,' obstime', obstime
        print*, ' it ion_time(it)', ion_time(it)
        print*, ' slat slon /deg',slat,slon
        print*, 'londp ',londp
        print*, 'dilat', dilat, 'dilon',dilon
      endif
                                      
c.. Get the indices for the "lower-left" corners of the spatial box
c   (where "lower" is first in the array but usually higher in latitude
c   if the entries are north to south; longitude always west to east)
                    
C      latindx = idint(dint((slat-ilat1)/dilat)) + 1
C      lonindx = idint(dint((slon-ilon1)/dilon)) + 1   
      latindx = idint(dint((slat-ilat1)/dilat)) + 1
      lonindx = idint(dint((londp-ilon1)/dilon)) + 1   
c      print*, ' tec1 latindx',latindx,'  lonindx',lonindx

c.. Interpolate to the requested lon/lat at times before and after the observation
      p = dmod(dabs(londp),dabs(dilon))/dabs(dilon)
       if( londp.lt.0 ) p = 1.d0 - p    

      q = dmod(dabs(slat),dabs(dilat))/dabs(dilat)
      if( slat.gt.0 ) q = 1.d0 -q

      tec1 = 
     .    (1.d0-p)*(1.d0-q)*ion_val(lonindx,latindx,it) +
     .         p*(1.d0-q)*ion_val(lonindx+1,latindx,it) +
     .         q*(1.d0-p)*ion_val(lonindx,latindx+1,it) +
     .           p*q*ion_val(lonindx+1,latindx +1,it)      
C Check not influenced by null values from ionex file 
C	9999 in file, 999.9 in ion_val
      tec1null=.false.
      tec2null=.false.
      if (ion_val(lonindx,latindx,it).gt.999 .or.
     . ion_val(lonindx+1,latindx,it).gt.999 .or.
     . ion_val(lonindx,latindx+1,it).gt.999 .or.
     . ion_val(lonindx+1,latindx +1,it).gt.999) then
         tec1null = .true. 
       call report_stat('WARNING','MODEL','iondel',' '
     .      , 'Null value(s) in IONEX file, tec1',0)
      endif
C--------------------------------------------------------------------------------------------------
C Try for tec two, recalculating the lonindx as is using the next map, so rotates backwards, actually.
       londp = slon + (obstime -ion_time(it+1))*360.d0
      if( londp.gt.360.d0 ) londp = londp - 360.d0
      if( londp.gt.180.d0 ) londp = londp - 360.d0
       if( londp.lt.-180.d0) londp = londp + 360.d0

      if( debug.ge.4 ) then       
C        print *,' obstime', obstime
C        print*, ' it ion_time(it)', ion_time(it)
C        print*, ' slat slon /deg',slat,slon
        print*, 'MODEL/iondel:londp for tec2 ',londp
        print*, 'dilat', dilat, 'dilon',dilon
      endif
                                      
c.. Get the indices for the "lower-left" corners of the spatial box
c   (where "lower" is first in the array but usually higher in latitude
c   if the entries are north to south; longitude always west to east)
                    
      latindx = idint(dint((slat-ilat1)/dilat)) + 1
      lonindx = idint(dint((londp-ilon1)/dilon)) + 1   


c.. Interpolate to the requested lon/lat at map times before and after the observation

      p = dmod(dabs(londp),dabs(dilon))/dabs(dilon)
      if( londp.lt.0 ) p = 1.d0 - p

      q = dmod(dabs(slat),dabs(dilat))/dabs(dilat)
      if( slat.gt.0 ) q = 1.d0 -q


c.. Interpolate to the right time
                          
c   If outside the map boundaries, then just use the rotated tec1:

      If (( obstime.gt.(ion_time(ntion)).or.
     .     obstime.lt.(ion_time(1)))) then
         tec = tec1   
         if(debug.ge.3) then                 
           Print*,'MODEL/iondel:tec is tec1'      
         endif  
C     Set to zero if affected by null values
         If (tec1null) then
	    tec = 0.d0
            call report_stat('WARNING','MODEL','iondel',' '
     .      , 'Null value(s) in IONEX file, tec set zero',0) 
         endif
      Else  
C     Calculate tec2
      tec2 = 
     .    (1.d0-p)*(1.d0-q)*ion_val(lonindx,latindx,it+1) +
     .         p*(1.d0-q)*ion_val(lonindx+1,latindx,it+1) +
     .         q*(1.d0-p)*ion_val(lonindx,latindx+1,it+1) +
     .          p*q*ion_val(lonindx+1,latindx+1,it+1)     

C Check not influenced by null values from ionex file
      if (ion_val(lonindx,latindx,it+1).gt.999 .or.
     . ion_val(lonindx+1,latindx,it+1).gt.999 .or.
     . ion_val(lonindx,latindx+1,it+1).gt.999 .or.
     . ion_val(lonindx+1,latindx +1,it+1).gt.999) then
         tec2null = .true. 
       call report_stat('WARNING','MODEL','iondel',' '
     .      , 'Null value(s) in IONEX file, tec2',0) 
      endif                                                                         

      if( debug.ge.3 ) then           
        print*, 'latindx',latindx,'  lonindx',lonindx
        print *,' p ' , p
        print *,' q ' , q
        print *,' tec1 ' , tec1, ' tec2 ', tec2
      endif

C Deal with any null values:
      If ((.not.(tec1null)).and.(.not.(tec2null))) then
C       Interpolate tec1 and tec2           
        dt = ion_time(it+1) - ion_time(it)                    
        tec = tec1 * (ion_time(it+1)-obstime) / dt +
     .       tec2 * (obstime-ion_time(it)) / dt

      Else if (tec1null.and.tec2null) then
C       Set tec to zero
	    tec = 0.d0
            call report_stat('WARNING','MODEL','iondel',' '
     .      , 'Null value(s) in IONEX file, tec set zero',0) 
      Else if (tec1null.and.(.not.tec2null)) then
         tec=tec2
      Else
         tec=tec1
      endif
        

       if(debug.ge.3) then
         Print*,'MODEL/iondel:tec is from tec1 and tec2' 
       endif 
      Endif
      if (debug.ge.3) then                                              
       Print*,'MODEL/iondel_pot ',prn, obstime, tec
      endif
      if( debug.ge.4 ) then
        print *,'slat ilat1 dilat latindx ',slat,ilat1,dilat,latindx
        print *,'slon ilon1 dilon lonindx ',slon,ilon1,dilon,lonindx
        print *,'obstime it ion_time1 ion_time2 dt tec '
     .        , obstime,it,ion_time(it),ion_time(it+1),dt,tec
      endif
C Having obtained interpolated VTEC, now obtain slant TEC
C      Work out zenith angle
       zenith=90.d0-(elev(ichan)*360.d0/(2.d0*pi))
C Debug
        if(debug.ge.2) then
       Print*,'MODEL\iondel zenith angle/deg', zenith 
       Print*,'MODEL\iondel az/deg', azim(ichan)*360.d0/(2*pi)
C       Print*,'MODEL\iondel tec in tecu',tec
c       az      : terrestrial azimuth from north - station to satellite    R*8(maxsat)
c       Print*, 'Model\iondel: radiM',radiM
        endif

C Subroutine to calculate slant TEC from the vertical TEC interpolated from the IONEX file
C using appropriate mapping function according to date.     
       Call slant_tec(zenith,radiM,shellht1,shellht2,MSLMhght,alpha,
     .  tec,STEC,iyr,doy,IRC2,debug)
C Debug
C       Print*, 'MODEL\iondel: Slant TEC in TECU = ',STEC

c----------------------------------------------------------------------------

C Obtain magnetic field vector
     
c     Call approprite magnetic field option        
                       
      if ( magfield .eq.'DIPOLE') then
         call mag_dipole( jdobs,latr,lonr,slat,slon,htion,mag_vec
     .                  , magsize, debug)
C         Print*, 'Dipole corrections'
      elseif ( magfield.eq.'IGRF10' .or. magfield.eq.'IGRF11'
     .    .or. magfield.eq.'IGRF12' .or. magfield.eq.'IGRF13' ) then
         call call_mag (datedec,slat,slon,mag_vec,magsize,debug,htion
     .                 ,magfield)
C         Print*,'IGRF corrections'
      else
       call report_stat('FATAL','MODEL','iondel',magfield
     .      , 'Unrecognized magnetic field model:',0)
      endif

      if(debug.ge.3) then   
        Print*, 'MODEL\iondel magvec: ',mag_vec
        Print*, 'MODEL\iondel magsize: ',magsize 
      endif                                                                     

c---------------------------------------------------------------------------------   

C Calculate the third-order (and second order) ionospheric correction
                         
      freq(1) = fl1(ichan)
      freq(2) = fl2(ichan)
      Call iono_corr( tec,STEC,vel_light,freq,mag_vec,unvec
     .              , phsdel3(1),phsdel3(2),phsdel2(1),phsdel2(2)
     .              , debug )
      if(debug.ge.2) then
        Print*,'MODEL\iondel L1 3rd order phase delay /m', phsdel3(1)
        Print*,'MODEL\iondel L2 3rd order phase delay /m', phsdel3(2)
        Print*,'MODEL\iondel L1 2nd order phase delay /m', phsdel2(1)
        Print*,'MODEL\iondel L2 2nd order phase delay /m', phsdel2(2)
      endif
C Calculate group delays
      Do j=1,2
      grpdel2(j)=-2.d0*phsdel2(j)
      grpdel3(j)=-3.d0*phsdel3(j)
      Enddo
        if(debug.ge.4) then
      Print*,'MODEL\iondel L1 3rd order group delay /m', grpdel3(1)
      Print*,'MODEL\iondel L2 3rd order group delay /m', grpdel3(2)
      Print*,'MODEL\iondel L1 2nd order group delay /m', grpdel2(1)
      Print*,'MODEL\iondel L2 2nd order group delay /m', grpdel2(2)
        endif
C Add 2nd and 3rd order corrections together
      phsdel(1)=phsdel2(1)+phsdel3(1)
      phsdel(2)=phsdel2(2)+phsdel3(2)
      grpdel(1)=grpdel2(1)+grpdel3(1)
      grpdel(2)=grpdel2(2)+grpdel3(2)
C        if(debug.ge.2) then
        if(debug.ge.3) then
      Print*,'Model\iondel phsdel:',phsdel
      Print*,'Model\iondel grpdel:',grpdel
        endif

C     convert from meters to seconds
      do i=1,2
        grpdel(i) = grpdel(i)/vel_light
        phsdel(i) = phsdel(i)/vel_light 
      enddo     
      if(debug.ge.4) then                       
         print *,'DEBUG grpdel phsdel cyc ',grpdel,phsdel
      endif
      return
      end
