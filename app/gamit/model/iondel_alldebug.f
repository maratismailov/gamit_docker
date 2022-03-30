      Subroutine IONDEL( ionsrc,jdobs, tobs, latr, lonr, height
     .                  , prn, sitpos, rhat, elev, azim, pnsmat
     .                  , freql1, freql2, grpdel, phsdel )
                            
c       Calculate the higher order ionospheric delay using one of several models.
c         R. King Apr 2007  
 
c*** Temporary: There are redundant inputs here to allow some checking of
c               calculations; later use only cartesian (Fritshe et al ) or
c               trig (Kedar et al)  rwk 070515
c
c      Input
c
c        ionsrc :  Source of TEC values (must be IONEX file in current code)
c        jdobs  :  PEP Julian Day of observation    
c        tobs   :  Seconds of day observation  
c        lonr   :  Longitude in radians
c        latr   :  Geodetic latitude in radians
c        height :  Site geodetic height (km)
c        prn    :  Satellite PRN number (for warning) 
c        sitpos :  Station coordinates (xyz, km)
c        rhat   :  Unit vector in direction of satellite    
c        elev   :  Elevation angle of satellite
c        azim   :  Azimuth of satellite  
c        pnsmat :  Rotation matrix from E-fixed to inertial
c        freq1, freq2 : Frequencies in Hz of the L1 and L2 signals    
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
      include '../includes/model.h'   

      integer*4 jdobs,doy,iyr,prn,latindx,lonindx,it,i
    
      real*4 lond,latd
          
      real*8 tobs,lonr,latr,height,sitpos(3),rhat(3),rhatm(3),elev,azim 
     .     , pnsmat(3,3),freql1,freql2,eps,b0(3),htion,re,rm,rmag,bg
     .     , nmax,radian,obstime,pi,vlight,p,q,dt,tec1,tec2,tec
     .     , theta,thetam,thetamp,phim,magmat(3,3),sbeta,cbeta
     .     , cdelta,sdelta,stheta,ctheta,sphi,cphi,sphim,cphim
     .     , emhat(3),nmhat(3),umhat(3),emhatg(3),nmhatg(3),umhatg(3)
     .     , ycof,zcof,azimm,lambda(2),dcap,geocn(3),rotmat(3,3)
     .     , grpdel3(2),phsdel3(2),grpdel(2),phsdel(2)
     .     , pnsmatt(3,3),rhate(3),twopi 
         
c     for debug only 
      integer*4 debug   
      real*8 sitposm(3),azim1(maxsat),elev1(maxsat)
     .     , nhat1(3),ehat1(3),uhat1(3),nhat(3),uhat(3)
     .     , tmpl1,tmpl2

  
      character*256 message  
      character*4 ionsrc

      logical first,found,start_warning,end_warning
                      
c       Functions                    
      real*8 dot



      data vlight/2.99792458D+08/,twopi/6.283185307179586d0/

                          
c     eps is a small number tolerance to allow use of the end values
c     currently set to be 90s if times are in days
      data eps/.001/


c     Initialize 'first' for setup call
      data first /.true./,start_warning/.false./,end_warning/.false./      

      save first, radian, thetam, start_warning, end_warning
              
      debug = 1
 
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
       
c..   At first call, get the geomagnetic coordinates of the station and
c      the magnetic field vector (assumed constant)
 
      if( first ) then       
c        print *,'IONDEL DEBUG ionsrc first  ',ionsrc,first



c---debug dump maps  
cd        print *,'nilon nilat ntion ',nilon,nilat,ntion
cd        write(*,'(a,f7.3,5(/,20f5.1))') 
cd     .       't1 lat1 lon ',ion_time(1),(ion_val(i,1,1),i=1,nilon) 
cd        write(*,'(/,a,f7.3,5(/,20f5.1))') 
cd     .       't5 lat5 lon ',ion_time(5),(ion_val(i,5,5),i=1,nilon)
cd        write(*,'(/,a,f7.3,5(/,20f5.1))') 
cd     .       'tlast latlast lon ',ion_time(ntion)
cd     .        ,(ion_val(i,nilat,ntion),i=1,nilon)
c-----------
        pi = 4.d0*datan(1.d0)   
c       degrees to radians
        radian = datan(1.d0)/45.d0
c       get unit vector of site in the Earth-centered geodetic frame
        rmag = dsqrt(sitpos(1)**2+sitpos(2)**2+sitpos(3)**2)
c       define the axes (E N U) of the geomagnetic frame at the site
        emhat(1) = 1.0d0 
        emhat(2) = 0.0d0 
        emhat(3) = 0.0d0
        nmhat(1) = 0.0d0 
        nmhat(2) = 1.0d0 
        nmhat(3) = 0.0d0
        umhat(1) = 0.0d0 
        umhat(2) = 0.0d0 
        umhat(3) = 1.0d0 
c       rotation matrix from geomagnetic to geodetic coordinates 
        sbeta = dsin( 291.d0 * radian )
        cbeta = dcos( 291.d0 * radian )
        sdelta = dsin( 11.5d0 * radian )
        cdelta = dcos( 11.5d0 * radian )  
c       Bassiri & Hajj Eqn 15 
        magmat(1,1) = cdelta*cbeta 
        magmat(1,2) = cdelta*sbeta
        magmat(1,3) = -sdelta
        magmat(2,1) = -sbeta
        magmat(2,2) =  cbeta
        magmat(2,3) = 0.d0
        magmat(3,1) = sdelta*cbeta
        magmat(3,2) =  sdelta*sbeta 
        magmat(3,3) = cdelta  
c        print *,'magmat ',magmat     
c       this for debug only: rotate the station to the magnetic frame
        call matmpy(magmat,sitpos,sitposm,3,3,1)
c        print *,'sitposm ',sitposm                
c       get the magnetic frame unit vectors in the local geodetic frame               
        call matmpy (magmat,emhat,emhatg,3,3,1) 
        call matmpy (magmat,nmhat,nmhatg,3,3,1) 
        call matmpy (magmat,umhat,umhatg,3,3,1)   
c        print *,'emhat nmhat umhat ',emhat,nmhat,umhat
c        print *,'emhatg nmhatg unhatg ',emhatg,nmhatg,umhatg
c       get the magnetic frame unit vectors in the E-centered geodetic frame
        call rotate_geod (emhatg, emhat,'NEU','XYZ',sitpos,geocn,rotmat) 
        call rotate_geod (nmhatg, nmhat,'NEU','XYZ',sitpos,geocn,rotmat) 
        call rotate_geod (umhatg, umhat,'NEU','XYZ',sitpos,geocn,rotmat) 
c        print *,'geocentric emhat nmhat umhat ',emhat,nmhat,umhat
c       magnetic co-latitude 
        stheta = dsin(90.d0*radian-latr)
        ctheta = dcos(90.d0*radian-latr)
        sphi = dsin(lonr)
        cphi = dcos(lonr)
        thetam = dacos(  sdelta*cbeta*stheta*cphi + 
     .           sdelta*sbeta*stheta*sphi+ cdelta*ctheta )
        dcap =cdelta*(cbeta*stheta*cphi+sbeta*stheta*sphi)-sdelta*ctheta
        phim = datan( (-sbeta*stheta*cphi+cbeta*stheta*sphi)/dcap )
        sphim = dsin(phim)
        cphim = dcos(phim)   
c        print *,'sdelta cbeta stheta cphi sbeta sphi cdelta ctheta '
c     .         , sdelta,cbeta,stheta,cphi,sbeta,sphi,cdelta,ctheta 
c        print *,' dcap thetam phim (deg/rad)'
c     .           , dcap,thetam,thetam/radian,phim,phim/radian
c       magnitude of B at the equator (Tesla)
        bg = 3.12d-5    
        re = 6378.13655d0 
        ycof = bg*(re/rmag)**3*dsin(thetam)                      
        zcof = -2.d0* bg*(re/rmag)**3*dcos(thetam)
        do i=1,3
          b0(i) = ycof*nmhat(i) + zcof*(umhat(i))
        enddo
c        print *,'END OF FIRST CALL  ycof zcof B ',ycof,zcof,b0
        first = .false.
        return
      endif

c     print *,'IONDEL DEBUG '
             
c   Height of ionosphere (km)

      htion = 400.d0
           
c.. Initialize return variables to zero
      phsdel(1) = 0.d0
      phsdel(2) = 0.d0
      grpdel(1) = 0.d0
      grpdel(2) = 0.d0
    
c.. Get the observation time

      call dayjul( jdobs,iyr,doy)         
      obstime = dfloat(doy) + tobs/86400.d0 
c      print *,'obstime ion_time1 ion_timelast '
c     .     ,obstime,ion_time(1),ion_time(ntion)
                                            

c.. See if the requested time is outside the available array of values

      if( obstime.lt.(ion_time(1)-eps) ) then
        if( .not.start_warning ) then
          write(message,'(a,f8.4,a,f8.4)') 'Obs epoch '
     .        ,obstime,' before start of IONEX file ',ion_time(1)  
          call report_stat('WARNING','MODEL','iondel',' ',message,0)
          start_warning = .true.
        endif
        return
      elseif( obstime.gt.(ion_time(ntion)+eps) ) then
        if( .not.end_warning ) then 
          write(message,'(a,f8.4,a,f8.4)') 'Obs epoch '
     .     ,obstime,' after end of IONEX file ',ion_time(ntion)  
          call report_stat('WARNING','MODEL','iondel',' ',message,0) 
          end_warning = .true.
        endif
        return
      endif
        
c      print *,'jdobs tobs ',jdobs,tobs   
c      print *,'latr lonr height sitpos ',latr,lonr,height,sitpos
c      print *,' prn rhat elev azim ',prn,rhat,elev/radian,azim/radian
c      print *,'  azim radians ',azim


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
      if( .not.found ) then
         write(message,'(a,f8.4,a,2f8.4,a)') 'Obs epoch '
     .     ,obstime,' not found in IONEX file ('
     .     ,ion_time(1),ion_time(ntion),')'
          call report_stat('WARNING','MODEL','iondel',' ',message,0) 
          return
      endif
        
c.. Get lat and lon in degrees to interpolate from the ion maps  
c   
c     Shift the longitude by the time since the first of the two bracketing
c     maps to compensate for the correlation between the ionsophere and the
c     Sun's position [Schaer, Gurtner, and Feltens, IONEX Version 1, Proceedings
c     of the IGS AC Workshop, Darmstadt, 1998]
         
      latd = latr/radian 
      lond = lonr/radian   
c      print *,'latr latd lonr lond ',latr,latd,lonr,lond
      lond = lonr/radian + (obstime -ion_time(it))*360.    
      if( lond.gt.360.d0 ) lond = lond - 360.d0
      if( lond.gt.180.d0 ) lond = lond - 360.d0
          
c      print *,' obstime it ion_time(it) londp '
c     .    ,obstime,it,ion_time(it),lond
                                             
c.. Get the indices for the "lower-left" corners of the spatial box
c   (where "lower" is first in the array but usually higher in latitude
c   if the entries are north to south; longitude always west to east)
                    
      latindx = int(aint((latd-ilat1)/dilat)) + 1
      lonindx = int(aint((lond-ilon1)/dilon)) + 1
c      print *,'latd ilat1 dilat latindx ',latd,ilat1,dilat,latindx
c      print *,'lond ilon1 dilon lonindx ',lond,ilon1,dilon,lonindx

c.. Interpolate to the requested lon/lat at times before and after the observation
         
      p = amod(abs(lond),abs(dilon))/abs(dilon)
      if( lond.lt.0 ) p = 1. - p
      q = amod(abs(latd),abs(dilat))/abs(dilat)
      if( latd.lt.0 ) q = 1. -q
      tec1 = 
     .    (1.-p)*(1.-q)*ion_val(lonindx,latindx,it) +
     .         p*(1.-q)*ion_val(lonindx+1,latindx,it) +
     .         q*(1.-p)*ion_val(lonindx,latindx+1,it) +
     .              p*q*ion_val(lonindx+1,latindx+1,it)                                                                               
      tec2 = 
     .    (1.-p)*(1.-q)*ion_val(lonindx,latindx,it+1) +
     .         p*(1.-q)*ion_val(lonindx+1,latindx,it+1) +
     .         q*(1.-p)*ion_val(lonindx,latindx+1,it+1) +
     .              p*q*ion_val(lonindx+1,latindx+1,it+1)                                                                               
           
c      print *,' p q tec1 tec2 ',p,q,tec1,tec2

c.. Interpolate to the right time
                                       
      dt = ion_time(it+1) - ion_time(it)                    
      tec = tec1 * (ion_time(it+1)-obstime) / dt +
     .       tec2 * (obstime-ion_time(it)) / dt
                       
c       print *,'obstime it ion_time1 ion_time2 dt tec '
c     .        , obstime,it,ion_time(it),ion_time(it+1),dt,tec


c.. Get the unit vector in the direction of the satellite in the geomagnetic system
                          
c     rotate from inertial to E-fixed     
c      print *,'PNSMAT ',pnsmat
      call transp(pnsmat,pnsmatt,3,3)        
c      print *,'PNSMAT-t ',pnsmatt
      call matmpy (pnsmatt,rhat,rhate,3,3,1)
c     then to geomagnetic
      call matmpy (magmat,rhate,rhatm,3,3,1)  
c      print *,' rhat rhate rhatm ',pnsmatt,rhat,rhate,rhatm
      
        
c ...  or get the azimuth in the geomagnetic system (Bassiri & Hajj) (elevation same as geodetic)
                              

c** For debug, do this first with vectors    
c       first redo the geodetic values
c      nhat(1) = 1.d0
c      nhat(2) = 0.d0
c      nhat(3) = 0.d0
c      uhat(1) = 0.d0
c      uhat(2) = 0.d0
c      uhat(3) = 1.d0 
c      call az_elev( .true.,nhat,uhat,pnsmat,sitpos,rhat
c     .            , twopi,1,azim1,elev1,nhat1,ehat1,uhat1 )
c      print *,'Recalc geod sitpos rhat azim1 elev1 '
c     .       , sitpos,rhat,azim1(1),elev1(1)
c      call az_elev( .true.,nmhat,umhat,pnsmat,sitposm,rhatm
c     .            , twopi,1,azim1,elev1,nhat1,ehat1,uhat1 )
c      print *,'Calc mag sitposm rhatm azim1, elev1 '
c     .            , sitposm,rhatm,azim1(1),elev1(1)
c      theta = 90.d0*radian - latr                              
      dcap = cdelta*(cbeta*stheta*cphi +
     .               sbeta*stheta*sphi ) - sdelta*ctheta
      phim = datan(  (-sbeta*stheta*cphi
     .                +cbeta*stheta*sphi)/dcap )
      azimm = azim + dacos( sphi*sphim*cdelta*cbeta +
     .                      cphi*cphim*cbeta +
     .                      sphi*cphim*sbeta -
     .                      cphi*sphim*cdelta*sbeta )
c      print *,'Trig theta thetam phi phim azim azimm '
c     .       ,  theta,thetam,lonr,phim,azim,azimm
                       
c.. Calculate the second order ionspheric group and phase delays
     
      lambda(1) = vlight/freql1  
      lambda(2) = vlight/freql2
  
c --------Kedar

      thetamp = thetam - htion/(re*dsin(elev))*dcos(azimm)*dcos(elev)
c      print *,'Kedar thetam thetamp ',thetam,thetamp  
      rm = re + htion
      do i=1,2
        grpdel(i) = 2.61d-18*(lambda(i)*re/rm)**3 *
     .      ( dsin(thetamp)*dcos(elev)*dcos(azimm) 
     .      - 2.d0*dcos(thetamp)*dsin(elev) ) 
     .          * tec*1.d16
        phsdel(i) = -0.5d0 * grpdel(i)
      enddo    
      tmpl1 = phsdel(1)
      tmpl2 = phsdel(2)
c      print *,'Kedar grpdel phsdel ',grpdel,phsdel

c ------ Fritsche et al formulation
                                       
      do i=1,2
        grpdel(i)=(7527.d0/vlight**2)*(lambda(i)**3)*dot( b0,rhatm)
c          1 TEC = 1.e16 electrons/m**2
     .     *tec*1.d16
        phsdel(i) = -0.5d0 * grpdel(i)
      enddo      

c alternative formulation, using B-vector and k-vector in inertial frame
c     note rotate_geod wants N E U (not E N U)
      bmneu(1) = bg*sthetam
      bmneu(2) = 0.d0
      bmneu(3) = -2.d0*bg*cthetam
c     rotate the local magnetic vector to E-center magnetic at the site
      call rotate_geod(bmneu,bmxyz,'NEU','XYZ',sitposm,geocn,rotmat) 
c     rotate the B-vector from the magnetic to the geodetic system
      call matmpy(bmxyz,bmxyzg,magmat)
c     rotate the B-vector from the geodetic system to the inertial system
      call matmpy(bmxyzg,bmxygi,pnsmat)
      grpdel2 =           


c.. Calculate the third order group and phase delays (Fritsche et al only)

c      nmax = 14.d12/3.27d18 * tec *1.d16
c      do i=1,2
c        grpdel3(i) = 2437.d0*((lambda(i)/vlight)**4)*nmax
c     .                    * tec*1.d16
c        phsdel3(i) = -grpdel3(i)/3.d0 
c      enddo
c      print *,'grpdel3 phsdel3 ',grpdel3,phsdel3

      if(debug.ge.1 ) then
        write(*,'(a,f6.2,a,i2,f5.1,2f8.2,4f6.1)') 
     .    'T PN TEC el az  L1(mm)  Kedar   Fritche: '
     .   , obstime,'  P ',prn,tec,elev/radian,azimm/radian
     .   , tmpl1*1.d3,tmpl2*1.d3,phsdel(1)*1.d3,phsdel(2)*1.d3
      endif

      return
      end
