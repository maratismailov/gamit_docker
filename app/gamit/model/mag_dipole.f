      Subroutine MAG_DIPOLE( jdobs,latr,lonr,slat,slon,htion,vecdipxyz 
     .                     , magsize,debug )

C   Calculate the Earth's magnetic field for 2nd/3rd order ionospheric  corrections using
c   a simple dipole (for comparison with better field from IGRF).  This routine activated
c   by undocumented 'DIPOLE' after ion model in the model batch file.  Routine created by
c   R. King from E. Petrie code put originally into iondel.f.  081107/090429/090504

      implicit none

      integer iyr,idoy,jdobs,debug

      real*8 azim,azimm,phim,cphim,sphim,cbeta,sbeta,cphi,sphi
     .     , cdelta,sdelta,thetam,ctheta,stheta,slat,thetamppt
     .     , phimppt,phim2,lonr,latr,slon
     .     , sthetappt,cthetappt,sphippt,cphippt,sphimppt,cphimppt
     .     , dcap,dcapppt,bodipole
      real*8 mag_vecdip(3),diprot(3,3),mag_vecdipxyz(3),vecdipxyz(3)
     .     , dipNEDtoXYZ(3,3),sazimm,cazimm,htion,radian,re,magsize

C     colat, longitude of magnetic dipole
      Real*8 magcolat,maglong
                  
c       degrees to radians and Earth radius
        radian = datan(1.d0)/45.d0
        re = 6378.13655d0    

C   3.12d-5 is strength of magnetic field at magnetic equator at Earth's surface.
C        Bassiri & Hajj (1993)
         bodipole = (3.12d-5)*(re/(re+htion))**3   

c  Set the orientation of the nominal magnetic field 
c       (Bassiri & Hajj, p.284)

C Calculate 
       if(debug.ge.3) then
          print *,'using magnetic dipole'
        endif
C Assume constant relationship between magnetic and geographic.
C        sbeta = dsin( 291.d0 * radian )
C        cbeta = dcos( 291.d0 * radian )
C        sdelta = dsin( 11.5d0 * radian )
C        cdelta = dcos( 11.5d0 * radian )  
C Use relationship for day 302 2003 from ION files
C        sbeta = dsin( 288.01d0 * radian )
C        cbeta = dcos( 288.01d0 * radian )
C        sdelta = dsin( 10.41d0 * radian )
C        cdelta = dcos( 10.41d0 * radian ) 

C Geographic-geomagnetic relationships from ION files:
        call dayjul( jdobs,iyr,idoy)
        if (iyr.eq.1998) then 
		magcolat=90.d0-79.41d0
		maglong=-71.64d0+360.d0
        elseif (iyr.eq.1999) then 
		magcolat=90.d0-79.45d0
		maglong=-71.71d0+360.d0
        elseif (iyr.eq.2000) then 
		magcolat=90.d0-79.48d0
		maglong=-71.78d0+360.d0
        elseif (iyr.eq.2001) then 
		magcolat=90.d0-79.52d0
		maglong=-71.85d0+360.d0
        elseif (iyr.eq.2002) then 
		magcolat=90.d0-79.55d0
		maglong=-71.92d0+360.d0
        elseif (iyr.eq.2003) then 
		magcolat=90.d0-79.59d0
		maglong=-71.59d0+360.d0
        elseif (iyr.eq.2004) then 
		magcolat=90.d0-79.62d0
		maglong=-71.06d0+360.d0
        elseif (iyr.eq.2005) then 
		magcolat=90.d0-79.77d0
		maglong=-71.80d0+360.d0
        elseif (iyr.eq.2006) then 
		magcolat=90.d0-79.81d0
		maglong=-71.84d0+360.d0
        elseif (iyr.eq.2007) then 
		magcolat=90.d0-79.85d0
		maglong=-71.88d0+360.d0
        elseif (iyr.eq.2008) then 
		magcolat=90.d0-79.89d0
		maglong=-71.92d0+360.d0
	else
C Set all others as 2007 value for now
		magcolat=90.d0-79.85d0
		maglong=-71.88d0+360.d0
	endif
C year day  lat		long
C 1995 365, 0.000000  0.000000 
C 1996 365, 0.000000  0.000000 
C 1997 365, 0.000000  0.000000 
C 1998 365, 79.410000  -71.640000 
C 1999 365, 79.450000  -71.710000 
C 2000 365, 79.480000  -71.780000 
C 2001 365, 79.520000  -71.850000 
C 2002 365, 79.550000  -71.920000 
C 2003 365, 79.590000  -71.990000 
C 2004 365, 79.620000  -71.060000 
C 2005 365, 79.770000  -71.800000 
C 2006 365, 79.810000  -71.840000 
C 2007 365, 79.850000  -71.880000 
        sbeta = dsin( maglong * radian )
        cbeta = dcos( maglong * radian )
        sdelta = dsin( magcolat * radian )
        cdelta = dcos( magcolat * radian ) 
            
      if( debug.ge.2 ) then
         print *,'iyr maglong magcolat ',iyr,maglong,magcolat
      endif

c  Get the geomagnetic colatitude and longitude of the site 
c       (Bassiri & Hajj Eqs 22 & 23)

        stheta = dsin(90.d0*radian-latr)
        ctheta = dcos(90.d0*radian-latr)
        sphi = dsin(lonr)
        cphi = dcos(lonr)
C Magnetic colatitude (of the site)
        thetam = dacos(  sdelta*cbeta*stheta*cphi + 
     .           sdelta*sbeta*stheta*sphi+ cdelta*ctheta )
        dcap =cdelta*(cbeta*stheta*cphi+sbeta*stheta*sphi)-sdelta*ctheta
C Magnetic longitude (of the site )
        phim = datan2( -sbeta*stheta*cphi+cbeta*stheta*sphi,dcap )
        sphim = dsin(phim)
        cphim = dcos(phim)   
        if( debug.ge.2 ) then
          print *,'sdelta cbeta stheta cphi sbeta sphi cdelta ctheta '
     .         , sdelta,cbeta,stheta,cphi,sbeta,sphi,cdelta,ctheta 
          print *,' dcap thetam phim (rad/deg) '
     .           , dcap,thetam,thetam/radian,phim,phim/radian
          print *,'latr, lonr, latr/deg, lonr/deg, radian'
     .           ,latr,lonr,latr/radian,lonr/radian,radian
      endif
      

C Obtain magnetic field vector

C    Method 2: Field calculated using adapted ppt colatitude & rotations
c    Get the geomagnetic colatitude and longitude of the pierce point 
c   (Bassiri & Hajj Eqs 22 & 23)
C 
C -------
C Debug - setup slat slon for test cases
C        slat = 79.59d0-90.d0
C        slon = -71.59d0+360.d0
C -----
          sthetappt = dsin(90.d0*radian-slat*radian)
          cthetappt = dcos(90.d0*radian-slat*radian)
          sphippt = dsin(slon*radian)
          cphippt = dcos(slon*radian)
C Magnetic colatitude (of the pierce point) in radians
          thetamppt = dacos(  sdelta*cbeta*sthetappt*cphippt + 
     .           sdelta*sbeta*sthetappt*sphippt+ cdelta*cthetappt )
          dcapppt =cdelta*(cbeta*sthetappt*cphippt+sbeta*sthetappt
     .    *sphippt)-sdelta*cthetappt
C Debug
C          Print*, 'dcapppt: ',dcapppt,'sphippt',sphippt
C          Print*, 'cdelta',cdelta,'cbeta',cbeta,'sthetappt',sthetappt
C          Print*, 'cphippt',cphippt,'sbeta',sbeta,'stheta',stheta
C          Print*, 'sdelta',sdelta,'cthetappt',cthetappt
C         thetamppt = 135.d0*radian
C Magnetic longitude (of the pierce point ) in radians
          phimppt = datan2( -sbeta*sthetappt*cphippt+cbeta*sthetappt*
     .           sphippt,dcapppt )
          sphimppt = dsin(phimppt)
          cphimppt = dcos(phimppt)   
C Debug 
C        phimppt = 180.d0*radian
         if (debug.ge.3) print*,'thetamppt',thetamppt,'phimppt',phimppt
C Magnetic vector in East North Up (magnetic frame)
C        mag_vecdip(1) = 0
C        mag_vecdip(2) = bodipole*dsin(thetamppt)
C        mag_vecdip(3) = -2.0d0*bodipole*dcos(thetamppt)

C Magnetic vector in North East Down (magnetic frame, at pierce point)
        mag_vecdip(1) = bodipole*dsin(thetamppt)
        mag_vecdip(2) = 0.d0
        mag_vecdip(3) = 2.0d0*bodipole*dcos(thetamppt)
        if(debug.ge.3) then        
           print *,'MODEL\mag_dipole bodipole thetamppt '
     .           ,bodipole,thetamppt
           Print*, 'MODEL\mag_dipole magvecdip: ',mag_vecdip
        endif
C Change to magnetic XYZ
C    now get rotation matrix
        Call vec_xyz((90.d0-(thetamppt/radian)), phimppt/radian, 
     .             dipNEDtoXYZ)
C     Rotate dipole vector to X Y Z 
        Call matmpy(dipNEDtoXYZ,mag_vecdip,mag_vecdipxyz,3,3,1)
        if(debug.ge.3) then
          Print*, 'MODEL\mag_dipole dipNEDtoXYZ: ',dipNEDtoXYZ
          Print*, 'MODEL\mag_dipole magvecdipxyz: ',mag_vecdipxyz
        endif
C Now need to rotate back to geographic coords

        diprot(1,1) = cbeta*cdelta
        diprot(1,2) = -sbeta
        diprot(1,3) = cbeta*sdelta
        diprot(2,1) = sbeta*cdelta
        diprot(2,2) = cbeta
        diprot(2,3) = sbeta*sdelta
        diprot(3,1) = -sdelta
        diprot(3,2) = 0.d0
        diprot(3,3) = cdelta

          Call matmpy(diprot,mag_vecdipxyz,vecdipxyz,3,3,1)  
          magsize=dsqrt(vecdipxyz(1)**2+vecdipxyz(2)**2+vecdipxyz(3)**2)
          if(debug.ge.3) then
            Print*, 'MODEL\mag_dipole diprot : ',diprot
          endif
          if(debug.ge.2) then
            Print*, 'MODEL\mag_dipole vecdipxyz: ',vecdipxyz
          endif               

C -------------------------------------------------------------------

      return
      end
