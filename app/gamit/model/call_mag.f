C Program to call mag.f and get an X,Y,Z and total magnetic field
C and convert the components to fit in a geocentric cartesian XYZ system

      Subroutine call_mag(DATE,XLT,XLN,mag_xyz,F,debug,htion,magfield)

C	In: 	DATE 	 	date in decimal years
C		XLT,XLN		geocentric latitude, longitude (decimal degrees) at which magnetic field is required
C               htion 		height of the ionospheric thin layer /km
C	Out: 	mag_xyz(3) 	the X, Y, Z, components of the magnetic field vector
C				where X,Y,X are the cartesian geocentric axes 
C Declare variables
      implicit none
      Double precision DATE,ALT,XLT,XLN,X,Y,Z,F,mag_ned(3),mag_xyz(3),
     . rot_mat(3,3), rot_mat2(3,3), mag_xyza(3),htion
      Character*30 FNM
      Character*20 NAME
      Character*6 magfield
      Integer*4 ITYPE, IOPT, IDM,debug
C Define variables
C Filename for output from mag.f - ' ' means output from mag.f should go to the screen
C      FNM='testigrf'
      FNM=' '
C			ITYPE - 1:geodetic coords (spheroidal earth)
C				2:geocentric coords (spherical earth)
      ITYPE=2
C			IOPT - 	1: want value(s) at locations and dates
C				2: want value yearly intervals at one loc
C				3: want values on a lat/long grid at 1 date
      IOPT=1
CC			IDM - lat long format 	1: degrees and minutes
C						2: decimal degrees
      IDM=2
C			DATE - time in decimal years AD - now input through MODEL
C      DATE=2007.d0

C			ALT - height/km if ITYPE=1
C			    - radial dist from earth centre/km if ITYPE=2
C htion now defined in iondel.f and passed through
      ALT=6371.0d0+htion
C 6371 km is for the radius of the spherical earth used for the piercepoint calculation
C 450 km is the value for the height of the pierce point above the sphere

C			XLT - Latitude/decimal degrees
C			XLN - Longitude/decimal degrees
C      XLT=30.5d0
C      XLN=0.0d0
C      XLT=-90.0d0
C      XLN=-180.0d0
C			NAME - SITEepochSAT
      Name='TESTtestingPRNno'
C 
C      Print*,'MODEL\call_mag: Call mag.f - date is: ',DATE
c      IGRF10 geomagnetic mode
      if ( magfield .eq.'IGRF10') then
       call mag (FNM,ITYPE,IOPT,IDM,DATE,ALT,XLT,XLN,NAME,X,Y,Z,F,
     .           debug)
      elseif ( magfield .eq.'IGRF11') then
c         IGRF11 geomagnetic mode
       call mag11 (FNM,ITYPE,IOPT,IDM,DATE,ALT,XLT,XLN,NAME,X,Y,Z,F,
     .             debug)         
      elseif ( magfield .eq.'IGRF12') then
c         IGRF12 geomagnetic mode
       call mag12 (FNM,ITYPE,IOPT,IDM,DATE,ALT,XLT,XLN,NAME,X,Y,Z,F,
     .             debug)
      elseif ( magfield .eq.'IGRF13') then
c         IGRF13 geomagnetic mode (post 2020)
       call mag13 (FNM,ITYPE,IOPT,IDM,DATE,ALT,XLT,XLN,NAME,X,Y,Z,F,
     .             debug)

      else
       call report_stat('FATAL','MODEL','call_mag.f',magfield
     .      , 'Unrecognized magnetic field model:',0)      
      endif
C
C X is north component,Y is east component, Z is vertical component (+ve down)
C F is total intensity.
C So the magnetic vectors are in a (usually) geocentric coordinate system
C equivalent to a conversion from WGS84 spheroid and datum, so not geomagnetic coordinate system
C Check result
C      Print*, 'MODEL\call_mag N,E,down,magnitude X: ',x,'   Y: ',y,
C     .'   Z: ',z,'    F: ',f
C Convert from the nanoTeslas supplied by the IGRF code to Teslas 
      mag_ned(1)=x*1.d-9
      mag_ned(2)=y*1.d-9
      mag_ned(3)=z*1.d-9
      F=F*1.d-9
      mag_xyz(1)=0.d0
      mag_xyz(2)=0.d0
      mag_xyz(3)=0.d0  
C Debug
      If (debug.ge.3) then
        Print*, 'mag_ned', mag_ned
C       Print*, 'mag_xyz', mag_xyz  
      endif  
C Convert components to geocentric Cartesian XYZ coordinate system
C Calculate rotation vector
      Call vec_xyz(XLT,XLN,rot_mat)
C Multiply NED components by 1 rotation matrix
      Call MATMPY(rot_mat,mag_ned,mag_xyz,3,3,1)
C Debug
        if(debug.ge.2) then
      Print*, 'MODEL\call_mag cartesian x:',mag_xyz(1),'y:',mag_xyz(2)
      Print*,'z:',mag_xyz(3),
     .'nf:',sqrt((mag_xyz(1)**2)+(mag_xyz(2)**2)+(mag_xyz(3)**2))
        endif
      End
